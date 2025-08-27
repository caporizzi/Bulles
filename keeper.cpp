#include <mpi.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

// -----------------------------
// Union–Find (Disjoint Set Union)
// -----------------------------
// Used in Connected Component Labeling (CCL) to manage equivalences between labels.
// - Each label starts as its own "set" (parent[i] = i).
// - find(x): returns the root representative of x, with path compression to flatten chains.
// - unite(a,b): merges the sets containing a and b (so both labels are treated as the same).
// This allows us to efficiently handle cases where two different provisional labels
// actually belong to the same connected component.
// -----------------------------

struct UF {
    std::vector<int> parent;
    explicit UF(size_t n = 0) : parent(n + 1) {
        for (size_t i = 0; i <= n; ++i) parent[i] = (int)i;
    }
    void reset(size_t n) {
        parent.resize(n + 1);
        for (size_t i = 0; i <= n; ++i) parent[i] = (int)i;
    }
    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a != b) parent[b] = a;
    }
};

// -----------------------------
// Connected Component Labeling (CCL)
// -----------------------------
// This function runs CCL on a *local sub-image* (the rows owned by one MPI rank).
//
// Process:
//  1) First pass (local, row-by-row):
//     - Scan each pixel left-to-right, top-to-bottom.
//     - If the pixel is background → label = 0.
//     - Otherwise, check its north and west neighbors (4-connectivity):
//         * Both neighbors 0 → assign a new temporary label.
//         * One neighbor labeled → copy that label.
//         * Both neighbors labeled (and different) → assign one label,
//           and record an equivalence between the two using Union–Find.
//     - This creates *provisional* labels that are locally consistent,
//       but not yet compact or globally unique.
//
//  
//
// Note:
//  - This function only handles the local sub-image.
//  - Later steps in main() will take care of merging components that
//    cross MPI rank boundaries and making the labels globally unique.
// -----------------------------
// ===== in localCCL: return UF so we don’t lose equivalences =====
int localCCL(const std::vector<int> &mask, int nx, int ny,
             std::vector<int> &labels, UF &uf) {
    labels.assign(nx * ny, 0);
    auto idx = [nx](int x, int y) { return y * nx + x; };

    uf.reset(1);
    int nextLabel = 1;

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            if (mask[idx(x, y)] == 0) continue;

            int north = (y > 0) ? labels[idx(x, y - 1)] : 0;
            int west  = (x > 0) ? labels[idx(x - 1, y)] : 0;

            if (north == 0 && west == 0) {
                labels[idx(x, y)] = nextLabel;
                if ((int)uf.parent.size() <= nextLabel) uf.parent.resize(nextLabel + 1);
                uf.parent[nextLabel] = nextLabel;
                ++nextLabel;
            } else if (north != 0 && west == 0) {
                labels[idx(x, y)] = north;
            } else if (north == 0 && west != 0) {
                labels[idx(x, y)] = west;
            } else {
                int a = std::min(north, west);
                int b = std::max(north, west);
                labels[idx(x, y)] = a;
                uf.unite(a, b);  // keep local equivalence
            }
        }
    }
    return nextLabel - 1;
}


// -----------------------------
// build_counts_displs
// -----------------------------
// This function computes how to split a 2D image among multiple MPI ranks
// for parallel processing using MPI_Scatterv and MPI_Gatherv.
//
// Inputs:
//   nx   - number of columns in the image
//   ny   - number of rows in the image
//   size - number of MPI ranks (processes)
//
// Outputs (by reference):
//   counts       - number of elements (pixels) assigned to each rank
//   displs       - starting index of each rank's block in the flattened array
//   rowsPerRank  - number of rows assigned to each rank
//
// How it works:
// 1) Divide rows evenly: each rank gets at least ny / size rows (base)
// 2) Distribute any remaining rows (ny % size) to the first few ranks
//    so that all rows are processed
// 3) Compute counts = rowsPerRank[r] * nx (number of pixels per rank)
// 4) Compute displs = starting offset in the 1D flattened array for each rank
//    This ensures MPI knows where each rank's data begins
//
// Example:
//   nx = 4, ny = 10, size = 3
//   base = 10 / 3 = 3, remainder = 1
//   ranks get rows: [4, 3, 3]
//   counts = [16, 12, 12], displs = [0, 16, 28]
//
// Why it's useful:
// - Ensures each MPI rank gets a contiguous block of rows to work on
// - Handles uneven row division automatically
// - Required for MPI_Scatterv/Gatherv to distribute and gather data correctly
// - Allows independent processing per rank and proper assembly of the full image

void build_counts_displs(int nx, int ny, int size,
                         std::vector<int> &counts, std::vector<int> &displs,
                         std::vector<int> &rowsPerRank) {
    counts.resize(size);
    displs.resize(size);
    rowsPerRank.resize(size);

    int base = ny / size;
    int rem  = ny % size;
    int offset = 0;
    for (int r = 0; r < size; ++r) {
        int rows = base + (r < rem ? 1 : 0);
        rowsPerRank[r] = rows;
        counts[r] = rows * nx;
        displs[r] = offset;
        offset += counts[r];
    }
}

// -----------------------------
// Main
// -----------------------------
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <input.dat> <output.dat>\n";
        }
        MPI_Finalize();
        return 1;
    }

    const std::string inFile  = argv[1];
    const std::string outFile = argv[2];

    int nx = 0, ny = 0, nz = 0;
    std::vector<float> fullData;  // only on rank 0

    // -----------------------------
    // 1) Rank 0 reads the whole image
    // -----------------------------
    //
    // Only rank 0 (the master) reads the input file to ensure a single, consistent copy.
    // Even though processing is parallel, the file exists as a single binary stream.
    //
    // Steps:
    // 1) Open the file and read the header: nx (cols), ny (rows), nz (slices, ignored).
    // 2) Allocate `fullData` for all pixels and read the image data.
    // 3) Other ranks will later receive their portions via MPI_Scatterv.
    //
    // Why:
    // - Simplifies I/O: avoids multiple ranks reading the same file.
    // - Enables safe parallel processing by scattering rows to each rank afterwards.

    if (rank == 0) {
        std::ifstream in(inFile, std::ios::binary);
        if (!in) {
            std::cerr << "Failed to open input file\n";
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
        in.read((char *)&nx, sizeof(nx));
        in.read((char *)&ny, sizeof(ny));
        in.read((char *)&nz, sizeof(nz));
        if (!in) {
            std::cerr << "Failed to read header\n";
            MPI_Abort(MPI_COMM_WORLD, 3);
        }
        fullData.resize((size_t)nx * ny);
        in.read((char *)fullData.data(), sizeof(float) * nx * ny);
        if (!in) {
            std::cerr << "Failed to read image data\n";
            MPI_Abort(MPI_COMM_WORLD, 4);
        }
    }

    // -----------------------------
    // 2) Broadcast dimensions
    // -----------------------------
    //
    // After rank 0 reads the image, all ranks need to know its size to process their part.
    //
    // 1) MPI_Bcast sends nx (columns) and ny (rows) from rank 0 to all other ranks.
    //    - This ensures every rank knows the dimensions of the image.
    //
    // 2) Use `build_counts_displs` to calculate:
    //    - `rowsPerRank`: how many rows each rank will process
    //    - `counts`: number of pixels per rank (rows * nx)
    //    - `displs`: starting index of each rank's rows in the flattened array
    //
    // 3) `localNy` stores the number of rows assigned to this rank.
    //    - If there are more ranks than rows, some ranks may get 0 rows (safe case).
    //
    // Purpose:
    // - This prepares the data for scattering the image rows across ranks for parallel processing.

    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // (nz is ignored for processing; we handle a single 2D slice)

    // Safety: if there are more ranks than rows, some ranks get 0 rows (works)
    std::vector<int> counts, displs, rowsPerRank;
    build_counts_displs(nx, ny, size, counts, displs, rowsPerRank);
    int localNy = rowsPerRank[rank];

    // -----------------------------
    // 3) Scatter rows to all ranks
    // -----------------------------
    //
    // Each rank receives only the rows it is responsible for processing.
    //
    // 1) Allocate `localData` to hold this rank's portion (localNy rows * nx columns).
    //
    // 2) If only one rank, just copy the full image (no scattering needed).
    //
    // 3) Otherwise, use `MPI_Scatterv` to distribute rows from rank 0 to all ranks:
    //    - Rank 0 provides the full image, counts, and displacements.
    //    - Other ranks provide nullptr for send parameters and receive only their portion.
    //
    // Purpose:
    // - This ensures each rank has its local sub-image for parallel processing.
    // - After this, each rank can independently threshold and label its rows without conflicts.

    std::vector<float> localData((size_t)localNy * nx, 0.0f);
    if (size == 1) {
        // trivial: already in fullData
        localData = fullData;
    } else {
        if (rank == 0) {
            MPI_Scatterv(fullData.data(), counts.data(), displs.data(), MPI_FLOAT,
                         localData.data(), counts[0], MPI_FLOAT, 0, MPI_COMM_WORLD);
        } else {
            MPI_Scatterv(nullptr, nullptr, nullptr, MPI_FLOAT,
                         localData.data(), counts[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
    }

    // -----------------------------
    // 4) Threshold to mask (foreground = value > 0)
    // -----------------------------
    //
    // Convert the local floating-point image into a binary mask for labeling.
    //
    // 1) Allocate `mask` with the same size as `localData`.
    // 2) Loop through each pixel:
    //    - If the pixel value > 0 → mark as foreground (1)
    //    - Otherwise → background (0)
    //
    // Purpose:
    // - Simplifies the image for Connected Component Labeling (CCL).
    // - Only the foreground pixels (1) will be labeled as connected components.

    std::vector<int> mask((size_t)localNy * nx, 0);
    for (size_t i = 0; i < mask.size(); ++i) mask[i] = (localData[i] > 0.0f) ? 1 : 0;

    // -----------------------------
    // 5) Local CCL on each rank
    // -----------------------------
    //
    // Each rank runs Connected Component Labeling (CCL) on its local sub-image.
    //
    // 1) `localLabels` will store the provisional labels for this rank's pixels.
    // 2) If the rank has rows (`localNy > 0`):
    //    - Call `localCCL` on the binary `mask` to assign temporary labels.
    //    - `localMaxTempLabel` records the largest label used locally.
    // 3) If the rank has no rows, initialize `localLabels` as empty.
    //
    // Purpose:
    // - Each rank labels connected foreground regions independently.
    // - Labels are still local and temporary; they will be made globally unique later.

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "================ DEBUG: step 5 ================" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);


    std::vector<int> localLabels;
    int localMaxTempLabel = 0;
    // ===== in main(), step 5 =====
    UF localUF;   // keep Union-Find
    localMaxTempLabel = localCCL(mask, nx, localNy, localLabels, localUF);


    std::cout << "[DEBUG][Rank " << rank << "] Starting local CCL: "
            << "localNy=" << localNy << ", nx=" << nx << std::endl;

    if (localNy > 0) {
        std::cout << "[DEBUG][Rank " << rank << "] Running localCCL on mask with "
                << (size_t)nx * localNy << " pixels" << std::endl;

        

        std::cout << "[DEBUG][Rank " << rank << "] localCCL finished: "
                << "localMaxTempLabel=" << localMaxTempLabel
                << ", localLabels.size()=" << localLabels.size() << std::endl;

        // Optional: print first few labels for sanity check
        if (!localLabels.empty()) {
            std::cout << "[DEBUG][Rank " << rank << "] First 10 local labels: ";
            for (size_t i = 0; i < std::min<size_t>(10, localLabels.size()); ++i) {
                std::cout << localLabels[i] << " ";
            }
            std::cout << std::endl;
        }

    } else {
        localLabels.assign(0, 0);
        std::cout << "[DEBUG][Rank " << rank << "] No rows assigned → localLabels empty, localMaxTempLabel=0" << std::endl;
    }


    // -----------------------------
    // 6) Give each rank a unique label range via EXSCAN
    // -----------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "================ DEBUG: step 6 ================" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int base = 0; // starting offset for this rank's label IDs

    std::cout << "[DEBUG][Rank " << rank << "] Before EXSCAN: localMaxTempLabel="
            << localMaxTempLabel << std::endl;

    MPI_Exscan(&localMaxTempLabel, &base, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Rank 0 has no previous ranks, set base = 0
    if (rank == 0) base = 0;

    std::cout << "[DEBUG][Rank " << rank << "] After EXSCAN: base=" << base
            << " (shift for local labels)" << std::endl;

    // Offset labels to be globally unique
    for (auto &v : localLabels) {
        if (v > 0) v += base;
    }


    std::cout << "[DEBUG][Rank " << rank << "] Applied base offset. "
            << "Adjusted label range: [1+" << base << " ... "
            << localMaxTempLabel + base << "]" << std::endl;

    // Optional: sanity check — show first few adjusted labels
    if (!localLabels.empty()) {
        std::cout << "[DEBUG][Rank " << rank << "] First 10 adjusted local labels: ";
        for (size_t i = 0; i < std::min<size_t>(10, localLabels.size()); ++i) {
            std::cout << localLabels[i] << " ";
        }
        std::cout << std::endl;
    }


    // -----------------------------
    // 7) Find cross-rank equivalences (only vertical neighbors between ranks)
    // -----------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "================ DEBUG: step 7 ================" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    

    std::vector<int> sendBottom(nx, 0), recvTop(nx, 0);

    // 1) Extract bottom row from local labels
    if (localNy > 0) {
        for (int x = 0; x < nx; ++x) 
            sendBottom[x] = localLabels[(localNy - 1) * nx + x];
    }

    std::cout << "[DEBUG][Rank " << rank << "] Prepared bottom row for send: ";
    for (int i = 0; i < std::min(10, nx); ++i) std::cout << sendBottom[i] << " ";
    if (nx > 10) std::cout << "...";
    std::cout << std::endl;

    // 2) Exchange bottom/top rows
    MPI_Status st{};
    MPI_Sendrecv(
        sendBottom.data(), nx, MPI_INT, (rank + 1 < size ? rank + 1 : MPI_PROC_NULL), 42,
        recvTop.data(),   nx, MPI_INT, (rank - 1 >= 0 ? rank - 1 : MPI_PROC_NULL), 42,
        MPI_COMM_WORLD, &st);

    std::cout << "[DEBUG][Rank " << rank << "] Received top row from rank-1: ";
    for (int i = 0; i < std::min(10, nx); ++i) std::cout << recvTop[i] << " ";
    if (nx > 10) std::cout << "...";
    std::cout << std::endl;

    // 3) Build equivalence pairs
    std::vector<int> eqPairs;

    // --- add inter-rank (border) pairs as before ---
    if (rank > 0 && localNy > 0) {
        for (int x = 0; x < nx; ++x) {
            int a = recvTop[x];
            int b = localLabels[x];
            if (a > 0 && b > 0 && a != b) {
                eqPairs.push_back(a);
                eqPairs.push_back(b);
            }
        }
    }

    // --- add intra-rank pairs from localUF ---
    if (localMaxTempLabel > 0) {
        for (int id = 1; id <= localMaxTempLabel; ++id) {
            int localRoot = localUF.find(id);   // ✅ safe: within [1..localMaxTempLabel]
            if (localRoot != id) {
                int A = id + base;              // shift both ends into global space
                int B = localRoot + base;
                eqPairs.push_back(A);
                eqPairs.push_back(B);
            }
        }
    }


    // -----------------------------
    // 8) Gather all equivalence pairs on rank 0
    // -----------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "================ DEBUG: step 8 ================" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int localPairCount = (int)(eqPairs.size() / 2); // pairs, not ints
    std::cout << "[DEBUG][Rank " << rank << "] Local equivalence pairs count: " 
            << localPairCount << std::endl;

    // Gather counts from all ranks
    std::vector<int> pairCounts(size, 0);
    MPI_Gather(&localPairCount, 1, MPI_INT, pairCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> allPairs;          // flattened pairs across ranks
    std::vector<int> pairDisplsInts;    // displs for ints (each pair has 2 ints)
    if (rank == 0) {
        int totalPairs = std::accumulate(pairCounts.begin(), pairCounts.end(), 0);
        std::cout << "[DEBUG][Rank 0] Total equivalence pairs from all ranks: " 
                << totalPairs << std::endl;

        allPairs.resize(totalPairs * 2, 0);
        pairDisplsInts.resize(size, 0);
        std::vector<int> pairCountsInts(size, 0);
        for (int r = 0, off = 0; r < size; ++r) {
            pairCountsInts[r] = pairCounts[r] * 2;
            pairDisplsInts[r] = off;
            off += pairCountsInts[r];
            std::cout << "[DEBUG][Rank 0] Rank " << r 
                    << " contributes " << pairCounts[r] << " pairs, displ=" 
                    << pairDisplsInts[r] << std::endl;
        }

        MPI_Gatherv(eqPairs.data(), localPairCount * 2, MPI_INT,
                    allPairs.data(), pairCountsInts.data(), pairDisplsInts.data(), MPI_INT,
                    0, MPI_COMM_WORLD);

        // Optional: preview the first few gathered pairs
        std::cout << "[DEBUG][Rank 0] First few gathered equivalence pairs: ";
        for (size_t i = 0; i + 1 < std::min<size_t>(allPairs.size(), 10); i += 2)
            std::cout << "(" << allPairs[i] << "," << allPairs[i + 1] << ") ";
        std::cout << std::endl;
    } else {
        MPI_Gatherv(eqPairs.data(), localPairCount * 2, MPI_INT,
                    nullptr, nullptr, nullptr, MPI_INT,
                    0, MPI_COMM_WORLD);
    }


    // -----------------------------
    // 9) Rank 0: Global Union-Find to merge all labels
    // -----------------------------
    //
    // After gathering all cross-rank equivalence pairs, rank 0 performs the final global merging.
    //
    // Steps:
    // 1) Compute totalLabelCapacity: sum of localMaxTempLabel across all ranks using MPI_Reduce.
    //    - This gives the maximum number of temporary labels that exist across all ranks.
    // 2) Initialize a Union-Find structure (guf) to manage equivalences for all labels globally.
    // 3) Merge cross-rank equivalences:
    //    - For each pair of labels (a, b) in allPairs, unite them in guf.
    //    - This ensures that labels that represent the same component across ranks are treated as one.
    // 4) Build a compact relabeling (globalMap):
    //    - Each unique root in the Union-Find gets a new consecutive label 1..K.
    //    - This creates a mapping from old (possibly scattered) label IDs to compact, globally consistent labels.
    //
    // Purpose:
    // - Ensure that connected components spanning multiple ranks are merged correctly.
    // - Produce a global mapping to relabel all local labels consistently across the entire image.
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "================ DEBUG: step 9 ================" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int totalLabelCapacity = 0;
    MPI_Reduce(&localMaxTempLabel, &totalLabelCapacity, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    std::vector<int> globalMap; // old -> new compact labels
    if (rank == 0) {
        std::cout << "[DEBUG][Rank 0] Total label capacity (sum of localMaxTempLabel): " 
              << totalLabelCapacity << std::endl;
        UF guf((size_t)totalLabelCapacity); // parents 0..totalLabelCapacity

        // Merge cross-rank equivalences
        // Merge cross-rank equivalences
        std::cout << "[DEBUG][Rank 0] Merging equivalence pairs..." << std::endl;
            for (size_t i = 0; i + 1 < allPairs.size(); i += 2) {
            int a = allPairs[i], b = allPairs[i + 1];
            if (a > 0 && b > 0) guf.unite(a, b);
        }
        // Build compact relabeling 1..K
        globalMap.assign(totalLabelCapacity + 1, 0);
        std::unordered_map<int, int> root2new;
        root2new.reserve((size_t)totalLabelCapacity);
        int next = 1;
        for (int id = 1; id <= totalLabelCapacity; ++id) {
            int r = guf.find(id);
            auto it = root2new.find(r);
            if (it == root2new.end()) {
            root2new[r] = next;
            globalMap[id] = next;
            std::cout << "[DEBUG][Rank 0] Assign new global label: old=" << id 
                      << ", root=" << r << ", new=" << next << std::endl;
            ++next;
            } else {
            globalMap[id] = it->second;
            std::cout << "[DEBUG][Rank 0] Reuse global label: old=" << id 
                      << ", root=" << r << ", new=" << it->second << std::endl;
            }
        }
        std::cout << "[DEBUG][Rank 0] Global mapping complete. Total unique labels: " << next - 1 << std::endl;
    }

    // -----------------------------
    // 10) Broadcast global mapping and remap local labels
    // -----------------------------
    //
    // After rank 0 has built the global label mapping (globalMap),
    // all other ranks need to update their local labels to match the global labeling.
    //
    // Steps:
    // 1) Broadcast the size of globalMap from rank 0 to all ranks using MPI_Bcast.
    //    - Each rank knows how big the mapping array is.
    // 2) Non-zero ranks resize their local copy of globalMap to receive the broadcasted data.
    // 3) Broadcast the globalMap array itself so that all ranks have the same mapping.
    // 4) Remap localLabels using globalMap:
    //    - Each local label > 0 is replaced by its compact, globally consistent label.
    //
    // Purpose:
    // - Ensure that each rank’s local labels correspond to the global labeling.
    // - After this step, all ranks share a consistent numbering for connected components.
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "================ DEBUG: step 10 ================" << std::endl;
    }   
    int mapSize = 0;
    if (rank == 0) mapSize = (int)globalMap.size();
    MPI_Bcast(&mapSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) globalMap.resize(mapSize, 0);
    if (mapSize > 0) {
        MPI_Bcast(globalMap.data(), mapSize, MPI_INT, 0, MPI_COMM_WORLD);
    }

    std::cout << "[DEBUG][Rank " << rank << "] Received globalMap size: " << mapSize << std::endl;

    // Optional: print first few entries of globalMap
    if (!globalMap.empty()) {
        std::cout << "[DEBUG][Rank " << rank << "] First 10 entries of globalMap: ";
        for (size_t i = 0; i < std::min<size_t>(10, globalMap.size()); ++i) {
            std::cout << globalMap[i] << " ";
        }
        std::cout << std::endl;
    }

    for (auto &v : localLabels) {
        if (v > 0) v = globalMap[v]; // compact final label
    }

    // Optional: print first few remapped local labels
    if (!localLabels.empty()) {
        std::cout << "[DEBUG][Rank " << rank << "] First 10 remapped local labels: ";
        for (size_t i = 0; i < std::min<size_t>(10, localLabels.size()); ++i) {
            std::cout << localLabels[i] << " ";
        }
        std::cout << std::endl;
    }

    // -----------------------------
    // 11) Gather final labels and write output on rank 0
    // -----------------------------
    //
    // After all ranks have updated their local labels to the global numbering:
    // 1) Prepare a container on rank 0 to hold the full image labels (fullLabels).
    // 2) If running on a single rank, just copy localLabels to fullLabels.
    // 3) If running on multiple ranks, gather all localLabels from each rank into fullLabels on rank 0
    //    using MPI_Gatherv, which respects different counts and displacements per rank.
    //
    // 4) On rank 0, write the full labeled image to a binary file:
    //    - Write the dimensions (nx, ny, nz_out = 1) first.
    //    - Convert integer labels to float for compatibility with the existing pipeline.
    //    - Write the label data to the output file.
    //
    // Purpose:
    // - Collect all distributed label data back to the master rank.
    // - Produce a single labeled image file with globally consistent component IDs.
    // - Finalizes the MPI environment.
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "================ DEBUG: step 11 ================" << std::endl;
    }

    std::vector<int> fullLabels;
    if (rank == 0) fullLabels.resize((size_t)nx * ny);

    if (size == 1) {
        fullLabels = localLabels;
        std::cout << "[DEBUG][Rank 0] Single rank: fullLabels copied from localLabels, size=" 
              << fullLabels.size() << std::endl;
    } else {
        std::cout << "[DEBUG][Rank " << rank << "] Sending " << counts[rank] 
              << " labels to rank 0" << std::endl;
        MPI_Gatherv(localLabels.data(), counts[rank], MPI_INT,
                    fullLabels.data(), counts.data(), displs.data(), MPI_INT,
                    0, MPI_COMM_WORLD);
        if (rank == 0) {
        std::cout << "[DEBUG][Rank 0] Received all labels from ranks" << std::endl;
        std::cout << "[DEBUG][Rank 0] First 20 labels in fullLabels: ";
        for (size_t i = 0; i < std::min<size_t>(20, fullLabels.size()); ++i)
            std::cout << fullLabels[i] << " ";
        std::cout << std::endl;
        }
    }
    if (rank == 0) {
    std::cout << "================ DEBUG: final labeled grid ================\n";
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            std::cout << fullLabels[y * nx + x] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "===========================================================\n";
    }


    if (rank == 0) {
        std::ofstream out(outFile, std::ios::binary);
        if (!out) {
            std::cerr << "Failed to open output file\n";
            MPI_Abort(MPI_COMM_WORLD, 5);
        }
        out.write((char *)&nx, sizeof(nx));
        out.write((char *)&ny, sizeof(ny));
        int nz_out = 1;
        out.write((char *)&nz_out, sizeof(nz_out));

        // Write as float for compatibility with your original pipeline
        std::vector<float> outData((size_t)nx * ny, 0.0f);
        for (size_t i = 0; i < outData.size(); ++i) outData[i] = (float)fullLabels[i];
        out.write((char *)outData.data(), outData.size() * sizeof(float));
        out.close();
    }


    MPI_Finalize();
    return 0;
}
