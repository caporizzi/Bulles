# Bubble Tracking and Characterization in Multiphase Flow Using HPC

## Goals 

The objective is to accurately detect, count, and characterize gas bubbles in multiphase flow simulations in order to improve the understanding of gas–liquid mixing phenomena. 
The algorithm is designed to identify and label bubbles, track their temporal evolution, and extract quantitative characteristics such as size. 
At the same time, the implementation must remain efficient and scalable, which is why we rely on distributed-memory parallelism with MPI to handle large-scale datasets on high-performance computing platforms. 
By combining robust detection methods with parallel computing techniques, the goal is to obtain a solution that is both computationally efficient and suitable for analyzing realistic two dimensional simulations.

## Conception

### Understanding and Drawing

To conceive the algorithm, the main source was "A parallel Eulerian interface tracking/Lagrangian point particle multi-scale coupling procedure" [1] where both the sequential and the parallel solution were given.
Here, this figure shows the sequential algorithm. The main idea is that it does a breath first search. Each time we find a node we look for the direct neighbor mark them then execute this on the neighbor. Once we can't find direct neighbors all the marked node are of the same id and we continue to iterate through the grid. 

![Figure 1 - Basic serial single block identification algorithm](https://i.imgur.com/dALZYxE.png) 

Here is a visual representation of the sequential algorithm labeling a node and their neighbor forming a block with the same label.

![enter image description here](https://i.imgur.com/VcUjxCI.png)

Here, this figure shows the parallel algorithm. 
Here the main idea is that the grid is divided in blocks. Each block is explored giving indexes to each node. Once all nodes are tagged we compare if node of different block are neighbors and update the node id.

![Figure 2 - Parallel, multi-block tag synchronization algorithm](https://i.imgur.com/8W6BvQj.png)

Here is a visual representation of the parallel algorithm. It checks if neighboring cell of different block exist and update the label.
![enter image description here](https://i.imgur.com/K34D2fm.png)



## Implementation

### Why this algorithm

The goal is to create something parallel friendly and simple based on the time given. 
Bubble identification can be formulated as a connected component labeling (CCL) problem on a binary grid. CCL is a natural choice because it provides a clear mapping from individual cells to bubble IDs and can be parallelized efficiently.

#### Equivalence 
Then, based on Discrete Math, Math based functions are always fast so went with Union based algorithm which will solve the equivalence part. 
Union-Find has nearly constant-time operations with path compression and union by rank, making it highly efficient for merging label sets both locally (intra-rank) and across process boundaries (inter-rank). 

Union-Find was chosen instead of BFS or DFS because it scales better for large grids. BFS and DFS both rely on traversing neighbors repeatedly, which can create redundant work when bubbles are large or spread across multiple ranks. DFS also suffers from recursion depth issues when components get big.
Union-Find avoids these issues because it works with a simple parent array.  
In practice this means it is faster and more parallel friendly than BFS or DFS, especially when equivalences must be resolved across MPI ranks.

#### Connectivity
Connectivity was restricted to 4-neighbors (up, down, left, right) instead of 8-neighbors. This choice avoids diagonal merging of bubbles, which can artificially connect distinct structures and bias size distribution. Moreover, 4-connectivity reduces ambiguity and simplifies communication at process boundaries, making it more aligned with physical interpretations in gas–liquid multiphase flows.



### About data parallelism

For simplicity purposes, the choice is to only divide my grid horizontally and to exclude vertical.
Each MPI rank gets rows. In this case it means that we only care about the communication from the rank above and below.

### Step-by-step code implementation


Rank 0 begins by reading the entire input file to ensure a single, consistent copy of the image. After that, the image dimensions are broadcast to all ranks using `MPI_Bcast`, so every rank knows the size of the data.

Next, we determine how many rows each rank will process using `build_counts_displs`, which calculates the number of rows per rank, the total number of elements per rank, and the starting offsets for each rank. The rows are then distributed to each rank using `MPI_Scatterv`, giving each rank its own portion of the image to work on independently, which avoids conflicts.

Once each rank has its local rows, a threshold is applied to convert the image into a binary mask with values 0 and 1, identifying background and foreground pixels.

At this stage, parallel processing begins. Each rank performs Connected Component Labeling (CCL) on its local sub-image. This assigns provisional labels to connected foreground pixels. The maximum label used by each rank is recorded.

To ensure globally unique labels, an offset is computed for each rank using `MPI_Exscan`. Each rank adds its offset to its local labels, so no two ranks share the same label ID.

Next, cross-rank equivalences are identified. Each rank sends its bottom row to the next rank and receives the top row from the previous rank using `MPI_Sendrecv`. Any labels that match across these rows are recorded as equivalences. In addition, a local Union-Find structure resolves intrarank equivalences between provisional labels.

All equivalence pairs (inter- and intra-rank) are gathered on rank 0 using `MPI_Gatherv`. Rank 0 then builds a global mapping (`globalMap`) that assigns compact, consecutive labels for all connected components across ranks.

This global mapping is broadcast to all ranks, which then remap their local labels accordingly. Finally, the fully labeled image is gathered on rank 0 using `MPI_Gatherv` and written to the output file.


## Verification
Here, my code pass through several test cases to prove it works as intended. 

### Testcase 1
From the Testcase 1, launched with 3 process, we expect the result to have 2 labels max.

The grid is initiated, divided in 3, ranging from rank 0 to rank 2. 
Each rank scan locally each value from left-to-right, top to bottom. 
If one local neighbor labeled, it copy that label.
If both neighbors are labeled, it copy the one on the left.

![enter image description here](https://i.imgur.com/Pg3di1r.png)

Then, each rank register their max label and we proceed to an offset.

Then each rank sends their bottom row to their rank+1, and receive the top rank from their rank-1.

![enter image description here](https://i.imgur.com/Qv9iexN.png)

Here comes the pair gathering. 
We get interrank pairs and intrarank pairs.

![enter image description here](https://i.imgur.com/e56o9LU.png)

Here comes the pair equivalence.
From the previous step, we know that all my 4 are 2. 
Also know that all my 3 are 4. This means that all my 3 which are 4, are in reality 2.
Then we know that all my 5 are 3, this means that all my 5 which are 3 are in reality 4 which are in reality 2.
So we get the global and final grid, resulting as 2 maximal global label as expected.

![enter image description here](https://i.imgur.com/8yiseSZ.png)

### Testcase 2

Same grid as Testcase 1, but with only 2 process, we expect the same result.

We get "First few gathered equivalence pairs: (3,4) (3,4) (2,4) (2,4)" 
The offset and the pair are what is intended after the global merge, we reach the expected result.

### Testcase 3

Here I ask a LLM to give me some testcase.

It outputs a test case which ressembles a snake, testing if the corner, rank communication and gathering works as intended.
It outputs a test case which looks like little boxe.
Both outputs were as intended.

![enter image description here](https://i.imgur.com/48Ds8fV.png) 

 ## Analysis

This output is processed by a Python program designed to count the number of bubbles in my domain and to measure the size of the gas regions.

Here is an example for the initial iteration. 
![enter image description here](https://i.imgur.com/z9mMI1O.png)

Here, I choose to view it as a percentage of gasses versus liquidity.
![enter image description here](https://i.imgur.com/bIF9Zr1.png)

Here, one of the last iteration in my simulation.

As expected most of the bubble have risen to the surface.
![enter image description here](https://i.imgur.com/DSmegJD.png)

![enter image description here](https://i.imgur.com/e6Nv6nR.png)

The domain is mostly liquid, with ~70% of cells as liquid and ~30% representing gas regions. 

The main bottleneck will likely come from the global equivalence resolution step. This stage currently relies on gathering equivalence pairs on rank 0, which introduces a synchronization point and limits parallel efficiency as the number of processes increases. 	

Since all the execution was done locally, performance metrics such as runtime scaling or memory consumption were excluded from the analysis. The main focus was on verifying correctness and ensuring that the algorithm works in a distributed-memory setting. However, the algorithm was designed with scalability in mind, and the next step is to run it on a supercomputer. With access to larger resources, it will be possible to measure speed-up, efficiency, and memory usage across hundreds or thousands of MPI ranks. This will provide a real evaluation of the parallel efficiency of the Union-Find based implementation and confirm its suitability for high-performance computing platforms.


## Limitations and Future Work

### Limitations 
While the current implementation successfully identifies and tracks bubbles in a parallel framework, several limitations remain. First, the algorithm only employs horizontal partitioning, which simplifies communication but limits scalability for very large domains.
As the grid size increases, the communication between MPI ranks also grows, particularly during the exchange of boundary rows and the global equivalence resolution. This may limit the efficiency gains from additional processes in large-scale simulations.
### Future work

One idea is to implement a fully parallel Union-Find with global synchronization, which would allow equivalence resolution to be done more efficiently without relying on a central process.

These improvements would not only extend the algorithm’s functionality but also could help increase speed-up and improve scalability on high-performance computing systems. It should be noted that in the current work, the primary focus was on correctness and parallel implementation of the algorithm, rather than on optimizing performance or achieving maximum speed-up.

### How to execute

Convert your `.vtk` into `.dat` using `./conv`

Execute the program `mpirun -np n ./track yourfile.dat yourresult.dat` where n is the number of process

## Conclusion

As expected, the bubbles move upward and merge together, which reduces the number of separate regions over time. 
This matches the physical behavior of multiphase flows.

## References 
[1] "A parallel Eulerian interface tracking/Lagrangian point particle multi-scale coupling procedure" M. Herrmann, 2009 


