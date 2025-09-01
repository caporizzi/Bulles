# Bubble Tracking and Characterization in Multiphase Flow Using HPC

## Goals 

We want to accurately identify, count and characterize bubbles for understanding gas-liquid mixing. We track bubble's evolution over time and extract characteristic such as size. We develop algorithms for bubble detection and leverage parallel computing techniques. 

## Conception

The conception of the algorithm was done in three phases:
1- Understanding & drawing both given algorithm
2- Reusing concept from previous courses 
3- Planning, coding

### Understanding and Drawing

To conceive the algorithm, the main source was "A parallel Eulerian interface tracking/Lagrangian point particle multi-scale coupling procedure" [1] where both the sequential and the parallel solution were given.
Here, Figure 1 shows the sequential algorithm. The main idea is that it does a breath first search. Each time we find a node we look for the direct neighbor mark them then execute this on the neighbor. Once we can't find direct neighbors all the marked node are of the same id and we continue to iterate through the grid. 

![Figure 1 - Basic serial single block identification algorithm](https://i.imgur.com/dALZYxE.png) 

Here is a visual representation of the sequential algorithm labeling a node and their neighbor forming a block with the same label.

![enter image description here](https://i.imgur.com/VcUjxCI.png)

Here, Figure 2 shows the parallel algorithm. 
Here the main idea is that the grid is divided in blocks. Each block is explored giving indexes to each node. Once all nodes are tagged we compare if node of different block are neighbors and update the node id.

![Figure 2 - Parallel, multi-block tag synchronization algorithm](https://i.imgur.com/8W6BvQj.png)

Here is a visual representation of the parallel algorithm. It checks if neighboring cell of different block exist and update the label.
![enter image description here](https://i.imgur.com/K34D2fm.png)



## Implementation

Here I discuss my algorithm choice

### Why this algorithm

The goal is to create something parallel friendly and simple based on the time given. 
The first choice is to implement a connected component labelling to create label. As given in the algorithm.
Then, based on Discrete Math, Math based functions are always fast so went with Union based algorithm which will solve the equivalence part. 
Plus, I followed the choice of the studied algorithm and stayed in a 4-connectivity to avoid diagonals connectivity.

### About data parallelism

For simplicity purposes, the choice is to only divide my grid horizontally and to exclude vertical.
Each MPI rank gets rows. In this case it means that we only care about the communication from the rank above and below.





## Verification
Here, my code pass through several test cases to prove it works as intended 

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

Here, the output is treated by a python program which aims to count the number of bubble inside my domain and also to see how big gaz are in my domain.

Here is an example for one iteration. My visit display this, I added the label.
![enter image description here](https://i.imgur.com/z9mMI1O.png)

Here, I choose to view it as a percentage of gasses versus liquidity.
![enter image description here](https://i.imgur.com/bIF9Zr1.png)
## References 
[1] "A parallel Eulerian interface tracking/Lagrangian point particle multi-scale coupling procedure" M. Herrmann, 2009 

## Figures
Figure 1 - Basic serial single block identification algorithm
Figure 2 - Parallel, multi-block tag synchronization algorithm


