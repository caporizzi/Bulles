# Bubble Tracking and Characterization in Multiphase Flow Using HPC

## Goals 

We want to accurately identify, count and characterize bubbles for understanding gas-liquid mixing. We track bubble's evolution over time and extract characteristic such as size. We develop algorithms for bubble detection and leverage parallel computing techniques. 

## Conception

To conceive the algorithm, the main source was "A parallel Eulerian interface tracking/Lagrangian point particle multi-scale coupling procedure" [1] where both the sequential and the parallel solution were given.
Here, Figure 1 shows the sequential algorithm. The main idea is that it does a breath first search. Each time we find a node we look for the direct neighbor mark them then execute this on the neighbor. Once we can't find direct neighbors all the marked node are of the same id and we continue to iterate through the grid. 
![Figure 1 - Basic serial single block identification algorithm](https://i.imgur.com/dALZYxE.png) 
Here, Figure 2 shows the parallel algorithm. 
Here the main idea is that the grid is divided in blocks. Each block is explored giving indexes to each node. Once all nodes are tagged we compare if node of different block are neighbors and update the node id.
![Figure 2 - Parallel, multi-block tag synchronization algorithm](https://i.imgur.com/8W6BvQj.png)

The conception was done in three phases:
1- Understanding & drawing both algorithm
2- Reusing concept from previous courses 
3- Planning, coding, testing and verifying



## References 
[1] "A parallel Eulerian interface tracking/Lagrangian point particle multi-scale coupling procedure" M. Herrmann, 2009 

## Figures
Figure 1 - Basic serial single block identification algorithm
Figure 2 - Parallel, multi-block tag synchronization algorithm


