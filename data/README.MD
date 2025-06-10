# Data for replication

Folder "random" contains 15 randomly-generated networks and folder "real life" contains 15 real-life networks. Each file contains 1+X+Y lines, where X and Y are the number of vertices and the number of edges in the network, respectively. It starts with two numbers indicating the number of vertices X and the number of edges Y. Then, it follows with X lines of two numbers. Note that we do not assign indices to the vertices. The index of a vertex is the same as the line number it appears on. The first number in each line is the weight of a vertex and the second number is its blocking cost. The last Y lines contain three numbers. The first and second numbers in each line are the head and the tail of an edge. Note that the edges are undirected. The third number is the blocking cost of an edge.

For the randomly-generated networks, we choose some networks from Buchanan et. al.'s paper cited below and generate random weights and blocking costs. 
A. Buchanan, J.S. Sung, S. Butenko, E.L. Pasiliao. An integer programming approach for fault-tolerant connected dominating sets. INFORMS Journal on Computing, 27(1):178-188, 2015
https://austinlbuchanan.github.io/files/CDS.zip

For the real-life instances, we chose some networks from Buchanan et. al.'s paper cited above and some networks from different categories of Network Repository cited below and generate random weights and blocking costs.

R. A. Rossi and N. K. Ahmed. The network data repository with interactive graph analytics and visualization. In AAAI, 2015. URL https://networkrepository.com.
