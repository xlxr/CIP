# Data for replication

The folder "random" contains 15 randomly generated networks, and the folder "real life" contains 16 real-life networks. Each file consists of 1 + X + Y lines, where X and Y are the number of vertices and the number of edges in the network, respectively. The file starts with two numbers indicating the number of vertices X and the number of edges Y. Then, it follows with X lines of one number. Note that we do not assign indices to the vertices. The index of a vertex is the same as the line number it appears on. The number in each line is the weight of a vertex. The last Y lines contain three numbers. The first and second numbers in each line are the head and the tail of an edge. Note that the edges are undirected. The third number is the blocking cost of an edge.

For the randomly generated networks, we selected some networks from Buchanan et al.'s paper (cited below) and assigned random weights and blocking costs. 
A. Buchanan, J.S. Sung, S. Butenko, E.L. Pasiliao. An integer programming approach for fault-tolerant connected dominating sets. INFORMS Journal on Computing, 27(1):178-188, 2015
https://austinlbuchanan.github.io/files/CDS.zip

For the real-life instances, we selected some networks from Buchanan et al.'s paper (as cited above), along with additional networks from various categories in the Network Repository (cited below), and assigned random weights and blocking costs.

R. A. Rossi and N. K. Ahmed. The network data repository with interactive graph analytics and visualization. In AAAI, 2015. URL https://networkrepository.com.

The Rome99 network can be downloaded from https://www.diag.uniroma1.it/challenge9/download.shtml.
