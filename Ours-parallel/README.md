# Maximum k-Plex Computation: Theory and Practice

This repository implements the maximum k-plex computation algorithm kPlexT proposed in our SIGMOD 2024 paper. If you are using the code, please cite our paper.
<pre>
Lijun Chang, Kai Yao:
Maximum k-Plex Computation: Theory and Practice.
Proc. ACM Manag. Data 2(1): 63:1-63:26 (2024)
</pre>

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "kPlexT", which corresponds to the kPlexT algorithm.

## Run the code

```sh
$ ./kPlexT -g {path_to_graph} -k {k_value}
```

An example of computing the exact maximum 3-plex for the dataset CA-GrQc is as follows
```sh
$ ./kPlexT -g datasets/CA-GrQc -k 3
```

## Data format
Two data formats are supported. The default data format is "edges.txt", which contains a list of undirected edges represented as vertex pairs. The first line contains two numbers n and m, representing the number of vertices and the number of undirected edges, respectively. Note that, the vertex ids must be between 0 and n-1.

The more time-efficient format is the binary format; to read the input graph from this format, please add "-b" when running the code. Each graph is represented by two binary files, b_adj.bin and b_degree.bin (e.g. datasets/CA-GrQc/b_adj.bin and datasets/CA-GrQc/b_degree.bin). More details of the data format can be found in [https://lijunchang.github.io/Cohesive_subgraph_book/datasets](https://lijunchang.github.io/Cohesive_subgraph_book/datasets)

## Get datasets
10th DIMACS graphs collection: https://networkrepository.com/dimacs10.php

Real-world graphs collection: http://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz
