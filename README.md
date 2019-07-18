# A Local Expansion and Conductance Minimizing C++ implementation

A C++ implementation of the Local Expansion and Conductance Minimizing (LECM) [1] algorithm. The LECM is an overlapping community detection algorithm for large-scale graphs proposed by Gao, Zhang and Zhang (2019) [1].

Since I am not a LECM author, I tried to made this implementation as close as possible to the algorithm description provided by the authors in [1].

### Disclaimer

> Note that I am not a LECM author, so this LECM version may has errors and/or discrepancies with the actual Gao, Zhang and Zhang [1] LECM algorithm.

## Prerequisites

* GNU Make

* gcc 4.8 compiler or an early version.

## Compile and run instructions

Go to the lecm folder and to compile type:

```sh
$ make
```

To run execute the lecm file:

```sh

$ ./lecm [flag] <value>

```

The flags are the parameters of the LECM algorithm (see [1] for a full description of each parameter):

| flag | description |
| --- | --- |
| -f | input graph path |
| -a | alpha value |
| -b | beta value |
| -c | chi value |
| -e | epsilon value |
| -t | theta value |

`-f` must to be specified. For the others flags, the algorithm can use default values following [1]:

`a = 0.99`, `b = 0.8`, `c = 0.5`, `e = 1e-4` and `t = 0.5`.

LECM algorithm reads file containing the graph's list of edges. In this file, edges are increasing ordered considering the source vertex and they are repeated even if the graph is undirected. In addition, the first file entry has the number of vertices. For example:
```
4
1	2
1	4
2	1
2	3
3	2
4	1
```
The example above represents an undirected graph with four vertices and three edges `(1, 2)`, `(1, 4)` and `(2, 3)`.

After execution, is produced a file `clustering.dat` at the `./output/` directory. This file contains the clustering generated by LECM following the structure below:
```
0 2 5
1
3 6 7 4
```
Means that the clustering has three clusters where each row represents a cluster and its memberships. In the example above, the clustering has three clusters in which vertices `0`, `2` and `5` are contained in the first cluster, vertex `1` belongs to the second one and vertices `3`, `6`, `7` and `4` are contained in the third cluster.

### Example:

```sh
$ ./lecm -f network_1000_yyyy.lfi -a 0.98 -b 0.75 -c 0.55 -e 1e-5 -t 0.6
```

## License

This project is licensed under the GNU General Public License - see the [LICENSE.md](LICENSE.md) file for details.

## References

**[\[1\] Y. Gao, H. Zhang and Y. Zhang. Overlapping community detection based on conductance optimization in large-scale networks. Physica A 522 (2019) 69-79.](https://https://www.sciencedirect.com/science/article/pii/S0378437119301487)**