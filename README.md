# How to Run UMKP
## Serial Execution
Go to Ours-serial folder and run: 
```
make
./UMKP <dataset> <k> [parameters]
```
The kplex is written in `kplex.txt` file. 

for example: 
```
./UMKP ../datasets/socfb-Duke14.bin 2
```

## Parallel Execution
Go to Ours-parallel folder and run:
```
make
export OMP_NUM_THREADS=N
./UMKP <dataset> <k>
```

for example: 
```
./UMKP ../datasets/socfb-Duke14.bin 5
```


## Configuration of Parameters

Several switches can be configured by providing optional parameters to the program. The following table summarizes usage of parameters. Values in bold are the default values. 


| Parameter |Explaination | Values |
|----------|-------|--------------------|
| -branching | Branching method selection | **Default (uses S-Br for K<=5, Pivot-Br otherwise)**, S-Br, R-Br, SR-Br, Pivot-Br, Binary-Br  |
| -bounding  | Upper bounding method to terminate a branch | **None**, S-Bound, R-Bound, SR-Bound |
| -UBR2 | Select whether to apply pruning of vertices based on upper bounding of current branch | **true**, false|
| -topCTCP | Select whether to apply CTCP at top level | **true**, false|
| -twoHopG | Select whether to form two hop graph at second top level | **true**, false|
| -dense | Select whether to find densest k-plex, also computes all k-plexes (output in `all_kplexes.txt` file) |true, **false** |

### Enable CTCP in the Recursion

To enable CTCP in the recursion, compile the program with `-D_SECOND_ORDER_PRUNING_` 
the switch is enabled in following make target

```
make CTCP
```
While enabling CTCP all other parameters above can be configured. 

# How to Run Baselines

## Running *kPlexT*

```
cd kPlexT
make
./kPlexT -g <dataset> -k <k> -b
```

## Running *Maple*

```
cd Maple
make
./Maple <dataset> <k> 
```

## Running *DiseMKP*

```
cd DiseMKP
make
./DiseMKP <dataset> <k> 
```

# Dataset Preparation
Go to datasets folder and convert network-repository or SNAP graphs to binary format. 

```
cd datasets
g++ toBin.cpp -o toBin
./toBin graph.txt
```

The DIMACS formatted graph should be provided with a `.dim` file, e.g. 
```
p edge 4 3
e 1 2
e 2 3
e 4 1
```

The SNAP formatted graph should be provided with a `.txt` file, e.g. 
```
1 2
2 3
3 1
```