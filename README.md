# MKP-Study
## Serial Execution
Go to Ours-serial folder and run: 
```
make
./UMKP <dataset> <k> [parameters]
```
The kplex is written in kplexes.txt file. 

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
## Dataset Preparation
Go to datasets folder and convert network-repository or SNAP graphs to binary format. 

```
cd datasets
g++ toBin.cpp -o toBin
./toBin graph.txt
```

## Configuration of Parameters

Several switches can be configured by providing optional parameters to the program. The following table summarizes usage of parameters.


| Parameter | Values |
|----------|--------------------|
| -branching | **Default (uses S-Br for K<=5, Pivot-Br otherwise)**, S-Br, R-Br, SR-Br, Pivot-Br, Binary-Br  |
| -bounding  | **None**, S-Bound, R-Bound, SR-Bound |
| -UBR2 | **true**, false|
| -topCTCP  | **true**, false|
| -twoHopG  | **true**, false|

# Enable CTCP in the Recursion

To enable CTCP in the recursion, compile the program with -D_SECOND_ORDER_PRUNING_, the switch is enabled in following make target

```
make CTCP
```

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

