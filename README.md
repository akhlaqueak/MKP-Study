# MKP-Study
## Serial Execution
Go to Ours-serial folder and run: 
```
make
./UMKP <dataset> <k>
```

for example: 
```
./UMKP ../datasets/socfb-Duke14.bin
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
Go to datasets folder and convert network-repository or SNAP graphs to binary format

```
./toBin graph.txt
```

## Create all versions
Go to datasets folder and convert network-repository or SNAP graphs to binary format

```
make versions
```