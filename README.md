# KDTree_C_DataMining ReadMe File
Data Mining & High Performance Computing
 9/23/2018

## File compilation

Windows
GNU bash, version 4.1.11(2)-release (x86_64-unknown-cygwin)
```
gcc main.c lib/utilities.c lib/kdTree.c lib/kdTreeSearch.c -o KDTree -std=gnu99
```
run ```KDTree.exe``` file in Windows

Ubuntu Linux 18.04 LTS
gcc version 7.3.0 (Ubuntu 7.3.0-16ubuntu3)
```
gcc main.c ./lib/utilities.c ./lib/kdTree.c lib/kdTreeSearch.c -std+gnu99 -o kdtree
```
run ```./kdtree``` in ubuntu


## main.c Documentation

main.c is used to generate the dataset that will then be sorted into a KD-Tree
- delacares and pre-initializes output values for the kdTree function
- generates random dataset data by using (double)(rand()%10)/10 to generate numbers between 0 and 1
- generates datasets containing n_data points, n_dimensions, and K^2 clusters(ex 2,4,8,16,....ect)
- passes generated dataset to the kdTree function for sorting into a KD-Tree



## kdTree.c Documentation

kdTree.c receives data from main.c for building the KD-Tree
- builds the entries of the KD-Tree based on bi-partition function
- calculates the max variance and centroid of each cluster that is passed through it max variance is computed using: var x = 1/n (sum (xi- |x|)^2) 
- cluster centroid is computed using: mean x = 1/n sum(x1+x2+...)
<<<<<<< HEAD
- clusters can b located by using cluster_start and cluster_size
=======
- clusters are divided to smaller clusters at the dimension that has largest value of variance. kdTree.c shifts the elements have larger value than value of centroid at that dimension up in the KD-Tree to make room for new clusters so variance is preserved between clusters 

- clusters can be located by using cluster_start and cluster_size

>>>>>>> 60d79419f19d37081981ed19174f90fd63a017f7
- calculates the cluster boundaries of each cluster in n-dimensions as to show an accurate physicial discription of each cluster
- uses the biPartition function to split each cluster into 2 subclusters until the number of clusters has reached k^2 clusters.

## kdTreeSearch.c Documentation

- kdTreeSearch take in inputs of a KdTree (dimension, number of cluster, data, clusters' information ....) and a query data point
- outputs are an integer indicating number of data points that are checked being return out of the function and a pointer of doubles that points to the closest data point
- kdTreeSearch function calculates distances from the query data point to each cluster
- kdTreeSearch function exhausted searches the closest cluster to find the data point in the cluster that closest to the query point
- function checks the distance from query to all clusters  to verify no cluster is closer than the data point found
- If there is cluster that closer to data point than the data point found than the kdTreeSearch function do exhausted search within the cluster again until there is no cluster that is closer to the query than the data point found

## utilities.c Documentation

utilities.c is a group of helper functions used to facilitate the coding process
