# KDTree_C_DataMining ReadMe File
# CS 4331 Data Mining & High Performance Computing
# 9/23/2018
# Huyen Vu
# Kurtis Jackson
# Samuel Eash



# File compilation#

# Windows
# GNU bash, version 4.1.11(2)-release (x86_64-unknown-cygwin)
# gcc main.c lib/utilities.c lib/kdtree.c -o KDTree -std=gnu99
# run KDTree.exe file in Windows

# Ubuntu Linux 18.04 LTS
# gcc version 7.3.0 (Ubuntu 7.3.0-16ubuntu3)
# gcc main.c ./lib/utilities.c ./lib/kdTree.c -std+gnu99 -o kdtree
# run ./kdtree in ubuntu


# main.c Documentation#

# main.c is used to generate the dataset that will then be sorted into a KD-Tree
# main.c delacares and pre-initializes output values for the kdTree function
# main.c generates random dataset data by using 1/rand() to generate numbers between 0 and 1
# main.c generates datasets containing n_data points, n_dimensions, and K^2 clusters(ex 2,4,8,16,....ect)
# main.c passes generated dataset to the kdTree function for sorting into a KD-Tree



# kdTree.c Documentation#

# kdTree.c receives data from main.c for building the KD-Tree
# kdTree.c builds the entries of the KD-Tree based on bi-partition function
# kdTree.c calculates the max variance and centroid of each cluster that is passed through it
# max variance is computed using: var x = 1/n (sum (xi- |x|)^2)
# cluster centroid is computed using: mean x = 1/n sum(x1+x2+...)

# kdTree.c shifts the elements of previous clusters in the KD-Tree to make room for new clusters so 
# variance is preserved between clusters

# kdTree.c calculates the cluster boundaries of each cluster in n-dimensions as to show an accurate 
# physicial discription of each cluster

# kdTree.c uses the biPartition function to split each cluster into 2 subclusters until the 
# number of clusters has reached k^2 clusters.



# utilities.c Documentation#

# utilities.c recieves data from kdTree.c and main.c
# utilities.c is used to allocate pointer variables for both kdTree.c and main.c
# print out starting data, cluster data, and final results of the KD-Tree construction
