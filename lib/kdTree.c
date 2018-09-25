#include <malloc.h>
#include <stdio.h>
#include "kdTree.h"
#include "utilities.h"

int *allocateTempIndices(int kk, int current_kk) {
    //invariant: index for next cluster will be previous index + kk/current_kk
    int *temp_indices = allocatePointerVariable(current_kk);
    temp_indices[0] = 0;
    for(int i = 1; i< current_kk; i++) {
        temp_indices[i] = temp_indices[i-1]+kk/current_kk;
    }
    return temp_indices;
}

int calculateClusterCentroids(int dim, double *data, int cluster_start, int cluster_size, double *cluster_centroid) {
//invariant: mean x = 1/n sum(x1+x2+...)
    for(int i=0; i<dim; i++) {
        cluster_centroid[i] = 0;
    }

    if(cluster_size == 0) {
        return 0;
    }

    for(int i = 0; i < cluster_size; i++) {
        struct dataPoint e = getElement(dim, cluster_start + i, data);
        for(int j = 0; j< dim; j++) {
            cluster_centroid[j]+= e.data[j];
        }
    }

    for(int i = 0; i < dim ; i++) {
        cluster_centroid[i] = cluster_centroid[i] / cluster_size;
    }
}


int findMaxVariance(int dim, double *data, int cluster_start, int cluster_size, double* cluster_centroid) {
//invariant: calculate variances of each dimension using formula var x = 1/n (sum (xi- |x|)^2)
    double *variances = (double*) malloc(dim * sizeof(double));
    for(int i = 0; i< dim; i++) {
        variances[i] = 0;
    }

    for(int i = 0; i< cluster_size; i++) {
        struct dataPoint element = getElement(dim, cluster_start + i, data);
        for (int j =0; j< dim; j++) {
            variances[j] += square(element.data[j] - cluster_centroid[j]);
        }
    }

    double max_variance = 0;
    int max_variance_index = 0;
    for(int i = 0; i< dim; i++) {
        variances[i] = variances[i] /cluster_size;
        if (max_variance < variances[i]) {
            max_variance = variances[i];
            max_variance_index = i;
        }
    }
    return max_variance_index;
}

void shiftElementUpward(int dim, double *data, int counter, int cluster_start, int elementIdx) {
    int swapElementIndex = cluster_start + counter;
    struct dataPoint topElement = getElement(dim, swapElementIndex, data);
    struct dataPoint element = getElement(dim, elementIdx, data);

    for(int i = 0; i< dim; i++){
        data[swapElementIndex * dim+i] = element.data[i];
        data[elementIdx * dim + i] = topElement.data[i];
    }
}

int calculateClusterBoundaries(int dim, double *data, int cluster_start, int cluster_size, double *cluster_bdry) {
    double tempMin[dim];
    double tempMax[dim];
    for (int i=0; i< dim; i++) {
        tempMin[i] = 1;
        tempMax[i] = 0;
    }

    for (int i =0; i< cluster_size; i++ ) {
        struct dataPoint element = getElement(dim, cluster_start + i, data);
        for(int j =0; j< dim; j++) {
            if(tempMin[j] > element.data[j]) {
                tempMin[j] = element.data[j];
            }

            if(tempMax[j] < element.data[j]) {
                tempMax[j] = element.data[j];
            }
        }
    }

    int cluster_bdry_idx = 0;
    for (int i =0; i< dim; i++) {
        cluster_bdry[cluster_bdry_idx] = tempMin[i];
        cluster_bdry_idx++;
        cluster_bdry[cluster_bdry_idx] = tempMax[i];
        cluster_bdry_idx++;
    }
}

void calculateNewOutputs(int dim, int kk, int upper_cluster_size, int current_kk, int idx, double *data, int *cluster_start, int *cluster_size,
                         double **cluster_centroid, double **cluster_bdry){
    //invariant: find new cluster_starts, cluster_sizes, cluster_bdrys for two new clusters based on cluster_indices
    int lower_cluster_size = cluster_size[idx] - upper_cluster_size;
    int lower_cluster_idx = idx + (kk/(current_kk * 2));

    //assign new cluster_starts
    cluster_start[lower_cluster_idx] = cluster_start[idx] + upper_cluster_size;

    //assign new cluster_sizes
    cluster_size[idx] = upper_cluster_size;
    cluster_size[lower_cluster_idx] = lower_cluster_size;

    //calculate cluster centroids
    calculateClusterCentroids(dim, data, cluster_start[idx], cluster_size[idx], cluster_centroid[idx]);
    calculateClusterCentroids(dim, data, cluster_start[lower_cluster_idx], cluster_size[lower_cluster_idx], cluster_centroid[lower_cluster_idx]);

    //calculate cluster boundaries:
    calculateClusterBoundaries(dim, data, cluster_start[idx], cluster_size[idx], cluster_bdry[idx]);
    calculateClusterBoundaries(dim, data,  cluster_start[lower_cluster_idx], cluster_size[lower_cluster_idx], cluster_bdry[lower_cluster_idx]);
}

int rearrangeData(int dim, int partition_dimension, double *data, int cluster_start, int cluster_size,const double* cluster_centroid) {
    //invariant: based on cluster_centroid and value of datapoint at dimension = partitionDimension, separate the cluster into 2 small clusters
    //return number of upper cluster
    int counter = 0;
    for(int i= 0; i< cluster_size; i++) {
        int idx = cluster_start + i;
        struct dataPoint element = getElement(dim, idx, data);
        if(cluster_centroid[partition_dimension] < element.data[partition_dimension]) {
            //printf("%f --- %f\n", cluster_centroid[partition_dimension], element.data[partition_dimension]);
            shiftElementUpward(dim, data, counter, cluster_start, idx);
            counter ++;
        }
    }

    return counter;
}

int biPartition(int dim, double *data, int cluster_start, int cluster_size, double* cluster_centroid) {
    int partitionDimension = findMaxVariance(dim, data, cluster_start, cluster_size, cluster_centroid);
    int upper_cluster_size = rearrangeData(dim, partitionDimension, data, cluster_start, cluster_size, cluster_centroid);
    return upper_cluster_size;
}

int biPartitionClusters(int dim, int kk, int current_kk, const int *cluster_indices, double *data, int *cluster_start, int *cluster_size,
                        double **cluster_bdry, double **cluster_centroid) {
    for( int i = 0; i < current_kk; i++) {
        int idx = cluster_indices[i];
        int upper_cluster_size = biPartition(dim, data, cluster_start[idx], cluster_size[idx], cluster_centroid[idx]);
        calculateNewOutputs(dim, kk, upper_cluster_size, current_kk, idx, data, cluster_start, cluster_size, cluster_centroid, cluster_bdry);
    }
}

int kdTree(int dim, int ndata, double *data, int kk, int *cluster_start, int *cluster_size,
           double **cluster_bdry, double **cluster_centroid) {
    int current_kk = 1;
    cluster_size[0] = ndata;
    cluster_start[0] = 0;
    calculateClusterCentroids(dim, data, cluster_start[0], ndata, cluster_centroid[0]);

    while(current_kk < kk) {
        int *cluster_indices = allocateTempIndices(kk, current_kk);
        biPartitionClusters(dim, kk, current_kk, cluster_indices, data, cluster_start, cluster_size, cluster_bdry, cluster_centroid);
        current_kk *= 2;
    }
}