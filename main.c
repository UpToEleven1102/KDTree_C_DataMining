#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


//calculate cluster boundary
//cluster centroid


struct dataPoint {
    double *data;
};

//Utility thingies
int *allocatePointerVariable(int array_size) {
    int *cluster_variable = (int *)malloc(array_size*sizeof(int));
    for(int i=0; i<array_size; i++) {
        cluster_variable[i]= -1;
    }
    return cluster_variable;
}

double **allocateDoublePointerVariable(int kk, int array_size) {
    double **cluster_variable = (double**)malloc(kk*sizeof(int *));
    for(int i=0; i < kk; i++) {
        cluster_variable[i] = (double*) malloc(sizeof(array_size*sizeof(double)));
    }
    return cluster_variable;
}

double square(double i) {
    return i*i;
}

struct dataPoint emptyDataPoint(int dim) {
    struct dataPoint out;
    out.data = (double*)malloc(dim * sizeof(double));
    return out;
};

struct dataPoint getElement(int dim, int index, const double *data) {
    struct dataPoint out = emptyDataPoint(dim);
    int n = 0;

    for(int i = index * dim; i < (index + 1) * dim; i++) {
        out.data[n] = data[i];
        n++;
    }
    return out;
}

void printArray(int dim, int ndata, double *data) {
    for (int i = 0; i < dim * ndata; i++) {
        printf("%f \n", data[i]);
    }
}

void printData(int dim, int ndata, double *data) {
    for(int i = 0; i< ndata; i++) {
        struct dataPoint e = getElement(dim, i, data);
        printf("%d-----\n", i);
        printArray(dim, 1, e.data);
    }
}

void printCluster(int dim, int idx, int cluster_start, int cluster_size, double *cluster_bdry, double *cluster_centroid,  double *data) {
    printf("%d, %d, %d, %d \n", dim, idx, cluster_start, cluster_size);

    printf("---------Data Points--------\n");
    for (int i = 0; i < cluster_size; i++) {
        struct dataPoint e = getElement(dim, cluster_start/dim + i, data);
        printf("%d-----\n", i);
        printArray(dim, 1, e.data);
    }

    printf("-------Cluster Boundaries-------\n");
    for(int i=0; i< dim; i++) {
        if(i%2 ==0){
            printf("Dimension %d - min %f \n", i, cluster_bdry[i*2-1]);
        } else {
            printf("Dimension %d - max %f \n", i, cluster_bdry[i*2]);
        }
    }

//    printf("-------Cluster Centroid---------\n");
//    for(int i=0; i<dim; i++){
//        printf("Dimension %d: %f", i, cluster_centroid[i]);
//    }
}

void printResult(int dim, int ndata, double *data, int kk, int *cluster_start, int *cluster_size, double **cluster_bdry, double **cluster_centroid){
    for (int i=0; i< kk; i++) {
        printf("------------------cluster %d ----------------- \n", i);
        printCluster(dim, i, cluster_start[i], cluster_size[i], cluster_bdry[i], cluster_centroid[i], data);
    }
}


//KDTree thingies
int *allocateTempIndices(int kk, int current_kk) {
    int *temp_indices = (int*)allocatePointerVariable(current_kk);
    temp_indices[0] = 0;
    for(int i = 1; i< current_kk; i++) {
        temp_indices[i] = temp_indices[i-1]+kk/current_kk;
    }
    return temp_indices;
}

int calculateClusterCentroids(int dim, double *data, int cluster_start, int cluster_size, double *cluster_centroid, double *cluster_bdry) {
//invariant: mean x = 1/n sum(x1+x2+...)
    for(int i = 0; i< cluster_size; i++) {
        int index = cluster_start/dim + i;
        struct dataPoint e = getElement(dim, index, data);
        for(int j = 0; j< dim; j++) {
            cluster_centroid[j]+= e.data[j];
        }
    }

    for(int i = 0; i< dim ; i++) {
        cluster_centroid[i] = cluster_centroid[i] /cluster_size;
    }

    printf("-----------------------cluster centroid---------------\n");
    for(int i =0; i<dim; i++) {
        printf("%f \n", cluster_centroid[i]);
    }
}


int findMaxVariant(int dim, double *data, int cluster_start, int cluster_size, double* cluster_centroid) {
//invariant: calculate variants of each dimension using formula var x = 1/n (sum (xi- |x|)^2)
    double *variants = (double*) malloc(dim * sizeof(double));
    for(int i = 0; i< dim; i++) {
        variants[i] = 0;
    }
    for(int i = 0; i< cluster_size; i++) {
        int index = cluster_start/dim + i;
        struct dataPoint element = getElement(dim, index, data);
        for (int j =0; j< dim; j++) {
            variants[j] += square(element.data[j] - cluster_centroid[j]);
        }
    }
    double max_variant = 0;
    int max_variant_index = 0;
    printf("------------------variants----------------- \n");
    for(int i = 0; i< dim; i++) {
        variants[i] = variants[i] /cluster_size;
        if (max_variant < variants[i]) {
            max_variant = variants[i];
            max_variant_index = i;
        }
        printf("%f    \n", variants[i]);
    }
    return max_variant_index;
}

void shiftElementUpward(int dim, double *data, int counter, int cluster_start, int elementIdx) {
    int swapElementIndex = cluster_start/dim + counter;
    struct dataPoint topElement = getElement(dim, swapElementIndex, data);
    struct dataPoint element = getElement(dim, elementIdx, data);

    for(int i = 0; i< dim; i++){
       data[swapElementIndex * dim+i] = element.data[i];
       data[elementIdx * dim + i] = topElement.data[i];
    }
}


int calculateClusterBoundaries(int dim, int cluster_idx, double *data, int cluster_start, int cluster_size, double **cluster_bdry) {
    double tempMin[dim];
    double tempMax[dim];
    for (int i=0; i< dim; i++) {
        tempMin[i] = 0;
        tempMax[i] = 0;
    }

    for (int i =0; i< cluster_size; i++ ) {
        struct dataPoint element = getElement(dim, cluster_start/dim + i, data);
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
        cluster_bdry[cluster_idx][cluster_bdry_idx] = tempMin[i];
        cluster_bdry_idx++;
        cluster_bdry[cluster_idx][cluster_bdry_idx] = tempMax[i];
        cluster_bdry_idx++;
    }
}

int calculateNewOutputs(int dim, int kk, int upper_cluster_size, int current_kk, int idx, double *data, int *cluster_start, int *cluster_size,
        double **cluster_bdry){
    //invariant: find new cluster_starts, cluster_sizes, cluster_bdrys for two new clusters based on cluster_indices
    int lower_cluster_size = cluster_size[idx] - upper_cluster_size;
    int lower_cluster_idx = idx + (kk/(current_kk * 2));

    //assign new cluster_starts
    cluster_start[lower_cluster_idx] = cluster_start[idx] + (upper_cluster_size*dim);

    //assign new cluster_sizes
    cluster_size[idx] = upper_cluster_size;
    cluster_size[lower_cluster_idx] = lower_cluster_size;

//    printf("-------------cluster_size---------------\n");
//    printf("upper cluster: %d\n", upper_cluster_size);
//    printf("lower cluster: %d\n", lower_cluster_size);
//
//    printf("-------------cluster_start---------------\n");
//    printf("lower cluster: %d\n", cluster_start[idx]);
//    printf("upper cluster: %d\n", cluster_start[lower_cluster_idx]);

    //calculate cluster boundaries:
    calculateClusterBoundaries(dim, idx, data, cluster_start[idx], cluster_size[idx], cluster_bdry);
    calculateClusterBoundaries(dim, lower_cluster_idx, data,  cluster_start[lower_cluster_idx], cluster_size[lower_cluster_idx], cluster_bdry);
}

int rearrangeData(int dim, int partition_dimension, double *data, int cluster_start, int cluster_size,const double* cluster_centroid) {
    //invariant: based on cluster_centroid and value of datapoint at dimension = partitionDimension, separate the cluster into 2 small clusters
    //return number of upper cluster
    int counter = 0;

    printf(" -----------partitionDimension--------------- \n");
    printf("%d \n", partition_dimension);

    for(int i= 0; i< cluster_size; i++) {
        int idx = cluster_start/dim + i;
        struct dataPoint element = getElement(dim, idx, data);
        if(cluster_centroid[partition_dimension] < element.data[partition_dimension]) {
            printf("%f --- %f\n", cluster_centroid[partition_dimension], element.data[partition_dimension]);
            shiftElementUpward(dim, data, counter, cluster_start, idx);
            counter ++;
        }
    }

    printf("counter: %d \n", counter);

    return counter;
}

int biPartition(int dim, double *data, int cluster_start, int cluster_size,
        double *cluster_bdry, double* cluster_centroid) {
    calculateClusterCentroids(dim, data, cluster_start, cluster_size, cluster_centroid, cluster_bdry);
    int partitionDimension = findMaxVariant(dim, data, cluster_start, cluster_size, cluster_centroid);
    int upper_cluster_size = rearrangeData(dim, partitionDimension, data, cluster_start, cluster_size, cluster_centroid);
    return upper_cluster_size;
}

int biPartitionClusters(int dim, int kk, int current_kk, int *cluster_indices, double *data, int *cluster_start, int *cluster_size,
        double **cluster_bdry, double **cluster_centroid) {
    for( int i = 0; i < current_kk; i++) {
        int idx = cluster_indices[i];
        int upper_cluster_size = biPartition(dim, data, cluster_start[idx], cluster_size[idx], cluster_bdry[idx],
                cluster_centroid[idx]);
        calculateNewOutputs(dim, kk, upper_cluster_size, current_kk, idx, data, cluster_start, cluster_size, cluster_bdry);
    }
}

int kdTree(int dim, int ndata, double *data, int kk, int *cluster_start, int *cluster_size,
           double **cluster_bdry, double **cluster_centroid) {

    int current_kk = 1;
    cluster_size[0] = ndata;
    cluster_start[0] = 0;

    while(current_kk < kk) {
        int *cluster_indices = allocateTempIndices(kk, current_kk);
        biPartitionClusters(dim, kk, current_kk, cluster_indices, data, cluster_start, cluster_size, cluster_bdry, cluster_centroid);
        current_kk *= 2;
    }
}

//generate dataSet
void generateDataPoint(int index, int dim, double *data) {
    for(int i=0; i<dim;i++) {
        data[index + i] = 1/(double)(rand()%20 +1);
    }
}

void generateData(int dim, int ndata, double *data) {
    for (int i=0; i< ndata*dim; i=i+dim) {
        generateDataPoint(i, dim, data);
    }
}

int main() {
    srand((unsigned)time(NULL));
    //declare input variables
    const int N_DATA = 10;
    const int DIM = 2;
    const int KK = 2;
    double *data;

    //declare output variables
    int *cluster_start, *cluster_size, *cluster_assign;
    cluster_start = allocatePointerVariable(KK);
    cluster_size = allocatePointerVariable(KK);
    cluster_assign = allocatePointerVariable(KK);

    double **cluster_bdry, **cluster_centroid;
    cluster_bdry = allocateDoublePointerVariable(KK, 2*DIM);
    cluster_centroid = allocateDoublePointerVariable(KK, DIM);

    //generate data set
    data = (double*)malloc(N_DATA*DIM*sizeof(double));
    generateData(DIM, N_DATA, data);
    printf("-----------------Starting Data Set---------------------\n");
    printData(DIM, N_DATA, data);

    //start KDTree
    kdTree(DIM, N_DATA, data, KK, cluster_start, cluster_size, cluster_bdry, cluster_centroid);

    printData(DIM, N_DATA, data);
    printResult(DIM, N_DATA, data,  KK, cluster_start, cluster_size, cluster_bdry, cluster_centroid);
}
