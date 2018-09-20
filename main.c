#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

struct dataPoint {
    double *data;
};

int *allocatePointerVariable(int array_size) {
    int *cluster_variable = (int *)malloc(array_size*sizeof(int));
    for(int i=0; i<array_size; i++) {
        cluster_variable[i]= -1;
    }
    return cluster_variable;
}

int **allocateDoublePointerVariable(int kk, int array_size) {
    int **cluster_variable = (int**)malloc(kk*sizeof(double*));
    for(int i=0; i < kk; i++) {
        cluster_variable[i] = (int*) malloc(sizeof(array_size*sizeof(int)));
    }
    return cluster_variable;
}

//Utility thingies
double square(double i) {
    return i*i;
}

struct dataPoint emptyDataPoint(int dim) {
    struct dataPoint out;
    out.data = (double*)malloc(dim * sizeof(double));
    return out;
};

struct dataPoint getElement(int dim, int index, double *data) {
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

//KDTree thingies

int calculateClusterBoundary() {

}

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
        int index = cluster_start + i;
        struct dataPoint element = getElement(dim, index, data);
        for (int j =0; j< dim; j++) {
            variants[j] += square(element.data[j] - cluster_centroid[j]);
        }
    }
    double max_variant = 0;
    int max_variant_index = 0;
    printf("------------------variants----------------- \n");
    for(int i = 0; i< dim; i++) {
        variants[i] /= cluster_size;
        if (max_variant < variants[i]) {
            max_variant = variants[i];
            max_variant_index = i;
        }
        printf("%f    \n", variants[i]);
    }
    return max_variant_index;
}

int biPartition(int dim, int cluster_idx, int current_kk, double *data, int cluster_start, int cluster_size,
        double *cluster_bdry, double* cluster_centroid, int *cluster_assign) {
    calculateClusterCentroids(dim, data, cluster_start, cluster_size, cluster_centroid, cluster_bdry);
    int partitionDimension = findMaxVariant(dim, data, cluster_start, cluster_size, cluster_centroid);

    printf(" -----------partitionDimension--------------- \n");
    printf("%d \n", partitionDimension);
}

int biPartitionClusters(int dim, int current_kk, int *cluster_indices, double *data, int *cluster_start, int *cluster_size,
        double **cluster_bdry, double **cluster_centroid, int *cluster_assign) {
    for( int i = 0; i < current_kk; i++) {
        int idx = cluster_indices[i];
        biPartition(dim, idx, current_kk, data, cluster_start[idx], cluster_size[idx], cluster_bdry[idx], cluster_centroid[idx], cluster_assign);
    }
}



int kdTree(int dim, int ndata, double *data, int kk, int *cluster_start, int *cluster_size,
           double **cluster_bdry, double **cluster_centroid, int *cluster_assign) {

    int current_kk = 1;
    cluster_size[0] = ndata;
    cluster_start[0] = 0;

//    while(current_kk < kk) {
        int *cluster_indices = allocateTempIndices(kk, current_kk);
        biPartitionClusters(dim, current_kk, cluster_indices, data, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);
        current_kk *= 2;
//    }
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
    const int N_DATA = 5;
    const int DIM = 4;
    const int KK = 4;
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
    printData(DIM, N_DATA, data);

    //start KDTree
    kdTree(DIM, N_DATA, data, KK, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);
    //printResult(DIM, N_DATA, data,  KK, cluster_start, cluster_size, cluster_assign, cluster_bdry, cluster_centroid);
}
