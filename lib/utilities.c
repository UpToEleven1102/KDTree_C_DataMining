#include "utilities.h"
#include <stdio.h>
#include <malloc.h>

int *allocatePointerVariable(int array_size) {
    int *cluster_variable = (int *)malloc(array_size*sizeof(int));
    for(int i=0; i<array_size; i++) {
        cluster_variable[i]= -1;
    }
    return cluster_variable;
}

double **allocateDoublePointerVariable(int kk, int array_size) {
    double **cluster_variable = (double**)malloc(kk*sizeof(double*));
    for(int i=0; i < kk; i++) {
        cluster_variable[i] = (double*) malloc(array_size*sizeof(double));
        for(int j=0; j< array_size; j++){
            cluster_variable[i][j] = 0;
        }
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

int printCluster(int dim, int cluster_start, int cluster_size, double *cluster_bdry, double *cluster_centroid,  double *data) {
    printf("Cluster Start: %d, Cluster Size: %d \n", cluster_start, cluster_size);

    if (cluster_size == 0) {
        return 0;
    }
    printf("---------Data Points--------\n");
    for (int i = 0; i < cluster_size; i++) {
        struct dataPoint e = getElement(dim, cluster_start + i, data);
        printf("%d-----\n", i);
        printArray(dim, 1, e.data);
    }

    printf("-------Cluster Boundaries-------\n");
    for(int i=0; i< 2*dim; i++) {
        if(i%2 ==0){
            printf("min %f \n", cluster_bdry[i]);
        } else {
            printf("max %f \n", cluster_bdry[i]);
        }
    }

    printf("-------Cluster Centroid---------\n");
    for(int i=0; i<dim; i++){
        printf("Dimension %d: %f\n", i, cluster_centroid[i]);
    }
}

void printResult(int dim, double *data, int kk, int *cluster_start, int *cluster_size, double **cluster_bdry, double **cluster_centroid){
    printf("\n\n\nFinal KD TREE: \n");
    for (int i=0; i< kk; i++) {
        printf("------------------cluster %d ----------------- \n", i);
        printCluster(dim, cluster_start[i], cluster_size[i], cluster_bdry[i], cluster_centroid[i], data);
    }
}
