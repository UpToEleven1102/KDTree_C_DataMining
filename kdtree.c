#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

struct dataPoint {
    double *data;
};

//Utility thingies


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

void printData(int dim, int i0, int in, double *data) {
    struct dataPoint ele;
    for(int i = i0; i< in; i++) {
        ele = getElement(dim, i, data);
        printf("%d --- \n", i);
        for(int j = 0; j < dim; j++) {
            printf("%f", ele.data[j]);
        }
    }
}

void printResult(DIM, N_DATA, data,  KK, cluster_start, cluster_size, cluster_assign, cluster_bdry, cluster_centroid) {

}

//KDTree thingies
int calculateClusterCentroids(int dim, double *data, int cluster_start, int cluster_size, int *cluster_centroid, int *cluster_bdry) {
//invariant: mean x = 1/n sum(x1+x2+...)
    for(int i = 0; i< cluster_size; i++) {
        int index = cluster_start/dim + i;
        struct dataPoint e = getElement(dim, index, data);
        printf("--------\n");
        printData(dim, index, index+1, e.data);
    }
}

int biPartition(int dim, int i0, int in, double *data, int cluster_start, int cluster_size,
                double *cluster_bdry, double *cluster_centroid, int *cluster_assign) {
//what is i0, in, cluster_assign
    printf("No problem");
    calculateClusterCentroids(dim, data, cluster_start, cluster_size, cluster_centroid, cluster_bdry);
}

int kdTree(int dim, int ndata, double *data, int kk, int *cluster_start, int *cluster_size,
           double **cluster_bdry, double **cluster_centroid, int *cluster_assign) {
    int current_kk = 1;
    cluster_size[0] = ndata;
    cluster_start[0] = 0;

    while(current_kk < kk) {
        for(int i = 0; i < current_kk; i++) {
           // biPartition(dim, i, current_kk, data, cluster_start[i], cluster_size[i], cluster_bdry[i], cluster_centroid[i], cluster_assign[i]);
        }
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

/**this function doesn't work. wth???**/
//got it, pointer in C doesn't work like in GO. piff
void allocatePointerVariable(int kk, int *cluster_variable) {
    cluster_variable = (int*) malloc(kk*sizeof(int));
    for(int i = 0; i<kk; i++) {
        cluster_variable[i] = 0;
        printf("%d", cluster_variable[i]);
    }
}

int **allocateDoublePointerVariable(int kk, int array_size) {
    int **cluster_variable = (int**)malloc(kk*sizeof(double*));
    for(int i=0; i < kk; i++) {
        cluster_variable[i] = (int*) malloc(sizeof(array_size*sizeof(int)));
    }
    return cluster_variable;
}

int main() {
    srand((unsigned)time(NULL));
    //declare input variables
    const int N_DATA = 20;
    const int DIM = 8;
    const int KK = 8;
    double *data;

    //declare output variables
    int *cluster_start = (int*) malloc(KK*sizeof(int));
    int *cluster_size = (int*) malloc(KK*sizeof(int));
    int *cluster_assign = (int*)malloc(KK*sizeof(int));
    //allocatePointerVariable(KK, cluster_assign);

    double **cluster_bdry, **cluster_centroid;
    cluster_bdry = allocateDoublePointerVariable(KK, 2*DIM);
    cluster_centroid = allocateDoublePointerVariable(KK, DIM);

    //generate data set
    data = (double*)malloc(N_DATA*DIM*sizeof(double));
    generateData(DIM, N_DATA, data);
    //printData(DIM, N_DATA, data);
    //printf("-----------------\n");

    //start KDTree
    //kdTree(DIM, N_DATA, data, KK, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);
    //printResult(DIM, N_DATA, data,  KK, cluster_start, cluster_size, cluster_assign, cluster_bdry, cluster_centroid);
}
