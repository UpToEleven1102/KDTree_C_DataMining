#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int calculateClusterCentroids() {

}

void generateClusterBoundary(int dim, int ndata, double **cluster_bdry) {
    cluster_bdry = (double**) malloc(ndata*sizeof(double));
    for (int i=0 ; i< ndata; i++) {
        cluster_bdry[i] = (double*) malloc(dim* sizeof(double));
    }
}

int kdTree(int dim, int ndata, double *data, int kk, int *cluster_start, int *cluster_size,
           double **cluster_bdry, double **cluster_centroid, int *cluster_assign) {
    int clusterNum = 1;
    while( clusterNum< kk ) {
        for(int i=0; i<clusterNum; i++) {
            biPartition(dim, )
        }
    }
}

int biPartition(int dim, int i0, int in, double *data, int cluster_start, int cluster_size,
                double *cluster_bdry, double *cluster_centroid, int *cluster_assign) {

}

//int biPartition(int dim, int i0, int in, double *data, int cluster_start[2], int cluster_size[2],
//                double *cluster_bdry[2], double *cluster_centroid[2], int *cluster_assign) {

//}


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

void printData(int dim, int ndata, double *data) {
    for(int i=0; i< dim * ndata; i++) {
        printf("%f \n", data[i]);
    }
}

void printResult(DIM, N_DATA, data,  KK, cluster_start, cluster_size, cluster_assign, cluster_bdry, cluster_centroid) {

}

int main() {
    srand((unsigned)time(NULL));
    //declare input variables
    const int N_DATA = 100;
    const int DIM = 8;
    const int KK = 8;
    double *data;

    //declare output variables
    int *cluster_start, *cluster_size, *cluster_assign;
    double **cluster_bdry, **cluster_centroid;

    //generate data set
    data = (double*)malloc(N_DATA*DIM*sizeof(double));
    generateData(DIM, N_DATA, data);
    printData(DIM, N_DATA, data);


    //start KDTree
    kdTree(DIM, N_DATA, data, KK, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);

    printResult(DIM, N_DATA, data,  KK, cluster_start, cluster_size, cluster_assign, cluster_bdry, cluster_centroid);
}
