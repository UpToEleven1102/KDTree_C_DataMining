#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "lib/utilities.h"
#include "lib/kdTree.h"
#include "lib/kdTreeSearch.h"

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
    printf("-----------------Starting Data Set---------------------\n");
    printData(DIM, N_DATA, data);

    //start KDTree
    kdTree(DIM, N_DATA, data, KK, cluster_start, cluster_size, cluster_bdry, cluster_centroid);
    printResult(DIM, data,  KK, cluster_start, cluster_size, cluster_bdry, cluster_centroid);

    //start Search
    double *query = (double*) malloc(DIM * sizeof(double));
    generateDataPoint(0, DIM, query);

    kdTreeSearch(DIM, KK, query, data, cluster_start, cluster_size, cluster_bdry, cluster_centroid);

}
