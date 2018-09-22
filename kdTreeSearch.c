//
// Created by Huyen VU on 9/22/2018.
//

#include <malloc.h>
#include <stdio.h>
#include "kdTreeSearch.h"

//implement search method: find the data point that is closest to the query data point
double distanceToCluster(int dim, const double *query, double *data, int cluster_start, int cluster_size, double *cluster_bdry, double *cluster_centroid) {
    // calculate distance of data point to the cluster using formula: 1/dim *(sum di^2) with di = {xi - Mi//x > Mi, 0 //mi<x<Mi, xi - mi //xi < mi}
    printf("testing");

}

double *kdTreeSearch(int dim,int kk, const double *query, double *data, int *cluster_start, int *cluster_size, double **cluster_bdry, double **cluster_centroid){
    //invariant: find the cluster that is closest to the query data point, exhausted search the cluster to find the closest data point.
    // Verify the closest data point with neighbor clusters
    double *cluster_distance = (double*) malloc(kk * sizeof(double));

    for (int i =0; i< kk; i++ ) {
        cluster_distance[i] = distanceToCluster(dim, query, data, cluster_start[i], cluster_size[i], cluster_bdry[i], cluster_centroid[i]);
    }


}