//
// Created by Huyen VU on 9/22/2018.
//

#include <malloc.h>
#include <stdio.h>
#include "kdTreeSearch.h"


//utilities:
void printSearchResults(int dim, int kk, const double *query, double *data, int *cluster_start, int *cluster_size, double **cluster_bdry, double **cluster_centroid, double *cluster_distance) {
    //TODO: for testing purpose, might be deleted afterwards
    printf("\n\n\n=============Search Result=============\n");
    printf("--------Query Data Point------\n");

    for(int i=0; i<dim; i++) {
        printf("%f\n", query[i]);
    }

    printf("--------Cluster Distance-------\n");
    for (int j = 0; j < kk; ++j) {
        printf("%f\n", cluster_distance[j]);
    }
}

double square(double n) {
    return n*n;
}

//implement search method: find the data point that is closest to the query data point
double distanceToDataPoint(int dim, int idx, double *query, const double *data){

    for (int i = 0; i < dim; ++i) {

    }
}

int findClosestDataPoints(int dim, int cluster_idx, const double *data, int cluster_start, int cluster_size, double *shortest_distance) {
    //invariant: exhausted search over the cluster set. push distance to the closest point into *shortest distance, return the index with the shortest distance
    int closest_datapoint_idx = 0;
    *shortest_distance = 99999;



    return closest_datapoint_idx;
}

double distanceToCluster(int dim, const double *query, double *cluster_bdry) {
    //invariant: calculate distance of data point to the cluster using formula: 1/dim *(sum di^2) with di = {xi - Mi//x > Mi, 0 //mi<x<Mi, xi - mi //xi < mi}
    double distance = 0;
    for (int i = 0; i < dim*2; i=i+2) {
        int j = i+1;
        double min = cluster_bdry[i];
        double max = cluster_bdry[j];
        double d;
        if(query[i/2] < min) {
            d = min - query[i/2];
        } else if(query[i/2] > max) {
            d = query[i/2] - max;
        } else {
            d = 0;
        }

        distance += square(d);
    }
    return distance/dim;
}

int findClosestCluster(int dim, int kk,  const double *query, double **cluster_bdry, double *cluster_distance) {
    //invariant: push the distance to each cluster to *cluster_distance, return the index of the cluster that closest to the query
    int minIndex =0;

    for (int i =0; i< kk; i++ ) {
        cluster_distance[i] = distanceToCluster(dim, query, cluster_bdry[i]);
    }

    for (int i = 0; i < kk; ++i) {
        if(cluster_distance[i]< cluster_distance[0]){
            minIndex = i;
        }
    }

    return minIndex;
}

double *kdTreeSearch(int dim,int kk, const double *query, double *data, int *cluster_start, int *cluster_size, double **cluster_bdry, double **cluster_centroid){
    //invariant: find the cluster that is closest to the query data point, exhausted search the cluster to find the closest data point.
    // Verify the closest data point with neighbor clusters
    double *cluster_distance = (double*) malloc(kk * sizeof(double));
    int closest_cluster_idx = findClosestCluster(dim, kk, query, cluster_bdry, cluster_distance);

    double *shortest_distance = (double*) malloc(sizeof(double));
    int closest_datapoint_idx = findClosestDataPoints(dim, closest_cluster_idx, data, cluster_start[closest_cluster_idx], cluster_size[closest_cluster_idx],
            shortest_distance);

    printf("debug %d, %f", closest_datapoint_idx, *shortest_distance);

    for(int i =0; i<kk; i++) {
        if(cluster_distance[i] < *shortest_distance && i!= closest_cluster_idx) {
            //search the cluster
        }
    }

    printSearchResults(dim, kk, query, data, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_distance);
}
