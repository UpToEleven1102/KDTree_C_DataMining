#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include "kdTreeSearch.h"
#include "utilities.h"

//utilities:
void printSearchResults(int dim, int kk, const double *query, double *data, int *cluster_start, int *cluster_size,
                        double **cluster_bdry, double *cluster_distance, int data_point_idx,
                        double *shortest_distance) {
    //TODO: for testing purpose, might be deleted afterwards
    printf("\n\n\n=============Search Result=============\n");
    printf("--------Query Data Point------\n");

    for (int i = 0; i < dim; i++) {
        printf("%f\n", query[i]);
    }

    printf("--------Cluster Distance-------\n");
    for (int j = 0; j < kk; ++j) {
        printf("%f\n", cluster_distance[j]);
    }

    printf("--------Closest Data Point------\n");
    struct dataPoint e = getElement(dim, data_point_idx, data);
    printf("Index: %d - Distance: %f \n", data_point_idx, *shortest_distance);
    printArray(dim, 1, e.data);
}

//implement search method: find the data point that is closest to the query data point
double distanceToDataPoint(int dim, int idx, const double *query, const double *data) {
    //invariant: return the distance to the data point
    double distance = 0;
    struct dataPoint element = getElement(dim, idx, data);
    for (int i = 0; i < dim; ++i) {
        distance += square(query[i] - element.data[i]);
    }

    return sqrt(distance);
}

int findClosestDataPoint(int dim, int cluster_idx, const double *query, const double *data, int cluster_start,
                         int cluster_size, double *shortest_distance, int *datapoint_count) {
    //invariant: exhausted search over the cluster set. push distance to the closest point into *shortest distance, return the index with the shortest distance
    int closest_datapoint_idx = 0;
    *shortest_distance = 99999;

    for (int i = 0; i < cluster_size; ++i) {
        int idx = cluster_start + i;
        double distance = distanceToDataPoint(dim, idx, query, data);
        if (distance < *shortest_distance) {
            *shortest_distance = distance;
            closest_datapoint_idx = idx;
        }
        *datapoint_count = *datapoint_count +1;
    }

    return closest_datapoint_idx;
}

double distanceToCluster(int dim, const double *query,const double *cluster_bdry) {
    //invariant: calculate distance of data point to the cluster using formula: 1/dim *(sum di^2) with di = {xi - Mi//x > Mi, 0 //mi<x<Mi, xi - mi //xi < mi}
    double distance = 0;
    double min , max;
    for (int i = 0; i < dim * 2; i = i + 2) {
        int j = i + 1;
        min = cluster_bdry[i];
        max = cluster_bdry[j];
        double d;
        if (query[i / 2] < min) {
            d = min - query[i / 2];
        } else if (query[i / 2] > max) {
            d = query[i / 2] - max;
        } else {
            d = 0;
        }

        distance += d*d;
    }
    return sqrt(distance);
}

int findClosestCluster(int dim, int kk, const double *query, double **cluster_bdry, double *cluster_distance) {
    //invariant: push the distance to each cluster to *cluster_distance, return the index of the cluster that closest to the query
    int minIndex = 0;

    for (int i = 0; i < kk; i++) {
        cluster_distance[i] = distanceToCluster(dim, query, cluster_bdry[i]);
    }

    for (int i = 0; i < kk; ++i) {
        if (cluster_distance[i] < cluster_distance[minIndex]) {
            minIndex = i;
        }
    }

    return minIndex;
}

int kdTreeSearch(int dim, int kk, const double *query, double *data, int *cluster_start, int *cluster_size,
                     double **cluster_bdry, double *result_ptr) {
    //invariant: find the cluster that is closest to the query data point, exhausted search the cluster to find the closest data point.
    // Verify the closest data point with neighbor clusters
    double *cluster_distance = (double *) malloc(kk * sizeof(double));
    int *datapoint_count = (int*)malloc(sizeof(int));
    *datapoint_count = 0;
    int closest_cluster_idx = findClosestCluster(dim, kk, query, cluster_bdry, cluster_distance);

    double *shortest_distance = (double *) malloc(sizeof(double));
    int closest_datapoint_idx = findClosestDataPoint(dim, closest_cluster_idx, query, data,
                                                     cluster_start[closest_cluster_idx],
                                                     cluster_size[closest_cluster_idx],
                                                     shortest_distance, datapoint_count);

    for (int i = 0; i < kk; i++) {
        if (cluster_distance[i] < *shortest_distance && i != closest_cluster_idx) {
            double *temp_distance = (double *) malloc(sizeof(double));
            int temp_idx = findClosestDataPoint(dim, i, query, data, cluster_start[i], cluster_size[i], temp_distance, datapoint_count);
            if (*temp_distance < *shortest_distance) {
                *shortest_distance = *temp_distance;
                closest_datapoint_idx = temp_idx;
            }
        }
    }

    struct dataPoint e = getElement(dim, closest_datapoint_idx, data);
    for (int i = 0; i < dim; ++i) {
        result_ptr[i] = e.data[i];
    }

    printSearchResults(dim, kk, query, data, cluster_start, cluster_size, cluster_bdry,
                       cluster_distance, closest_datapoint_idx, shortest_distance);
    return *datapoint_count;
}
