#ifndef KD_TREE_UTILITIES_H
#define KD_TREE_UTILITIES_H
struct dataPoint {
    double *data;
};


int *allocatePointerVariable(int array_size);
double **allocateDoublePointerVariable(int kk, int array_size);
double square(double i);
struct dataPoint emptyDataPoint(int dim);
struct dataPoint getElement(int dim, int index, const double *data);
void printArray(int dim, int ndata, double *data);
void printData(int dim, int ndata, double *data);
int printCluster(int dim, int cluster_start, int cluster_size, double *cluster_bdry, double *cluster_centroid,  double *data);
void printResult(int dim, double *data, int kk, int *cluster_start, int *cluster_size, double **cluster_bdry, double **cluster_centroid);
#endif //KD_TREE_UTILITIES_H
