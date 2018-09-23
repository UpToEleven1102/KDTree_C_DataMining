#ifndef KD_TREE_KDTREESEARCH_H
#define KD_TREE_KDTREESEARCH_H
double *kdTreeSearch(int dim,int kk, const double *query, double *data, int *cluster_start, int *cluster_size, double **cluster_bdry, double **cluster_centroid);
#endif //KD_TREE_KDTREESEARCH_H
