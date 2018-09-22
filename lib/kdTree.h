//
// Created by Huyen VU on 9/22/2018.
//

#ifndef KD_TREE_KDTREE_H
#define KD_TREE_KDTREE_H

int kdTree(int dim, int ndata, double *data, int kk, int *cluster_start, int *cluster_size,
           double **cluster_bdry, double **cluster_centroid);

#endif //KD_TREE_KDTREE_H
