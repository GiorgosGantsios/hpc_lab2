/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         seq_kmeans.c  (sequential version)                        */
/*   Description:  Implementation of simple k-means clustering algorithm     */
/*                 This program takes an array of N data objects, each with  */
/*                 M coordinates and performs a k-means clustering given a   */
/*                 user-provided value of the number of clusters (K). The    */
/*                 clustering results are saved in 2 arrays:                 */
/*                 1. a returned array of size [K][N] indicating the center  */
/*                    coordinates of K clusters                              */
/*                 2. membership[N] stores the cluster center ids, each      */
/*                    corresponding to the cluster a data object is assigned */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department, Northwestern University                        */
/*            email: wkliao@ece.northwestern.edu                             */
/*                                                                           */
/*   Copyright (C) 2005, Northwestern University                             */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "kmeans.h"


/*----< euclid_dist_2() >----------------------------------------------------*/
/* square of Euclid distance between two multi-dimensional points            */
__inline static
float euclid_dist_2(int    numdims,  /* no. dimensions */
                    float *coord1,   /* [numdims] */
                    float *coord2)   /* [numdims] */
{
    int i;
    float ans=0.0;

    // #pragma omp parallel for schedule(static) reduction(+:ans)
    for (i=0; i<numdims; i++){
        //printf("numthreads: %d\n", omp_get_num_threads());
        float diff = coord1[i] - coord2[i];
        ans += diff * diff;
    }

    return(ans);
}

/*----< find_nearest_cluster() >---------------------------------------------*/
__inline static
int find_nearest_cluster(int     numClusters, /* no. clusters */
                         int     numCoords,   /* no. coordinates */
                         float  *object,      /* [numCoords] */
                         float **clusters)    /* [numClusters][numCoords] */
{
    int   index, i;
    float dist, min_dist;

    /* find the cluster id that has min distance to object */
    index    = 0;
    min_dist = euclid_dist_2(numCoords, object, clusters[0]);

    for (i=1; i<numClusters; i++) {
        dist = euclid_dist_2(numCoords, object, clusters[i]);
        /* no need square root */
        if (dist < min_dist) { /* find the min and its array index */
            min_dist = dist;
            index    = i;
        }
    }
    return(index);
}

/*----< seq_kmeans() >-------------------------------------------------------*/
/* return an array of cluster centers of size [numClusters][numCoords]       */
int seq_kmeans(float **objects,      /* in: [numObjs][numCoords] */
               int     numCoords,    /* no. features */
               int     numObjs,      /* no. objects */
               int     numClusters,  /* no. clusters */
               float   threshold,    /* % objects change membership */
               int    *membership,   /* out: [numObjs] */
               float **clusters)     /* out: [numClusters][numCoords] */

{
    int      i, j, index, loop=0;
    int     *newClusterSize; /* [numClusters]: no. objects assigned in each
                                new cluster */
    float    delta;          /* % of objects change their clusters */
    float  **newClusters;    /* [numClusters][numCoords] */

    /* initialize membership[] */
    for (i=0; i<numObjs; i++) membership[i] = -1;

    /* need to initialize newClusterSize and newClusters[0] to all 0 */
    newClusterSize = (int*) calloc(numClusters, sizeof(int));
    assert(newClusterSize != NULL);

    newClusters    = (float**) malloc(numClusters *            sizeof(float*));
    assert(newClusters != NULL);
    newClusters[0] = (float*)  calloc(numClusters * numCoords, sizeof(float));
    assert(newClusters[0] != NULL);

    for (i=1; i<numClusters; i++)
        newClusters[i] = newClusters[i-1] + numCoords;

do {
    delta = 0.0;
    // #pragma omp parallel for private(index, j) reduction(+:delta) reduction(+:newClusterSize[:numClusters])
    for (i = 0; i < numObjs; i++) {
        /* find the nearest cluster center */
        index = find_nearest_cluster(numClusters, numCoords, objects[i], clusters);
        /* if membership changes, increase delta by 1 */
        if(membership[i] != index){
            delta += 1;
            /* assign the membership to object i */
            membership[i] = index;
        }
        
        /* Update cluster centers using atomic operations for thread safety */
        // #pragma omp critical
        newClusterSize[index]++;

        #pragma omp parallel for num_threads(2) schedule(static)
        for (j = 0; j < numCoords; j++) {
            newClusters[index][j] += objects[i][j];
        }
    }

    /* Average the sum and replace old cluster centers with newClusters */
    // #pragma omp parallel for num_threads(4)
    for (i = 0; i < numClusters; i++) {
        int newClusterSizeTemp = newClusterSize[i];
        #pragma omp parallel for num_threads(2) private(newClusterSizeTemp)
        for (j = 0; j < numCoords; j++) {
            // printf("%d\n",omp_get_thread_num());
            if (newClusterSizeTemp > 0)
                clusters[i][j] = newClusters[i][j] / newClusterSizeTemp;
            newClusters[i][j] = 0.0;   /* set back to 0 */
        }
        newClusterSize[i] = 0;   /* set back to 0 */
    }

    delta /= numObjs;
} while (delta > threshold && loop++ < 500);


    free(newClusters[0]);
    free(newClusters);
    free(newClusterSize);

    return 1;
}

