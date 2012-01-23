/* bs.c
 *
 * Mark A. Ruzon
 * Stanford Vision Laboratory
 * 8 September 1998
 *
 * This code implements the basic Binary Splitting algorithm described in 
 * the 1991 IEEE Trans. on Sig. Proc. article "Color Quantization of Images"
 * by Michael Orchard and Charles Bouman, pp. 2677-90.
 *
 * The input is a 1-D array of data points.  The #define'd constant DIM 
 * holds the number of dimensions.  The output is a set of cluster centers
 * and is also a 1-D array.  The number of clusters is also returned.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bs.h"

#define MAX_ITERATIONS 50
#define MIN_LAMBDA 100.0             /* Do not go below 3.0 */
#define TRUE 1
#define FALSE 0
#define DEBUG 0
#define BS_UNDERFLOW 1E-4

/* An internal data structure */

typedef struct cluster_data {
  Coord R[DIM][DIM];               /* Summation matrix for covariance */
  Coord m[DIM];                    /* summation matrix for mean */
  Coord N;                         /* Total weight of cluster */
  Coord lambda;                    /* Largest eigenvalue of covariance */
  Coord e[DIM];                    /* Eigenvector corresponding to lambda */
} Cluster_Data;


void ComputeEigenStuff(Cluster_Data *ptr)
{
  int j, k, iterations = 0, iszero = 1;
  Coord cov[DIM][DIM], q[DIM], x[DIM], xold[DIM];
  float ratio, oldratio, max;

  /* Compute local data used for eigenvalues */
  for (j = 0; j < DIM; j++) {
    for (k = 0; k < DIM; k++)
      if (k >= j)
	cov[j][k] = ptr->R[j][k] - ptr->m[j] * ptr->m[k] / ptr->N;
      else
	cov[j][k] = cov[k][j];
    q[j] = ptr->m[j] / ptr->N;
    xold[j] = 1.0 / (float) (j + 1);
  }

  if (DEBUG)
    for (j = 0; j < DIM; j++) {
      for (k = 0; k < DIM; k++)
	fprintf(stderr,"%10.4f ",cov[j][k]);
      fprintf(stderr," | %10.4f\n",q[j]);
    }

  /* Find out if the diagonal is zero (cluster of similar pixels) */
  for (j = 0; j < DIM; j++)
    if (abs(cov[j][j]) > BS_UNDERFLOW) {
      iszero = 0;
      break;
    }
  if (iszero) {
    ptr->lambda = 0.0;
    return;
  }

  /* Divide by the element with the largest absolute value */
  max = 0.0;
  for (j = 0; j < DIM; j++) 
    for (k = 0; k < DIM; k++)
      if (fabs(cov[j][k]) > fabs(max))
	max = cov[j][k];
  for (j = 0; j < DIM; j++) 
    for (k = 0; k < DIM; k++)
      cov[j][k] /= max;
   
  /* Use power method to compute eigenvector corresponding to largest 
     eigenvalue */
  oldratio = 1E10;
  while (iterations < MAX_ITERATIONS) {
    for (j = 0; j < DIM; j++) {
      x[j] = 0.0;
      for (k = 0; k < DIM; k++)
	x[j] += cov[j][k] * xold[k];
    }
    /* Set ratio to maximum of ratios of the components */
    ratio = (xold[0] == 0) ? 0 : x[0] / xold[0];
    for (j = 1; j < DIM; j++)
      if (xold[j] && x[j] / xold[j] > ratio)
	ratio = x[j] / xold[j];
    /* Check successive ratios */
    if (ratio / oldratio > 0.9999 && ratio / oldratio < 1.0001) 
      break;
    else {
      iterations++;
      for (j = 0; j < DIM; j++)
	xold[j] = x[j];
      oldratio = ratio;
    }
  }

  ptr->lambda = ratio * max;
  for (j = 0; j < DIM; j++)
    ptr->e[j] = x[j];

  if (DEBUG)
    if (iterations == MAX_ITERATIONS)
      fprintf(stderr,"Warning: Maximum iterations reached in power method; ratio = %.2f\n",ratio * max);
    else
      fprintf(stderr,"%d iterations\n",iterations);
}


/* This routine calculates summations of second-order statistics of the 
 * pixels with the correct label and stores the results in the Cluster_Data
 * structure. 
 */
void CalculateClusterData(Cluster_Data *ptr, Coord *points, int npoints, 
			  int *index, int label)
{
  int i, j, k;

  /* Compute data to be put in the data structure */
  for (i = 0; i < npoints; i++)
    if (index[i] == label) {
      for (j = 0; j < DIM; j++) {
	for (k = 0; k <= j; k++)
	  ptr->R[j][k] += points[i * DIM + j] * points[i * DIM + k];
	ptr->m[j] += points[i * DIM + j];
      }
      ptr->N++;
    }

  /* Copy the lower half of the matrix to the upper half */
  for (j = 0; j < DIM; j++) 
    for (k = j + 1; k < DIM; k++)
      ptr->R[j][k] = ptr->R[k][j];

  ComputeEigenStuff(ptr);
}


/* This routine creates a new cluster in the space marked by child.
 */
int Split_Cluster(Cluster_Data *parent, Coord *points, int npoints, 
		   int *index, int parent_label, int child_label)
{
  int i, j, parent_notempty = FALSE, child_notempty = FALSE;
  float sum, threshold = 0.0;

  /* First compute the threshold (function of eigenvector and mean) */
  for (i = 0; i < DIM; i++)
    threshold += parent->m[i] * parent->e[i];
  threshold /= parent->N;

  /* Compare each point with the parent's label with the threshold */
  for (i = 0; i < npoints; i++) 
    if (index[i] == parent_label) {
      sum = 0.0;
      for (j = 0; j < DIM; j++)
	sum += points[i * DIM + j] * parent->e[j];
      if (sum > threshold) {
	index[i] = child_label;
	child_notempty = TRUE;
      } else
	parent_notempty = TRUE;
    }

  return(child_notempty && parent_notempty);    
}


/* Recompute the parent's auxilliary information using the child that just
 * split off. 
 */
void CalcluateComplementaryClusterData(Cluster_Data *parent, 
				       Cluster_Data *child)
{
  int i, j;

  parent->N -= child->N;

  if (DEBUG)
    fprintf(stderr,"Cluster of %3.0f pixels split into %3.0f and %3.0f\n",
	    parent->N+child->N, parent->N,child->N);

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++)
      parent->R[i][j] -= child->R[i][j];
    parent->m[i] -= child->m[i];
  }

  ComputeEigenStuff(parent);
}


void bs(Coord *points, int npoints, int maxclusters, int *nclusters, 
	Coord **clusters, int *index)
{
  int maxlambda_index, i, j, good_split, index_alloc = 0;
  float maxlambda;
  Cluster_Data *cd;

  if (DEBUG)
    mexPrintf("Number of points: %d, max: %d, MAX: %d\n",npoints,
	      MAXCLUSTERS,maxclusters);

  if (npoints < maxclusters)
    maxclusters = npoints;
  if (maxclusters > MAXCLUSTERS)
    maxclusters = MAXCLUSTERS;

  /* Index says what cluster a pixel belongs to; if it is already allocated,
   * the user wants the index returned; otherwise, we must allocate it 
   * ourselves.
   */
  if (!index) {
    index_alloc = 1;
    index = (int *)calloc(npoints, sizeof(int));
    if (!index) { /* Still */
      fprintf(stderr,"Error in binary split: index could not be allocated\n");
      return;
    }
  }

  /* cd holds auxilliary data about each cluster. */
  cd = (Cluster_Data *)calloc(maxclusters, sizeof(Cluster_Data));
  if (!cd) {
    fprintf(stderr,"Error in binary split: cd could not be allocated\n");
    return;
  }

  CalculateClusterData(&cd[0], points, npoints, index, 0);
  if (DEBUG)
    fprintf(stderr,"Total points: %d / %f\n",npoints,cd[0].N);
  *nclusters = 1;
  while (*nclusters < maxclusters) {  /* Split */

    maxlambda = -100000.0;
    for (i = 0; i < *nclusters; i++)
      if (cd[i].lambda > maxlambda) {
	maxlambda = cd[i].lambda;
	maxlambda_index = i;
      }

    if (DEBUG) {
      fprintf(stderr,"Lambda maximum #%d: %10.4f\n",maxlambda_index,maxlambda);
      fprintf(stderr,"New cluster #%d\n",*nclusters);
    }

    /* Check to make sure there's something to split */
    if (maxlambda > MIN_LAMBDA)
      good_split = Split_Cluster(&cd[maxlambda_index], points, npoints, 
				 index, maxlambda_index, *nclusters);
    else
      break;

    if (good_split) {
      CalculateClusterData(&cd[*nclusters], points, npoints, index, 
			   *nclusters);
      CalcluateComplementaryClusterData(&cd[maxlambda_index], &cd[*nclusters]);
      ++*nclusters;
    } else
      break;
  }

  /* clusters holds the output */
  (*clusters) = (Coord *)calloc(*nclusters * DIM, sizeof(Coord));
  if (!*clusters) {
    fprintf(stderr,"Error in binary split: clusters could not be allocated\n");
    return;
  }

  /* Transfer the data to the output data structure */
  for (i = 0; i < *nclusters; i++) {
    for (j = 0; j < DIM; j++)
      (*clusters)[i * DIM + j] = cd[i].m[j] / cd[i].N;
    if (DEBUG)
      fprintf(stderr,"(%.1f,%.1f) ",cd[i].N,cd[i].lambda);
    }
  if (DEBUG)
    fprintf(stderr,"\nNumber of clusters: %d\n",*nclusters);

  free(cd);
  if (index_alloc)
    free(index);
}
