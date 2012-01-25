/***************************************************************************
 * The Compass Operator -- an algorithm to find edges and corners in color *
 *                   images by representing regions with distributions     *
 *                                                                         *
 * Mark A. Ruzon                                                           *
 * Stanford Vision Laboratory                                              *
 * October, 1997                                                           *
 * Copyright 1997-1999.  All rights reserved.                              *
 *                                                                         *
 * Details -- The compass operator uses a circular window centered at a    *
 *   junction where 4 pixel squares meet.  The needle is a diameter at a   *
 *   given orientation.  The color distributions of the two semicircles    *
 *   are computed, and the distance between them is measured using the     *
 *   Earth Mover's Distance (EMD).  The maximum EMD over all orientations  *
 *   gives the edge strength and orientation at that point.  Uncertainty   *
 *   and abnormality are also measured.                                    *
 *                                                                         *
 ***************************************************************************/

#include <string.h> // for memset
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ruby.h>
#include "compass.h"
#include "bs.h"
#include "emd.h"

#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define PI 3.14159265358979323846
#define MAXRADIUS 64
#define PERCEPTUAL_THRESH 11.0
#define GAMMA 14.0
#define BLACK 40
#define SAMPLE_PERCENTAGE 1.0
#define WEDGEWT 1.0
#define DEBUG 0
#define PARTIAL_MATCH 1
#define PRINT_LOC 0
#define PRINT_WEDGES 0
#define PRINT_GAUSSIAN 0
#define PRINT_MASK 0
#define PRINT_MASKSUM 0
#define PRINT_SAMPLING 0
#define PRINT_CLUSTERS 0
#define PRINT_COST 0
#define PRINT_HIST 0

/* This global is for the EMD cost matrix */
double cost[MAXCLUSTERS][MAXCLUSTERS];

void RGB2Lab(float, float, float, float *, float *, float *);

/**************************************************/
/***     Level 3 Functions                      ***/
/**************************************************/

/* This routine implements the MATLAB routine CArea.m, which integrates a
 * circle of radius r centered at the origin from Xlow to Xhigh in Quadrant
 * I.  Since this function is only interested in square pixels, the Y value 
 * is included so that the area below Y is ignored.
 */
double CArea(double Xhigh, double Xlow, double Y, double r)
{
  return(0.5 * ((Xhigh * sqrt(r*r - Xhigh * Xhigh) + r*r * asin(Xhigh / r)) -
	 (Xlow * sqrt(r*r - Xlow * Xlow) + r*r * asin(Xlow / r))) -
	 Y * (Xhigh - Xlow));
}


/**************************************************/
/***     Level 2 Functions                      ***/
/**************************************************/

/* This routine abstracts a little of the complexity of placing pixels and
 * their weights into the unclustered array 
 */
void PlacePoint(float *L, float *a, float *b, int row, int col, int nrows, 
		float *point_array, float *weight_array, float weight)
{
  int index = col * nrows + row;

  *point_array = L[index];
  *(point_array+1) = a[index];
  *(point_array+2) = b[index];
  *weight_array = weight;
}


/* Same thing here */
void AddWeight(float *L, float *a, float *b, double weight, int row,
	       int col, int nrows, float data[][DIM], double *hist,
	       float com[][DIM], float sum[], int nclusters)
{
  int index = col * nrows + row, min_index, i;
  double dist, Ldist, adist, bdist, min_dist = 100 * 240 * 240;

  min_index = 0;
  for (i = 0; i < nclusters; i++) {
    Ldist = L[index] - data[i][0];
    adist = a[index] - data[i][1];
    bdist = b[index] - data[i][2];
    dist = Ldist * Ldist + adist * adist + bdist * bdist;
    if (dist < min_dist) {
      min_dist = dist;
      min_index = i;
    }
  }

  hist[min_index] += weight;
  sum[min_index] += weight;
  com[min_index][0] += weight * L[index];
  com[min_index][1] += weight * a[index];
  com[min_index][2] += weight * b[index];

}


/* Distance function for the EMD */
float dist(feature_t *F1, feature_t *F2) 
{ 
  return((float)cost[*F1][*F2]);
}


/* Order statistics function, from Cormen, Leiserson, and Rivest, Ch. 8 & 10 */
float RandomizedSelect(float *A, int p, int r, int i)
{
  int k, q, m, n;
  float temp, x;

  /* Page 187 */
  if (p == r)
    return(A[p]);

  /* Page 162 */
  q = rand() % (r - p + 1) + p; 
  temp = A[p];
  A[p] = A[q];
  A[q] = temp;

  /* Page 154 */
  x = A[p];
  m = p;
  n = r;
  while (1) {
    while (A[n] > x)
      n--;
    while (A[m] < x)
      m++;
    if (m < n) {
      temp = A[m];
      A[m] = A[n];
      A[n] = temp;
      n--;
      m++;
    } else
      break;
  }

  /* Page 187 */
  k = n - p + 1;
  if (i <= k)
    return(RandomizedSelect(A,p,n,i));
  else
    return(RandomizedSelect(A,n+1,r,i-k));
}

/* This is a C implementation of the MATLAB file MakeQtrMask2.m.  r is the
 * radius of the circle, and n is the number of wedges.  The result is an
 * array of length ceil(r) * ceil(r) * n where the entries state how much
 * area of each wedge of radius r (1/n of the quarter-circle) is inside 
 * each pixel.  The pixels are in Quadrant I, listed in column-major order 
 * (just as we would get by calling a MATLAB routine).  This is admittedly
 * a mess, as there are at least 27 different cases for a circle and a radial
 * line and a square to interact.
 */
double *MakeQtrMask(double r, int nwedges)
{
  int i, j, x, y, R = ceil(r);
  int ULC, ULL, ULR, URC, URL, URH, LRC, LRL, LRH, LLC, LLL, LLH;
  int InCircle, NoLine, LowLine, HighLine, LowIntersect, HighIntersect;
  double *mask, mlow, mhigh, lowangle, highangle;
  double CA, BA, AA, LA, BXC, BXL, BXH, TXC, TXL, TXH, LYC, LYL, LYH;
  double RYC, RYL, RYH, AAC, BAC, AAN, BAN, XLC, YLC, XHC, YHC;

  mask = (double *)calloc(R * R * nwedges, sizeof(double));
  
  /* Iterate over the lower left hand corner of each pixel */
  for (x = 0; x <= R - 1; x++)
    for (y = R - 1; y >= 0; y--) {
      /* Start by computing the pixel's area in the circle */
      if (x * x + y * y >= r * r) /* Pixel entirely outside circle */
	continue;
      if ((x+1) * (x+1) + (y+1) * (y+1) <= r * r) { /* Pixel entirely inside */
	CA = 1.0;
	InCircle = 1;
	URC = 1;
      } else { /* Tricky part; circle intersects pixel */
	URC = 0;
	ULC = x * x + (y+1) * (y+1) <= r * r;
	LRC = (x+1) * (x+1) + y * y <= r * r;
	BXC = sqrt(r * r - y * y);
	TXC = sqrt(r * r - (y+1) * (y+1));
	if (!ULC && !LRC)
	  CA = CArea(BXC, x, y, r);
	else if (ULC && !LRC)
	  CA = CArea(BXC, TXC, y, r) + TXC - x;
	else if (!ULC && LRC)
	  CA = CArea(x + 1, x, y, r);
	else /* if (ULC && LRC) */
	  CA = CArea(x + 1, TXC, y, r) + TXC - x;
	InCircle = 0;  /* Therefore, it must be on the border */
      }

      /* Check through each wedge */
      for (i = 0; i < nwedges; i++) {
	/* Compute area above lower radial line of wedge */
	lowangle = i * PI / (2 * nwedges);
	mlow = tan(lowangle);
	TXL = (y+1)/mlow;
	BXL = y/mlow;
	if (TXL <= x)
	  AA = 0.0;
	else if (i == 0 || BXL >= x+1)
	  AA = 1.0;
	else {
	  LLL = BXL > x;
	  URL = TXL > x+1;
	  LYL = mlow * x;
	  RYL = mlow * (x + 1);
	  if (LLL && URL)
	    AA = 1 - 0.5 * (RYL-y) * (x+1-BXL);
	  else if (!LLL && URL)
	    AA = 0.5 * (2*(y+1)-LYL-RYL);
	  else if (LLL && !URL)
	    AA = 0.5 * (TXL+BXL-2*x);
	  else
	    AA = 0.5 * (y+1-LYL) * (TXL-x);
	}
	LowLine = AA < 1.0 && AA > 0.0;

	/* Compute area below upper radial line of wedge */
	/* The cases are reversed from the lower line cases */
	highangle = (i+1) * PI / (2 * nwedges);
	mhigh = tan(highangle);
	TXH = (y+1)/mhigh;
	BXH = y/mhigh;
	RYH = mhigh*(x+1);
	LYH = mhigh*x;
	if (i == nwedges-1 || TXH <= x)
	  BA = 1.0;
	else if (BXH >= x+1)
	  BA = 0.0;
	else {
	  LLH = BXH < x;
	  URH = TXH < x+1;
	  if (LLH && URH)
	    BA = 1 - 0.5 * (y+1-LYH) * (TXH-x);
	  else if (!LLH && URH)
	    BA = 1 - 0.5 * (BXH+TXH-2*x);
	  else if (LLH && !URH)
	    BA = 0.5 * (LYH+RYH-2*y);
	  else /* if (!LLH && !URH) */
	    BA = 0.5 * (RYH-y) * (x+1-BXH);
	}
	HighLine = BA < 1.0 && BA > 0.0;
	LA = BA + AA - 1.0;
	if (LA == 0.0) /* Pixel not in wedge */
	  continue;
	NoLine = LA == 1.0;

	/* Finish the cases we know about so far */
	if (InCircle) {
	  mask[i * R*R + x * R + R - 1 - y] = LA;
	  continue;
	} else if (NoLine) {
	  mask[i * R*R + x * R + R - 1 - y] = CA;
	  continue;
	}
	
	/* We can now assert (~InCircle && (HighLine || LowLine)) */
	/* But this does not ensure the circular arc intersects the line */
	LYC = sqrt(r * r - x * x);
	RYC = sqrt(r * r - (x+1) * (x+1));
	LowIntersect = LowLine &&
	  ((!ULC && !LRC && ((LLL && BXL < BXC) || (!LLL && LYL < LYC))) ||
	   (!ULC && LRC) || (ULC && !LRC) ||
	   (ULC && LRC && ((!URL && TXL >= TXC) || (URL && RYL >= RYC))));

	HighIntersect = HighLine &&
	  ((!ULC && !LRC && ((!LLH && BXH < BXC) || (LLH && LYH < LYC))) ||
	   (!ULC && LRC) || (ULC && !LRC) ||
	   (ULC && LRC && ((URH && TXH >= TXC) || (!URH && RYH >= RYC))));

	/* Recompute BA and AA (now BAC and AAC) given the intersection */
	if (LowIntersect) {
	  XLC = cos(lowangle) * r;
	  YLC = sin(lowangle) * r;
	  if (!LRC && LLL)
	    AAC = CA - 0.5 * (XLC - BXL) * (YLC - y) - CArea(BXC, XLC, y, r);
	  else if (!LRC && !LLL)
	    AAC = CA - 0.5 * (XLC - x) * (YLC + LYL - 2 * y) -
	      CArea(BXC, XLC, y, r);
	  else if (LRC && LLL)
	    AAC = CArea(XLC, x, y, r) - 0.5 * (YLC - y) * (XLC - BXL);
	  else /* if (LRC && !LLL) */
	    AAC = CA - CArea(x+1, XLC, y, r) - 
	      0.5 * (YLC + LYL - 2 * y) * (XLC - x);
	}
	  
	if (HighIntersect) {
	  XHC = cos(highangle) * r;
	  YHC = sin(highangle) * r;
	  if (!LRC && !LLH)
	    BAC = 0.5 * (XHC - BXH) * (YHC - y) + CArea(BXC, XHC, y, r);
	  else if (!LRC && LLH)
	    BAC = 0.5 * (XHC - x) * (YHC + LYH - 2 * y) + 
	      CArea(BXC, XHC, y, r);
	  else if (LRC && LLH)
	    BAC = CArea(x+1, XHC, y, r) + 
	      0.5 * (YHC + LYH - 2 * y) * (XHC - x);
	  else /* if (LRC && !LLH) */
	    BAC = CArea(x+1, XHC, y, r) + 0.5 * (YHC - y) * (XHC - BXH);
	}
	  
	/* Compute area for a few more cases */
	if (LowIntersect && !HighLine) {
	  mask[i * R*R + x * R + R - 1 - y] = AAC;
	  continue;
	} else if (HighIntersect && !LowLine) {
	  mask[i * R*R + x * R + R - 1 - y] = BAC;
	  continue;
	} else if (HighIntersect && LowIntersect) {  
	  mask[i * R*R + x * R + R - 1 - y] = AAC + BAC - CA;
	  continue;
	}
	
	/* Here we can assert (~InCircle && (HighLine || LowLine) &&
	 * !LowIntersect && !HighIntersect).  There are still many 
	 * possible answers.  Start by computing BAN and AAN (N for No
	 * Intersection)
	 */
	if (LowLine && !LowIntersect) {
	  if (!ULC && !LLL)
	    AAN = 0;
	  else if (!LRC && LLL)
	    AAN = CA;
	  else if (LRC && URL && LLL)
	    AAN = CA - 0.5 * (RYL - y) * (x+1 - BXL);
	  else if (ULC && URL && !LLL)
	    AAN = CA - 0.5 * (RYL + LYL - 2 * y);
	  else /* if (ULC && !URL) */
	    AAN = AA;
	}

	if (HighLine && !HighIntersect) {
	  if (!ULC && LLH)
	    BAN = CA;
	  else if (!LRC && !LLH)
	    BAN = 0;
	  else if (LRC && !URH && !LLH)
	    BAN = BA;
	  else if (ULC && !URH && LLH)
	    BAN = 0.5 * (RYL + LYL - 2 * y);
	  else if (ULC && URH)
	    BAN = CA + BA - 1;
	}

	if (LowLine && !LowIntersect && HighLine && !HighIntersect) {
	  mask[i * R*R + x * R + R - 1 - y] = AAN + BAN - CA;
	  continue;
	} else if (LowIntersect && HighLine && !HighIntersect) {
	  mask[i * R*R + x * R + R - 1 - y] = AAC + BAN - CA;
	  continue;
	} else if (LowLine && !LowIntersect && HighIntersect) {
	  mask[i * R*R + x * R + R - 1 - y] = AAN + BAC - CA;
	  continue;
	} else if (LowLine && !LowIntersect) {
	  mask[i * R*R + x * R + R - 1 - y] = AAN;
	  continue;
	} else if (HighLine && !HighIntersect) {
	  mask[i * R*R + x * R + R - 1 - y] = BAN;
	  continue;
	} else {
	  fprintf(stderr,"Big nasty horrible bug just happened\n");
	  mask[i * R*R + x * R + R - 1 - y] = 0.0;
	}
      }
    }

  if (PRINT_WEDGES)
    for (i = 0; i < nwedges; i++) {
      fprintf(stderr,"Wedge %d\n",i);
      for (y = 0; y < R; y++) {
	for (x = 0; x < R; x++)
	  fprintf(stderr,"%.4f ",mask[i * R*R + x * R + y]);
	fprintf(stderr,"\n");
      }
    }

  return(mask);      
}

	  
/**************************************************/
/***     Level 1 Functions                      ***/
/**************************************************/

/* This routine converts the image data from the previously unknown type
 * to the internal representation, going through RGB2Lab if necessary.
 * The image is also cropped according to the given dimensions.
 */
void ConvertImage(VALUE imgdata, int imgrows, int imgcols, enum imgtype type,
		 double *dimensions, int maxscale, int *rows, int *cols,
		 float **L, float **a, float **b)
{
  int sheetelts, npixels, index1, index2, i, j;

  *rows = dimensions[2] - dimensions[0] + 2 * maxscale;
  *cols = dimensions[3] - dimensions[1] + 2 * maxscale;
  sheetelts = *rows * *cols;
  npixels = imgrows * imgcols;
  *L = (ImgElt *)malloc(sheetelts * sizeof(ImgElt));
  *a = (ImgElt *)malloc(sheetelts * sizeof(ImgElt));
  *b = (ImgElt *)malloc(sheetelts * sizeof(ImgElt));
  if (!*L || !*a || !*b)
    rb_raise(rb_eStandardError, "Image could not be allocated (out of memory)");

  VALUE const * const imgdataptr = RARRAY_PTR(imgdata);
  if (type == LabImg) /* imgdata is double * */

    for (i = dimensions[0] - maxscale; i < dimensions[2] + maxscale; i++)
      for (j = dimensions[1] - maxscale; j < dimensions[3] + maxscale; j++) {
	index1 = (j-dimensions[1]+maxscale)*(*rows)+(i-dimensions[0]+maxscale);
	index2 = j * imgrows + i;
	(*L)[index1] = NUM2DBL(imgdataptr[index2]);
	(*a)[index1] = NUM2DBL(imgdataptr[index2+npixels]);
	(*b)[index1] = NUM2DBL(imgdataptr[index2+2*npixels]);
      }

  else /* imgdata is unsigned char * */

    for (i = dimensions[0] - maxscale; i < dimensions[2] + maxscale; i++)
      for (j = dimensions[1] - maxscale; j < dimensions[3] + maxscale; j++) {
	index1 = (j-dimensions[1]+maxscale)*(*rows)+(i-dimensions[0]+maxscale);
	index2 = j * imgrows + i;
	RGB2Lab((float)((unsigned char)NUM2INT(imgdataptr[index2])),
		(float)((unsigned char)NUM2INT(imgdataptr[index2+npixels])),
		(float)((unsigned char)NUM2INT(imgdataptr[index2+2*npixels])),
		&(*L)[index1], &(*a)[index1], &(*b)[index1]);
      }
}


/* This routine creates the mask that is one-quarter of the size of the circle
 * and nwedges deep.  It also creates a mask sum over all the wedges and
 * computes the number of pixels with non-zero weights.
 */
int CreateMask(double sigma, int wedges, int radius, double **mask, 
	       double **sum)
{
  double *gauss, r;
  int npoints = 0, i, j, k;

  *mask = MakeQtrMask(sigma * 3.0, wedges);
  if (!*mask)
    rb_raise(rb_eStandardError, "mask could not be allocated (out of memory)");

  /* Create a normalized Gaussian mask in Quadrant I */
  /* MATLAB matrices (and hence this mask) are in column order */
  /* ISOTROPIC */
  /* Note that there is no normalizing constant in front */
  gauss = (double *)malloc(radius*radius * sizeof(double));
  for (i = 0; i < radius; i++)
    for (j = 0; j < radius; j++) {
      r = sqrt((radius-i-0.5) * (radius-i-0.5) + (j+0.5) * (j+0.5));
      gauss[j*radius+i] = r * exp(-(r*r)/(2*sigma*sigma));  /* RAYLEIGH */
      /* gauss[j*radius+i] = exp(-(r*r)/(2*sigma*sigma)); GAUSSIAN */
      /* gauss[j*radius+i] = 1 - r / radius;  LINEAR */
    }

  if (PRINT_GAUSSIAN) {
    printf("Gaussian\n[");
    for (j = 0; j < radius; j++) {
      for (k = 0; k < radius; k++)
	printf("%.3f ",gauss[k*radius + j]);
      printf("\n");
    }
    printf("];\n");
  }

  *sum = (double *)malloc(radius*radius * sizeof(double));
  if (!*sum)
    rb_raise(rb_eStandardError, "sum could not be allocated (out of memory)");
  for (i = 0; i < radius * radius * wedges; i++) {
    (*mask)[i] *= gauss[i % (radius * radius)];
    (*sum)[i % (radius * radius)] += (*mask)[i];
  }

  if (PRINT_MASK) {
    printf("Mask\n");
    for (i = 0; i < wedges; i++) {
      printf("[");
      for (j = 0; j < radius; j++) {
	for (k = 0; k < radius; k++)
	  printf("%.3f ",(*mask)[i*radius*radius + k*radius + j]);
	printf("\n");
      }
      printf("];\n");
    }
  }

  if (PRINT_MASKSUM) {
    printf("Mask Sum\n[");
    for (j = 0; j < radius; j++) {
      for (k = 0; k < radius; k++)
	printf("%.3f ",(*sum)[k*radius + j]);
      printf("\n");
    }
    printf("];\n");
  }

  for (i = 0; i < radius * radius; i++)
    npoints += ((*sum)[i] > 0); /* Count up non-zero pixels */
  return(4 * npoints);
}


/* Note (1 Nov 99): some day soon, index should be returned for quicker 
 * building of the histograms.
 */
int ClusterPoints(float *L, float *a, float *b, double *masksum, 
		  int radius, int totalpoints, int rows, int r, int c, 
		  double totalweight, float maxweight, int maxclusters, 
		  float cluster[][DIM])
{
  float *temp, *temp2, *points, *weights;
  int i, j, nclusters, npoints, simple, *index = NULL; 
  float cutoff = 0.0, mprob, ratio, mweight;
  
  points = (float *)malloc(DIM * totalpoints * sizeof(float));
  if (!points)
    rb_raise(rb_eStandardError, "points could not be allocated (out of memory)");
  temp = points; /* temp will walk through the data array */

  weights = (float *)malloc(totalpoints * sizeof(float));
  if (!weights)
    rb_raise(rb_eStandardError, "weights could not be allocated (out of memory)");
  temp2 = weights; /* same here */

  if (PRINT_SAMPLING)
    printf("Total weight: %f (%.2f) Max: %f  Desired: %f\n", totalweight,
	  totalweight / totalpoints, masksum[radius-1],
	  ((float) totalpoints * SAMPLE_PERCENTAGE));

  /* To get the proper sampling percentage, we have two options.  The simple
   * one is to multiply all the mask weights (which we treat as probabilities)
   * by a constant.  However, if this constant is too high, the maximum weight
   * will go above 1.0.  Therefore, we must use a more complex option.
   */
  if (maxweight * totalpoints * SAMPLE_PERCENTAGE / totalweight > 1.0) {
    simple = 0;
    ratio = (1.0 - SAMPLE_PERCENTAGE) / (1.0 - totalweight / totalpoints);
  } else {
    simple = 1;
    ratio = (float)totalpoints * SAMPLE_PERCENTAGE / totalweight;
  }
  
  /* Determine if the window is large enough for sampling, and if so,
   * sample each pixel in the mask with a Bernoulli trial.  Otherwise,
   * use all the pixels for clustering.
   */
  npoints = 0;
  for (i = 0; i < radius; i++)
    for (j = 0; j < radius; j++)
      if (masksum[i*radius+j] > 0) {
	mweight = masksum[i*radius+j];
	PlacePoint(L,a,b,r-radius+j,c+i,rows,temp,temp2,mweight);
	temp += DIM;
	temp2++;
	npoints++;
	PlacePoint(L,a,b,r-i-1,c-radius+j,rows,temp,temp2,mweight);
	temp += DIM;
	temp2++;
	npoints++;
	PlacePoint(L,a,b,r+radius-j-1,c-i-1,rows,temp,temp2,mweight);
	temp += DIM;
	temp2++;
	npoints++;
	PlacePoint(L,a,b,r+i,c+radius-j-1,rows,temp,temp2,mweight);
	temp += DIM;
	temp2++;
	npoints++;
      }

  if (PRINT_SAMPLING)
    printf("Simple=%d Ratio: %.3f  Points: %d of %d  Pixel percentage: %.2f\n", simple, ratio, npoints, totalpoints, (double) npoints / (double) totalpoints);

  /* Cluster using binary split algorithm */
  /* Changed 18 May 99 to include indexing */
  bs(points, npoints, maxclusters, &nclusters, &temp, index);

  /* Copy output from temp to cluster */
  for (i = 0; i < nclusters; i++)
    for (j = 0; j < DIM; j++)
      cluster[i][j] = temp[i * DIM + j];

  if (PRINT_CLUSTERS)
    for (i = 0; i < nclusters; i++)
      printf("%3d: %6.2f %6.2f %6.2f\n",i,cluster[i][0],cluster[i][1],
		cluster[i][2]);

  free(points);
  free(weights);
  free(temp);
  free(index);      /* This must be passed back for use in CreateWedgeHist */

  return(nclusters);
}


void CreateWedgeHistograms(float *L, float *a, float *b, double *mask, int r, 
			   int c, int radius, int nwedges, int rows, 
			   float cluster[][DIM], int nclusters, double *hist)
{
  int i, j, k, index;
  float com[MAXCLUSTERS][DIM], sum[MAXCLUSTERS];

  memset(hist, 0, nwedges * 4 * MAXCLUSTERS * sizeof(double));
  memset(com, 0, MAXCLUSTERS * DIM * sizeof(float));
  memset(sum, 0, MAXCLUSTERS * sizeof(float));

  for (k = 0; k < nwedges; k++) {
    for (i = 0; i < radius; i++)
      for (j = 0; j < radius; j++) {
	index = k*radius*radius+i*radius+j;
	if (mask[index] > 0) {
	  AddWeight(L, a, b, mask[index], r-radius+j, c+i, rows, cluster,
		    hist+k*MAXCLUSTERS, com, sum, nclusters);
	  AddWeight(L, a, b, mask[index], r-i-1, c-radius+j, rows, cluster,
		    hist+(nwedges+k)*MAXCLUSTERS, com, sum, nclusters);
	  AddWeight(L, a, b, mask[index], r+radius-j-1, c-i-1, rows, cluster,
		    hist+(2*nwedges+k)*MAXCLUSTERS, com, sum, nclusters);
	  AddWeight(L, a, b, mask[index], r+i, c+radius-j-1, rows, cluster,
		    hist+(3*nwedges+k)*MAXCLUSTERS, com, sum, nclusters);
	}
      }
  }

}


void ComputeCostMatrix(float cluster[][DIM], int nclusters)
{
  int i, j;
  float dist, temp_dist;

  /* cost has to be a global variable because of the way the EMD code is */
  memset(cost, 0, MAXCLUSTERS*MAXCLUSTERS*sizeof(double));

  for (i = 0; i < nclusters; i++) {
    for (j = 0; j < nclusters; j++) {
      if (i > j)
	cost[i][j] = cost[j][i];
      else if (i == j || cluster[i][0] == BLACK && cluster[j][0] == BLACK)
	cost[i][j] = 0.0;
      else {
	dist = 0;
	temp_dist = (cluster[i][0] - cluster[j][0]);
	dist += temp_dist * temp_dist; 
	temp_dist = (cluster[i][1] - cluster[j][1]);
	dist += temp_dist * temp_dist;
	temp_dist = (cluster[i][2] - cluster[j][2]);
	dist += temp_dist * temp_dist;

	/* dist = sqrt(dist) / PERCEPTUAL_THRESH; 
	   cost[i][j] = (dist >= 1) ? 1 : dist;*/

        cost[i][j] = 1 - exp(-sqrt(dist) / GAMMA);
        /* cost[i][j] = sqrt(dist) / 100; */
      }
      if (PRINT_COST)
	printf("%.2f ",cost[i][j]);
      }
    if (PRINT_COST)
      printf("\n");
  }
  if (PRINT_COST)
    printf("\n");
}


void ComputeOutputParameters(float *work, int nwedges, int nori, int index, 
			     int pagesize, int plot, double *str, double *ab, 
			     double *unc, double *ori)
{
  int i, strindex = 0, abindex = 0, l_uncert, r_uncert, responses;
  int startindex = 0;
  double angle, maxEMD = 0.0, minEMD = 1.0, unc_thresh, maxangle, tmax;
  double a, b, c, m, b1, b2, strength, orientation, dangle, unc_length;
  double wedgesize, maxEMDori;

  wedgesize = 90.0 / nwedges;
  maxangle = (double)(nori * wedgesize); /* Always 180 or 360 */

  /* Compute Minimum and Maximum EMD values */
  for (i = 0; i < nori; i++) {
    if (work[i] > maxEMD) {
      maxEMD = work[i];
      strindex = i;
      maxEMDori = strindex * wedgesize;
    }
    if (work[i] < minEMD) {
      minEMD = work[i];
      abindex = i;
    }
  }

  /* The strength and orientation of an edge (or a corner) lie not at the
   * maximum EMD value but rather at the vertex of the parabola that runs
   * through the maximum and the two points on either side.  The first
   * computation of orientation assumes the maximum is the y-intercept.  
   * After computing the strength, we adjust the orientation.  
   */
  a = work[strindex];
  b = work[(strindex+nori-1) % nori];
  c = work[(strindex+nori+1) % nori];
  if (b + c != 2 * a) {
    orientation = (wedgesize / 2) * ((b - c) / (b + c - 2 * a));
    strength = a + ((c - b) / (2 * wedgesize)) * orientation +
      ((b + c - 2 * a) / (2 * wedgesize * wedgesize)) * orientation * 
      orientation;
    orientation = fmod(maxEMDori + orientation + 360, maxangle);
  } else { /* Uncertainty abounds */
    strength = a;
    orientation = maxEMDori;
  }

  /* Assuming no compass plots, there is only one value for str and ab */
  if (!plot) /* User did not ask for compass plots, just the max */
    str[index] = strength;
  if (ab)
    ab[index] = minEMD;

  /* Compute Uncertainties and Orientations */
  unc_thresh = 0.95 * maxEMD;

  /* First, find a value below the threshold, or declare total uncertainty */
  startindex = 0;
  while (startindex < nori && work[startindex] > unc_thresh)
    startindex++;
  if (startindex == nori) { /* Total uncertainty */
    if (unc)
      unc[index] = maxangle;
    if (ori)
      for (responses = 0; responses < MAXRESPONSES; responses++)
	ori[responses * pagesize + index] = -1;
    return;
  }

  /* Now scan until we find values above the threshold and record them */
  /* Do this as many times as needed */
  responses = 0;
  i = startindex;
  while (responses < MAXRESPONSES) {
    while (work[(i + nori) % nori] <= unc_thresh && i != startindex + nori)
      i++;
    l_uncert = i;
    if (i == startindex + nori) /* Finished all orientations */
      break;
    while (work[(i + nori) % nori] > unc_thresh && i != startindex + nori)
      i++;
    r_uncert = i - 1;

    /* Uncertainty exists only if three or more consecutive values are all
     * above unc_thresh (meaning r_uncert - l_uncert is 2 or greater).  
     * Otherwise, the strength and orientation estimation routine above can 
     * handle it.
     */

    if (unc) {
      if (r_uncert - l_uncert > 1) {
	unc[responses * pagesize + index] = (r_uncert - l_uncert) * wedgesize;
	if (unc[responses * pagesize + index] > maxangle)
	  unc[responses * pagesize + index] = maxangle;
      } else
	unc[responses * pagesize + index] = 0;
    }
  
    /* Compute Orientation 
     * Our orientation estimate can only be used when we are near the  
     * maximum; other high responses must use a more roughly quantized value.
     */
    if (ori) {
      angle = (r_uncert + l_uncert) * wedgesize / 2;
      unc_length = r_uncert - l_uncert;
      while (angle < 0)
	angle += maxangle;
      while (angle >= maxangle)
	angle -= maxangle;
      dangle = MIN(abs(maxEMDori - angle), 
		   abs(maxEMDori - (angle + maxangle)));
      if (unc_length < 2 && dangle <= unc_length * wedgesize / 2 ) { 
	/* Use orientation estimate */
	ori[responses * pagesize + index] = orientation;
      } else
	ori[responses * pagesize + index] = angle;
    }

    responses++;
  }

  if (ori) /* Fill the remaining values with -1's */
    while (responses < MAXRESPONSES) {
      ori[responses * pagesize + index] = -1;
      responses++;
    }
}


void Compass(VALUE imgdata, int imgrows, int imgcols, enum imgtype type,
	     double *sigmas, int numsigmas, int maxradius, double *spacing, 
	     double *dimensions, double *angles, int numangles, int nwedges, 
	     double *maxclusters, int plot, double **strength, 
	     double **abnormality, double **orientation, 
	     double **uncertainty)
{
  int subimgrows, subimgcols, npoints, outputrows, outputcols, cellindex, mc;
  int h, i, j, k, r, c, index, s, space, nclusters, nori, anglewedges, masksz;
  double *mask, *masksum, *hist, *str, *ab = NULL, *unc = NULL, *ori = NULL;
  double radius, sum, sigma;
  float cluster[MAXCLUSTERS][DIM], hist1[MAXCLUSTERS], hist2[MAXCLUSTERS];
  float hist1norm[MAXCLUSTERS], hist2norm[MAXCLUSTERS];
  float *L, *a, *b, *work, maxwt = 0.0, inmass, outmass;
  feature_t f[MAXCLUSTERS];
  signature_t s1, s2;
  float temp;

  /* Set up all parameters and memory that is initialized once only */
  ConvertImage(imgdata, imgrows, imgcols, type, dimensions, maxradius, 
	       &subimgrows, &subimgcols, &L, &a, &b);
  /* 4 * nwedges is the maximum number of orientations */
  work = (float *)malloc(4 * nwedges * sizeof(float));
  if (!work)
    rb_raise(rb_eStandardError, "work could not be allocated (out of memory)");
  hist = (double *)malloc(nwedges * 4 * MAXCLUSTERS * sizeof(double));
  if (!hist)
    rb_raise(rb_eStandardError, "hist could not be allocated (out of memory)");
  for (i = 0; i < MAXCLUSTERS; i++)
    f[i] = i;
  s1.Features = f;
  s2.Features = f;
  s1.Weights = hist1norm;
  s2.Weights = hist2norm;
  srand(clock());

  /* Loop over every sigma value */
  for (s = 0; s < numsigmas; s++) {

    /* Retrieve data for this sigma */
    /* Rewritten 31 July 1999 to allow floating-point values for sigma */
    sigma = sigmas[s];
    masksz = (int)ceil(3 * sigma);  
    space = (int)spacing[s];
    mc = (int)maxclusters[s];

    npoints = CreateMask(sigma, nwedges, masksz, &mask, &masksum);

    /* Compute max weight of entire circle and sum of weights over one wedge 
     * of the circle.
     */
    sum = 0.0;
    for (i = 0; i < masksz * masksz; i++) { /* sums over 1/4 of circle */
      sum += masksum[i];
      if (masksum[i] > maxwt)
	maxwt = masksum[i];
    }
    sum /= nwedges;

    /* Loop over every desired image position */
    for (r = maxradius; r <= subimgrows - maxradius; r += space)
      for (c = maxradius; c <= subimgcols - maxradius; c += space) {
	if (PRINT_LOC)
	  printf("Sigma %.2f, Row %3d, Column %3d\n", sigma,
		    (int)dimensions[0]+r-maxradius, 
		    (int)dimensions[1]+c-maxradius);

	nclusters = ClusterPoints(L, a, b, masksum, masksz, npoints, 
				  subimgrows, r, c, sum * nwedges * 4, maxwt, 
				  mc, cluster);

	CreateWedgeHistograms(L, a, b, mask, r, c, masksz, nwedges, 
			      subimgrows, cluster, nclusters, hist);

	ComputeCostMatrix(cluster, nclusters);

	/* Loop over each desired angle */
	for (h = 0; h < numangles; h++) {

	  /* Compute desired data and cell pointers for each angle */
	  anglewedges = angles[h] * nwedges / 90.0;
	  inmass = sum * ((anglewedges <= 2) ? anglewedges * WEDGEWT :
			  WEDGEWT * 2 + anglewedges - 2);
	  outmass = sum * (WEDGEWT * 2 + 4 * nwedges - anglewedges - 2); 
	  nori = (angles[h] == 180) ? 2 * nwedges : 4 * nwedges;
	  cellindex = s * numangles + h;
	  str = strength[cellindex];
	  outputrows = subimgrows;
	  if (plot)
	    outputcols = subimgcols / nori;
	  else
	    outputcols = subimgcols;
	  if (abnormality)
	    ab = abnormality[cellindex];
	  if (orientation)
	    ori = orientation[cellindex];
	  if (uncertainty)
	    unc = uncertainty[cellindex];

	  /* Compute initial histogram sums */
	  for (i = 0; i < nclusters; i++) {
	    hist1[i] = 0.0;
	    for (j = 0; j < anglewedges; j++)
	      if (j == 0 || j == anglewedges - 1)
		hist1[i] += hist[j * MAXCLUSTERS + i] * WEDGEWT;
	      else
		hist1[i] += hist[j * MAXCLUSTERS + i];

	    hist2[i] = 0.0;
	    for (j = anglewedges; j < 4 * nwedges; j++)
	      if (j == anglewedges || j == 4 * nwedges - 1)
		hist2[i] += hist[j * MAXCLUSTERS + i] * WEDGEWT;
	      else
		hist2[i] += hist[j * MAXCLUSTERS + i];
	  }

	  /* Loop over every orientation */
	  for (i = 0; i < nori; i++) {

	    /* Normalize the histograms */
	    for (j = 0; j < nclusters; j++) {

	      hist1norm[j] = hist1[j] / inmass;
	      if (PARTIAL_MATCH) /* Normalize both by same amount */
		hist2norm[j] = hist2[j] / inmass;
	      else  /* Normalize both by the number of wedges */
		hist2norm[j] = hist2[j] / outmass;
	    }

	    if (PRINT_HIST) {
	      for (j = 0; j < nclusters; j++)
		printf("%.2f ",hist1norm[j]);
	      printf("\n");
	      for (j = 0; j < nclusters; j++)
		printf("%.2f ",hist2norm[j]);
	      printf("\n\n");
	    }
	  
	    /* Compute EMD */
	    s1.n = nclusters;
	    s2.n = nclusters;
	    work[i] = emd(&s1, &s2, dist, NULL, NULL);
	    if (work[i] > 1)
	      work[i] = 1.0;
	    else if (work[i] < 0)
	      work[i] = 0.0;
	    if (plot)
	      str[i * outputrows * outputcols + (c-maxradius) / space *
		 outputrows + (r - maxradius) / space] = work[i];

	    if (DEBUG)
	      printf("%.4f ",work[i]);
	  
	    /* Update the histograms except for the last iteration */
	    if (i < nori - 1)
	      for (j = 0; j < MAXCLUSTERS; j++) {
		hist1[j] += -hist[i*MAXCLUSTERS+j] * WEDGEWT +
		  -hist[((i+1)%(4*nwedges))*MAXCLUSTERS+j] * (1 - WEDGEWT) +
		  hist[((i+anglewedges)%(4*nwedges))*MAXCLUSTERS+j] * WEDGEWT +
		  hist[((i+anglewedges-1)%(4*nwedges))*MAXCLUSTERS+j] * 
		      (1 - WEDGEWT);
		hist2[j] += -hist[((i+anglewedges)%(4*nwedges))*MAXCLUSTERS+j]
		  * WEDGEWT +
		  -hist[((i+anglewedges+1)%(4*nwedges))*MAXCLUSTERS+j] *
		  (1 - WEDGEWT) +
		  hist[i*MAXCLUSTERS+j] * WEDGEWT +
		  hist[((i+4*nwedges-1)%(4*nwedges))*MAXCLUSTERS+j] * 
		  (1 - WEDGEWT);
	      }
	  }
	  if (DEBUG)
	    printf("\n");


	  ComputeOutputParameters(work, nwedges, nori, 
              (c - maxradius) / space * outputrows + (r - maxradius) / space, 
				  outputrows * outputcols, plot,
				  str, ab, unc, ori);
    
 
	} /* Next angle */

      } /* Next image position */

    free(mask);
    free(masksum);
  } /* Next sigma */

  free(work);
  free(hist);
  free(L);
  free(a);
  free(b);
}
