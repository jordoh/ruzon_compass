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
 * 2 August 1999 -- Modified for use with greyscale images.  No VQ is      *
 *   performed, and the EMD is implemented directly.                       *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "greycompass.h"

#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ROUND(x) floor((x) + 0.5)
#define PI 3.14159265358979323846
#define MAXRADIUS 64
#define MAXBRIGHTNESS 255.0
#define MAXCLUSTERS 51
#define PERCEPTUAL_THRESH 11.0
#define GAMMA 14.0
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

typedef struct bin_struct {
  float wsum;
  float weight;
  float value;
} Bin;

/* This global defines the distance between two intensities as a function of
 * their Euclidean distance.  Since intensity is measured under Lab as a 
 * number from 0-100, the distance array is in this range, measured in steps
 * of 0.5 units.  The function is 1 - exp(-x/GAMMA), where GAMMA = 14.0. */
/*float cost[201] = {
         0,   0.0351,   0.0689,   0.1016,   0.1331,   0.1635,   0.1929,
    0.2212,   0.2485,   0.2749,   0.3003,   0.3249,   0.3486,   0.3714,
    0.3935,   0.4147,   0.4353,   0.4551,   0.4742,   0.4927,   0.5105,
    0.5276,   0.5442,   0.5602,   0.5756,   0.5905,   0.6049,   0.6187,
    0.6321,   0.6450,   0.6575,   0.6695,   0.6811,   0.6923,   0.7031,
    0.7135,   0.7235,   0.7332,   0.7426,   0.7516,   0.7603,   0.7688,
    0.7769,   0.7847,   0.7923,   0.7995,   0.8066,   0.8134,   0.8199,
    0.8262,   0.8323,   0.8382,   0.8439,   0.8494,   0.8546,   0.8597,
    0.8647,   0.8694,   0.8740,   0.8784,   0.8827,   0.8868,   0.8908,
    0.8946,   0.8983,   0.9019,   0.9053,   0.9086,   0.9118,   0.9149,
    0.9179,   0.9208,   0.9236,   0.9263,   0.9288,   0.9313,   0.9337,
    0.9361,   0.9383,   0.9405,   0.9426,   0.9446,   0.9465,   0.9484,
    0.9502,   0.9520,   0.9536,   0.9553,   0.9568,   0.9584,   0.9598,
    0.9612,   0.9626,   0.9639,   0.9652,   0.9664,   0.9676,   0.9687,
    0.9698,   0.9709,   0.9719,   0.9729,   0.9738,   0.9747,   0.9756,
    0.9765,   0.9773,   0.9781,   0.9789,   0.9796,   0.9803,   0.9810,
    0.9817,   0.9823,   0.9829,   0.9835,   0.9841,   0.9847,   0.9852,
    0.9857,   0.9862,   0.9867,   0.9872,   0.9876,   0.9881,   0.9885,
    0.9889,   0.9893,   0.9897,   0.9900,   0.9904,   0.9907,   0.9910,
    0.9913,   0.9917,   0.9919,   0.9922,   0.9925,   0.9928,   0.9930,
    0.9933,   0.9935,   0.9937,   0.9939,   0.9942,   0.9944,   0.9946,
    0.9948,   0.9949,   0.9951,   0.9953,   0.9955,   0.9956,   0.9958,
    0.9959,   0.9961,   0.9962,   0.9963,   0.9965,   0.9966,   0.9967,
    0.9968,   0.9969,   0.9970,   0.9971,   0.9972,   0.9973,   0.9974,
    0.9975,   0.9976,   0.9977,   0.9978,   0.9979,   0.9979,   0.9980,
    0.9981,   0.9981,   0.9982,   0.9983,   0.9983,   0.9984,   0.9984,
    0.9985,   0.9985,   0.9986,   0.9986,   0.9987,   0.9987,   0.9988,
    0.9988,   0.9989,   0.9989,   0.9989,   0.9990,   0.9990,   0.9991,
    0.9991,   0.9991,   0.9992,   0.9992,   0.9992}; */

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
 * their weights into the histograms
 */
void AddWeight(float *intensity, double weight, int row, int col, int nrows, 
	       Bin *hist)
{
  int index = col * nrows + row, min_index, i;

  min_index = floor(intensity[index] * (MAXCLUSTERS - 1) / MAXBRIGHTNESS);
  hist[min_index].weight += weight;
  hist[min_index].wsum += weight * intensity[index];
  /* hist[min_index].value is not computed */
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
void ConvertImage(void *imgdata, int imgrows, int imgcols, enum imgtype type,
		 double *dimensions, int maxscale, int *rows, int *cols,
		 float **intensity)
{
  int sheetelts, index1, index2, i, j;
  float dummy1, dummy2;

  *rows = dimensions[2] - dimensions[0] + 2 * maxscale;
  *cols = dimensions[3] - dimensions[1] + 2 * maxscale;
  sheetelts = *rows * *cols;
  *intensity = (ImgElt *)calloc(sheetelts, sizeof(ImgElt));
  if (!*intensity)
    Error("Image could not be allocated (out of memory)");

  if (type == LabImg) /* imgdata is double * */

    for (i = dimensions[0] - maxscale; i < dimensions[2] + maxscale; i++)
      for (j = dimensions[1] - maxscale; j < dimensions[3] + maxscale; j++) {
	index1 = (j-dimensions[1]+maxscale)*(*rows)+(i-dimensions[0]+maxscale);
	index2 = j * imgrows + i;
	(*intensity)[index1] = ((double *)imgdata)[index2];
      }

  else /* imgdata is unsigned char * */

    for (i = dimensions[0] - maxscale; i < dimensions[2] + maxscale; i++)
      for (j = dimensions[1] - maxscale; j < dimensions[3] + maxscale; j++) {
	index1 = (j-dimensions[1]+maxscale)*(*rows)+(i-dimensions[0]+maxscale);
	index2 = j * imgrows + i;
	RGB2Lab((float)((unsigned char *)imgdata)[index2],
		(float)((unsigned char *)imgdata)[index2],
		(float)((unsigned char *)imgdata)[index2],
		&(*intensity)[index1], &dummy1, &dummy2);
	/*(*intensity)[index1] = (float)((unsigned char *)imgdata)[index2];*/
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
    Error("mask could not be allocated (out of memory)");

  /* Create a normalized Gaussian mask in Quadrant I */
  /* MATLAB matrices (and hence this mask) are in column order */
  /* ISOTROPIC */
  /* Note that there is no normalizing constant in front */
  gauss = (double *)calloc(radius*radius,sizeof(double));
  for (i = 0; i < radius; i++)
    for (j = 0; j < radius; j++) {
      r = sqrt((radius-i-0.5) * (radius-i-0.5) + (j+0.5) * (j+0.5));
      gauss[j*radius+i] = r * exp(-(r*r)/(2*sigma*sigma));  /* RAYLEIGH */
      /* gauss[j*radius+i] = exp(-(r*r)/(2*sigma*sigma)); GAUSSIAN */
      /* gauss[j*radius+i] = 1 - r / radius;  LINEAR */
    }

  if (PRINT_GAUSSIAN) {
    fprintf(stderr,"Gaussian\n[");
    for (j = 0; j < radius; j++) {
      for (k = 0; k < radius; k++)
	fprintf(stderr,"%.3f ",gauss[k*radius + j]);
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"];\n");
  }

  *sum = (double *)calloc(radius*radius,sizeof(double));
  if (!*sum)
    Error("sum could not be allocated (out of memory)");
  for (i = 0; i < radius * radius * wedges; i++) {
    (*mask)[i] *= gauss[i % (radius * radius)];
    (*sum)[i % (radius * radius)] += (*mask)[i];
  }

  if (PRINT_MASK) {
    fprintf(stderr,"Mask\n");
    for (i = 0; i < wedges; i++) {
      fprintf(stderr,"[");
      for (j = 0; j < radius; j++) {
	for (k = 0; k < radius; k++)
	  fprintf(stderr,"%.3f ",(*mask)[i*radius*radius + k*radius + j]);
	fprintf(stderr,"\n");
      }
      fprintf(stderr,"];\n");
    }
  }

  if (PRINT_MASKSUM) {
    fprintf(stderr,"Mask Sum\n[");
    for (j = 0; j < radius; j++) {
      for (k = 0; k < radius; k++)
	fprintf(stderr,"%.3f ",(*sum)[k*radius + j]);
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"];\n");
  }

  for (i = 0; i < radius * radius; i++)
    npoints += ((*sum)[i] > 0); /* Count up non-zero pixels */
  return(4 * npoints);
}


void CreateWedgeHistograms(float *intensity, double *mask, int r, 
			   int c, int radius, int nwedges, int rows, 
			   Bin *hist)
{
  int i, j, k, index;

  memset(hist, 0, nwedges * 4 * MAXCLUSTERS * sizeof(Bin));

  for (k = 0; k < nwedges; k++) {
    for (i = 0; i < radius; i++)
      for (j = 0; j < radius; j++) {
	index = k*radius*radius+i*radius+j;
	if (mask[index] > 0) {
	  AddWeight(intensity, mask[index], r-radius+j, c+i, rows,
		    hist+k*MAXCLUSTERS);
	  AddWeight(intensity, mask[index], r-i-1, c-radius+j, rows, 
		    hist+(nwedges+k)*MAXCLUSTERS);
	  AddWeight(intensity, mask[index], r+radius-j-1, c-i-1, rows, 
		    hist+(2*nwedges+k)*MAXCLUSTERS);
	  AddWeight(intensity, mask[index], r+i, c+radius-j-1, rows, 
		    hist+(3*nwedges+k)*MAXCLUSTERS);
	}
      }
  }
}


/* This implementation of the Earth Mover's Distance (EMD) is for one-
 * dimensional datasets where the distance between clusters can be easily
 * computed on the fly.
 */
float GreyEMD(Bin *dirt, Bin *hole)
{
  float work = 0.0, leftoverdirt = 0.0, leftoverhole = 0.0;
  float dirtamt, holeamt, massmoved;
  int i = -1, j = -1;

  /* We exit from inner loops, so this one must run forever */
  while (1) {

    /* Compute the amount of mass in the lowest numbered bin that hasn't 
     * been moved yet from the piles of dirt */
    if (leftoverdirt == 0.0) {
      i++;
      while(dirt[i].weight == 0.0 && i < MAXCLUSTERS)
	i++;
      if (i == MAXCLUSTERS)
	return(work);
      else
	dirtamt = dirt[i].weight;
    } else
      dirtamt = leftoverdirt;

    /* Do the same for the holes */
    if (leftoverhole == 0.0) {
      j++;
      while(hole[j].weight == 0.0 && j < MAXCLUSTERS)
	j++;
      if (j == MAXCLUSTERS)
	return(work);
      else
	holeamt = hole[j].weight;
    } else
      holeamt = leftoverhole;

    /* Compute the work done moving the smaller amount of mass and decide
     * how much is left over in each bin. */
    massmoved = MIN(dirtamt, holeamt);
    /*work += massmoved * cost[(int)ROUND(2 * abs(hole[j].value - 
      dirt[i].value))];*/
    work += massmoved * (1 - exp(-abs(hole[j].value - dirt[i].value) 
      / GAMMA));
    leftoverdirt = dirtamt - massmoved;
    leftoverhole = holeamt - massmoved;
  }
}

void ComputeOutputParameters(float *work, int nwedges, int nori, int index, 
			     int pagesize, int sheets, double *str, 
			     double *ab, double *unc, double *ori)
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
  if (sheets == 1) /* User did not ask for compass plots */
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


void Compass(void *imgdata, int imgrows, int imgcols, enum imgtype type,
	     double *sigmas, int numsigmas, int maxradius, double *spacing, 
	     double *dimensions, double *angles, int numangles, int nwedges, 
	     Matrix *strength, Matrix *abnormality, Matrix *orientation, 
	     Matrix *uncertainty)
{
  int subimgrows, subimgcols, npoints, outputrows, outputcols, cellindex;
  int h, i, j, k, r, c, index, s, space, nori, anglewedges, masksz;
  double *mask, *masksum, *str, *ab = NULL, *unc = NULL, *ori = NULL;
  double radius, sum, sigma;
  Bin *hist, hist1[MAXCLUSTERS], hist2[MAXCLUSTERS];
  Bin hist1norm[MAXCLUSTERS], hist2norm[MAXCLUSTERS];
  float *L, *work, maxwt = 0.0, inmass, outmass;
  float temp;

  /* Set up all parameters and memory that is initialized once only */
  ConvertImage(imgdata, imgrows, imgcols, type, dimensions, maxradius, 
	       &subimgrows, &subimgcols, &L);
  /* 4 * nwedges is the maximum number of orientations */
  work = (float *)calloc(4 * nwedges, sizeof(float));
  if (!work)
    Error("work could not be allocated (out of memory)");
  hist = (Bin *)calloc(nwedges * 4 * MAXCLUSTERS,sizeof(Bin));
  if (!hist)
    Error("hist could not be allocated (out of memory)");
  srand(clock());

  /* Loop over every sigma value */
  for (s = 0; s < numsigmas; s++) {

    /* Retrieve data for this sigma */
    /* Rewritten 31 July 1999 to allow floating-point values for sigma */
    sigma = sigmas[s];
    masksz = (int)ceil(3 * sigma);  
    space = (int)spacing[s];

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
	  fprintf(stderr,"Sigma %.2f, Row %3d, Column %3d\n", sigma,
		    (int)dimensions[0]+r-maxradius, 
		    (int)dimensions[1]+c-maxradius);

	CreateWedgeHistograms(L, mask, r, c, masksz, nwedges, subimgrows, 
			      hist);

	/* Loop over each desired angle */
	for (h = 0; h < numangles; h++) {

	  /* Compute desired data and cell pointers for each angle */
	  anglewedges = angles[h] * nwedges / 90.0;
	  inmass = sum * ((anglewedges <= 2) ? anglewedges * WEDGEWT :
			  WEDGEWT * 2 + anglewedges - 2);
	  outmass = sum * (WEDGEWT * 2 + 4 * nwedges - anglewedges - 2); 
	  nori = (angles[h] == 180) ? 2 * nwedges : 4 * nwedges;
	  cellindex = s * numangles + h;
	  str = strength[cellindex].ptr;
	  outputrows = strength[cellindex].rows;
	  outputcols = strength[cellindex].cols;
	  if (abnormality)
	    ab = abnormality[cellindex].ptr;
	  if (orientation)
	    ori = orientation[cellindex].ptr;
	  if (uncertainty)
	    unc = uncertainty[cellindex].ptr;

	  /* Compute initial histogram sums */
	  for (i = 0; i < MAXCLUSTERS; i++) {
	    hist1[i].weight = 0.0;
	    hist1[i].wsum = 0.0;
	    for (j = 0; j < anglewedges; j++)
	      if (j == 0 || j == anglewedges - 1) {
		hist1[i].weight += hist[j * MAXCLUSTERS + i].weight * WEDGEWT;
		hist1[i].wsum += hist[j * MAXCLUSTERS + i].wsum * WEDGEWT;
	      } else {
		hist1[i].weight += hist[j * MAXCLUSTERS + i].weight;
		hist1[i].wsum += hist[j * MAXCLUSTERS + i].wsum;
	      }

	    hist2[i].weight = 0.0;
	    hist2[i].wsum = 0.0;
	    for (j = anglewedges; j < 4 * nwedges; j++)
	      if (j == anglewedges || j == 4 * nwedges - 1) {
		hist2[i].weight += hist[j * MAXCLUSTERS + i].weight * WEDGEWT;
		hist2[i].wsum += hist[j * MAXCLUSTERS + i].wsum * WEDGEWT;
	      } else {
		hist2[i].weight += hist[j * MAXCLUSTERS + i].weight;
		hist2[i].wsum += hist[j * MAXCLUSTERS + i].wsum;
	      }
	  }

	  /* Loop over every orientation */
	  for (i = 0; i < nori; i++) {

	    /* Normalize the histograms */
	    for (j = 0; j < MAXCLUSTERS; j++) {
	      hist1norm[j].value = hist1[j].wsum / hist1[j].weight;
	      hist1norm[j].weight = hist1[j].weight / inmass;
	      hist2norm[j].value = hist2[j].wsum / hist2[j].weight;

	      if (PARTIAL_MATCH) /* Normalize both by same amount */
		hist2norm[j].weight = hist2[j].weight / inmass;
	      else /* Normalize both by the number of wedges */
		hist2norm[j].weight = hist2[j].weight / outmass;
	    }

	    if (PRINT_HIST) {
	      for (j = 0; j < MAXCLUSTERS; j++)
		fprintf(stderr,"%.2f ",hist1norm[j].weight);
	      fprintf(stderr,"\n");
	      for (j = 0; j < MAXCLUSTERS; j++)
		fprintf(stderr,"%.2f ",hist2norm[j].weight);
	      fprintf(stderr,"\n\n");
	    }
	  
	    /* Compute EMD */
	    work[i] = GreyEMD(hist1norm, hist2norm);
	    if (work[i] > 1)
	      work[i] = 1.0;
	    else if (work[i] < 0)
	      work[i] = 0.0;
	    if (strength[cellindex].sheets > 1)
	      str[i * outputrows * outputcols + (c-maxradius) / space *
		 outputrows + (r - maxradius) / space] = work[i];

	    if (DEBUG)
	      fprintf(stderr,"%.4f ",work[i]);
	  
	    /* Update the histograms except for the last iteration */
	    if (i < nori - 1)
	      for (j = 0; j < MAXCLUSTERS; j++) {
		hist1[j].weight += -hist[i*MAXCLUSTERS+j].weight * WEDGEWT +
		  -hist[((i+1)%(4*nwedges))*MAXCLUSTERS+j].weight * 
		  (1 - WEDGEWT) +
		  hist[((i+anglewedges)%(4*nwedges))*MAXCLUSTERS+j].weight * 
		  WEDGEWT +
		  hist[((i+anglewedges-1)%(4*nwedges))*MAXCLUSTERS+j].weight * 
		  (1 - WEDGEWT);
		hist1[j].wsum += -hist[i*MAXCLUSTERS+j].wsum * WEDGEWT +
		  -hist[((i+1)%(4*nwedges))*MAXCLUSTERS+j].wsum * 
		  (1 - WEDGEWT) +
		  hist[((i+anglewedges)%(4*nwedges))*MAXCLUSTERS+j].wsum * 
		  WEDGEWT +
		  hist[((i+anglewedges-1)%(4*nwedges))*MAXCLUSTERS+j].wsum * 
		  (1 - WEDGEWT);
		hist2[j].weight += 
		  -hist[((i+anglewedges)%(4*nwedges))*MAXCLUSTERS+j].weight
		  * WEDGEWT +
		  -hist[((i+anglewedges+1)%(4*nwedges))*MAXCLUSTERS+j].weight *
		  (1 - WEDGEWT) +
		  hist[i*MAXCLUSTERS+j].weight * WEDGEWT +
		  hist[((i+4*nwedges-1)%(4*nwedges))*MAXCLUSTERS+j].weight * 
		  (1 - WEDGEWT);
		hist2[j].wsum += 
		  -hist[((i+anglewedges)%(4*nwedges))*MAXCLUSTERS+j].wsum
		  * WEDGEWT +
		  -hist[((i+anglewedges+1)%(4*nwedges))*MAXCLUSTERS+j].wsum *
		  (1 - WEDGEWT) +
		  hist[i*MAXCLUSTERS+j].wsum * WEDGEWT +
		  hist[((i+4*nwedges-1)%(4*nwedges))*MAXCLUSTERS+j].wsum * 
		  (1 - WEDGEWT);
	      }
	  }
	  if (DEBUG)
	    fprintf(stderr,"\n");

	  ComputeOutputParameters(work, nwedges, nori, 
              (c - maxradius) / space * outputrows + (r - maxradius) / space, 
				  outputrows * outputcols, 
				  strength[cellindex].sheets,
				  str, ab, unc, ori);

 
	} /* Next angle */

      } /* Next image position */

    free(mask);
    free(masksum);
  } /* Next sigma */

  free(work);
  free(hist);
  free(L);
}
