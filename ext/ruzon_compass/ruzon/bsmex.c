#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mex.h>
#include "matrix.h"
#include "bs.h"

#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define DEBUG 0
#define ROWS 5

enum imgtype {RGBImg, LabImg};

void RGB2Lab(float, float, float, float *, float *, float *);
void Lab2RGB(float, float, float, float *, float *, float *);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int min_index, i, j, pt, imgpixels, npoints, nclusters, maxclusters, qindex;
  const int *dims;
  enum imgtype type;
  double *Labdata, *answer, *mask, *map, dist, tdist;
  unsigned char *RGBdata, *quantized;
  Coord *clusters, *points;
  float L, a, b;
  int *index;

  /* Error checking */
  if (nrhs < 1 || nrhs > 3)
    mexErrMsgTxt("Incorrect number of arguments supplied");

  if (nlhs > 3)
    mexErrMsgTxt("No more than 3 outputs");

  if (mxGetNumberOfDimensions(prhs[0]) != 3)
    mexErrMsgTxt("Image should be M x N x 3");

  if (!mxIsUint8(prhs[0]) && !mxIsDouble(prhs[0]))
    mexErrMsgTxt("Image data should be uint8 or double");

  if (nrhs >= 2 && !(mxGetM(prhs[1]) == 1 && mxGetN(prhs[1]) == 1) &&
      !(mxGetM(prhs[1]) == mxGetM(prhs[0]) && 
	mxGetN(prhs[1]) == mxGetN(prhs[0]) / 3))
    mexErrMsgTxt("Mask must be a scalar or same height and width as image");

  if (nrhs >= 3 && !(mxGetM(prhs[2]) == 1 && mxGetN(prhs[2]) == 1))
    mexErrMsgTxt("Maximum number of clusters must be a scalar");

  /* Collect meta-data about each argument */
  imgpixels = mxGetM(prhs[0]) * mxGetN(prhs[0]) / 3; /* mxGetN counts 2-d */
  if (mxIsUint8(prhs[0])) {
    type = RGBImg;
    RGBdata = (unsigned char *)mxGetData(prhs[0]);
    Labdata = NULL;
  } else {
    type = LabImg;
    Labdata = mxGetData(prhs[0]);
    RGBdata = NULL;
  }

  if (nrhs >= 2 && (mxGetM(prhs[1]) == mxGetM(prhs[0]) && 
		    mxGetN(prhs[1]) == mxGetN(prhs[0]) / 3))
    mask = mxGetData(prhs[1]);
  else
    mask = NULL;

  if (nrhs >= 3)
    maxclusters = mxGetScalar(prhs[2]);
  else
    maxclusters = MAXCLUSTERS;

  /* Determine size of data */
  if (nrhs < 2 || mxGetM(prhs[1]) == 1 && mxGetN(prhs[1]) == 1)
    npoints = imgpixels;
  else {
    npoints = 0;
    for (i = 0; i < imgpixels; i++)
      npoints += mask[i] != 0;
    if (npoints == 0) { /* Mask is empty */
      plhs[0] = mxCreateDoubleMatrix(3, 0, mxREAL);
      return;
    }
  }

  /* Allocate output arguments */

  plhs[0] = mxCreateDoubleMatrix(ROWS, maxclusters, mxREAL);
  answer = mxGetPr(plhs[0]);

  if (nlhs > 1) { /* Allocate quantized image map */
    dims = mxGetDimensions(prhs[0]);
    plhs[1] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
    quantized = (unsigned char *)mxGetPr(plhs[1]);
  } else
    quantized = NULL;

  if (nlhs > 2) { 
    dims = mxGetDimensions(prhs[0]);
    plhs[2] = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    map = mxGetPr(plhs[2]);
  } else
    map = NULL;

  /* Set up data structure for vector quantization */
  points = (Coord *)mxCalloc(npoints * DIM, sizeof(Coord));
  index = (int *)mxCalloc(npoints, sizeof(int));
  pt = 0;
  for (i = 0; i < imgpixels; i++)
    if (npoints == imgpixels || mask[i] == 1) {
      if (type == RGBImg) {
	RGB2Lab((float)RGBdata[i], (float)RGBdata[i+imgpixels],
		(float)RGBdata[i+2*imgpixels], &L, &a, &b);
	points[pt * DIM] = (Coord) L;
	points[pt * DIM + 1] = (Coord) a;
	points[pt * DIM + 2] = (Coord) b;
      } else {
	points[pt * DIM] = (Coord) Labdata[i];
	points[pt * DIM + 1] = (Coord) Labdata[i+imgpixels];
	points[pt * DIM + 2] = (Coord) Labdata[i+2*imgpixels];
      }
      /*points[pt * DIM + 3] = floor(i / dims[0]) * 0.0;*/
      /*points[pt * DIM + 4] = (i - floor(i / dims[0]) * dims[0]) * 0.0;*/
      pt++;
    }
	
  bs(points, npoints, maxclusters, &nclusters, &clusters, index);

  /* Copy clusters to output array */
  mxSetN(plhs[0], nclusters);
  for (i = 0; i < nclusters; i++) {
    answer[i * ROWS] = (double) clusters[i * DIM];
    answer[i * ROWS + 1] = (double) clusters[i * DIM + 1];
    answer[i * ROWS + 2] = (double) clusters[i * DIM + 2];
  }

  /* Map pixels to clusters */
  qindex = -1;
  for (i = 0; i < npoints; i++) {
    if (DEBUG)
      mexPrintf("%d) %.2f %.2f %.2f %.2f %.2f | ", i + 1, points[i * DIM], 
		points[i * DIM + 1], points[i * DIM + 2],
		points[i * DIM + 3], points[i * DIM + 4]);

    /* Increment the counter of the number of pixels */
    answer[index[i] * ROWS + 3]++;

    /* Compute the variance */
    dist = 0.0;
    tdist = points[i * DIM] - answer[index[i] * ROWS];
    dist += tdist * tdist;
    tdist = points[i * DIM + 1] - answer[index[i] * ROWS + 1];
    dist += tdist * tdist;
    tdist = points[i * DIM + 2] - answer[index[i] * ROWS + 2];
    dist += tdist * tdist;
    answer[index[i] * ROWS + 4] += dist;

    if (quantized) {
      qindex++;
      while (mask && mask[qindex] == 0)
	qindex++;
      Lab2RGB((Coord)clusters[index[i] * DIM], 
	      (Coord)clusters[index[i] * DIM + 1],
	      (Coord)clusters[index[i] * DIM + 2], &L, &a, &b);
      quantized[qindex] = (unsigned char) L;
      quantized[qindex + imgpixels] = (unsigned char) a;
      quantized[qindex + 2 * imgpixels] = (unsigned char) b;
      if (map)
	map[qindex] = index[i];
    }
  }

  /* Compute the final variances */
  for (i = 0; i < nclusters; i++)
    if (answer[i * ROWS + 3] > 0)
      answer[i * ROWS + 4] /= answer[i * ROWS + 3]; /* Biased estimate */

  mxFree(points);
  mxFree(index);
  free(clusters);
}




