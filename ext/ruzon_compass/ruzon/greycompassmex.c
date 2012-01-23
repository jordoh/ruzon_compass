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
 * This file contains the MATLAB wrapper for the compass operator.         *
 *                                                                         *
 * 2 August 1999 -- extended for use with greyscale images                 *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mex.h>
#include "matrix.h"
#include "greycompass.h"

#define STRENGTH 0
#define ORIENTATION 1
#define ABNORMALITY 2
#define UNCERTAINTY 3

void Error(char *text)
{
  mexErrMsgTxt(text);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int imgrows, imgcols, numsigmas, numangles, nwedges, maxradius;
  int i, j, k, anglewedges, index, plot = 0;
  enum imgtype type;
  double *sigmas, *spacing, *dimensions, *angles, wedgeangle;
  double maxsigma = -1.0;
  Matrix *strength = NULL, *orientation = NULL, *uncertainty = NULL;
  Matrix *abnormality = NULL;
  void *imgdata;
  int dims[3];
  char buf[20] = "xxxxxxxxxxxxxxxxxxx\0";

  /* Before error checking, remove the optional string argument at the end */
  if (nrhs > 0 && mxIsChar(prhs[nrhs-1])) {
    if (mxGetString(prhs[nrhs-1], buf, 20) == 0 && buf[0] == 'p')
      plot = 1;
    nrhs--; /* Don't want to process it numerically later on */
  }

  /* Error checking (duh) */
  if (nrhs < 2 || nrhs > 6)
    mexErrMsgTxt("Incorrect number of arguments supplied");

  if (nlhs < 1 || nlhs > 4)
    mexErrMsgTxt("Between 1 and 4 output matrices should be supplied");

  if (mxGetNumberOfDimensions(prhs[0]) != 2)
    mexErrMsgTxt("Image should be M x N (x 1)");

  if (!mxIsUint8(prhs[0]) && !mxIsDouble(prhs[0]))
    mexErrMsgTxt("Image data should be uint8 or double");

  if (mxGetM(prhs[1]) != 1)
    mexErrMsgTxt("Sigmas must be a scalar or a row vector");

  if (mxGetN(prhs[1]) == 0)
    mexErrMsgTxt("Empty sigma vector is not allowed");

  if (nrhs >= 3 && mxGetM(prhs[2]) != 1)
    mexErrMsgTxt("Spacing must be a scalar or a row vector");

  if (nrhs >= 3 && mxGetN(prhs[2]) != 1 && mxGetN(prhs[2]) != mxGetN(prhs[1]))
    mexErrMsgTxt("Spacing must be a scalar or same length as scale vector");

  if (nrhs >= 4 && mxGetM(prhs[3]) != 1)
    mexErrMsgTxt("Image sub-dimensions must be a scalar or a row vector");

  if (nrhs >= 4 && (mxGetN(prhs[3]) == 3 || mxGetN(prhs[3]) > 4))
    mexErrMsgTxt("Image sub-dimension vector must have 1, 2, or 4 entries");

  if (nrhs >= 5 && mxGetM(prhs[4]) != 1)
    mexErrMsgTxt("Angles must be a scalar or a row vector");

  /* Collect meta-data about each argument */
  /* Argument #0: the image */
  imgrows = mxGetM(prhs[0]);
  imgcols = mxGetN(prhs[0]);
  if (mxIsUint8(prhs[0]))
    type = RGBImg;
  else
    type = LabImg;
  imgdata = mxGetData(prhs[0]);

  /* Argument #1: one or more standard deviation values */
  /* Modified 31 July 1999 so that standard deviations (sigmas) rather than
   * radii are specified; also, sigmas can be floating-point values */
  numsigmas = mxGetN(prhs[1]);
  sigmas = mxGetPr(prhs[1]);
  for (i = 0; i < numsigmas; i++)
    if (sigmas[i] > maxsigma)
      maxsigma = sigmas[i];
  maxradius = ceil(3 * maxsigma);
  if (maxradius * 2 > imgrows || maxradius * 2 > imgcols)
    mexErrMsgTxt("Image is too small for maximum scale chosen");

  /* Argument #2: the spacing for each application of the operator */
  if (nrhs >= 3 && mxGetN(prhs[2]) > 1)
    spacing = mxGetPr(prhs[2]);
  else { 
    spacing = (double *)mxCalloc(numsigmas,sizeof(double));
    for (i = 0; i < numsigmas; i++)
      spacing[i] = (nrhs >= 3) ? mxGetScalar(prhs[2]) : 1;
  }

  /* Argument #3: Dimensions of the sub-image to find edges over */
  if (nrhs >= 4 && mxGetN(prhs[3]) == 4) {
    dimensions = mxGetPr(prhs[3]); /* Dimensions specified */
    if (dimensions[0] < maxradius) {
      mexWarnMsgTxt("Dimension entry 1 below minimum row");
      dimensions[0] = maxradius;
    }
    if (dimensions[1] < maxradius) {
      mexWarnMsgTxt("Dimension entry 2 below minimum column");
      dimensions[1] = maxradius;
    }
    if (dimensions[2] > imgrows - maxradius) {
      mexWarnMsgTxt("Dimension entry 3 above maximum row");
      dimensions[2] = imgrows - maxradius;
    }
    if (dimensions[3] > imgcols - maxradius) {
      mexWarnMsgTxt("Dimension entry 4 above maximum column");
      dimensions[3] = imgcols - maxradius;
    }
    if (dimensions[2] < dimensions[0] || dimensions[3] < dimensions[1])
      mexErrMsgTxt("Dimensions are [TR LC BR RC] with BR >= TR and RC >= LC");
  } else {
    dimensions = (double *)mxCalloc(4,sizeof(double));
    if (nrhs < 4 || mxGetN(prhs[3]) == 1) { /* Whole image */
      dimensions[0] = maxradius;
      dimensions[1] = maxradius;
      dimensions[2] = imgrows - maxradius;
      dimensions[3] = imgcols - maxradius;
    } else { /* One point */
      dimensions[0] = mxGetPr(prhs[3])[0];
      dimensions[1] = mxGetPr(prhs[3])[1];
      dimensions[2] = mxGetPr(prhs[3])[0];
      dimensions[3] = mxGetPr(prhs[3])[1];
      if (dimensions[0] < maxradius || dimensions[1] < maxradius ||
	  dimensions[2] > imgrows - maxradius || 
          dimensions[3] > imgcols - maxradius)
	mexErrMsgTxt("Point is inoperable by maximum scale");
    }
  }

  /* Argument #4: Number of angles (alpha values in ICCV99) to compute */
  if (nrhs >= 5 && mxGetN(prhs[4]) >= 1) {
    numangles = mxGetN(prhs[4]);
    angles = mxGetPr(prhs[4]);
    for (i = 0; i < numangles; i++)
      if (angles[i] <= 0 || angles[i] > 180) {
	mexErrMsgTxt("All angles must be between 0 and 180");
      }
  } else {
    numangles = 1;
    angles = (double *)mxCalloc(numangles,sizeof(double));
    angles[0] = 180;
  }

  /* Argument #5: the number of wedges in one quarter of the circle */
  nwedges = (nrhs >= 6) ? mxGetScalar(prhs[5]) : 6;
  if (nwedges * 2 > MAXWEDGES)
    mexErrMsgTxt("Too many wedges");
  
  wedgeangle = 90.0 / nwedges;
  for (i = 0; i < numangles; i++)
    if (angles[i]/wedgeangle != floor(angles[i]/wedgeangle))
      mexErrMsgTxt("Angles chosen not compatible with number of wedges");

  /* Allocate output arguments */
  for (i = 0; i < nlhs; i++) {
    if (i == STRENGTH && plot == 1) { /* Allocate compass plots */
      if (numsigmas == 1 && numangles == 1) {
	dims[0] = (int)(dimensions[2] - dimensions[0])/(int)spacing[0] + 1;
	dims[1] = (int)(dimensions[3] - dimensions[1])/(int)spacing[0] + 1;
	dims[2] = (angles[0] == 180) ? 2 * nwedges : 4 * nwedges; 
	plhs[STRENGTH] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      } else { /* Create a list of compass plots */
	plhs[STRENGTH] = mxCreateCellMatrix(numangles, numsigmas);
	for (j = 0; j < numsigmas; j++) {
	  dims[0] = (int)(dimensions[2] - dimensions[0])/(int)spacing[j] + 1;
	  dims[1] = (int)(dimensions[3] - dimensions[1])/(int)spacing[j] + 1;
	  for (k = 0; k < numangles; k++) {
	    dims[2] = (angles[k] == 180) ? 2 * nwedges : 4 * nwedges;
	    mxSetCell(plhs[STRENGTH], j * numangles + k, 
		      mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
	  }
	}
      }
    } else {
      /* Allow multiple pages for uncertainty and orientation */
      dims[2] = (i == STRENGTH || i == ABNORMALITY) ? 1 : MAXRESPONSES;
      if (numsigmas == 1 && numangles == 1) { /* One matrix per argument */
	dims[0] = (int)(dimensions[2] - dimensions[0])/(int)spacing[0] + 1;
	dims[1] = (int)(dimensions[3] - dimensions[1])/(int)spacing[0] + 1;
	plhs[i] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      } else { /* Create a list for each argument */
	plhs[i] = mxCreateCellMatrix(numangles, numsigmas);
	for (j = 0; j < numsigmas; j++) {
	  dims[0] = (int)(dimensions[2] - dimensions[0])/(int)spacing[j] + 1;
	  dims[1] = (int)(dimensions[3] - dimensions[1])/(int)spacing[j] + 1;
	  for (k = 0; k < numangles; k++) 
	    mxSetCell(plhs[i], j * numangles + k, mxCreateNumericArray(3,
                      dims, mxDOUBLE_CLASS, mxREAL));
	}
      }
    }
  }
  switch (nlhs) {  /* This has no breaks; be careful when modifying */
  case 4: 
    uncertainty = (Matrix *)mxCalloc(numsigmas*numangles,sizeof(Matrix));
    if (numsigmas == 1 && numangles == 1) {
      uncertainty[0].rows = mxGetM(plhs[UNCERTAINTY]);
      uncertainty[0].cols = mxGetN(plhs[UNCERTAINTY]) / MAXRESPONSES;
      uncertainty[0].sheets = MAXRESPONSES;
      uncertainty[0].ptr = mxGetPr(plhs[UNCERTAINTY]);
    } else
      for (i = 0; i < numsigmas * numangles; i++) {
	uncertainty[i].rows = mxGetM(mxGetCell(plhs[UNCERTAINTY],i));
	uncertainty[i].cols = mxGetN(mxGetCell(plhs[UNCERTAINTY],i)) / 
	  MAXRESPONSES;
	uncertainty[i].sheets = MAXRESPONSES;
	uncertainty[i].ptr = mxGetPr(mxGetCell(plhs[UNCERTAINTY],i));
      }
  case 3: 
    abnormality = (Matrix *)mxCalloc(numsigmas*numangles,sizeof(Matrix));
    if (numsigmas == 1 && numangles == 1) {
      abnormality[0].rows = mxGetM(plhs[ABNORMALITY]);
      abnormality[0].cols = mxGetN(plhs[ABNORMALITY]);
      abnormality[0].sheets = 1;
      abnormality[0].ptr = mxGetPr(plhs[ABNORMALITY]);
    } else
      for (i = 0; i < numsigmas * numangles; i++) {
	abnormality[i].rows = mxGetM(mxGetCell(plhs[ABNORMALITY],i));
	abnormality[i].cols = mxGetN(mxGetCell(plhs[ABNORMALITY],i));
	abnormality[i].sheets = 1;
	abnormality[i].ptr = mxGetPr(mxGetCell(plhs[ABNORMALITY],i));
      }
  case 2: 
    orientation = (Matrix *)mxCalloc(numsigmas*numangles,sizeof(Matrix));
    if (numsigmas == 1 && numangles == 1) {
      orientation[0].rows = mxGetM(plhs[ORIENTATION]);
      orientation[0].cols = mxGetN(plhs[ORIENTATION]) / MAXRESPONSES;
      orientation[0].sheets = MAXRESPONSES;
      orientation[0].ptr = mxGetPr(plhs[ORIENTATION]);
    } else
      for (i = 0; i < numsigmas * numangles; i++) {
	orientation[i].rows = mxGetM(mxGetCell(plhs[ORIENTATION],i));
	orientation[i].cols = mxGetN(mxGetCell(plhs[ORIENTATION],i)) / 
	  MAXRESPONSES;
	orientation[i].sheets = 1;
	orientation[i].ptr = mxGetPr(mxGetCell(plhs[ORIENTATION],i));
      }
  case 1: 
    strength = (Matrix *)mxCalloc(numsigmas*numangles,sizeof(Matrix));
    if (numsigmas == 1 && numangles == 1) {
      strength[0].rows = mxGetM(plhs[STRENGTH]);
      strength[0].sheets = plot ? ((angles[0] == 180) ? 2 * nwedges : 
				   4 * nwedges) : 1;
      strength[0].cols = mxGetN(plhs[STRENGTH]) / strength[0].sheets;
      strength[0].ptr = mxGetPr(plhs[STRENGTH]);
    } else
      for (i = 0; i < numsigmas; i++)
	for (j = 0; j < numangles; j++) {
	  if (numangles == 1) /* Still a 1-D cell array */
	    index = i;
	  else
	    index = i * numsigmas + j;
	strength[index].rows = mxGetM(mxGetCell(plhs[STRENGTH],index));
	strength[index].sheets = plot ? ((angles[j] == 180) ? 2 * nwedges : 
					 4 * nwedges) : 1;
	strength[index].cols = mxGetN(mxGetCell(plhs[STRENGTH],index)) / 
	  strength[index].sheets;
	strength[index].ptr = mxGetPr(mxGetCell(plhs[STRENGTH],index));
      }
  }

  Compass(imgdata, imgrows, imgcols, type, sigmas, numsigmas, maxradius, 
	  spacing, dimensions, angles, numangles, nwedges, strength, 
	  abnormality, orientation, uncertainty);
}
