#include <ruby.h>
#include <math.h>
#include "bs.h"
#include "compass.h"

static VALUE ruzon_compass(VALUE self, VALUE data, VALUE width, VALUE height, VALUE sigmas, VALUE angles) {
  Check_Type(data, T_ARRAY);
  Check_Type(sigmas, T_ARRAY);
  if (angles != Qnil) {
    Check_Type(angles, T_ARRAY);
  }

  int i, j;

  const int rows = NUM2INT(height);
  const int cols = NUM2INT(width);
  const int imageDataLength = RARRAY_LEN(data);
  VALUE const * const imageDataPtr = RARRAY_PTR(data);
  
  /* Argument #0: the image */
  const enum imgtype imageDataType = TYPE(imageDataPtr[0]) == T_FLOAT ? LabImg : RGBImg;

  /* Argument #1: one or more standard deviation values */
  /* Modified 31 July 1999 so that standard deviations (sigmas) rather than
   * radii are specified; also, sigmas can be floating-point values */
  const int numSigmas = RARRAY_LEN(sigmas);
  VALUE const * const sigmasPtr = RARRAY_PTR(sigmas);
  
  double* sigmasData = (double*)malloc(sizeof(double) * numSigmas);
  for (i = 0; i < numSigmas; ++i) {
    sigmasData[i] = NUM2DBL(sigmasPtr[i]);
  }

  double maxSigma = -1.0;
  for (i = 0; i < numSigmas; ++i)
    if (sigmasData[i] > maxSigma)
        maxSigma = sigmasData[i];

  double maxRadius = ceil(3.0 * maxSigma);
  if (maxSigma * 2 > rows || maxSigma * 2 > cols) {
    rb_raise(rb_eStandardError, "Image is too small for maximum scale chosen");
  }

  /* Argument #2: the spacing for each application of the operator */
  double* spacing = (double*)malloc(sizeof(double) * numSigmas);
  for (i = 0; i < numSigmas; ++i) {
    spacing[i] = 1.0;
  }

  /* Argument #3: Dimensions of the sub-image to find edges over */
  double* dimensions = (double*)malloc(sizeof(double) * 4);
  dimensions[0] = maxRadius;
  dimensions[1] = maxRadius;
  dimensions[2] = rows - maxRadius;
  dimensions[3] = cols - maxRadius;

  /* Argument #4: Number of angles (alpha values in ICCV99) to compute */
  double* anglesData;
  int numAngles;
  if (angles != Qnil) {
    numAngles = RARRAY_LEN(angles);
    anglesData = (double*)malloc(sizeof(double) * numAngles);
    
    VALUE const * const anglesPtr = RARRAY_PTR(angles);
    for (i = 0; i < numAngles; i++) {
      anglesData[i] = NUM2DBL(anglesPtr[i]);
      if (anglesData[i] <= 0 || anglesData[i] > 180) {
        rb_raise(rb_eStandardError, "All angles must be between 0 and 180");
      }
    }
  } else {
    numAngles = 1;
    anglesData = (double*)malloc(sizeof(double) * numAngles);
    
    anglesData[0] = 180.0;
  }

  /* Argument #5: the number of wedges in one quarter of the circle */
  const int numWedges = 6;
  if (numWedges * 2 > MAXWEDGES) {
    rb_raise(rb_eStandardError, "Too many wedges");
  }
  
  const double wedgeAngle = 90.0 / numWedges;
  for (i = 0; i < numAngles; i++) {
    if (anglesData[i]/wedgeAngle != floor(anglesData[i]/wedgeAngle)) {
      rb_raise(rb_eStandardError, "Angles chosen not compatible with number of wedges");
    }
  }

  /* Argument #6: the maximum number of clusters in a color signature */
  double* maxClusters = (double*)malloc(sizeof(double) * numSigmas);
  for (i = 0; i < numSigmas; i++) {
    maxClusters[i] = 10;
    if (maxClusters[i] > MAXCLUSTERS) {
  	  rb_raise(rb_eStandardError, "Maximum number of clusters greater than allowed");
  	}
  }

  /* Allocate output arguments */ 
  int strengthHeight = rows;//(int)(dimensions[2] - dimensions[0])/(int)spacing[0] + 1;
  int strengthWidth = cols;//(int)(dimensions[3] - dimensions[1])/(int)spacing[0] + 1;
  int strengthLength = strengthWidth * strengthHeight;
  double** strength = (double**)malloc(sizeof(double*) * numSigmas * numAngles);
  strength[0] = (double*)malloc(sizeof(double) * strengthLength);

  Compass(
    /* VALUE imgdata        */ data,
    /* int imgrows          */ rows,
    /* int imgcols          */ cols,
    /* enum imgtype type    */ imageDataType,
    /* double* sigmas       */ sigmasData, 
    /* int numsigmas        */ numSigmas, 
    /* int maxradius        */ maxRadius, 
    /* double* spacing      */ spacing,
    /* double* dimensions   */ dimensions, 
    /* double* angles       */ anglesData, 
    /* int numangles        */ numAngles, 
    /* int nwedges          */ numWedges,
    /* double* maxclusters  */ maxClusters, 
    /* int plot             */ 0, 
    /* double** strength    */ strength, 
    /* double** abnormality */ NULL, 
    /* double** orientation */ NULL, 
    /* double** uncertainty */ NULL
  );

  VALUE result = rb_ary_new2(numSigmas * numAngles);
  for (i = 0; i < numSigmas * numAngles; ++i) {
    double* strengthImage = strength[i];
    
    VALUE strengthResult = rb_ary_new2(strengthLength);;
    for (j = 0; j < strengthLength; ++j) {
      rb_ary_store(strengthResult, j, DBL2NUM(strengthImage[j]));
    }
    free(strengthImage);
    
    rb_ary_push(result, strengthResult);
  }
  free(strength);
  
  free(maxClusters);
  free(anglesData);
  free(dimensions);
  free(spacing);
  free(sigmasData);
  
  return result;
}

void Init_ruzon_compass(void) {
  VALUE target = rb_define_module("RuzonCompass");
  rb_define_singleton_method(target, "compass", ruzon_compass, 5);
}
