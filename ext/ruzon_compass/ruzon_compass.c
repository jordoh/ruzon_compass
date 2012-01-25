#include <ruby.h>
#include <math.h>
#include "compass.h"

static VALUE ruzon_compass(VALUE self, VALUE data, VALUE width, VALUE height, VALUE sigmas) {
  Check_Type(data, T_ARRAY);
  Check_Type(sigmas, T_ARRAY);

  int i, j;

  int rows = NUM2INT(height);
  int cols = NUM2INT(width);

  int image_data_length = RARRAY_LEN(data);
//  unsigned char* img_data = (unsigned char*)malloc(sizeof(unsigned char) * image_data_length);
//  for (i = 0; i < image_data_length; ++i)
//      img_data[i] = NUM2INT(RARRAY_PTR(data)[i]);
  double* img_data = (double*)malloc(sizeof(double) * image_data_length);
  for (i = 0; i < image_data_length; ++i)
      img_data[i] = NUM2DBL(RARRAY_PTR(data)[i]);

  /* Argument #1: one or more standard deviation values */
  /* Modified 31 July 1999 so that standard deviations (sigmas) rather than
   * radii are specified; also, sigmas can be floating-point values */
  int num_sigmas = (int)RARRAY_LEN(sigmas);
  double* sigmas_data = (double*)malloc(sizeof(double) * num_sigmas);
  for (i = 0; i < num_sigmas; ++i)
    sigmas_data[i] = rb_num2dbl(RARRAY_PTR(sigmas)[i]);

  double max_sigma = -1.0;
  for (i = 0; i < num_sigmas; ++i)
    if (sigmas_data[i] > max_sigma)
        max_sigma = sigmas_data[i];

  double max_radius = ceil(3.0 * max_sigma);
  //if (max_sigma * 2 > rows || max_sigma * 2 > cols)
  //    mexErrMsgTxt("Image is too small for maximum scale chosen");

  /* Argument #2: the spacing for each application of the operator */
  double* spacing = (double*)malloc(sizeof(double) * num_sigmas);
  for (i = 0; i < num_sigmas; ++i)
    spacing[i] = 1.0;

  /* Argument #3: Dimensions of the sub-image to find edges over */
  double* dimensions = (double*)malloc(sizeof(double) * 4);
  dimensions[0] = max_radius;
  dimensions[1] = max_radius;
  dimensions[2] = rows - max_radius;
  dimensions[3] = cols - max_radius;

  /* Argument #4: Number of angles (alpha values in ICCV99) to compute */
  int num_angles = 1;
  double* angles = (double*)malloc(sizeof(double) * num_angles);
  angles[0] = 180;


  /* Argument #5: the number of wedges in one quarter of the circle */
  int num_wedges = 6;

  /* Argument #6: the maximum number of clusters in a color signature */
  double* max_clusters = (double*)malloc(sizeof(double) * num_sigmas);
  for (i = 0; i < num_sigmas; i++) {
    max_clusters[i] = 10;
    //if (maxclusters[i] > MAXCLUSTERS)
  	//mexErrMsgTxt("Maximum number of clusters greater than allowed");
  }

  /* Allocate output arguments */ 
  int stength_height = rows;//(int)(dimensions[2] - dimensions[0])/(int)spacing[0] + 1;
  int stength_width = cols;//(int)(dimensions[3] - dimensions[1])/(int)spacing[0] + 1;
  int strength_image_length = stength_width * stength_height;
  double** strength = (double**)malloc(sizeof(double*) * num_sigmas * num_angles);
  strength[0] = (double*)malloc(sizeof(double) * strength_image_length);

//    void Compass(void *imgdata, int imgrows, int imgcols, enum imgtype type,
//    	     double *sigmas, int numsigmas, int maxradius, double *spacing,
//    	     double *dimensions, double *angles, int numangles, int nwedges,
//    	     double *maxclusters, int plot, mxArray **strength,
//    	     mxArray **abnormality, mxArray **orientation,
//    	     mxArray **uncertainty)
  Compass(
    /* void *imgdata          */ img_data,
    /* int imgrows             */ rows,
    /* int imgcols              */ cols,
    /* enum imgtype type */ LabImg,//RGBImg,
    sigmas_data, num_sigmas, max_radius, spacing,
    dimensions, angles, num_angles, num_wedges,
    max_clusters, 0, strength, NULL, NULL, NULL
  );

  VALUE result = rb_ary_new();
  for (i=0; i<num_sigmas * num_angles; ++i) {
    double* strength_image = strength[i];
    
    for (j=0; j<strength_image_length; ++j) {
      rb_ary_push(result, DBL2NUM(strength_image[j]));
    }
    
    free(strength_image);
  }
  free(strength);

  free(max_clusters);
  free(angles);
  free(dimensions);
  free(spacing);
  free(sigmas_data);
  free(img_data);
  
  return result;
}

/* ruby calls this to load the extension */
void Init_ruzon_compass(void) {
  VALUE target = rb_define_module("RuzonCompass");
  rb_define_singleton_method(target, "compass", ruzon_compass, 4);
}