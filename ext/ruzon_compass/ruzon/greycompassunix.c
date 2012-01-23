/***************************************************************************
 * The Greyscale Compass Operator -- an algorithm to find edges in         *
 *   greyscale images by representing regions with distributions           *
 *                                                                         *
 * Mark A. Ruzon                                                           *
 * Stanford Vision Laboratory                                              *
 * October, 1997                                                           *
 * Copyright 1997-1999.  All rights reserved.                              *
 *                                                                         *
 * Details -- The compass operator uses a circular window centered at a    *
 *   junction where 4 pixel squares meet.  The needle is a diameter at a   *
 *   given orientation.  The intensity distributions of the two            *
 *   semicircles are computed, and the distance between them is measured   *
 *   using the Earth Mover's Distance (EMD).  The maximum EMD over all     *
 *   orientations gives the edge strength and orientation at that point.   *
 *   Uncertainty and abnormality are also measured.                        *
 *                                                                         *
 * This file contains the Unix wrapper for the compass operator.           *
 *                                                                         *
 * 2 August 1999 -- extended for use with greyscale images                 *
 *                                                                         *
 * 7 October 1999 -- creation of the Unix wrapper.  Most of the            *
 *   functionality of the MATLAB version (e.g. subimages, spacing, plot,   *
 *   maxclusters, multiple sigmas, abnormality, uncertainty, etc.)         *
 *   is being removed.                                                     *
 *                                                                         *
 * 21 October 1999 -- allowed only compass operator and NMS routine, or    *
 *                    only thresholding, or both to be computed.           *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "image.h"
#include "greycompass.h"

#define MAXVAL 255
#define MAXLINE 2000
#define PI 3.14159265358979323846
#define MODEC 0x2
#define MODET 0x1

void Error(char *text)
{
  fprintf(stderr,text);
  exit(1);
}


void PrintHelp(char *programname)
{
  fprintf(stderr,"Function: %s <filename> [options]\n", programname);
  fprintf(stderr,"   <filename> is an input PGM or double (mode T) file\n");
  fprintf(stderr,"   and options include:\n");
  fprintf(stderr,"   -o <outfile>, the edge map or edge strength (mode C) file (default stdout)\n");
  fprintf(stderr,"   -s <sigma>, the standard deviation of the Gaussian (default 1.0)\n");
  fprintf(stderr,"   -m C|T, the mode, compass operator or thresholding (default: both)\n"); 
  fprintf(stderr,"   -l <low>, the low threshold for edge detection (default 0.5)\n");
  fprintf(stderr,"   -h <high>, the high threshold for edge detection (default 0.7)\n");
  exit(1);
}


void ParseArguments(int argc, char *argv[], it_image **input, FILE **output,
		    double sigma[], double *low, double *high, int *mode)
{
  FILE *FP;
  int curarg = 2;
  char modechar;

  /* Filename is assumed to be in argv[1] */
  if (argc < 2)
    PrintHelp(argv[0]);

  FP = fopen(argv[1],"r");
  if (!FP)
    Error("Could not open image file for reading\n");

  while (curarg < argc) {
    if (argv[curarg][0] == '-') {  /* New argument has appeared */
      switch(argv[curarg][1]) {
      case 'o': /* Output */
	if (curarg < argc) {
	  *output = i_open_file(argv[curarg+1],-1,"w");
	  if (!*output)
	    Error("Output file could not be opened\n");
	} else
	  Error("No output file specified\n");
	break;	
      case 's': /* Sigma */
	if (curarg < argc) {
	  sigma[0] = strtod(argv[curarg+1],NULL);
	  if (sigma[0] <= 0.0)
	    Error("Sigma must be positive\n");
	} else
	  Error("No sigma value specified\n");
	break;
      case 'l': /* Low */
	if (curarg < argc) {
	  *low = strtod(argv[curarg+1],NULL);
	  if (*low < 0.0 || *low > 1.0)
	    Error("Low threshold must lie in [0,1]\n");
	} else
	  Error("No low threshold specified\n");
	break;
      case 'h': /* High */
	if (curarg < argc) {
	  *high = strtod(argv[curarg+1],NULL);
	  if (*high < 0.0 || *high > 1.0)
	    Error("High threshold must lie in [0,1]\n");
	} else
	  Error("No high threshold specified\n");
	break;
      case 'm': /* Mode */
	if (curarg < argc) {
	  modechar = argv[curarg+1][0];
	  if (modechar == 'C' || modechar == 'c')
	    *mode = MODEC;
	  else if (modechar == 'T' || modechar == 't')
	    *mode = MODET;
	  else
	    Error("Mode must be 'C' or 'T'\n");
	} else
	  Error("No mode specified\n");
	break;	
      default: /* Help */
	fprintf(stderr,"Unrecognized command option\n");
	PrintHelp(argv[0]);
	break;
      }
      curarg += 2;
    }
  }
  if (*high < *low)
    Error("High threshold cannot be below low threshold\n");

  if (*mode & MODEC)
    *input = i_read_image_file(FP, IT_BYTE, IM_CONTIG);
  else
    *input = i_read_image_file(FP, IT_DOUBLE, IM_CONTIG);
  if (!*input)
    Error("Bad input file format for given mode\n");

}


/*  NonMaximalSuppression
 *
 *  The standard Canny non-maximal suppression algorithm; looks in a 3x3
 *  neighborhood to determine if the pixel is a maximum in the direction
 *  perpendicular to that of the edge orientation.
 */
Matrix *NonMaximalSuppression(Matrix *strength, Matrix *orientation)
{
  int x, y, i, rows = strength->rows, cols = strength->cols;
  int pixels = strength->rows * strength->cols, maximum;
  double str1, str2;              /* interpolated edge strength */
  double a1, a2, b1, b2, c1, c2;  /* nearest pixels' edge strength */
  float ux, uy;                   /* weights of a, b, c, and str */
  double ori, str;       /* strength and orientation at center */
  Matrix *newstrength;

  /* Newstrength holds the NMS'ed strength values */
  newstrength = (Matrix *)malloc(sizeof(Matrix));
  newstrength->rows = strength->rows;
  newstrength->sheets = 1;
  newstrength->cols = strength->cols;
  newstrength->ptr = (double *)calloc(pixels, sizeof(double));

  /* For each pixel (except those on the boundary), check if 
     * the edge strength is a local (3x3 neighbourhood) maximum in 
     * the orientation perpendicular to the edge 
     */
  for (x = 1; x < rows - 1; x++)
    for (y = 1; y < cols - 1; y++) {

      str = strength->ptr[y * rows + x];
      maximum = 0;

      if (str == 0.0)
	continue;

      for (i = 0; i < MAXRESPONSES; i++) {
	ori = orientation->ptr[i * pixels + y * rows + x];

	if (ori == -1)
	  break;
	
	ux = fabs(cos(ori*PI/180));
	uy = fabs(sin(ori*PI/180));
	b1 = strength->ptr[y * rows + (x-1)];  /* Pixel above */
	b2 = strength->ptr[y * rows + (x+1)];  /* Pixel below */	  

	if (ori > 90) {
	  a1 = strength->ptr[(y+1) * rows + (x-1)];
	  c1 = strength->ptr[(y+1) * rows + x];
	  a2 = strength->ptr[(y-1) * rows + (x+1)];
	  c2 = strength->ptr[(y-1) * rows + x];
	} else {
	  a1 = strength->ptr[(y-1) * rows + (x-1)];
	  c1 = strength->ptr[(y-1) * rows + x];
	  a2 = strength->ptr[(y+1) * rows + (x+1)];
	  c2 = strength->ptr[(y+1) * rows + x];
	}
		
	str1 = ux*uy*a1 + ux*(1-uy)*b1 + (1-ux)*uy*c1 + (1-ux)*(1-uy)*str;
	str2 = ux*uy*a2 + ux*(1-uy)*b2 + (1-ux)*uy*c2 + (1-ux)*(1-uy)*str;
	    
	if (str > str1 && str >= str2) 
	  maximum = 1;
      }

      if (maximum)
	newstrength->ptr[y * rows + x] = str;
    }
  return(newstrength);
}


it_image *Hysteresis(Matrix *strength, double low, double high, int rows,
		     int cols, int radius)
{
  int x, y, start, end, nbr[2][8] = {{ 1,  1,  1,  0, -1, -1, -1, 0},
				     { 1,  0, -1, -1, -1,  0,  1, 1}};
  int edge, i, xp, yp, X[MAXLINE], Y[MAXLINE];
  it_image *img;
  
  img = i_create_image(cols, rows, IT_BYTE, IM_CONTIG);
  if (!img)
    Error("Could not create output image\n");

  /* First, mark all pixels above the high and below the low */	
  /* We don't check the borders because we didn't in NonMaximalSuppression */
  for (x = 0; x < strength->rows; x++)
    for (y = 0; y < strength->cols; y++) {
      if (x == 0 || y == 0 || x == strength->rows - 1 || 
	  y == strength->cols - 1)
	strength->ptr[y * strength->rows + x] = 0.0;
      if (strength->ptr[y * strength->rows + x] >= high) {
	im_byte_value(img, y + radius, x + radius) = MAXVAL;
	strength->ptr[y * strength->rows + x] = 0.0;
      } else if (strength->ptr[y * strength->rows + x] < low)
	strength->ptr[y * strength->rows + x] = 0.0;
    }

  /* The only pixels left with non-zero strength fall between the two
   * thresholds.  Group adjacent pixels together, and if they touch a
   * marked edge, they all become marked.
   */
  for (x = 1; x < strength->rows - 1; x++)
    for (y = 1; y < strength->cols - 1; y++)
      if (strength->ptr[y * strength->rows + x] > 0.0) {
	strength->ptr[y * strength->rows + x] = 0.0;
	start = 0;
	edge = 0;
	X[start] = x;
	Y[start] = y;
	end = 1;
	while (start != end) {
	  for (i = 0; i < 8; i++) {
	    xp = X[start] + nbr[0][i];
	    yp = Y[start] + nbr[1][i];
	    edge = edge || 
	      (im_byte_value(img, yp + radius, xp + radius) == MAXVAL);
	    if (strength->ptr[yp * strength->rows + xp] > 0.0) {
	      strength->ptr[yp * strength->rows + xp] = 0.0;
	      X[end] = xp;
	      Y[end] = yp;
	      if (++end == MAXLINE)
		Error("Edge Buffer Full");
	    }
	  }
	  start++;
	}
	
	if (edge)
	  for (i = 0; i < end; i++)
	    im_byte_value(img, Y[i] + radius, X[i] + radius) = MAXVAL;
      }
  return(img);
}


int main(int argc, char *argv[])
{
  int imgrows, imgcols, numsigmas = 1, numangles = 1, nwedges = 6, radius;
  int i, j, k, anglewedges, index;
  enum imgtype type = RGBImg;
  double sigma[1] = {1.0}, spacing[1] = {1.0}, dimensions[4];
  double angles[1] = {180}, wedgeangle;
  double maxsigma = -1.0, *q;
  Matrix *strength, *orientation, *uncertainty = NULL, *abnormality = NULL;
  Matrix *NMS;
  unsigned char *imgdata, *p;
  int dims[3];

  double low = 0.5, high = 0.7;
  it_image *input, *output;
  FILE *FP = stdout;
  int mode = MODEC | MODET;

  ParseArguments(argc, argv, &input, &FP, sigma, &low, &high, &mode);

  /* Collect meta-data about the input */
  imgrows = input->height;
  imgcols = input->width;

  /* The it_image data structure represents images as an array of pointers
   * to rows of pixel values.  Since the compass operator was originally
   * developed in the context of MATLAB, we must convert this structure to
   * a one-dimensional array in column-major order.
   */
  if (mode & MODEC) {
    imgdata = (unsigned char *)calloc(imgrows*imgcols, sizeof(unsigned char));
    p = imgdata;
    for (i = 0; i < imgcols; i++)
      for (j = 0; j < imgrows; j++)
	*p++ = im_byte_value(input,i,j);
    
    radius = ceil(3 * sigma[0]);
    if (radius * 2 > imgrows || radius * 2 > imgcols)
      Error("Image is too small for sigma chosen");

    dimensions[0] = radius;
    dimensions[1] = radius;
    dimensions[2] = imgrows - radius;
    dimensions[3] = imgcols - radius;

    /* Allocate output arguments */
    orientation = (Matrix *)malloc(sizeof(Matrix));
    orientation->rows = dimensions[2] - dimensions[0] + 1;
    orientation->cols = dimensions[3] - dimensions[1] + 1;
    orientation->sheets = MAXRESPONSES;
    orientation->ptr = (double *)calloc(orientation->rows*orientation->cols*
					MAXRESPONSES, sizeof(double));
    strength = (Matrix *)malloc(sizeof(Matrix));
    strength->rows = orientation->rows;
    strength->sheets = 1;
    strength->cols = orientation->cols;
    strength->ptr = (double *)calloc(strength->rows*strength->cols,
				     sizeof(double));

    /* Execute the compass operator */
    Compass((void *)imgdata, imgrows, imgcols, type, sigma, numsigmas, 
	    radius, spacing, dimensions, angles, numangles, nwedges, 
	    strength, abnormality, orientation, uncertainty);

    /* Perform non-maximal supression */
    NMS = NonMaximalSuppression(strength, orientation);
    free(strength->ptr);
    free(strength);
    free(orientation->ptr);
    free(orientation);

    /* If thresholding is not specified, write the output file */
    if (!(mode & MODET)) {
      /* Output file should be as big as original image */
      output = i_create_image(imgcols, imgrows, IT_DOUBLE, IM_CONTIG);
      if (!output)
	Error("Could not create output image\n");
      q = NMS->ptr;
      for (i = 0; i < NMS->cols; i++)
	for (j = 0; j < NMS->rows; j++)
	  im_double_value(output,i+radius,j+radius) = *q++;
    }
  }
     
  /* Perform hysteresis thresholding with edge following */
  if (mode & MODET) {
    /* If the compass operator is not specified, read the input file */
    if (!(mode & MODEC)) {
      NMS = (Matrix *)malloc(sizeof(Matrix));
      NMS->rows = imgrows;
      NMS->sheets = 1;
      NMS->cols = imgcols;
      NMS->ptr = (double *)calloc(imgrows*imgcols, sizeof(double));
      q = NMS->ptr;
      for (i = 0; i < imgcols; i++)
	for (j = 0; j < imgrows; j++)
	  *q++ = im_double_value(input,i,j);
      /* In addition, there is no radius value anymore */
      radius = 0;
    }

    output = Hysteresis(NMS, low, high, imgrows, imgcols, radius);
  }

  if (mode & MODET)
    i_write_image_file(FP, output, IT_BYTE);
  else
    i_write_image_file(FP, output, IT_DOUBLE);

  fclose(FP);
  free(NMS->ptr);
  free(NMS);
  i_destroy_image(input);
  i_destroy_image(output);
}

