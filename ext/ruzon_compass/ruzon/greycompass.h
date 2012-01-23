#define MAXWEDGES 30
#define MAXRESPONSES 2

typedef float ImgElt;
typedef struct mymatrix {
  int rows;               /* Number of rows in each matrix */
  int cols;               /* Number of columns in each matrix */
  int sheets;             /* Number of sheets in each matrix (tensor) */
  double *ptr;            /* Pointer to each matrix */
} Matrix;

enum imgtype {RGBImg, LabImg};

void GreyCompass(void *, int, int, enum imgtype, double *, int, int, double *, 
		 double *, double *, int, int, Matrix, Matrix, Matrix, Matrix);

void Error(char *); /* To be defined by user */
