#define MAXWEDGES 30
#define MAXRESPONSES 3

typedef float ImgElt;

enum imgtype {RGBImg, LabImg};

void Compass(void *, int, int, enum imgtype, double *, int, int, double *, 
	     double *, double *, int, int, double *, int, double **, 
	     double **, double **, double **);
