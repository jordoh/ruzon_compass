/****************************************************************************/
/*                                                                          */
/* Program:   image.h                                                       */
/* Author:    Simon A.J. Winder                                             */
/* Date:      Wed Feb 19 21:42:29 1992                                      */
/* Function:  Image Library Public Header File                              */
/* Copyright (C) 1994 Simon A.J. Winder                                     */
/*                                                                          */
/****************************************************************************/

typedef unsigned char byte;

/* supported image field value types */
typedef byte it_bit;
typedef byte it_byte;
typedef long it_long;
typedef float it_float;
typedef double it_double;
typedef struct {byte r,g,b;} it_rgb;
typedef struct {float Re,Im;} it_complex;
typedef struct {float rad,ang;} it_polar;

/* general image structure */
typedef struct {
  int width,height;
  int valid_x,valid_y;
  int valid_width,valid_height;
  double min_value,max_value;
  int mode,type;
  void **field;
} it_image;

/* mode of image memory allocation */
#define IM_CONTIG   0
#define IM_FRAGMENT 1

/* types of supported images */
#define IT_NONE     0
#define IT_BIT      (1<<0)
#define IT_BYTE     (1<<1)
#define IT_LONG     (1<<2)
#define IT_FLOAT    (1<<7)
#define IT_DOUBLE   (1<<3)
#define IT_RGB      (1<<4)
#define IT_COMPLEX  (1<<5)
#define IT_POLAR    (1<<6)
#define IT_ANY      (~0)

/* xy coordinate structure */
typedef struct {
  int x;
  int y;
} it_point;

/* main image function prototypes */
void      i_error(char *);
it_image  *i_create_image(int,int,int,int);
it_image  *i_create_any_image(int,int,int,int);
void      i_destroy_image(it_image *);
void      **i_create_image_field(it_image *,int,int,int,int);
void      i_destroy_image_field(it_image *);
int       i_get_bit_value(it_image *,int,int);
void      i_put_bit_value(it_image *,int,int,int);
void      i_set_valid_region(it_image *,int,int,int,int);
void      i_set_min_max_values(it_image *,double,double);

/* image macros */
#define im_get_bit_value(p,x,y)                                   \
  (((((it_bit **)((p)->field))[y][(x)>>3])&(0x80>>((x)&7)))?1:0)
#define im_put_bit_value(p,x,y,v)                                 \
  ( (v)                                                           \
  ?((((it_bit **)((p)->field))[y][(x)>>3])|=(0x80>>((x)&7)))      \
  :((((it_bit **)((p)->field))[y][(x)>>3])&=(~(0x80>>((x)&7)))) )

#define im_byte_value(p,x,y)    (((it_byte **)((p)->field))[y][x])
#define im_long_value(p,x,y)    (((it_long **)((p)->field))[y][x])
#define im_float_value(p,x,y)   (((it_float **)((p)->field))[y][x])
#define im_double_value(p,x,y)  (((it_double **)((p)->field))[y][x])
#define im_rgb_value(p,x,y)     (((it_rgb **)((p)->field))[y][x])
#define im_complex_value(p,x,y) (((it_complex **)((p)->field))[y][x])
#define im_polar_value(p,x,y)   (((it_polar **)((p)->field))[y][x])

#define im_bit_row(p,y)         (((it_bit **)((p)->field))[y])
#define im_byte_row(p,y)        (((it_byte **)((p)->field))[y])
#define im_long_row(p,y)        (((it_long **)((p)->field))[y])
#define im_float_row(p,y)       (((it_float **)((p)->field))[y])
#define im_double_row(p,y)      (((it_double **)((p)->field))[y])
#define im_rgb_row(p,y)         (((it_rgb **)((p)->field))[y])
#define im_complex_row(p,y)     (((it_complex **)((p)->field))[y])
#define im_polar_row(p,y)       (((it_polar **)((p)->field))[y])

#define im_bit_field(p)         ((it_bit **)((p)->field))
#define im_byte_field(p)        ((it_byte **)((p)->field))
#define im_long_field(p)        ((it_long **)((p)->field))
#define im_float_field(p)       ((it_float **)((p)->field))
#define im_double_field(p)      ((it_double **)((p)->field))
#define im_rgb_field(p)         ((it_rgb **)((p)->field))
#define im_complex_field(p)     ((it_complex **)((p)->field))
#define im_polar_field(p)       ((it_polar **)((p)->field))

#define im_rgb_equal(p,q)    ((p)->r==(q)->r&&(p)->g==(q)->g&&(p)->b==(q)->b)
#define im_set_rgb(p,rr,gg,bb)  {(p)->r=(rr);(p)->g=(gg);(p)->b=(bb);}
#define im_luminance(p)      ((p)->r*0.299+(p)->g*0.587+(p)->b*0.114)

/* HISTOGRAM */
/* histogram entry structure */
typedef struct _it_hist_entry {
  it_rgb colour;
  int index;
  long count;
  struct _it_hist_entry *next_col;
} it_hist_entry;

#define IM_HIST_HASHSIZE 101

/* histogram root structure */
typedef struct {
  int maxcolours;
  int numcolours;
  it_hist_entry *hash_col[IM_HIST_HASHSIZE];
  it_hist_entry *hist_list;
} it_hist;
  
/* histogram function prototypes */
it_hist        *i_create_histogram(int);
void           i_destroy_histogram(it_hist *);
it_hist_entry  *i_add_histogram_entry(it_hist *,it_rgb *);
it_hist_entry  *i_lookup_histogram_index(it_hist *,int);
it_hist_entry  *i_lookup_histogram_colour(it_hist *,it_rgb *);

/* histogram macro */
#define im_histogram_colour(p)  (&((p)->colour))

/* COLOURMAP */
/* colourmap entry structure */
typedef struct _it_cmap_entry {
  it_rgb colour;
  struct _it_cmap_entry *next_col;
} it_cmap_entry;

#define IM_CMAP_HASHSIZE 101

/* colourmap root structure */
typedef struct {
  int maxindex;
  it_cmap_entry *hash_col[IM_CMAP_HASHSIZE];
  it_cmap_entry *cmap_list;
} it_cmap;

/* colormap function prototypes */
it_cmap     *i_create_colourmap(int);
void        i_destroy_colourmap(it_cmap *);
void        i_set_colourmap_entry(it_cmap *,int,it_rgb *);
int         i_lookup_colourmap_colour(it_cmap *,it_rgb *);
it_rgb      *i_lookup_colourmap_index(it_cmap *,int);

/* IMAGE FILES */
/* binary and ascii flags */
#define IF_ASCII  0
#define IF_BINARY 1

#define     im_write_separator(fp) /* was: fprintf((fp),"*\n")*/

/* image file function prototypes */
void        i_write_image_file(FILE *,it_image *,int);
it_image    *i_read_image_file(FILE *,int,int);
int         i_save_image(char *,int,it_image *,int);
it_image    *i_load_image(char *,int,int,int);
FILE        *i_open_file(char *,int,char *);
char        *i_parse_filename(char *,int);

/* image type conversion function prototypes */
int         i_byte_to_float(it_image *,it_image *);
int         i_float_to_byte(it_image *,it_image *);
int         i_byte_to_double(it_image *,it_image *);
int         i_double_to_byte(it_image *,it_image *);
int         i_rgb_to_float(it_image *,it_image *,it_image *,it_image *);

/* fourier transform types */
#define FFFT  1
#define IFFT -1

/* image tool function prototypes */
it_image    *i_create_operator(int,int,double *);
int         i_smooth(it_image *,it_image *,double);
int         i_convolve(it_image *,it_image *,it_image *);
int         i_fourier_transform_1d(it_image *,int);
int         i_fourier_transform_2d(it_image *,int);
/* Version 1.0 (Oct 1994) */
