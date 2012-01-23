/****************************************************************************/
/*                                                                          */
/* Program:   i_main.c                                                      */
/* Author:    Simon A.J. Winder                                             */
/* Date:      Wed Feb 19 21:41:56 1992                                      */
/* Function:  Main image manipulation functions.                            */
/* Copyright (C) 1994 Simon A.J. Winder                                     */
/*                                                                          */
/****************************************************************************/

#include "lib.h"

/*#define debug(m) DEBUG("<main> " m)*/
#define debug(m)

/*****************************************************************************/

void i_error(char *e)
{
  fprintf(stderr,"Error: %s.\n",e);
  exit(1);
}

it_image *i_create_image(int width,int height,int type,int mode)
{
  int size;
  it_image *ptr;

  switch(type)
    {
    case IT_BIT:
      size=sizeof(it_bit);
      if((ptr=i_create_any_image((width+7)>>3,height,size,mode))!=NULL)
	{
	  ptr->type=type;
	  ptr->width=width;
	  ptr->valid_width=width;
	}
      return(ptr);
    case IT_BYTE:
      size=sizeof(it_byte);
      break;
    case IT_LONG:
      size=sizeof(it_long);
      break;
    case IT_DOUBLE:
      size=sizeof(it_double);
      break;
    case IT_FLOAT:
      size=sizeof(it_float);
      break;
    case IT_RGB:
      size=sizeof(it_rgb);
      break;
    case IT_COMPLEX:
      size=sizeof(it_complex);
      break;
    case IT_POLAR:
      size=sizeof(it_polar);
      break;
    default:
      return(NULL);
    }
  if((ptr=i_create_any_image(width,height,size,mode))!=NULL)
    ptr->type=type;
  return(ptr);
}

it_image *i_create_any_image(int width,int height,int size,int mode)
{
  it_image *ptr;
  void **field;

  if((ptr=(it_image *) calloc(1,sizeof(it_image)))==NULL)
    return(NULL);
  if((field=i_create_image_field(ptr,width,height,size,mode))==NULL)
    {
      free((char *) ptr);
      ptr=NULL;
    }
  return(ptr);
}

void i_destroy_image(it_image *ptr)
{
  if(ptr!=NULL)
    {
      i_destroy_image_field(ptr);
      free((char *) ptr);
    }
}

void **i_create_image_field(it_image *ptr,int width,int height,int size,
			    int mode)
{
  void **field;
  void *row;
  int i;

  if(width<1 || height <1) return(NULL);
  if((field=(void **) calloc(height+1,sizeof(void *)))==NULL)
    return(NULL);
  ptr->width=width;
  ptr->height=height;
  ptr->mode=mode;
  ptr->type=IT_NONE;
  ptr->valid_x=0;
  ptr->valid_y=0;
  ptr->valid_width=width;
  ptr->valid_height=height;
  ptr->max_value=0.0;
  ptr->min_value=0.0;
  ptr->field=field;
  if(mode==IM_CONTIG)
    {
      if((row=(void *) calloc(width*height,size))==NULL)
	{
	  free((char *) field);
	  return(NULL);
	}
      for(i=0;i<height;i++) field[i]=(void *)(((char *)row)+i*width*size);
    } else {
      for(i=0;i<height;i++)
	{
	  if((row=(void *) calloc(width,size))==NULL)
	    {
	      i_destroy_image_field(ptr);
	      return(NULL);
	    }
	  field[i]=row;
	}
    }
  return(field);
}

void i_destroy_image_field(it_image *ptr)
{
  void **field;
  int i;

  if(ptr!=NULL)
    {
      field=ptr->field;
      if(field!=NULL)
	{
	  if(ptr->mode==IM_CONTIG)
	    free((char *) field[0]);
	  else
	    for(i=0;field[i]!=NULL;i++)
	      free((char *) field[i]);
	  free((char *) field);
	  ptr->field=NULL;	
	}
    }
}

int i_get_bit_value(it_image *ptr,int x,int y)
{
  return(im_get_bit_value(ptr,x,y));
}

void i_put_bit_value(it_image *ptr,int x,int y,int val)
{
  im_put_bit_value(ptr,x,y,val);
}

void i_set_valid_region(it_image *ptr,int x,int y,int width,int height)
{
  /* Chack all the sizes are valid */
  if(x<0)
    x=0;
  else if(x>=ptr->width)
    x=ptr->width-1;
  if(y<0)
    y=0;
  else if(y>=ptr->height)
    y=ptr->height-1;
  if(x+width>ptr->width)
    width=ptr->width-x;
  else if(width<0)
    width=0;
  if(y+height>ptr->height)
    height=ptr->height-y;
  else if(height<0)
    height=0;

  /* Now set the values */
  ptr->valid_x=x;
  ptr->valid_y=y;
  ptr->valid_width=width;
  ptr->valid_height=height;
}

void i_set_min_max_values(it_image *ptr,double min,double max)
{
  ptr->min_value=min;
  ptr->max_value=max;
}
/* Version 1.0 (Oct 1994) */
