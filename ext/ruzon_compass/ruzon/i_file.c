/****************************************************************************/
/*                                                                          */
/* Program:   i_file.c                                                      */
/* Author:    Simon A.J. Winder                                             */
/* Date:      Wed Feb 26 11:52:12 1992                                      */
/* Function:  Image File Load and Save Routines                             */
/* Copyright (C) 1994 Simon A.J. Winder                                     */
/*                                                                          */
/****************************************************************************/

#include "lib.h"
#include <pwd.h>

/*#define debug(m) DEBUG("<file> " m)*/
#define debug(m)

void ii_write_size(FILE *,it_image *);
void ii_write_header(FILE *,it_image *);
void ii_read_header(FILE *,it_image *);
int ii_read_char(FILE *);
char *ii_read_field(FILE *);
double ii_get_scale(FILE *);
void ii_write_block(FILE *,it_image *,int);
int ii_read_block(FILE *,it_image *,int);
char *ii_parse_filename_pwd(char *);
char *ii_parse_filename_ext(char *file);

/****************************************************************************/

it_image *i_load_image(char *file,int seq_num,int type,int mode)
{
  FILE *fp;
  it_image *image;

  /* Open the file */
  if((fp=i_open_file(file,seq_num,"r"))==NULL)
    return(NULL);

  /* Load the image */
  image=i_read_image_file(fp,type,mode);
  fclose(fp);
  return(image);
}

int i_save_image(char *file,int seq_num,it_image *image,int mode)
{
  FILE *fp;
  int status=0;

  /* Open the file */
  if((fp=i_open_file(file,seq_num,"w"))==NULL)
    return(1);

  /* Save the image */
  i_write_image_file(fp,image,mode);
  if(ferror(fp))
    status=1;
  fclose(fp);
  return(status);
}

FILE *i_open_file(char *file,int seq_num,char *mode)
{
  char *name;
  FILE *fp;

  name=i_parse_filename(file,seq_num);
  if(name==NULL)
    return(NULL);
  if((fp=fopen(name,mode))==NULL)
    {
      free(name);
      return(NULL);
    }
  free(name);
  return(fp);
}

char *i_parse_filename(char *file,int seq_num)
{
  char seq_string[4],*ptr1,*ptr2,*ext,*name;

  /* First expand out any directories */
  name=ii_parse_filename_pwd(file);
  if(seq_num>=0)
    {
      /* Add sequence number to filename before extension */

      /* Find the extension part (if any) */
      ext=ii_parse_filename_ext(name);
      sprintf(seq_string,"%.3d",seq_num%1000);
      file=(char *) malloc((strlen(name)+4)*sizeof(char));
      if(file==NULL)
	return(NULL);

      /* Copy name into file inserting sequence number before ext */
      ptr1=file;
      for(ptr2=name;ptr2!=ext;ptr2++)
	*ptr1++ = *ptr2;
      strcpy(ptr1,seq_string);
      strcat(file,ext);
      free(name);
      name=file;
    }
  return(name);
}

void i_write_image_file(FILE *fp,it_image *image,int mode)
{
  int x,y,width,height;
  int count=0;
  long out;
  it_rgb *rgbptr;

  if(image==NULL) return;
  width=image->width;
  height=image->height;
  switch(image->type)
    {
    case IT_BIT:
      if(mode==IF_ASCII)
	{
	  fprintf(fp,"P1\n");
	  ii_write_size(fp,image);
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      {
		if(im_get_bit_value(image,x,y))
		  fputc('1',fp);
		else
		  fputc('0',fp);
		if(++count>68)
		  {
		    fputc('\n',fp);
		    count=0;
		  }
	      }
	  if(count) fputc('\n',fp);
	}
      else
	{
	  fprintf(fp,"P4\n");
	  ii_write_size(fp,image);
	  for(y=0;y<height;y++)
	    {
	      count=0;
	      out=0L;
	      for(x=0;x<width;x++)
		{
		  out<<=1;
		  out|=im_get_bit_value(image,x,y);
		  if(++count==8)
		    {
		      fputc((char) out,fp);
		      count=0;
		      out=0L;
		    }
		}
	      if(count) fputc((char) (out<<(8-count)),fp);
	    }
	}
      break;

    case IT_BYTE:
      if(mode==IF_ASCII)
	{
	  fprintf(fp,"P2\n");
	  ii_write_size(fp,image);
	  fprintf(fp,"255\n");
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      {
		if(count) fputc(' ',fp);
		fprintf(fp,"%d",im_byte_value(image,x,y));
		if(++count==16)
		  {
		    count=0;
		    fputc('\n',fp);
		  }
	      }
	  if(count) fputc('\n',fp);
	}
      else
	{
	  fprintf(fp,"P5\n");
	  ii_write_size(fp,image);
	  fprintf(fp,"255\n");
	  ii_write_block(fp,image,sizeof(it_byte));
	}
      break;

    case IT_LONG:
      if(mode==IF_ASCII)
	{
	  fprintf(fp,"ILibLong\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      {
		if(count) fputc(' ',fp);
		fprintf(fp,"%ld",im_long_value(image,x,y));
		if(++count==5)
		  {
		    count=0;
		    fputc('\n',fp);
		  }
	      }
	  if(count) fputc('\n',fp);
	}
      else
	{
	  fprintf(fp,"ILibRawLong\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  ii_write_block(fp,image,sizeof(it_long));
	}
      break;

    case IT_FLOAT:
      if(mode==IF_ASCII)
	{
	  fprintf(fp,"ILibFloat\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      fprintf(fp,"%.8g\n",im_float_value(image,x,y));
	}
      else
	{
	  fprintf(fp,"ILibRawFloat\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  ii_write_block(fp,image,sizeof(it_float));
	}
      break;

    case IT_DOUBLE:
      if(mode==IF_ASCII)
	{
	  fprintf(fp,"ILibDouble\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      fprintf(fp,"%.12g\n",im_double_value(image,x,y));
	}
      else
	{
	  fprintf(fp,"ILibRawDouble\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  ii_write_block(fp,image,sizeof(it_double));
	}
      break;

    case IT_RGB:
      if(mode==IF_ASCII)
	{
	  fprintf(fp,"P3\n");
	  ii_write_size(fp,image);
	  fprintf(fp,"255\n");
	  for(y=0;y<height;y++)
	    {
	      rgbptr=im_rgb_row(image,y);
	      for(x=0;x<width;x++)
		{
		  if(count) fputc(' ',fp);
		  fprintf(fp,"%d ",rgbptr->r);
		  fprintf(fp,"%d ",rgbptr->g);
		  fprintf(fp,"%d",rgbptr->b);
		  rgbptr++;
		  if(++count==5)
		    {
		      count=0;
		      fputc('\n',fp);
		    }
		}
	    }
	  if(count) fputc('\n',fp);
	}
      else
	{
	  fprintf(fp,"P6\n");
	  ii_write_size(fp,image);
	  fprintf(fp,"255\n");

	  /* Note that we can't use ii_write_block here because it will */
	  /* not write the structures in an architecture independent way */
	  for(y=0;y<height;y++)
	    {
	      rgbptr=im_rgb_row(image,y);
	      for(x=0;x<width;x++)
		{
		  fputc(rgbptr->r,fp);
		  fputc(rgbptr->g,fp);
		  fputc(rgbptr->b,fp);
		  rgbptr++;
		}
	    }
	}
      break;

    case IT_COMPLEX:
      if(mode==IF_ASCII)
	{
	  fprintf(fp,"ILibComplex\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      {
		fprintf(fp,"%.8g",im_complex_value(image,x,y).Re);
		fprintf(fp," %.8g\n",im_complex_value(image,x,y).Im);
	      }
	}
      else
	{
	  fprintf(fp,"ILibRawComplex\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  ii_write_block(fp,image,sizeof(it_complex));
	}
      break;

    case IT_POLAR:
      if(mode==IF_ASCII)
	{
	  fprintf(fp,"ILibPolar\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      {
		fprintf(fp,"%.8g",im_polar_value(image,x,y).rad);
		fprintf(fp," %.8g\n",im_polar_value(image,x,y).ang);
	      }
	}
      else
	{
	  fprintf(fp,"ILibRawPolar\n");
	  ii_write_size(fp,image);
	  ii_write_header(fp,image);
	  ii_write_block(fp,image,sizeof(it_polar));
	}
      break;

    default:
      fprintf(stderr,"Fault: unknown image type (i_write_image_file).\n");
    }
}

it_image *i_read_image_file(FILE *fp,int types,int mode)
{
  static char *magic[]={
    "P1","P4",
    "P2","P5",
    "ILibLong","ILibRawLong",
    "ILibFloat","ILibRawFloat",
    "ILibDouble","ILibRawDouble",
    "P3","P6",
    "ILibComplex","ILibRawComplex",
    "ILibPolar","ILibRawPolar",
  };
  static int mask[]={
    IT_BIT,IT_BIT,
    IT_BYTE,IT_BYTE,
    IT_LONG,IT_LONG,
    IT_FLOAT,IT_FLOAT,
    IT_DOUBLE,IT_DOUBLE,
    IT_RGB,IT_RGB,
    IT_COMPLEX,IT_COMPLEX,
    IT_POLAR,IT_POLAR,
  };
  it_image *ptr;
  int type,width,height,x,y,c,i;
  char *p;
  double scale;
  it_rgb *rgbptr;

  /* Make gcc happy */
  ptr=NULL;

  /* find out what the magic string is and compare */
  if((p=ii_read_field(fp))==NULL) return(NULL);
  for(type=0;type<16;type++)
    if(strcmp(p,magic[type])==0) break;

  /* is it the one we wanted ? */
  if(type==16) return(NULL);
  if(!(mask[type]&types)) return(NULL);
  debug("valid image type");

  /* find out the size */
  if((p=ii_read_field(fp))==NULL) return(NULL);

  width=atoi(p);
  if((p=ii_read_field(fp))==NULL) return(NULL);
  height=atoi(p);
  if(width<1 || height<1) return(NULL);
  debug("vaild image size");
  /*fprintf(stderr,"Debug: (%dx%d)\n",width,height);*/


  /* load the rest of the image */
  switch(type)
    {
    case 0:  /* ascii it_bit (pbm)      */
    case 1:  /* binary it_bit (pbmraw)  */
      if((ptr=i_create_image(width,height,IT_BIT,mode))==NULL)
	return(NULL);
      if(type==0)
	{
	  debug("loading it_bit (ascii)");
	  for(x=y=0;y<height && (p=ii_read_field(fp))!=NULL;)
	    for(;*p=='1' || *p=='0';p++)
	      {
		im_put_bit_value(ptr,x,y,*p-'0');
		if(++x==width)
		  {
		    x=0;
		    if(++y==height) break;
		  }
	      }
	} else {
	  debug("loading it_bit (binary)");
	  for(x=y=0;y<height && (c=fgetc(fp))!=EOF;)
	    for(i=7;i>=0;i--)
	      {
		im_put_bit_value(ptr,x,y,(c>>i)&1);
		if(++x==width)
		  {
		    x=0;
		    y++;
		    break;
		  }
	      }
	}
      break;
      
    case 2:  /* ascii it_byte (pgm)     */
    case 3:  /* binary it_byte (pgmraw) */
      if((scale=ii_get_scale(fp))==0.0) return(NULL);
      if((ptr=i_create_image(width,height,IT_BYTE,mode))==NULL)
	return(NULL);
      if(type==2)
	{
	  debug("loading it_byte (ascii)");
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      if((p=ii_read_field(fp))!=NULL)
		im_byte_value(ptr,x,y)=(it_byte) (atof(p)*scale);
	} else {
	  debug("loading it_byte (binary)");
	  if(scale!=1.0)
	    {
	      for(y=0;y<height;y++)
		for(x=0;x<width;x++)
		  if((c=fgetc(fp))!=EOF)
		    im_byte_value(ptr,x,y)=(it_byte) (((double) c)*scale);
	    } else {
	      if(ii_read_block(fp,ptr,sizeof(it_byte)))
		return(NULL);
	    }
	}
      break;
      
    case 4:  /* ascii it_long           */
    case 5:  /* binary it_long          */
      if((ptr=i_create_image(width,height,IT_LONG,mode))==NULL)
	return(NULL);
      ii_read_header(fp,ptr);
      if(type==4)
	{
	  debug("loading it_long (ascii)");
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      if((p=ii_read_field(fp))!=NULL)
		im_long_value(ptr,x,y)=atol(p);
	} else {
	  debug("loading it_long (binary)");
	  if(ii_read_block(fp,ptr,sizeof(it_long)))
	    return(NULL);
	}
      break;

    case 6:  /* ascii it_float          */
    case 7:  /* binary it_float         */
      if((ptr=i_create_image(width,height,IT_FLOAT,mode))==NULL)
	return(NULL);
      ii_read_header(fp,ptr);
      if(type==6)
	{
	  debug("loading it_float (ascii)");
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      if((p=ii_read_field(fp))!=NULL)
		im_float_value(ptr,x,y)=(it_float) atof(p);
	} else {
	  debug("loading it_float (binary)");
	  if(ii_read_block(fp,ptr,sizeof(it_float)))
	    return(NULL);
	}
      break;
      
    case 8:  /* ascii it_double         */
    case 9:  /* binary it_double        */
      if((ptr=i_create_image(width,height,IT_DOUBLE,mode))==NULL)
	return(NULL);
      ii_read_header(fp,ptr);
      if(type==8)
	{
	  debug("loading it_double (ascii)");
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      if((p=ii_read_field(fp))!=NULL)
		im_double_value(ptr,x,y)=atof(p);
	} else {
	  debug("loading it_double (binary)");
	  if(ii_read_block(fp,ptr,sizeof(it_double)))
	    return(NULL);
	}
      break;
      
    case 10:  /* ascii it_rgb (ppm)      */
    case 11:  /* binary it_rgb (ppmraw)  */
      if((scale=ii_get_scale(fp))==0.0) return(NULL);
      if((ptr=i_create_image(width,height,IT_RGB,mode))==NULL)
	return(NULL);
      if(type==10)
	{
	  debug("loading it_rgb (ascii)");
	  for(y=0;y<height;y++)
	    {
	      rgbptr=im_rgb_row(ptr,y);
	      for(x=0;x<width;x++)
		{
		  if((p=ii_read_field(fp))!=NULL)
		    rgbptr->r=(byte) (atof(p)*scale);
		  if((p=ii_read_field(fp))!=NULL)
		    rgbptr->g=(byte) (atof(p)*scale);
		  if((p=ii_read_field(fp))!=NULL)
		    rgbptr->b=(byte) (atof(p)*scale);
		  rgbptr++;
		}
	    }
	} else {
	  debug("loading it_rgb (binary)");
	  for(y=0;y<height;y++)
	    {
	      rgbptr=im_rgb_row(ptr,y);
	      for(x=0;x<width;x++)
		if((c=fgetc(fp))!=EOF)
		  {
		    rgbptr->r=(byte) (((double) c)*scale);
		    rgbptr->g=(byte) (((double) fgetc(fp))*scale);
		    rgbptr->b=(byte) (((double) fgetc(fp))*scale);
		    rgbptr++;
		  }
	    }
	}
      break;
      
    case 12:  /* ascii it_complex        */
    case 13:  /* binary it_complex       */
      if((ptr=i_create_image(width,height,IT_COMPLEX,mode))==NULL)
	return(NULL);
      ii_read_header(fp,ptr);
      if(type==12)
	{
	  debug("loading it_complex (ascii)");
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      {
		if((p=ii_read_field(fp))!=NULL)
		  im_complex_value(ptr,x,y).Re=atof(p);
		if((p=ii_read_field(fp))!=NULL)
		  im_complex_value(ptr,x,y).Im=atof(p);
	      }
	} else {
	  debug("loading it_complex (binary)");
	  if(ii_read_block(fp,ptr,sizeof(it_complex)))
	    return(NULL);
	}
      break;
      
    case 14: /* ascii it_polar          */
    case 15: /* binary it_polar         */
      if((ptr=i_create_image(width,height,IT_POLAR,mode))==NULL)
	return(NULL);
      ii_read_header(fp,ptr);
      if(type==14)
	{
	  debug("loading it_polar (ascii)");
	  for(y=0;y<height;y++)
	    for(x=0;x<width;x++)
	      {
		if((p=ii_read_field(fp))!=NULL)
		  im_polar_value(ptr,x,y).rad=atof(p);
		if((p=ii_read_field(fp))!=NULL)
		  im_polar_value(ptr,x,y).ang=atof(p);
	      }
	} else {
	  debug("loading it_polar (binary)");
	  if(ii_read_block(fp,ptr,sizeof(it_polar)))
	    return(NULL);
	}
      break;
      
    default:
      /* never reached */
      break;
    }
  debug("waiting for end");
  while((c=fgetc(fp))!=EOF && (c=='*'||isspace(c)));
  if(c!=EOF)
    ungetc(c,fp);
  return(ptr);
}

void ii_write_size(FILE *fp,it_image *ptr)
{
  fprintf(fp,"# Image Library S.A.J. Winder (c)1992\n");
  fprintf(fp,"%d %d\n",ptr->width,ptr->height);
}

void ii_write_header(FILE *fp,it_image *ptr)
{
  fprintf(fp,"%d %d %d %d\n%.12g %.12g\n",ptr->valid_x,ptr->valid_y,
	  ptr->valid_width,ptr->valid_height,ptr->min_value,ptr->max_value);
}

void ii_read_header(FILE *fp,it_image *ptr)
{
  char *p;

  if((p=ii_read_field(fp))!=NULL)
    ptr->valid_x=atoi(p);
  if((p=ii_read_field(fp))!=NULL)
    ptr->valid_y=atoi(p);
  if((p=ii_read_field(fp))!=NULL)
    ptr->valid_width=atoi(p);
  if((p=ii_read_field(fp))!=NULL)
    ptr->valid_height=atoi(p);
  if((p=ii_read_field(fp))!=NULL)
    ptr->min_value=atof(p);
  if((p=ii_read_field(fp))!=NULL)
    ptr->max_value=atof(p);
}

int ii_read_char(FILE *fp)
{
  int c;
  
  if((c=fgetc(fp))==EOF)
    return(c);
  if(c=='#')
    {
      do {
	if((c=fgetc(fp))==EOF)
	  return(c);
      } while(c!='\n'&&c!='\r');
    }

  return(c);
}

char *ii_read_field(FILE *fp)
{
  static char buf[71];
  int c,i;

  do {
    if((c=ii_read_char(fp))==EOF)
      return(NULL);
  } while(isspace(c));
  for(i=0;i<70&&!isspace(c);i++)
    {
      buf[i]=c;
      if((c=ii_read_char(fp))==EOF)
	break;
    }
  buf[i]='\0';
  return(buf);
}

double ii_get_scale(FILE *fp)
{
  char *p;
  double scale=0.0;

  if((p=ii_read_field(fp))==NULL) return(scale);
  if((scale=(double) atol(p))<=0.0) return(0.0);
  return(255.0/scale);
}

void ii_write_block(FILE *fp,it_image *image,int size)
{
  char *ptr;
  int y,height;

  height=image->height;
  for(y=0;y<height;y++)
    {
      ptr=((char **)(image->field))[y];
      fwrite(ptr,size,image->width,fp);
    }
}

int ii_read_block(FILE *fp,it_image *image,int size)
{
  char *ptr;
  int y,height,flag;

  flag=0;
  height=image->height;
  for(y=0;y<height && flag!=1;y++)
    {
      ptr=((char **)(image->field))[y];
      if(fread(ptr,size,image->width,fp)!=image->width)
	flag=1;
    }
  return(flag);
}

char *ii_parse_filename_pwd(char *file)
{
  char *ptr,*file_copy,*dir,user[100];
  int n;
  struct passwd *pwd_ptr;

  if(*file=='~')
    {
      for(n=0,ptr=file+1;*ptr!='/' && *ptr!='\0';ptr++)
	if(n<99)
	  user[n++]=*ptr;
      user[n]='\0';
      if(n==0)
	{
	  /* Must be just ~ or ~/whatever, so get HOME */
	  dir=getenv("HOME");
	  if(dir==NULL)
	    dir=".";
	} else {
	  /* Must be username so look it up in passwd file */
	  pwd_ptr=getpwnam(user);
	  if(pwd_ptr==NULL)
	    dir=user;
	  else
	    dir=pwd_ptr->pw_dir;
	}
      /* Add path to rest of filename */
      file_copy=(char *) malloc((strlen(dir)+strlen(ptr)+1)*sizeof(char));
      if(file_copy!=NULL)
	{
	  strcpy(file_copy,dir);
	  strcat(file_copy,ptr);
	}
    } else {
      /* Make a copy we can free */
      file_copy=strdup(file);
    }    
  /* Could return NULL */
  return(file_copy);
}

char *ii_parse_filename_ext(char *file)
{
  char *end,*dot;
  int cbd,cad,state;

  /* End of filename pointer (default return value) */
  end=file+strlen(file);

  /* Pointer gets set to last dot */
  /* cbd---chars before dot, cad---chars after dot */
  dot=NULL;
  cbd=cad=state=0;

  for(;*file!='\0';file++)
    {
      if(*file=='/')
	/* Reset on a slash */
	cbd=cad=state=0;
      else
	switch(state)
	  {
	  case 0:
	    /* Before a dot */
	    if(*file=='.')
	      {
		/* Switch state if dot found */
		dot=file;
		state=1;
	      }
	    else
	      /* Increment before dot count */
	      cbd++;
	    break;
	  case 1:
	  default:
	    /* After a dot */
	    if(*file=='.')
	      {
		/* Adjust things if a new dot found */
		dot=file;
		cbd+=1+cad;
		cad=0;
	      }
	    else
	      /* Increment after dot count */
	      cad++;
	    break;
	  }
    }
  if(cbd==0 || cad==0 || dot==NULL) dot=end;
  return(dot);
}
/* Version 1.0 (Oct 1994) */
