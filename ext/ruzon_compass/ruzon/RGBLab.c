/* RGBLab.c */

/* Convert between RGB and CIE-Lab color spaces
 * Uses ITU-R recommendation BT.709 with D65 as reference white.
 * Yossi Rubner 2/24/98
 *
 * Mark Ruzon 3/18/99 -- Added thresholds to truncate parts of the space
 * where changing a value has no visible effect on the color.
 */

#include <math.h>
#include <stdio.h>

#define BLACK 20.0
#define YELLOW 70.0

void RGB2Lab(float R, float G, float B,
	     float *L, float *a, float *b)
{
  float X, Y, Z, fX, fY, fZ;

  X = 0.412453*R + 0.357580*G + 0.180423*B;
  Y = 0.212671*R + 0.715160*G + 0.072169*B;
  Z = 0.019334*R + 0.119193*G + 0.950227*B;

  X /= (255 * 0.950456);
  Y /=  255;
  Z /= (255 * 1.088754);

  if (Y > 0.008856)
    {
      fY = pow(Y, 1.0/3.0);
      *L = 116.0*fY - 16.0;
    }
  else
    {
      fY = 7.787*Y + 16.0/116.0;
      *L = 903.3*Y;
    }

  if (X > 0.008856)
      fX = pow(X, 1.0/3.0);
  else
      fX = 7.787*X + 16.0/116.0;

  if (Z > 0.008856)
      fZ = pow(Z, 1.0/3.0);
  else
      fZ = 7.787*Z + 16.0/116.0;

  *a = 500.0*(fX - fY);
  *b = 200.0*(fY - fZ);

  if (*L < BLACK) {
    *a *= exp((*L - BLACK) / (BLACK / 4));
    *b *= exp((*L - BLACK) / (BLACK / 4));
    *L = BLACK;
  }
  if (*b > YELLOW)
    *b = YELLOW;

  /*printf("RGB=(%d,%d,%d) ==> Lab(%d,%d,%d)\n",(int)R,(int)G,(int)B,
    (int)*L,(int)*a,(int)*b); */
}


void Lab2RGB(float L, float a, float b, 
	     float *R, float *G, float *B)
{
  float X, Y, Z, fX, fY, fZ;
  float RR, GG, BB;

  fY = pow((L + 16.0) / 116.0, 3.0);
  if (fY < 0.008856)
      fY = L / 903.3;
  Y = fY;

  if (fY > 0.008856)
    fY = pow(fY, 1.0/3.0);
  else
    fY = 7.787 * fY + 16.0/116.0;

  fX = a / 500.0 + fY;      
  if (fX > 0.206893)
      X = pow(fX, 3.0);
  else
      X = (fX - 16.0/116.0) / 7.787;
 
  fZ = fY - b /200.0;      
  if (fZ > 0.206893)
      Z = pow(fZ, 3.0);
  else
      Z = (fZ - 16.0/116.0) / 7.787;

  X *= (0.950456 * 255);
  Y *=             255;
  Z *= (1.088754 * 255);

  RR =  3.240479*X - 1.537150*Y - 0.498535*Z;
  GG = -0.969256*X + 1.875992*Y + 0.041556*Z;
  BB =  0.055648*X - 0.204043*Y + 1.057311*Z;

  *R = (float)(RR < 0 ? 0 : RR > 255 ? 255 : RR);
  *G = (float)(GG < 0 ? 0 : GG > 255 ? 255 : GG);
  *B = (float)(BB < 0 ? 0 : BB > 255 ? 255 : BB);

  /*printf("Lab=(%f,%f,%f) ==> RGB(%f,%f,%f)\n",L,a,b,*R,*G,*B);*/
}
