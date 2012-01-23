/* bs.h
 *
 * Mark A. Ruzon
 * Stanford Vision Laboratory
 * 8 September 1998
 *
 * Header file for bs.c
 * 
 */

#define DIM 3
#define MAXCLUSTERS 30

typedef float Coord;

void bs(Coord *, int, int, int *, Coord **, int *);
