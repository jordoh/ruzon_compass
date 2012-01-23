/****************************************************************************/
/*                                                                          */
/* Program:   lib.h                                                         */
/* Author:    Simon A.J. Winder                                             */
/* Date:      Tue Feb  2 17:17:40 1993                                      */
/* Function:  Private image library header file                             */
/* Copyright (C) 1994 Simon A.J. Winder                                     */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "image.h"

/* library macros */
#define DEBUG(message) fprintf(stderr,"Debug: " message ". (Ilib)\n")
/*#define MEM_ERROR      i_error("out of memory (ILib)")*/

/* Version 1.0 (Oct 1994) */
