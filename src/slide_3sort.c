#define _FILE_OFFSET_BITS  64 // large file support
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mkl.h>
#include <R.h>
#include "slide.h"
#define DBL_CHUNK 1000

int main (int argc, char** argv)
{
  MKL_Set_Num_Threads(1); // avoid multi-threading
  char* orgf;
  int i,j,size,prevsize;
  FILE *fp;
  double chunk[DBL_CHUNK],*maxstats,*stats;

  if (argc < 3) print_usage();
  orgf = argv[1];

  // # sampling consistency check
  prevsize = -1;
  for (i = 2;i < argc;i++) {
    fp = readfile(argv[i]);
    size = 0;
    do {
      j = fread(chunk,sizeof(double),DBL_CHUNK,fp);
      size += j;
    } while (j == DBL_CHUNK);
    if (prevsize!=-1 && prevsize!=size) {
      fprintf(stderr,"Number of sampling has to be the same for all result files\n");exit(-1);
    }
    prevsize = size;
    fclose(fp);
  }

  // read result files
  maxstats = (double*)malloc(sizeof(double)*size);
  stats = (double*)malloc(sizeof(double)*size);
  for (j = 0;j < size;j++) maxstats[j] = 0.;
  for (i = 2;i < argc;i++) {
    fp = readfile(argv[i]);
    fread(stats,sizeof(double),size,fp);
    fclose(fp);
    for (j = 0;j < size;j++) maxstats[j] = max(maxstats[j], stats[j]);
  }

  // sort
  qsort(maxstats, size, sizeof(double), compare_double);

  // write
  fp = writefile(orgf);
  fwrite(maxstats,sizeof(double),size,fp);
  fclose(fp);

  return 0;
}

void print_usage() {
  fprintf(stderr,"usage:\n\
   ./slide_3sort [sorted file] [max stat file (chr1)] [max stat file (chr2)] ...\n");
  exit(-1);
}
int compare_double_rev( const void *s, const void *t) {
  double a = *((double*)s);
  double b = *((double*)t);
  return (a > b)? -1: a < b? 1:0;
}
int compare_double( const void *s, const void *t) {
  double a = *((double*)s);
  double b = *((double*)t);
  return (a > b)? 1: a < b? -1:0;
}
