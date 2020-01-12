#define _FILE_OFFSET_BITS  64 // large file support
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mkl.h>
#include <R.h>
#include "slide.h"
#define ONEWRITE 1000
int main (int argc, char** argv)
{
  MKL_Set_Num_Threads(1); // avoid multi-threading
  int ntags, szwin, nrand, *nsec;
  int i,j,l,seed,nloop,nsets,set,loop;
  char *inf, *outf;;
  char buf[1000], *line;
  FILE *fp;
  double r, m, f;
  double *convec, *constd, **scale, **thres, *rnds, *maxnull, *rgen, *rconstd;
  VSLStreamStatePtr stream;
  int ione = 1;

  if (argc != 5) print_usage();
  inf = argv[1];
  outf = argv[2];
  nrand = atoi(argv[3]);
  seed = atoi(argv[4]);

  // count # of SNPs
  ntags = 0;
  line = (char*)malloc(20000001*sizeof(char));
  sprintf(buf,"%s.txt_polySNPs",inf);
  fp = readfile(buf);
  while(fgets(line,20000000,fp)!=NULL) ntags++;
  fclose(fp);
  // read window size
  sprintf(buf,"%s.txt_wsize",inf);
  fp = readfile(buf);
  fscanf(fp, "%d", &szwin);
  fclose(fp);
  // alloc memory
  convec = (double*)malloc(sizeof(double)*ntags*szwin);
  constd = (double*)malloc(sizeof(double)*ntags);
  thres = (double**)malloc(sizeof(double*)*2*ntags);
  scale = (double**)malloc(sizeof(double*)*2*ntags);
  nsec = (int*)malloc(sizeof(int)*2*ntags);
  // read binaries
  sprintf(buf,"%s.bin_convec",inf);
  fp = readfile(buf);
  fread(convec,sizeof(double),ntags*szwin,fp);
  fclose(fp);
  sprintf(buf,"%s.bin_constd",inf);
  fp = readfile(buf);
  fread(constd,sizeof(double),ntags,fp);
  fclose(fp);
  sprintf(buf,"%s.txt_nsec",inf);
  fp = readfile(buf);
  for (i=0;i< 2*ntags;i++) fscanf(fp,"%d", &nsec[i]);
  fclose(fp);
  sprintf(buf,"%s.bin_scale",inf);
  fp = readfile(buf);
  for (i=0;i < 2*ntags;i++) {
    thres[i] = (double*)malloc(sizeof(double)*nsec[i]);
    scale[i] = (double*)malloc(sizeof(double)*nsec[i]);
    fread(&thres[i][0],sizeof(double),nsec[i],fp);
    fread(&scale[i][0],sizeof(double),nsec[i],fp);
  }
  fclose(fp);
  // prepare for sampling
  vslNewStream(&stream, VSL_BRNG_MT19937, seed);
  rnds = (double*)malloc(sizeof(double)*ntags);
  rgen = (double*)malloc(sizeof(double)*ntags);
  rconstd = (double*)malloc(sizeof(double)*ntags);
  maxnull = (double*)malloc(sizeof(double)*ONEWRITE*3);

  fp = writefile(outf);
  nsets = max(1, nrand/ONEWRITE);
  for(set = 0;set < nsets;set++) {
    if (set == nsets-1) {
      nloop = nrand - set*ONEWRITE;
    } else {
      nloop = ONEWRITE;
    }
    if (set > 0) fprintf(stderr, "Sampling %d\n", set*ONEWRITE);
    for(loop=0; loop < nloop; ++loop) {
      m = 0;
      vdRngGaussian(VSL_METHOD_DGAUSSIAN_BOXMULLER, stream, ntags, rgen, 0., 1.);
      vdMul(ntags, rgen, constd, rconstd);
      for(j=0; j < ntags; ++j) {
	r = rconstd[j];
	if ( j < szwin ) {
	  r += ddot(&j,&convec[j*szwin+szwin-j],&ione,&rnds[0],&ione);
	} else {
	  r += ddot(&szwin,&convec[j*szwin],&ione,&rnds[j-szwin],&ione);
	}
	rnds[j] = r;
	f = r;
	if (f <= 0) { // negative case
	  if (nsec[2*j] > 0) {
	    if (f >= thres[2*j][0]) { // if 1st section.
	      f *= scale[2*j][0];
	    } else {
	      for(l = 1;l < nsec[2*j];l++) {
		if (f > thres[2*j][l]) {
		  f *= scale[2*j][l-1] + (scale[2*j][l]-scale[2*j][l-1])/(thres[2*j][l]-thres[2*j][l-1])*(f-thres[2*j][l-1]); // linear interpolation
		  break;
		}
	      }
	      if (l == nsec[2*j]) f *= scale[2*j][l-1]; // over the final section
	    }
	  }
	} else { // positive case
	  if (nsec[2*j+1] > 0) {
	    if (f <= thres[2*j+1][0]) { // if 1st section
	      f *= scale[2*j+1][0];
	    } else {
	      for(l = 1;l < nsec[2*j+1];l++) {
		if (f < thres[2*j+1][l]) {
		  f *= scale[2*j+1][l-1] + (scale[2*j+1][l]-scale[2*j+1][l-1])/(thres[2*j+1][l]-thres[2*j+1][l-1])*(f-thres[2*j+1][l-1]); // linear interpolation
		  break;
		}
	      }
	      if (l == nsec[2*j+1]) f *= scale[2*j+1][l-1]; // over the final section
	    }
	  }
	}
	if (fabs(f) > m ) m = fabs(f);
      } //for(j)
      maxnull[loop] = m;
    } //for(loop)
    fwrite(maxnull,sizeof(double),nloop,fp);
/*     for(i=0;i<nloop;i++)printf("%.10lf\n",maxnull[i]); */
  }//for(set)
  fclose(fp);
  fprintf(stdout, ">>> Sampling result file %s has been successfully created.\n",outf);
  free(convec);
  free(constd);
  free(thres);
  free(scale);
  free(rnds);
  free(rgen);
  free(rconstd);
  free(maxnull);
  return 0;
}

void print_usage() {
  fprintf(stderr,"usage:\n\
   ./slide_2run [prep file] [max stat file] [#sampling] [seed]\n\
   note: [prep file] has to be the same name you specified when running slide_1prep.\n");
  exit(-1);
}

