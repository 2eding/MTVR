#define _FILE_OFFSET_BITS  64 // large file support
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mkl.h>
#include "slide.h"
#include <R.h>
#define DBL_CHUNK 1000
#define EPS 1E-100

int main (int argc, char** argv)
{
  MKL_Set_Num_Threads(1); // avoid multi-threading
  char *orgf, *outf, *mapf, *dataf;
  int i,j,size,nrow,ndata;
  FILE *fp;
  double chunk[DBL_CHUNK],*maxstats,x;
  double *data, *result1, *result2, dat, resul1, resul2, neff;
  char *name, *line, *token;

  enum MODE {P, T, IP, IT} mode;
  // receive arguments
  if (argc < 3) print_usage();
  if (!strcmp(argv[1], "-p")) { mode = P; }
  else if (!strcmp(argv[1], "-t")) { mode = T; }
  else if (!strcmp(argv[1], "-ip")) { mode = IP; }
  else if (!strcmp(argv[1], "-it")) { mode = IT; }
  else { print_usage(); }
  orgf = argv[2];
  // read orgf
  fp = readfile(orgf);
  size = 0;
  do {
    j = fread(chunk,sizeof(double),DBL_CHUNK,fp);
    size += j;
  } while (j == DBL_CHUNK);
  fclose(fp);
  if (size == 0) {fprintf(stderr,"orgainzed file empty\n");exit(-1);}
  maxstats = (double*)malloc(sizeof(double)*size);
  fp = readfile(orgf);
  fread(maxstats,sizeof(double),size,fp);
  fclose(fp);

  // by modes
  if (mode == P || mode == T) {// batch modes
    if (argc < 5) print_usage();
    dataf = argv[3];
    // read data
    fp = readfile(dataf);
    ndata = 0;
    while(fscanf(fp,"%lf",&x) != EOF) ndata++;
    fclose(fp);
    data = (double*)malloc(sizeof(double)*ndata);
    result1 = (double*)malloc(sizeof(double)*ndata);
    result2 = (double*)malloc(sizeof(double)*ndata);
    fp = readfile(dataf);
    for (i = 0;i < ndata;++i) fscanf(fp,"%lf",&data[i]);
    fclose(fp);
    outf = argv[4];
    //    for (i = 0;i < 10;i++) fprintf(stderr,"%.10lf\n",data[i]);
    // map file
    if (mode == P && argc == 6) {
      mapf = argv[5];
      name = (char*)malloc(sizeof(char)*ndata*100);
      line = (char*)malloc(20000001*sizeof(char));
      nrow = 0;
      fp = readfile(mapf);
      while(fgets(line,20000000,fp)!=NULL){
	nrow++;
	token = strtok(line," \t\n");
	if(token==NULL || token[0]=='#') nrow--;
      }
      fclose(fp);
      if (nrow != ndata) {fprintf(stderr,"# of pointwise p-values != # of SNPs in map file\n");exit(-1);}
      fp = readfile(mapf);
      for (i = 0;i < nrow;i++) {
	do {
	  fgets(line,20000000,fp);
	  token = strtok(line," \t\n");
	} while (token==NULL || token[0]=='#');
	sscanf(line,"%s",&name[100*i]);
	if (sscanf(line,"%s",&name[100*i]) == EOF)
	  {fprintf(stderr,"map file has to have second column (rsid)\n");exit(-1);}
      }
      fclose(fp);
    }
    fp = writefile(outf);
    fprintf(fp, "#Result based on %d sampling\n", size);
    if (mode == P) {
      for (i = 0;i < ndata;i++) correct_pvalue(maxstats,size,data[i],&result1[i],&result2[i]);
      fprintf(fp, "#SNP_id\tPointwise-P\tCorrected-P\tApprox.STDEV(%% of P)\tEffective#ofTests\tNote\n");
    } else if (mode == T) {
      for (i = 0;i < ndata;i++) per_marker_threshold(maxstats,size,data[i],&result1[i],&result2[i]);
      fprintf(fp, "#Signif-threshold\tPer-marker-thres\tApprox.STDEV(%% of P)\tEffective#ofTests\tNote\n");
    }
    for (i = 0;i < ndata;i++) {
      if (mode == P) {
	if (argc == 6) { fprintf(fp, "%s", &name[100*i]); }
	else { fprintf(fp, "SNP%d", i+1); }
      }
      fprintf(fp,"\t%.8le",data[i]);
      fprintf(fp,"\t%.8le",result1[i]);
      fprintf(fp,"\t%.8le",result2[i]);
      neff = (mode == P)?result1[i]/data[i]:data[i]/result1[i];
      if (result1[i] < EPS || result1[i] >= 1.-EPS) {
	fprintf(fp, "(%.5lf%%)", 0.);
	fprintf(fp, "\t----");
      } else {
	fprintf(fp, "(%.5lf%%)", result2[i]/result1[i]*100);
	fprintf(fp, "\t%.0lf",neff);
      }
      if (result1[i] < EPS) {
	fprintf(fp, "\tNeed more sampling");
      } else if (result1[i] >= 1.-EPS) {
	fprintf(fp, "\t");
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  } else { // interactive modes. (IP & IT)
    if (argc != 4) print_usage();
    dat = atof(argv[3]);
    fprintf(stdout, "#Result based on %d sampling\n", size);
    if (mode == IP) {
      correct_pvalue(maxstats,size,dat,&resul1,&resul2);
      printf("Pointwise P: %.8le\n", dat);
      printf("Corrected P: %.8le\n", resul1);
    } else if (mode == IT) {
      per_marker_threshold(maxstats,size,dat,&resul1,&resul2);
      printf("Signif threshold: %.8le\n", dat);
      printf("Per-marker thres: %.8le\n", resul1);
    }
    printf("STDEV(%% of P): %.8le", resul2);
    neff = (mode == IP)?resul1/dat:dat/resul1;
    if (resul1 < EPS || resul1 >= 1.-EPS) {
      fprintf(stdout, "(%.5lf%%)\n", 0.);
      fprintf(stdout, "# of Effective tests: ----\n");
    } else {
      fprintf(stdout, "(%.5lf%%)\n", resul2/resul1*100);
      fprintf(stdout, "# of Effective tests: %.0lf\n",neff);
    }
    if (resul1 < EPS) {
      fprintf(stdout, "Note: Need more sampling\n");
    } else if (resul1 >= 1.-EPS) {
      fprintf(stdout, "\n");
    }
  }
  return 0;
}

void correct_pvalue(double* maxstats, int size, double data, double* result1, double* result2) {
  int lfin, rfin, mid, x;
  double ltailp, negz, absz;
  ltailp = data/2.;
  vdCdfNormInv(1,&ltailp,&negz);
  absz = (-1.)*negz;
  lfin = 0;
  rfin = size-1;
  do {
    mid = (lfin+rfin)/2;
    if (absz < maxstats[mid]) {
      rfin = mid;
    } else if (absz > maxstats[mid]) {
      lfin = mid;
    } else if (fabs(absz - maxstats[mid])<EPS) {
      break;
    }
  } while (rfin - lfin > 2);
  x = mid;

  if (absz <= maxstats[x]) {
    while(x >= 0 && absz <= maxstats[x]) --x;
  } else if (absz > maxstats[x]) {
    while(x < size && absz > maxstats[x]) ++x;
    if (x < size) --x;
  }

  if (x == size) { *result1 = 0.; }
  else { *result1 = ((double)size - (x+1))/size; }
  *result2 = sqrt((*result1)*(1-(*result1))/size);
  return;
}

void per_marker_threshold(double* maxstats, int size, double data, double* result1, double* result2) {
  double negz, ltailp;
  negz = (-1.)*maxstats[size - (int)(size*data)];
  vdCdfNorm(1,&negz,&ltailp);
  *result1 = 2.*ltailp;
  *result2 = sqrt((*result1)*(1-(*result1))/size);
  return;
}

void print_usage() {
  fprintf(stderr,"usage:\n\
   P-value correction:\n\
      ./slide_4correct -p [sorted file] [pointwise-p file] [final output file] [mapfile(optional)]\n\
   Per-marker threshold:\n\
      ./slide_4correct -t [sorted file] [per-marker thres file] [final output file]\n\
   P-value correction (interactive):\n\
      ./slide_4correct -ip [sorted file] [pointwise-p]\n\
   Per-marker threshold (interactive):\n\
      ./slide_4correct -it [sorted file] [per-marker thres]\n\
");
  exit(-1);
}
