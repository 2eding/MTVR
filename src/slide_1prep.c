#define _FILE_OFFSET_BITS  64 // large file support
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>
#include "slide.h"
#include <R.h>
#define MAXNSEC 1000 // maximum number of sections for scaling factor
#define EIGEN_EPS 1E-10 // epsilon to add to the matrix diagonal
#define DBL_EPS_HARD 1E-20 // hard epsilon
#define DBL_EPS_SOFT 1E-10 // soft epsilon
#define P_THRES .05 // starting point of scaling factor computation
#define BONFFACTOR 10 // conservative estimate of bonferroni factor
#define DECREASE_RATIO .1 // control between-sections of scaling factor part
#define DOSE_EFF .5 // Armitage test
struct tables { double stat, prob; };

int main (int argc, char** argv)
{
  MKL_Set_Num_Threads(1); // avoid multi-threading
  int i,j,k,l,n11,n12,n21,n22,c_p,c_n,r0,r1,r2,g0,g1,g2;
  int nsnps, ncont, ncase, nindv, nmin, szwin, ws;
  int *a0, *a1, *a2, ntags, *tags, *connsec;
  char *bins, *ucpi, *ucpj;
  int lbound, rbound;
  char *binf, *outf, buf[1000], *line, *token;
  double f, *rs, *rsq, *convec, *constd, *scale, *thres, **conthres, **conscale, *denoms;
  double *vecx, *vecy, effntags;
  double cutpv, numer, max_st, *h1, *logfact, con, logfactcon, subt, stat, prob, statconst, cutstat, midp;
  double statcutconst_p, statcutconst_n, subcon, nextp, prevprob, maxprob;
  double negzcut, ltailp, *maxchi2, *negzs, *pointwiseps;
  double x = DOSE_EFF, xx = x*x ; // dose effect
  int ione = 1, info, nrow, ncol;
  FILE* fp;
  struct tables *t_p, *t_n;
  enum MODE {H, G, QH, QG, C} mode;
  // receive arguments
  if (argc < 2) print_usage();
  if (!strcmp(argv[1], "-H")) { mode = H; }
  else if (!strcmp(argv[1], "-G")) { mode = G; }
  else if (!strcmp(argv[1], "-QH")) { mode = QH; }
  else if (!strcmp(argv[1], "-QG")) { mode = QG; }
  else if (!strcmp(argv[1], "-C")) { mode = C; }
  else { print_usage(); }
  if (mode == H || mode == G) {
    if (argc != 7) print_usage();
    binf = argv[2];
    ncont = atoi(argv[3]); // # controls
    ncase = atoi(argv[4]); // # cases
    nindv = ncont + ncase;
    nmin = min(ncont, ncase);
    szwin = atoi(argv[5]); // window size
    outf = argv[6];
  } else if (mode == QH || mode == QG) {
    if (argc != 6) print_usage();
    binf = argv[2];
    nindv = atoi(argv[3]);
    szwin = atoi(argv[4]); // window size
    outf = argv[5];
    ncase = ncont = nmin = 0;
  } else if (mode == C) {
    if (argc != 5) print_usage();
    binf = argv[2];
    szwin = atoi(argv[3]);
    outf = argv[4];
    ncase = ncont = nmin = nindv = 0;
  }
  //count # of row and column (and sanity check)
  //  fprintf(stderr,">>> Counting SNPs\n");
  nrow = 0;
  line = (char*)malloc(20000001*sizeof(char));
  fp = readfile(binf);
  while(fgets(line,20000000,fp)!=NULL){
    nrow++;
    token = strtok(line," \t\n");
    if(token==NULL || token[0]=='#') nrow--;
  }
  if (nrow==0) { fprintf(stderr,"File empty\n");exit(-1); }
  fclose(fp);
  nsnps = nrow;
  if (mode == C) {
    nsnps++;
    ntags = nsnps;
  }
  fprintf(stderr,">>> Number of SNPs: %d\n", nsnps);

  // memory alloc
  bins = (char*)malloc(sizeof(char)*(szwin+1)*nindv);
  tags = (int*)malloc(sizeof(int)*nsnps);
  a0 = (int*)malloc(sizeof(int)*nsnps);
  a1 = (int*)malloc(sizeof(int)*nsnps);
  a2 = (int*)malloc(sizeof(int)*nsnps);
  h1 = (double*)malloc(sizeof(double)*nsnps);
  denoms = (double*)malloc(sizeof(double)*nsnps);
  rs = (double*)calloc(sizeof(double),szwin*szwin);
  rsq = (double*)malloc(sizeof(double)*(szwin)*(szwin));//denotes square matrix of r.. not r^2.
  vecx = (double*)malloc(sizeof(double)*(szwin));
  vecy = (double*)malloc(sizeof(double)*(szwin));
  constd = (double*)malloc(sizeof(double)*nsnps);
  convec = (double*)malloc(sizeof(double)*nsnps*szwin);
  conscale = (double**)malloc(sizeof(double*)*2*nsnps);
  conthres = (double**)malloc(sizeof(double*)*2*nsnps);
  scale = (double*)malloc(sizeof(double)*MAXNSEC);
  thres = (double*)malloc(sizeof(double)*MAXNSEC);
  connsec = (int*)malloc(sizeof(int)*2*nsnps);
  t_p = (struct tables*)malloc((nmin+2)*(nmin+1)/2 * sizeof(struct tables));
  t_n = (struct tables*)malloc((nmin+2)*(nmin+1)/2 * sizeof(struct tables));
  logfact = (double*)malloc(sizeof(double)*(nindv+1));
  negzs = (double*)malloc(sizeof(double)*nsnps);
  pointwiseps = (double*)malloc(sizeof(double)*nsnps);
  maxchi2 = (double*)malloc(sizeof(double)*nsnps);

  //read binary file (1st time) to obtain allele counts
  if (mode != C) {
    fprintf(stderr,">>> Sanity check of the matrix\n");
    ntags = nsnps;
    for(i=0;i<nsnps;i++) tags[i] = i;
    fp = readfile(binf);
    for(i=0, k=0; i < nsnps; ++i) {
      if (i>1&&i%10000==0) fprintf(stderr,"SNP %d checked\n",i);
      do {
	fgets(line,20000000,fp);
	token = strtok(line," \t\n");
      } while (token==NULL || token[0]=='#');
      ncol = 0;
      bins[ncol++] = token[0];
      while(token = strtok(NULL, " \t\n")) {
	if (ncol >= nindv) {fprintf(stderr, "# of column has to be exactly # of individuals)\n");exit(-1);}
	bins[ncol++] = token[0];
      }
      if (ncol < nindv) {fprintf(stderr, "# of column has to be exactly # of individuals)\n");exit(-1);}
      //    for (j=0; j < nindv; ++j) fscanf(fp, "%d", &bins[j]); // 0:homo major,  1:hetero,  2:homo minor
      a0[i] = a1[i] = a2[i] = r1 = r2 = 0;
      for (j=0; j < ncont; ++j) {
	if (bins[j] == '0') { ++a0[i]; }
	else if (bins[j] == '1') { ++a1[i]; }// allele 1 count
	else if (bins[j] == '2') { ++a2[i]; } // allele 2 count
	else { fprintf(stderr, "Only 0/1/2 are allowed in data\n");exit(-1); }
      }
      for (j=ncont; j < nindv; ++j) {
	if (bins[j] == '0') { ++a0[i]; }
	else if (bins[j] == '1') { ++a1[i]; ++r1; }// allele 1 count
	else if (bins[j] == '2') { ++a2[i]; ++r2; } // allele 2 count
	else { fprintf(stderr, "Only 0/1/2 are allowed in data\n");exit(-1); }
      }
      h1[i] = x*a1[i] + a2[i]; //  haplotypic allele count/2 (average allele dosage)
      denoms[i] = sqrt((double)nindv*(xx*a1[i]+a2[i]) - h1[i]*h1[i]);
      if (mode == H || mode == G) { // compute each SNP's statistic
	if (denoms[i] < DBL_EPS_HARD) {
	  negzs[i] = 0.;
	} else {
	  numer = ((double)nindv*(x*r1+r2)-ncase*h1[i]);
	  negzs[i] = (-1.)*sqrt(nindv*numer*numer/ncase/ncont/denoms[i]/denoms[i]);
	}
      }
      if (a0[i] == nindv || a1[i] == nindv || a2[i] == nindv) { //filter non-polymorphic sites
	for (l=k; l<ntags; ++l) { tags[l] = tags[l+1]; }
	--ntags;
      } else { ++k; }
      if ((mode == H || mode == QH) && a2[i] > 0) {
	fprintf(stderr, "For haplotype mode, '2' is not allowed in data\n");exit(-1);
      }
    }
    fclose(fp);
  }

  // scaling factor part. only work for H and G mode.
  if (mode == H || mode == G) {
    fprintf(stderr,">>> Calculating scaling factors\n");
    /* Step1. Compute the maximum statistic each SNP can achieve. */
    for(i=0; i < ntags; ++i) {
      r2 = min(a2[tags[i]], ncase); //extreme case 1
      r1 = min(a1[tags[i]], ncase - r2);
      numer = ((double)nindv*(x*r1+r2)-ncase*(x*a1[tags[i]]+a2[tags[i]]));
      max_st = nindv*numer*numer/ncase/ncont/denoms[tags[i]]/denoms[tags[i]];
      r0 = min(a0[tags[i]], ncase); //extreme case 2
      r1 = min(a1[tags[i]], ncase - r0);
      r2 = ncase - r0 - r1;
      numer = ((double)nindv*(x*r1+r2)-ncase*(x*a1[tags[i]]+a2[tags[i]]));
      max_st = max(max_st, nindv*numer*numer/ncase/ncont/denoms[tags[i]]/denoms[tags[i]]);
      maxchi2[tags[i]] = max_st;
    }
    /* Step2. Compute the scaling factor for each SNP. */
    effntags = ntags/BONFFACTOR;
    for(i=1, logfact[0]=0;i <= nindv;++i) logfact[i] = logfact[i-1] + log(i);
    logfactcon = logfact[ncase]+logfact[ncont]-logfact[nindv];
    for(i=0;i < ntags;++i) { // for each SNP
      if (i>1&&i%10000==0) fprintf(stderr,"SNP %d calculated\n",i);
      g0 = a0[tags[i]];
      g1 = a1[tags[i]];
      g2 = a2[tags[i]];
      //      fprintf(stderr,"%d %d %d %d %d\n",i, tags[i], g0, g1, g2);
      statconst = denoms[tags[i]]*sqrt(ncase*ncont)/nindv/sqrt(nindv);
      subt = (double)ncase*(x*g1+g2)/nindv;
      ltailp = P_THRES/effntags/2;
      vdCdfNormInv(1,&ltailp,&negzcut);
      cutstat = (-1.)*negzcut;
      statcutconst_p = (-1.)*negzcut * statconst + subt;
      statcutconst_n = negzcut * statconst + subt;
      //list all possible tables
      con = logfact[g0]+logfact[g1]+logfact[g2]+logfactcon;
      c_p = c_n = 0;
      for (r1=max(0,ncase-g2-g0);r1<=min(g1,ncase);++r1) {
	subcon = con-logfact[r1]-logfact[g1-r1];
	lbound = max(0,ncase-r1-g0);
	rbound = min(g2,ncase-r1);
	r2 = (int)rint((double)g2*(ncase-r1)/(g0+g2));
	r0 = ncase-r1-r2;
	maxprob = exp(subcon-logfact[r0]-logfact[g0-r0]-logfact[r2]-logfact[g2-r2]);//maximum probability
	if (maxprob >= DBL_EPS_HARD) {
	  prevprob = 0.;
	  for (r2=min(rbound,(int)ceil(statcutconst_n-x*r1));r2>=lbound;--r2) {
	    r0 = ncase-r1-r2;
	    stat = x*r1+r2 - subt;
	    prob = exp(subcon-logfact[r0]-logfact[g0-r0]-logfact[r2]-logfact[g2-r2]);
	    if (prob >= DBL_EPS_HARD) {
	      t_n[c_n].stat = stat;
	      t_n[c_n].prob = prob;
	      ++c_n;
	    } else if (prob < prevprob) break; //if prob is decreasing & below epsilon, stop
	    prevprob = prob;
	  }
	  prevprob = 0.;
	  for (r2=max(lbound,(int)floor(statcutconst_p-x*r1));r2<=rbound;++r2) {
	    r0 = ncase-r1-r2;
	    stat = x*r1+r2 - subt;
	    prob = exp(subcon-logfact[r0]-logfact[g0-r0]-logfact[r2]-logfact[g2-r2]);
	    if (prob >= DBL_EPS_HARD) {
	      t_p[c_p].stat = stat;
	      t_p[c_p].prob = prob;
	      ++c_p;
	    } else if (prob < prevprob) break; //if prob is decreasing & below epsilon, stop
	    prevprob = prob;
	  }
	}
      }
      qsort(t_n, c_n, sizeof(struct tables), compare_by_stat); // sort the tables (-)
      qsort(t_p, c_p, sizeof(struct tables), compare_by_stat); // sort the tables (+)
      for (j=c_p-2;j >= 0;--j) t_p[j].prob += t_p[j+1].prob; // accumulate probs to get C.D.F.
      for (j=1;j < c_n;++j) t_n[j].prob += t_n[j-1].prob;
      // compute scaling factor for each sections (-)
      cutpv = P_THRES;
      ltailp = cutpv/effntags/2;
      vdCdfNormInv(1,&ltailp,&negzcut);
      cutstat = negzcut*statconst;
      l = 0;
      for (j=c_n-1;j >= 0 && l < MAXNSEC-1;--j) {
	if (t_n[j].stat < cutstat) {
	  for (k = j-1;k >= 0 && t_n[j].stat-t_n[k].stat < DBL_EPS_SOFT;--k);//thru the same p-values.
	  nextp = (k == -1)?0.:t_n[k].prob;
	  midp = (t_n[j].prob+nextp)/2; // mid-p value
	  thres[l] = t_n[j].stat / statconst;
	  vdCdfNormInv(1,&midp,&negzcut);
	  scale[l] = thres[l] / negzcut;
	  ++l;
	  while(t_n[j].stat < cutstat) {
	    cutpv = cutpv * DECREASE_RATIO;
	    ltailp = cutpv/effntags/2;
	    vdCdfNormInv(1,&ltailp,&negzcut);
	    cutstat = negzcut*statconst;
	  }
	}
      }
      thres[l] = (-1.)*sqrt(maxchi2[i]);
      scale[l] = (l == 0)?0:scale[l-1];
      l++;
      thres[l] = thres[l-1]-.01;
      scale[l] = DBL_EPS_SOFT;
      l++;
      connsec[2*i] = l;
      conthres[2*i] = (double*)malloc(sizeof(double)*(l));
      conscale[2*i] = (double*)malloc(sizeof(double)*(l));
      memcpy(conthres[2*i],thres,sizeof(double)*l);
      memcpy(conscale[2*i],scale,sizeof(double)*l);
      // compute scaling factor for each sections (+)
      cutpv = P_THRES;
      ltailp = cutpv/effntags/2;
      vdCdfNormInv(1,&ltailp,&negzcut);
      cutstat = (-1.)*negzcut*statconst;
      l = 0;
      for (j=0;j < c_p && l < MAXNSEC-1;++j) {
	if (t_p[j].stat > cutstat) {
	  for (k = j+1;k < c_p && t_p[k].stat-t_p[j].stat < DBL_EPS_SOFT;++k);
	  nextp = (k == c_p)?0.:t_p[k].prob;
	  midp = (t_p[j].prob+nextp)/2; // mid-p value
	  thres[l] = t_p[j].stat / statconst;
	  vdCdfNormInv(1,&midp,&negzcut);
	  scale[l] = (-1.) * thres[l] / negzcut;
	  ++l;
	  while(t_p[j].stat > cutstat) {
	    cutpv = cutpv * DECREASE_RATIO;
	    ltailp = cutpv/effntags/2;
	    vdCdfNormInv(1,&ltailp,&negzcut);
	    cutstat = (-1.)*negzcut*statconst;
	  }
	}
      }
      thres[l] = sqrt(maxchi2[i]);
      scale[l] = (l == 0)?0:scale[l-1];
      l++;
      thres[l] = thres[l-1]+.01;
      scale[l] = DBL_EPS_SOFT;
      l++;
      connsec[2*i+1] = l;
      conthres[2*i+1] = (double*)malloc(sizeof(double)*(l));
      conscale[2*i+1] = (double*)malloc(sizeof(double)*(l));
      memcpy(conthres[2*i+1],thres,sizeof(double)*l);
      memcpy(conscale[2*i+1],scale,sizeof(double)*l);
    } // for(each snp)
  } else if (mode == QH || mode == QG || mode == C) {
    for(i=0;i < ntags;++i) {  // unit scaling.
      l = 0;
      thres[l] = -100;
      scale[l] = 1;
      l++;
      connsec[2*i] = l;
      conthres[2*i] = (double*)malloc(sizeof(double)*(l));
      conscale[2*i] = (double*)malloc(sizeof(double)*(l));
      memcpy(conthres[2*i],thres,sizeof(double)*l);
      memcpy(conscale[2*i],scale,sizeof(double)*l);
      l = 0;
      thres[l] = 100;
      scale[l] = 1;
      l++;
      connsec[2*i+1] = l;
      conthres[2*i+1] = (double*)malloc(sizeof(double)*(l));
      conscale[2*i+1] = (double*)malloc(sizeof(double)*(l));
      memcpy(conthres[2*i+1],thres,sizeof(double)*l);
      memcpy(conscale[2*i+1],scale,sizeof(double)*l);
    }
  }// end of scaling factor part

  //Computing, finally, convec and constd.
  fprintf(stderr,">>> Computing conditional distributions\n");
  constd[0] = 1.0;
  fp = readfile(binf); // read binary file (2nd time)
  if (mode != C) {
    for(l=0;l <= tags[0];++l) {
      do {
	fgets(line,20000000,fp);
	token = strtok(line," \t\n");
      } while (token==NULL || token[0]=='#');
      j = 0;
      bins[j++] = token[0];
      while(token = strtok(NULL, " \t\n"))
	bins[j++] = token[0];
    }
  }
  for(i=1; i < ntags ; ++i) {// 'i' is of interest
    if (i>1&&i%10000==0) fprintf(stderr,"SNP %d processed\n",i);
    ws = min(szwin, i);//effective window lookup size
    // compute "r" for one more line
    if (mode != C) {
      for(l=tags[i-1];l < tags[i];++l) {
	do {
	  fgets(line,20000000,fp);
	  token = strtok(line," \t\n");
	} while (token==NULL || token[0]=='#');
	j = 0;
	bins[(i%(szwin+1))*nindv+j] = token[0];
	j++;
	while(token = strtok(NULL, " \t\n")) {
	  bins[(i%(szwin+1))*nindv+j] = token[0];
	  j++;
	}
      }
      ucpi = &bins[(i%(szwin+1))*nindv];
      for(j=szwin-ws; j < szwin; ++j) {
	ucpj = &bins[((i-szwin+j)%(szwin+1))*nindv];
	n11 = n12 = n21 = n22 = 0;
	for( k=0; k < nindv; ++k)
	  if (ucpi[k] == '1') {if (ucpj[k] == '1') ++n11; else if (ucpj[k] == '2') ++n12;}
	  else if (ucpi[k] == '2') {if (ucpj[k] == '1') ++n21; else if (ucpj[k] == '2') ++n22;}
	rs[(i%szwin)*szwin+j] = ((double)nindv*(xx*n11+x*(n12+n21)+n22) - h1[tags[i]]*h1[tags[i-szwin+j]])/denoms[tags[i]]/denoms[tags[i-szwin+j]];
      }
    } else if (mode == C) {
      for(j=0;j < szwin; ++j) {
	fscanf(fp, "%lf",  &rs[(i%szwin)*szwin+j]);
      }
    }
    for(j=0; j < ws; ++j) {
      rsq[j*ws+j] = 1.0+EIGEN_EPS;
      for(k=0; k < j; ++k) { rsq[j*ws+k] = rs[((i-ws+j)%szwin)*szwin+(szwin-j+k)]; }
      vecx[j] = vecy[j] = rs[(i%szwin)*szwin+szwin-ws+j];
    }
    dpotrf("U",&ws,rsq,&ws,&info);
    dpotrs("U",&ws,&ione,rsq,&ws,vecy,&ws,&info);
    f = ddot(&ws, vecx, &ione, vecy, &ione);
    memset(&convec[i*szwin],0,sizeof(double)*(szwin-ws));
    memcpy(&convec[i*szwin+szwin-ws],vecy,sizeof(double)*ws);
    constd[i] = sqrt(max(1.+EIGEN_EPS-f,0));
  }
  fclose(fp);

  // Write the data.
  sprintf(buf,"%s.txt_polySNPs",outf);
  fp = writefile(buf);
  for(i=0;i < ntags;++i) fprintf(fp, "%d\n", tags[i]);
  fclose(fp);

  sprintf(buf,"%s.bin_constd",outf);
  fp = writefile(buf);
  fwrite(constd,sizeof(double),ntags,fp);
  fclose(fp);

  sprintf(buf,"%s.bin_convec",outf);
  fp = writefile(buf);
  fwrite(convec,sizeof(double),ntags*szwin,fp);
  fclose(fp);

  /// print purpose
/*   sprintf(buf,"%s.constd.txt",outf); */
/*   fp = writefile(buf); */
/*   for(i=0;i < ntags;++i) fprintf(fp, "%.3lf\n", constd[i]); */
/*   fclose(fp); */
/*   sprintf(buf,"%s.convec.txt",outf); */
/*   fp = writefile(buf); */
/*   for(i=0;i < ntags;++i) { */
/*     for(k=0; k < szwin;k++) { fprintf(fp, "%.3lf ", convec[i*szwin+k]); } */
/*     fprintf(fp, "\n"); */
/*   } */
/*   fclose(fp); */

  sprintf(buf,"%s.txt_nsec",outf);
  fp = writefile(buf);
  for(i=0;i < 2*ntags;++i) fprintf(fp, "%d\n", connsec[i]);
  fclose(fp);

  sprintf(buf,"%s.bin_scale",outf);
  fp = writefile(buf);
  for(i=0;i < 2*ntags;++i) {
    fwrite(&conthres[i][0],sizeof(double),connsec[i],fp);
    fwrite(&conscale[i][0],sizeof(double),connsec[i],fp);
  }
  fclose(fp);

  sprintf(buf,"%s.txt_wsize",outf);
  fp = writefile(buf);
  fprintf(fp, "%d\n", szwin);
  fclose(fp);

  if (mode == H || mode == G) { // each SNP's statistic
    //    vmlSetMode(VML_HA);
    vdCdfNorm(nsnps,negzs,pointwiseps);
    sprintf(buf,"%s.txt_pointwisep",outf);
    fp = writefile(buf);
    //    for(i=0;i < nsnps;++i) fprintf(fp, "%.10lf\t%.10lf\n", negzs[i]*negzs[i], pointwiseps[i]*2);
    for(i=0;i < nsnps;++i) {
      fprintf(fp, "%.10le\n", pointwiseps[i]*2);
      //      fprintf(stdout, "%.10le %.10le\n", negzs[i], pointwiseps[i]*2);
    }
    fclose(fp);
  }

  sprintf(buf,"%s.txt_log",outf);
  fp = writefile(buf);
  fprintf(fp,">>> Command invoked:");
  for (i = 0;i < argc;i++) fprintf(fp," %s",argv[i]); fprintf(fp,"\n");
  fprintf(fp,">>> Input file %s has been successfully pre-processed\n", binf);
  fprintf(fp,">>> Number of SNPs: %d\n", nsnps);
  fprintf(fp,">>> Number of removed non-polymorphic SNPs: %d\n", nsnps-ntags);
  fprintf(fp,">>> Number of polymorphic SNPs processed: %d\n", ntags);
  fprintf(fp,">>> The following files were created.\n");
  fprintf(fp,">>> %s.bin_convec\t:\tconditional distribution information\n",outf);
  fprintf(fp,">>> %s.bin_constd\t:\tconditional distribution information\n",outf);
  fprintf(fp,">>> %s.bin_scale\t:\tscaling factor information\n",outf);
  fprintf(fp,">>> %s.txt_nsec\t:\tscaling factor information\n",outf);
  fprintf(fp,">>> %s.txt_polySNPs\t:\tindex of polymorphic SNPs\n",outf);
  fprintf(fp,">>> %s.txt_wsize\t:\tstore window size\n",outf);
  if (mode == H || mode == G) fprintf(fp,">>> %s.txt_pointwisep\t:\tpointwise p-values of every SNP\n",outf);
  fprintf(fp,">>> %s.txt_log\t:\tlog file of this preprocess\n",outf);
  fclose(fp);

  fprintf(stdout,">>> Input file %s has been successfully pre-processed\n", binf);
  fprintf(stdout,">>> Number of removed non-polymorphic SNPs: %d\n", nsnps-ntags);
  fprintf(stdout,">>> Number of polymorphic SNPs processed: %d\n", ntags);
  fprintf(stdout,">>> The following files were created.\n");
  fprintf(stdout,">>> %s.bin_convec\t:\tconditional distribution information\n",outf);
  fprintf(stdout,">>> %s.bin_constd\t:\tconditional distribution information\n",outf);
  fprintf(stdout,">>> %s.bin_scale\t:\tscaling factor information\n",outf);
  fprintf(stdout,">>> %s.txt_nsec\t:\tscaling factor information\n",outf);
  fprintf(stdout,">>> %s.txt_polySNPs\t:\tindex of polymorphic SNPs\n",outf);
  fprintf(stdout,">>> %s.txt_wsize\t:\tstore window size\n",outf);
  if (mode == H || mode == G) fprintf(stdout,">>> %s.txt_pointwisep\t:\tpointwise p-values of every SNP\n",outf);
  fprintf(stdout,">>> %s.txt_log\t:\tlog file of this preprocess\n",outf);
  return 0;
}

int compare_by_stat(const void *a, const void *b) {
  double temp = ((struct tables*)a)->stat - ((struct tables*)b)->stat;
  return (temp > 0.)?(1):((temp < 0.)?(-1):0); }

void print_usage() {
  fprintf(stderr,"usage:\n\
  Case/control study\n\
      ./slide_1prep [-G|-H] [input file] [#control] [#case] [window size] [prep file]\n\
  Quantitative traits\n\
      ./slide_1prep [-QG|-QH] [input file] [#individuals] [window size] [prep file]\n\
  For pre-computed band covariance matrix\n\
     ./slide_1prep [-C] [input file] [window size] [prep file]\n\
");
  exit(-1);
}

