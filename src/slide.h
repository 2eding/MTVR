//slide.h

int compare_by_stat(const void *a, const void *b);
int compare_double( const void *s, const void *t);
int compare_double_rev( const void *s, const void *t);
void print_usage();

// for analysis
void correct_pvalue(double*, int, double, double*, double*);
void per_marker_threshold(double*, int, double, double*, double*);

// file io
FILE* readfile(char* filename);
FILE* writefile(char* filename);
FILE* readfile(char* filename) {
  FILE* fp;
  if ( (fp = fopen(filename,"r")) == NULL ) { fprintf(stderr,"Cannot read file %s\n",filename);exit(-1); }
  return(fp);
}
FILE* writefile(char* filename) {
  FILE* fp;
  if ( (fp = fopen(filename,"w")) == NULL ) { fprintf(stderr,"Cannot write to file %s\n",filename);exit(-1); }
  return(fp);
}
//
#define min(a,b) (((a)<(b)) ? (a) : (b))
#define max(a,b) (((a)>(b)) ? (a) : (b))
