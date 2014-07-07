#include "linalg.h"

long double *rsFirstEigenvector(const double** A, const long n, const long maxIterations, const double precision, const BOOL verbose)
{
    // Power iteration method
    // http://en.wikipedia.org/wiki/Power_iteration
    
    long double *tmp = malloc(n*sizeof(long double));
        
    // initialize seed eigenvector
    long double *b = malloc(n*sizeof(long double));
    for (long i=0; i<n; i=i+1) {
        b[i] = 1.0;
    }
    
    BOOL converged = FALSE;
    
    // run power iteration
    for (long iteration=1L; iteration<maxIterations; iteration=iteration+1) {
                
        // calculate the matrix-by-vector product tmp=Ab
        free(tmp);
        tmp = rsMatrixByVectorProduct(A, b, n, n);
        
        // calculate the vector norm(euclidean)
        const long double norm = rsEuclideanNorm(tmp, n);
        long double invNorm = ((long double)1) / norm;
        
        // normalize it for the next iteration
        rsScaleVector(tmp, n, invNorm);
        
        // check if it has converged
        rsVectorSwap(b, tmp, n);
        rsVectorSub(tmp, b, n);
        long double mean = fabsl(rsVectorMean(tmp, n));
        
        if (verbose) {
            fprintf(stdout, "Iteration %ld/%ld Mean Difference: %.10Lf Target: %.10f\n", iteration, maxIterations, mean, precision);
        }
        
        if ( mean <= precision ) {
            break;
        }
    }
    
    free(tmp);
    
    return b;
}

/*
 * Takes in a symmetric matrix and checks if all
 * eigenvalues are positive. If not it will 'fix'
 * those that are negative by setting them to the
 * average eigenvalue and then reconstructs the
 * original matrix. Finally, the number of altered
 * eigenvalues is returned.
 */
long rsMakePositiveDefiniteSymmetric(gsl_matrix* A)
{
    assert(A->size1 == A->size2);
	
	long nEVAltered   = 0;
    const long   n    = A->size1;
	const double zero = 1e-10;
	
	gsl_vector* eigenvalues  = gsl_vector_alloc(n);
    gsl_matrix* eigenvectors = gsl_matrix_alloc(n, n);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(n);

	// perform eigenvalue decomposition
    gsl_eigen_symmv(A, eigenvalues, eigenvectors, workspace);
	gsl_eigen_gensymmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_VAL_ASC);
	
	// compute mean eigenvalue(of those that are retained)
	double norm       = 0.0;
	long   normValues = 0;
	for ( long i=0; i<n; i=i+1 ) {
		const double v = gsl_vector_get(eigenvalues, i);
		if ( v > zero ) {
			norm       = norm + v;
			normValues = normValues + 1;
		}
	}
	double meanValue = norm / normValues;
	
	// check for negative eigenvalues
	for ( long i=0; i<n; i=i+1 ) {
		
		if ( gsl_vector_get(eigenvalues, i) > zero ) {
			continue;
		}
		
		gsl_vector_set(eigenvalues, i, meanValue);
		nEVAltered = nEVAltered + 1;
	}
	
	// if eigenvalues were changed reconstruct matrix
	if ( nEVAltered > 0 ) {
		// create diagional eigenvalue matrix
		gsl_matrix* eigenvaluesMatrix   = gsl_matrix_alloc(n, n);
	    gsl_vector_view eigenvaluesDiag = gsl_matrix_diagonal(eigenvaluesMatrix);
	    gsl_matrix_set_all(eigenvaluesMatrix, 0.0);
	    gsl_vector_memcpy(&eigenvaluesDiag.vector, eigenvalues);
	
		// reconstruct matrix
		gsl_matrix* tmp = gsl_matrix_alloc(n, n); 
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigenvectors, eigenvaluesMatrix, 0.0, tmp);
		gsl_blas_dgemm(CblasNoTrans,   CblasTrans, 1.0,          tmp,      eigenvectors, 0.0,   A);
		
	    gsl_matrix_free(tmp);
	    gsl_matrix_free(eigenvaluesMatrix);
	}

    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(eigenvectors);
    gsl_vector_free(eigenvalues);

	return nEVAltered;
}

double **d2matrix(int yh, int xh)
{
    long int j;
    long int nrow = yh+1;
    long int ncol = xh+1;
    double **t;
    
    
    /** allocate pointers to ydim */
    t=(double **) malloc((size_t)((nrow)*sizeof(double*)));
    if (!t) RSIOERR("d2matrix: allocation failure");
    
    /** allocate pointers for zdim */
    t[0]=(double *) malloc((size_t)((ncol*nrow)*sizeof(double)));
    if (!t[0]) RSIOERR("d2matrix: allocation failure");
    
    /** point everything to the data blob */
    for(j=1L;j<nrow;j++) t[j]=t[j-1]+ncol;
    
    return t;
}

/*
 * Computes the matrix by vector product
 * y = A*x
 * n..width of the matrix and vector length
 * m..height of the matrix
 * The matrix is assumed to be in the following format A[m][n]
 */
long double *rsMatrixByVectorProduct(const double **A, const long double *x, const long n, const long m) {
    long double *y = malloc(m*sizeof(long double));
    long i,j;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,j) shared(y)
    {
        #pragma omp for schedule(guided)
        for (i=0; i<m; i=i+1L) {
            long double product = 0.0;
            
            for (j=0; j<n; j=j+1) {
                product += (long double)A[i][j] * (long double)x[j];
            }
            
            y[i] = (double)product;
        }
    }
    
    return y;
}

long double rsEuclideanNorm(const long double *x, const long n)
{
    long double tss = 0.0;
    
    for (long i=0; i<n; i=i+1) {
        tss += x[i];
    }
    
    return sqrtl(tss);
}

void rsScaleVector(long double *x, const long n, const long double factor)
{
    long i;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i) shared(x)
    {
        #pragma omp for schedule(guided)
        for (i=0L; i<n; i++) {
            x[i] *= factor;
        }
    }
}

void rsVectorSub(long double *x, const long double *y, const long n)
{
    long i;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i) shared(x,y)
    {
        #pragma omp for schedule(guided)
        for (i=0L; i<n; i=i+1L) {
            x[i] -= y[i];
        }
    }
}

void rsVectorSwap(long double *x, long double *y, const long n)
{
    long double tmp;
    long i;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,tmp) shared(x,y)
    {
        #pragma omp for schedule(guided)
        for (i=0L; i<n; i++) {
            tmp  = x[i];
            x[i] = y[i];
            y[i] = tmp;
        }
    }
}

BOOL rsVectorContains(const long *x, const long n, const long element) {
	int i;
	for ( i=0; i<n; i=i+1 ) {
		if ( x[i] == element ) {
			return TRUE;
		}
	}
	
	return FALSE;
}

void rsMatrixConversion(double **A, const long m, const long n, const int mode)
{
    long i,j;
    
    #pragma omp parallel num_threads(rsGetThreadsNum()) private(i,j) shared(A)
    {
        #pragma omp for schedule(guided)
        for (i=0; i<m; i=i+1L) {
            if ( mode == RSMATRIXCONVERSION_ABSOLUTE ) {
                for (j=0; j<n; j=j+1) {
                    A[i][j] = fabs(A[i][j]);
                }
            } else if ( mode == RSMATRIXCONVERSION_POSITIVE ) {
                for (j=0; j<n; j=j+1) {
                    if ( A[i][j] < 0.0 ) {
                        A[i][j] = 0.0;
                    }
                }
            } else if ( mode == RSMATRIXCONVERSION_NEGATIVE ) {
                for (j=0; j<n; j=j+1) {
                    if ( A[i][j] >= 0.0 ) {
                        A[i][j] = 0.0;
                    } else {
                        A[i][j] = fabs(A[i][j]);
                    }
                }
            } else if ( mode == RSMATRIXCONVERSION_SCALED ) {
                // This is used to shift z-transformed correlations so that they're all positive.
                // Therefore the smallest possible(as implemented) correlation coefficient is computed.
                const double epsilon = 0.0000000001;
                const double shiftingConstant = fabs(log(epsilon / (2.0 - epsilon))) + epsilon;
                
                for (j=0; j<n; j=j+1) {
                    A[i][j] = A[i][j] + shiftingConstant;
                }
            }
        }
    }
}

long double rsVectorMean(const long double *x, const long n)
{
    long double mean = 0.0;
    
    for (long i=0L; i<n; i=i+1L) {
        mean += x[i];
    }
    
    return (double) (mean/(unsigned long)n);
}

void rs_matrix_fprintf(FILE *stream, const double **A, const long m, const long n, const char* fmt)
{
    for (long i=0; i<m; i=i+1) {
        for (long j=0; j<n; j=j+1) {
            fprintf(stream, fmt, A[i][j]);
            fprintf(stream, " ");
        }
        fprintf(stream, "\n");
    }
}

void rs_vector_fprintf(FILE *stream, const double *x, const long n, const char* fmt)
{
    for (long j=0; j<n; j=j+1) {
        fprintf(stream, fmt, x[j]);
        fprintf(stream, " ");
    }
    fprintf(stream, "\n");
}

void rs_vector_fprintfl(FILE *stream, const long double *x, const long n, const char* fmt)
{
        for (long j=0; j<n; j=j+1) {
            fprintf(stream, fmt, x[j]);
            fprintf(stream, " ");
        }
        fprintf(stream, "\n");
}

void rs_gsl_vector_fprintf(FILE *stream, gsl_vector *v, const char* fmt)
{
        for (long j=0; j<v->size; j=j+1) {
            fprintf(stream, fmt, gsl_vector_get(v, j));
            fprintf(stream, " ");
        }
        fprintf(stream, "\n");
}

int rs_gsl_matrix_fprintf(FILE *stream,gsl_matrix *m,char *fmt)
{
    size_t rows=m->size1;
    size_t cols=m->size2;
    size_t row,col,ml;
    int fill;
    char buf[100];
    gsl_vector *maxlen;
    
    maxlen=gsl_vector_alloc(cols);
    for (col=0;col<cols;++col) {
        ml=0;
        for (row=0;row<rows;++row) {
            sprintf(buf,fmt,gsl_matrix_get(m,row,col));
            if (strlen(buf)>ml)
                ml=strlen(buf);
        }
        gsl_vector_set(maxlen,col,ml);
    }
    
    for (row=0;row<rows;++row) {
        for (col=0;col<cols;++col) {
            sprintf(buf,fmt,gsl_matrix_get(m,row,col));
            fprintf(stream,"%s",buf);
            fill=gsl_vector_get(maxlen,col)+1-strlen(buf);
            while (--fill>=0)
                fprintf(stream," ");
        }
        fprintf(stream,"\n");
    }
    gsl_vector_free(maxlen);
    return 0;
}

BOOL rsSaveMatrix(const char *filename, const double** A, const long m, const long n)
{
    FILE *file;
    file = fopen(filename, "wb");
    BOOL success = TRUE;
    const size_t blockSize = sizeof(double);
    
    for (long row=0L; row<m && success; row=row+1) {
        success = success && fwrite(A[row], blockSize, n, file) == n;
    }
    fclose(file);
    
    return success;
}

BOOL rsLoadMatrix(const char *filename, double** A, const long m, const long n)
{
    FILE *file;
    file = fopen(filename, "r");
    BOOL success = TRUE;
    const size_t blockSize = sizeof(double);
    
    for (long row=0L; row<m && success; row=row+1) {
        success = success && fread(A[row], blockSize, n, file) == n;
    }
    fclose(file);
    
    return success;
}
