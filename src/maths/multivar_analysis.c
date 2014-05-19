#include "multivar_analysis.h"

// Originally taken from: https://gist.github.com/microo8/4065693
struct rsPCAResult rsGenericPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose)
{
    /*
    @param data - matrix of data vectors, MxN matrix, each column is a data vector, M - dimension, N - data vector count
    @param minVariance - min percentage of variance that should be retained
    @param nComponents - number of components that will be returned. ignored if less than 1
    */
	
    assert(data != NULL);
    long i;
    long j;
    long rows = data->size1;
    long cols = data->size2;
	int L = 0;
	struct rsPCAResult result;

	if ( verbose ) {
		fprintf(stdout, "Running PCA on a %ldx%ld matrix\n", rows, cols);
	}

    gsl_vector* mean = gsl_vector_alloc(rows);
	gsl_matrix* mean_substracted_data = gsl_matrix_alloc(rows, cols);
	gsl_matrix_memcpy(mean_substracted_data, data);
	
	gsl_matrix* covariance_matrix;
	gsl_vector* eigenvalues;
	gsl_matrix* eigenvectors;
	gsl_eigen_symmv_workspace* workspace;
	gsl_matrix_view L_eigenvectors;
	gsl_vector_view L_eigenvalues;
 
	#pragma omp parallel num_threads(rsGetThreadsNum()) shared(cols,rows,result,nComponents,verbose,mean,mean_substracted_data,covariance_matrix,eigenvalues,eigenvectors,workspace,L,L_eigenvectors,L_eigenvalues) private(i,j)
	{	
		// compute means
        #pragma omp for schedule(guided)
		for(i = 0; i < rows; i++) {
        	gsl_vector_set(mean, i, gsl_stats_mean(data->data + i * cols, 1, cols));
    	}
 
    	// get mean-substracted data into matrix mean_substracted_data.
		#pragma omp for schedule(guided)
    	for(i = 0; i < cols; i++) {
        	gsl_vector_view mean_substracted_point_view = gsl_matrix_column(mean_substracted_data, i);
        	gsl_vector_sub(&mean_substracted_point_view.vector, mean);
    	}

		#pragma omp sections
		{
			#pragma omp section
			{
		    	gsl_vector_free(mean);
 
		    	// Compute Covariance matrix
		    	covariance_matrix = gsl_matrix_alloc(rows, rows);
		    	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(cols - 1), mean_substracted_data, mean_substracted_data, 0.0, covariance_matrix);
		    	gsl_matrix_free(mean_substracted_data);
 
		    	// Get eigenvectors, sort by eigenvalue.
		    	eigenvalues = gsl_vector_alloc(rows);
		    	eigenvectors = gsl_matrix_alloc(rows, rows);
		    	workspace = gsl_eigen_symmv_alloc(rows);
		    	gsl_eigen_symmv(covariance_matrix, eigenvalues, eigenvectors, workspace);
		    	gsl_eigen_symmv_free(workspace);
		    	gsl_matrix_free(covariance_matrix);
 
		    	// Sort the eigenvectors
		    	gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);

				if ( verbose ) {
					fprintf(stdout, "\nPCAs(eigenvectors):\n");
		
					for ( i = 0; i<rows; i=i+1 ) {
						for ( j = 0; j<rows; j=j+1) {
		            		fprintf(stdout, "% 5.10f\t", gsl_matrix_get(eigenvectors, i, j));
						}
						fprintf(stdout, "\n");
		        	}
				}

				// Determine how many components to keep
				double explainedVariance = 0.0;
				double totalVariance = 0.0;
				for (unsigned int i=0; i<rows; i=i+1) {
					totalVariance = totalVariance + gsl_vector_get(eigenvalues, i);
				}
	
				if ( verbose ) {
					fprintf(stdout, "\nTotal variance: %.10f\n", totalVariance);
					fprintf(stdout, "\nEigenvalues:\n");
					for ( j = 0; j<rows; j=j+1) {
			        	fprintf(stdout, "%.10f\t", gsl_vector_get(eigenvalues, j)/totalVariance);
					}
					fprintf(stdout, "\n\n");
					fprintf(stdout, "Determining the number of components to keep:\n");
				}
		
				do {
					explainedVariance = explainedVariance + gsl_vector_get(eigenvalues, L);
		
					if ( verbose )
						fprintf(stdout, "% 3d %.10f\n", L, (explainedVariance/totalVariance));
			
					L = L + 1;
				} while( L < nComponents || ((explainedVariance/totalVariance) < minVariance && nComponents < 1) );
	
				if ( verbose )
					fprintf(stdout, "\n");
	
 	
			    L_eigenvectors = gsl_matrix_submatrix(eigenvectors, 0, 0, rows, L);
			    L_eigenvalues  = gsl_vector_subvector(eigenvalues, 0, L);
			}
		}

		#pragma omp sections
		{
			// Project the original dataset
			#pragma omp section
			{
    			result.transformed = gsl_matrix_alloc(L, cols); // transformed is a n LxN matrix, each column is the original data vector with reduced dimension from M to L
	    		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &L_eigenvectors.matrix, data, 0.0, result.transformed);
			}

			// Extract the first L eigenvectors
			#pragma omp section
			{
	 			result.eigenvectors = gsl_matrix_alloc(L, rows);
				gsl_matrix_transpose_memcpy(result.eigenvectors, &(L_eigenvectors.matrix));
			}
	
			// Extract the first L eigenvalues
			#pragma omp section
			{
				result.eigenvalues = gsl_vector_alloc(L);
				gsl_vector_memcpy(result.eigenvalues, &(L_eigenvalues.vector));
			}
		}
	}
	
	// Also store all eigenvalues
	result.eigenvalues_all = eigenvalues;
	
	gsl_matrix_free(eigenvectors);

	return result;
}

struct rsPCAResult rsPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose)
{
	return rsGenericPCA(data, minVariance, nComponents, verbose);
}

struct rsPCAResult rsTPCA(const gsl_matrix* data, double minVariance, int nComponents, BOOL verbose)
{
	gsl_matrix* dataT = gsl_matrix_alloc(data->size2, data->size1);
	gsl_matrix_transpose_memcpy(dataT, data);
	
	struct rsPCAResult result = rsGenericPCA(dataT, minVariance, nComponents, verbose);
	gsl_matrix_free(dataT);
	
	return result;
}

void rsPCAResultFree(struct rsPCAResult result)
{
	gsl_matrix_free(result.transformed);
	gsl_matrix_free(result.eigenvectors);
	gsl_vector_free(result.eigenvalues);
	gsl_vector_free(result.eigenvalues_all);
}

rsCSPResult *rsCSP(const gsl_matrix* A, const gsl_matrix* B, int nComponents, BOOL verbose)
{
	
    assert(A != NULL);
    assert(B != NULL);
	assert(A->size1 == B->size1);

    long i, j;
	int L = 2 * nComponents;
    long rows  = A->size1;
    long colsA = A->size2;
    long colsB = B->size2;
	rsCSPResult *result = rsMalloc(sizeof(rsCSPResult));
	
    gsl_vector* meanA = gsl_vector_alloc(rows);
    gsl_vector* meanB = gsl_vector_alloc(rows);
    gsl_matrix* demeanedA;
    gsl_matrix* demeanedB;
	gsl_matrix* covA;
	gsl_matrix* covB;
	gsl_vector* eigenvalues;
    gsl_matrix* eigenvectors;
    gsl_matrix* L_eigenvectors = gsl_matrix_alloc(rows, L);
    gsl_vector* L_eigenvalues  = gsl_vector_alloc(L);

	if ( verbose ) {
		fprintf(stdout, "Running Common Spatial Patterns Analysis on a %ldx%ld and %ldx%ld matrix\n", rows, colsA, rows, colsB);
	}
 	
	#pragma omp parallel num_threads(rsGetThreadsNum()) shared(colsA,colsB,rows,L,result,nComponents,verbose) private(i,j)
	{	
		// compute means
        #pragma omp for schedule(guided)
       	for(i = 0; i < rows; i=i+1) {
			double mean = 0.0;
			for(j=0; j<colsA; j=j+1) {
				mean = mean + gsl_matrix_get(A,i,j) / colsA;
			}
			gsl_vector_set(meanA, i, mean);
			
			mean = 0.0;
			for(j=0; j<colsB; j=j+1) {
				mean = mean + gsl_matrix_get(B,i,j) / colsB;
			}
        	gsl_vector_set(meanB, i, mean);
		}
 
	    // Get mean-substracted data into matrix demeanedA and demeanedB.
		#pragma omp sections
		{
			#pragma omp section
			{
				demeanedA = gsl_matrix_alloc(rows, colsA);
				gsl_matrix_memcpy(demeanedA, A);
			}
		
			#pragma omp section
			{
				demeanedB = gsl_matrix_alloc(rows, colsB);
	    		gsl_matrix_memcpy(demeanedB, B);
			}
		}
		
        #pragma omp for schedule(guided)
    	for(i = 0; i < colsA; i++) {
        	gsl_vector_view demeanedPointViewA = gsl_matrix_column(demeanedA, i);
        	gsl_vector_sub(&demeanedPointViewA.vector, meanA);
		}
		
        #pragma omp for schedule(guided)
    	for(i = 0; i < colsB; i++) {
        	gsl_vector_view demeanedPointViewB = gsl_matrix_column(demeanedB, i);
        	gsl_vector_sub(&demeanedPointViewB.vector, meanB);
		}
		
 
    	// Compute Covariance matrices
		#pragma omp sections
		{
			#pragma omp section
			{
			    gsl_vector_free(meanA);
				covA = gsl_matrix_alloc(rows, rows);
				gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(colsA - 1), demeanedA, demeanedA, 0.0, covA);
	    		gsl_matrix_free(demeanedA);
			}
	
			#pragma omp section
			{
			    gsl_vector_free(meanB);
				covB = gsl_matrix_alloc(rows, rows);
	    		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(colsB - 1), demeanedB, demeanedB, 0.0, covB);
	    		gsl_matrix_free(demeanedB);
				long evChanged = rsMakePositiveDefiniteSymmetric(covB);
				
				if ( evChanged > 0 ) {
					fprintf(stderr, "Fixed %ld eigenvalue(s) in the covariance matrix of the second mask to ensure it is positive-definite. Consider increasing its size!\n", evChanged);
				}
			}
		}
		
		// Get eigenvectors, sort by eigenvalue.
		#pragma omp sections
		{
			#pragma omp section
			{
				gsl_matrix_add(covA, covB); // covA = covA + covB

			    eigenvalues  = gsl_vector_alloc(rows);
			    eigenvectors = gsl_matrix_alloc(rows, rows);
			    gsl_eigen_gensymmv_workspace* workspace = gsl_eigen_gensymmv_alloc(rows);
			    gsl_eigen_gensymmv(covA, covB, eigenvalues, eigenvectors, workspace);
			    gsl_eigen_gensymmv_free(workspace);
			    gsl_matrix_free(covA);
			    gsl_matrix_free(covB);
 
			    gsl_eigen_gensymmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);
			}
		}
	
		// Extract the first and last nComponents eigenvectors
        #pragma omp for schedule(guided)
		for ( i = 0; i<nComponents; i=i+1 ) {
			const unsigned i2 = L-i-1;
			gsl_vector_view viewCol  = gsl_matrix_column(eigenvectors, i);
			gsl_vector_view viewCol2 = gsl_matrix_column(eigenvectors, rows-i-1);
			gsl_matrix_set_col(L_eigenvectors, i,  &(viewCol.vector));
			gsl_matrix_set_col(L_eigenvectors, i2, &(viewCol2.vector));
			gsl_vector_set(L_eigenvalues, i, gsl_vector_get(eigenvalues, i));
			gsl_vector_set(L_eigenvalues, i2, gsl_vector_get(eigenvalues, rows-i-1));
		}
 	
		// Prepare return value
		#pragma omp sections
		{
			// Project the original dataset
			#pragma omp section
			{
				result->transformedA = gsl_matrix_alloc(L, colsA);
			    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L_eigenvectors, A, 0.0, result->transformedA);
			}
		
			#pragma omp section
			{
	    		result->transformedB = gsl_matrix_alloc(L, colsB);
				gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L_eigenvectors, B, 0.0, result->transformedB);
			}
		
			// Store eigenvectors
			#pragma omp section
			{
				result->eigenvectors = gsl_matrix_alloc(L, rows);
				gsl_matrix_transpose_memcpy(result->eigenvectors, L_eigenvectors);
			}
		
			// Store eigenvalues
			#pragma omp section
			{
				result->eigenvalues     = L_eigenvalues;
				result->eigenvalues_all = eigenvalues;
			}
		}
	}

	return result;
}

void rsCSPResultFree(rsCSPResult *result)
{
	gsl_matrix_free(result->transformedA);
	gsl_matrix_free(result->transformedB);
	gsl_matrix_free(result->eigenvectors);
	gsl_vector_free(result->eigenvalues);
	gsl_vector_free(result->eigenvalues_all);
	rsFree(result);
}