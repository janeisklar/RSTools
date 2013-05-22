/*******************************************************************
 *
 * rstimecourse.c
 *
 * Extracts the time course from a 4D-Nifti by supplying either a voxel or a binary mask
 *
 * Usage: rstimecourse [-m <mask> -a mean] <volume> [<X> <Y> <Z>]
 *
 * 
 * André Hoffmann
 *******************************************************************/


#include <stdio.h>
#include <strings.h>

#include <nifti1.h>
#include <fslio.h>

int show_help( void )
{
   printf(
      "rstimecourse: Given a 4D-Nifti, this tool extracts the time course\n"
      "              for a single voxel or the meaned average of a region\n"
      "              specified by a binary mask\n"
      "\n"
      "    basic usage: rstimecourse [-m <mask> [-a <algorithm>]] [-p <X> <Y> <Z>] -input <volume>\n"
      "\n"
      "    options:     -help           : show this help\n"
      "                 -mask <mask>    : a mask specifying the ROI\n"
      "                 -a <algorithm>  : the algorithm used to aggregate the data within\n"
      "                                   a ROI, e.g. mean\n"
      "                 -p <X> <Y> <Z>  : speficies a voxel using nifti coordinates(0-based) from\n"
      "                                   which the timecourse is to be extracted\n"
      "\n");
   return 0;
}

int main(int argc, char * argv[])
{
	FSLIO *fslio;
	void *buffer;
	size_t buffsize;
	
	char *inputpath = NULL;
	char *maskpath = NULL;
	
	int x=-1, y=-1, z=-1, t=0;
	short xDim, yDim, zDim, vDim;
	short pixtype;
	size_t dt;
	
	int ac;

	if( argc < 2 ) return show_help();   /* typing '-help' is sooo much work */

	/* parse parameters */
	for( ac = 1; ac < argc; ac++ ) {
		if( ! strncmp(argv[ac], "-h", 2) ) {
			return show_help();
		} else if ( ! strcmp(argv[ac], "-input") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -input\n");
				return 1;
			}
			inputpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-m") ) {
			if( ++ac >= argc ) {
				fprintf(stderr, "** missing argument for -m\n");
				return 1;
			}
			maskpath = argv[ac];  /* no string copy, just pointer assignment */
		} else if ( ! strcmp(argv[ac], "-p") ) {
			if( ac+3 >= argc ) {
				fprintf(stderr, "** missing argument for -p, 3 coordinates must be supplied!\n");
				return 1;
			}
			ac++;
			x = atoi(argv[ac]);
			ac++;
			y = atoi(argv[ac]);
			ac++;
			z = atoi(argv[ac]);
		} else {
			fprintf(stderr, "\nError, unrecognized command %s\n",argv[ac]);
		}
	}
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified!\n");
		return 1;
	}
	
	if ( maskpath == NULL && (x<0||y<0||z<0) ) {
		fprintf(stderr, "Either a binary mask or a voxel coordinate must be specified!\n");
		return 1;
	}
	
	fprintf(stdout, "Input file: %s\n", inputpath);
	fprintf(stdout, "Mask file: %s\n", maskpath);
	fprintf(stdout, "Voxel: %d %d %d\n", x, y, z);
	
	/* open nifti dataset header */
	fslio = FslReadHeader(inputpath);
	if (fslio == NULL) {
		fprintf(stderr, "\nError, could not read header info for %s.\n",inputpath);
		return 1;
	}

	/* check inputs */
	if ( (x<0) || (x>=fslio->niftiptr->nx) ) {
		fprintf(stderr, "\nError: x index (%d) out of range [0..%d]\n",x,fslio->niftiptr->nx-1);
		return 1;
	}
	if ( (y<0) || (y>=fslio->niftiptr->ny) ) {
		fprintf(stderr, "\nError: y index (%d) out of range [0..%d]\n",y,fslio->niftiptr->ny-1);
		return 1;
	}
	if ( (z<0) || (z>=fslio->niftiptr->nz) ) {
		fprintf(stderr, "\nError: z index (%d) out of range [0..%d]\n",z,fslio->niftiptr->nz-1);
		return 1;
	}
	
	/* determine dimensions */
	FslGetDim(fslio, &xDim, &yDim, &zDim, &vDim);
	
	fprintf(stdout, "Dim: %d %d %d (%d Volumes)\n", xDim, yDim, zDim, vDim);
	
	/* determine datatype and initalize buffer */
	dt = FslGetDataType(fslio, &pixtype);
	
	fprintf(stdout, "Dt: %ld Pixtype: %d\n", dt, pixtype);

	
	buffsize = vDim*dt/8;
	
	fprintf(stdout, "Buffsize: %ld\n", buffsize);
	
    buffer = malloc(buffsize);
	
	/* read out timecourse */
	FslReadTimeSeries(fslio, buffer, x, y, z, vDim);
	
	for (t=0; t<vDim; t++) {
		int v = (char*)buffer + (t*dt);
		fprintf(stdout, "%d: %d\n", t, v);
	}
	
	/* clear buffer */
    free(buffer);
	
	return 0;
}
