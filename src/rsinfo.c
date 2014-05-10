//
//  rsroi.c
//  rstools
//
//  Created by André Hoffmann on 6/7/13.
//
//
#include <stdio.h>
#include <nifti1.h>
#include <fslio.h>
#include "src/nifti/rsniftiutils.h"
#include "src/maths/rsmathutils.h"

int show_help( void )
{
    printf(
	   RSTOOLS_VERSION_LABEL "\n"
       "rsinfo: Given a Nifti that has been created by any of the rstools,\n"
       "        this tool reads out the header information it entails.\n"
       "\n"
    );
    
    printf(
       "basic usage:  rsinfo <file>\n"
       "\n"
    );
    
    return 0;
}

int main(int argc, char * argv[])
{
    FSLIO *fslio;
    
	if( argc != 2 ) {
		return show_help();
	}
	
	// get filename
	char *inputpath = argv[1];
	
	if ( inputpath == NULL ) {
		fprintf(stderr, "No input volume specified!\n");
		return 1;
	}
	
	// open file
    fslio = FslOpen(inputpath, "rb");
    if (fslio == NULL) {
        fprintf(stderr, "\nError, could not read header info for %s.\n",inputpath);
        return 1;
    }
    
	// check extensions
	nifti_image *nim = fslio->niftiptr;
   
	if( nim->num_ext <= 0 || nim->ext_list == NULL ){
	    fprintf(stderr, "File does not contain any RSTools header information.\n");
		return 0;
	}
   	
	// find extension
	nifti1_extension *ext = nim->ext_list;
	
	int c;
	for ( c = 0; c < nim->num_ext; c++ ){
		if ( ext->ecode == NIFTI_ECODE_COMMENT && ext->edata != NULL ) { 
			break;
		}
	    ext++;
	}
	
	if ( c >= nim->num_ext ) {
		fprintf(stderr, "File does not contain any RSTools header information.\n");
		return 0;
	}
	
	// read extension
	int size = ext->esize;
	char data[size+1];
	strncpy(data, ext->edata, size);
    FslClose(fslio);

	// print info
	fprintf(stdout, "Header information:\n%s\n", data);

	return 0;
}