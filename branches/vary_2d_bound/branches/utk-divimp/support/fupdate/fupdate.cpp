/*  FUPDATE: This program is a support utility program for the DIVIMP
             monte carlo impurity transport code. The purpose of this
             program is to simplify updating input files that are 
             associated with earlier DIVIMP versions to a format that
             is compatible with the current version. It uses as a reference
             any file template that correctly works with the version of 
             DIVIMP at which the update is targetted. 

             The program requires three inputs - all are file names.

             1) The name of the reference input file
             2) The name of the old file that will be updated
             3) The name of the file where the results will be stored

             The format of the command line is:

             fupdate <filename#1> <filename#2> <filename#3>

             If an incorrect number of arguments is specified then 
             the code issues an error-message and briefly describes the
             format of the expected arguments.

             David Elder,   May 2, 1997                                   */


#include <iostream.h>
#include <strstream.h>
#include <fstream.h> 
#include <String.h>
#include <stdlib.h>
#include <math.h>


#include "Local_Types.h"
#include "DataItem.h"
#include "DataFile.h"


void main(int argc, char *argv[])

	{ 
	DataFile new_file,old_file;
        DataFile* output_file;
        String filename;

        cout << "START COMPARISON:\n";

	// Check to see if there are sufficient arguments

        if (argc!=4) 

	  {
	    cout << "FUPDATE ERROR!: Improper Argument Count = " << argc-1 << endl;
	    cout << " - Three file name arguments must be supplied" << endl;
	    cout << " - <template file> <file to update> <name for output file>" << endl;
	    cout << " - e.g. fupdate template5.05.d5i jet719ae.d4i jet719af.d5i " << endl;
	    return(1);
	  };

        // Echo file names to output 

        cout << "Template    File (Input) : " << argv[1] << endl;
        cout << "Old Input   File (Input) : " << argv[2] << endl;
        cout << "New Updated File (Output): " << argv[3] << endl;

        // Read in and analyse the reference file

        new_file.read(argv[1]);
        new_file.split_file();

        // Read in and analyse the old file to be updated

        old_file.read(argv[2]);
        old_file.split_file();

        // Update the old file to the reference file-> result is in output_file

        output_file = new_file.update(old_file);

        // Print out the updated file  

        output_file->write(argv[3]);

        cout << "FINISHED UPDATE" << endl ;
}





