#!/home/david/mambaforge/envs/omfit/bin/python3

import os
import sys
from oedge_input import OMFIToedgeInput

outfile_init = 'init_code.f90'
outfile_read = 'read_code.f90'
outfile_doc  = 'unstruc_doc.txt'


#
# Code processes a file with input specifications into unstructured input format
# using divrd and basic default values
#
# Sort tags?
#
#

def main(argv):
    nargs = len(argv)
    indata = argv[0]
    incode = argv[1]
    #print("test:",nargs)

    print("Processing code file:",incode)
    print("Processing data file:",indata)

    outread = 'read-unstruc.f90'
    outinit = 'init-unstruc.f90'
    indocs = ['mod_unstructured_input.f90','mod_sol23_input.f90','mod_sol22_input.f90','ero_interface.f90']
    outdoc  = 'docs-unstruc.txt'
    
    #fic=open(incode,"r")
    #fid=open(indata,"r")
    
    # load data file
    #fc=fic.readlines()
    #fd=fid.readlines()

    # output files
    # output read code = orc
    # output initialization code = oic
    orc=[]
    oic=[]

    
    
    # check to see if specified case exists
    # assume that if the file exists it is an OEDGE input file
    if os.path.isfile(indata):

        casename, filetype = os.path.splitext(indata.split(os.sep)[-1])

        print('Reading:',casename,filetype)
        # load the case and record where it is being stored in the tree
        case = OMFIToedgeInput(indata, casename=casename)
        case.load_code(incode)

    oic=[]

    
    
    # check to see if specified case exists
    # assume that if the file exists it is an OEDGE input file
    if os.path.isfile(indata):

        casename, filetype = os.path.splitext(indata.split(os.sep)[-1])

        print('Reading:',casename,filetype)
        # load the case and record where it is being stored in the tree
        case = OMFIToedgeInput(indata, casename=casename)
        case.load_code(incode)
        case.load_doc(indocs)
        
        case.create_read_code(outread)

        case.create_init_code(outinit)

        case.create_doc(outdoc)
        
        #cnt = 0
        #for item in list(case.keys()):
        #    cnt = cnt+1
        #    print("Cnt:",cnt,"TAG:",item)

        #cnt = 0    
        #for item in list(case.keys()):
        #    cnt = cnt+1
        #    if isinstance(case[item],str):
        #        print("Cnt:",cnt,"TAG:",case[item],isinstance(case[item],dict))
        #    else:
        #        print("Cnt:",cnt,"TAG:",case[item].keys(),isinstance(case[item],dict))

            

if __name__ == "__main__":
    main(sys.argv[1:])
    
