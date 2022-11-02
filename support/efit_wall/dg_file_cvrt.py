#!/usr/bin/python3

import os
import sys
import re


def main(argv):
    nargs = len(argv)

    if (nargs != 2):
        print("Script requires two arguments: <filename> <scaling factor> ",nargs)
        print("Arguments: ",argv[0],argv[1])

        quit()
        
    infilename = argv[0]
    cvrtf = argv[1]
    cvrt = float(cvrtf)
    outfilename= infilename + ".meters"
    
    print("Processing:",infilename)
    
    fi=open(infilename,"r")
    
    # load data file
    data_in=fi.readlines()
    data_out=[]
    
    cnt = 0
    
    for line in data_in:
        line_data = re.split(' ,',line.lstrip())
        #print("line_data:",line_data)
        if (len(line_data) == 2):
            cnt = cnt + 1
            data_out.append(f'{float(line_data[0].strip())*cvrt:.5f}' + " " +f'{float(line_data[1].strip())*cvrt:.5f}'+"\n") 
        elif (len(line_data) == 3):
            cnt = cnt + 1
            data_out.append(f'{float(line_data[0].strip())*cvrt:.5f}' + " " +f'{float(line_data[2].strip())*cvrt:.5f}'+"\n") 

    fi.close()

    data_out.insert(0, "NEUTRAL WALL: "+str(cnt)+"  0   1.0 \n")
    
    outfilename = infilename+'.out'
    fo=open(outfilename,"w")
    fo.writelines(data_out)
    fo.close()


if __name__ == "__main__":
    main(sys.argv[1:])
