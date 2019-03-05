#!/usr/bin/python2.7

import os
import sys



def main(argv):
    nargs = len(argv)
    infilename = argv[0]

    
    fi=open(infilename,"r")
    
    # load data file
    data_in=fi.readlines()
    data_out=[]
    
    cnt = 0
    insert_cnt = 0
    function_start = 0
    function_header = ''
    include_name = ''
    
    in_function = False
    for line in data_in:
        if ((line.lower().lstrip().find('subroutine')==0) or \
            (line.lower().lstrip().find('function')==0) or \
            (line.lower().lstrip().find('program')==0) or \
            ((line.lower().lstrip().find('function')>0) and \
             (line.lower().lstrip().find('real')==0 or \
              line.lower().lstrip().find('integer')==0 or \
              line.lower().lstrip().find('logical')==0 or \
              line.lower().lstrip().find('character')==0))):
            # found a subroutine/function/program
            if (insert_cnt != 0 and insert_cnt < function_start):
                print "WARNING: Function/Subrotuine without implicit declaration:"+function_header

            function_header = line    
            data_out.append(line)
            cnt = cnt + 1
            function_start = cnt
            in_function = True

        elif (in_function and line.lower().lstrip().find('implicit')==0):
            data_out.append(line)
            insert_cnt = cnt
            cnt = cnt+1
        
        elif (in_function and line.lower().lstrip().find('include')==0):
            # include statement found
            # included file name should be between '' or "" immediately after include
            include_name = ''
            line_data = line.lower().lstrip().split("'")
            if (len(line_data) > 2):
                include_name = line_data[1]
            else:
                line_data = line.lower().lstrip().split('"')
                if (len(line_data) > 2):
                    include_name = line_data[1]
                    

            if (include_name != ''):
                # found an include name
                module_line = '      use mod_'+include_name.strip()+'\n'
                data_out.insert(insert_cnt,module_line)
                cnt = cnt+1
                new_line = 'c'+line[1:]
                data_out.append(new_line)
                cnt=cnt+1
            else:
                print "WARNING: Include statement found with no name available:"+line.strip()

        else:
            data_out.append(line)
            cnt=cnt+1


    fi.close()

    outfilename = infilename+'.out'
    fo=open(outfilename,"w")
    fo.writelines(data_out)
    fo.close()


if __name__ == "__main__":
    main(sys.argv[1:])
