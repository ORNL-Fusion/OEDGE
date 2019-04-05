#!/usr/bin/python2.7

import os
import sys



def main(argv):
    nargs = len(argv)
    infilename = argv[0]

    print "Processing:",infilename

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
            # verify by checking that the character immediately before 'function' is a space
            pos = line.lower().lstrip().find('function')
            if (pos <= 0 or (pos > 0 and line.lower().lstrip()[pos-1]==' ')):
                # found a subroutine/function/program
                if (insert_cnt < function_start):
                    print "WARNING: Function/Subroutine without implicit declaration:"+function_header

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
                    
            if (insert_cnt < function_start):
                print "WARNING: Function/Subroutine without implicit declaration:"+function_header

            if (include_name != ''):
                # found an include name
                module_line = '      use mod_'+include_name.strip()+'\n'
                data_out.insert(insert_cnt,module_line)
                # move insert point down so that the modules appear in matching order to includes
                insert_cnt = insert_cnt + 1
                cnt = cnt+1
                new_line = 'c'+line[1:]
                data_out.append(new_line)
                cnt=cnt+1
            else:
                print "WARNING: Include statement found with no name available:"+line.strip()

        else:
            data_out.append(line)
            cnt=cnt+1


    # catch the case where the last function in a module doesn't have an implicit none statement
    if (insert_cnt < function_start):
        print "WARNING: Function/Subroutine without implicit declaration:"+function_header

    fi.close()

    outfilename = infilename+'.out'
    fo=open(outfilename,"w")
    fo.writelines(data_out)
    fo.close()


if __name__ == "__main__":
    main(sys.argv[1:])
