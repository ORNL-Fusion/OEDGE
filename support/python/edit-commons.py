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

    module_name = 'mod_'+infilename.lower()
    outfilename = module_name+'.f90'
    
    cnt = 0
    insert_cnt = 0
    function_start = 0
    function_header = ''
    include_name = ''

    next_line = None
    last_line = None


    # line types - save, common, parameter, declarations
    # common - comment out or delete
    # save - comment out or delete
    # parameter - replicate
    # declarations - change to public, replicate, convert arrays to dynamic
    # include - change to use and put at the beginning before implicit none

    data1 = []
    comment_list=[]
    
    for line in data_in:
        line = line.lower()
        # go through and eliminate continuation lines creating a new set of lines
        # 1- no character in column 1 or 6 - not a comment or continuation - newline
        # 2- c in column 1 - comment
        # 3 character in column 6 - continuation - add to previous line
        if (last_line is None):
            last_line = line.strip()
        elif (len(line)>1):
            if (line[0]!=' '):
                # comment line found
                comment_list.append(line)
            elif (len(line)>6):
                if (line[5]!=' '):
                    # continuation found
                    last_line = last_line + line[6:].strip()
                else:
                    # found a new statement
                    if (len(comment_list)>0):
                        data1.extend(comment_list) 
                    last_line = last_line+'\n'
                    data1.append(last_line)
                    last_line = line.strip()
                    comment_list = []

        else:
            # include blank lines
            data1.append(line)

            
                    
    # add the last_line
    if (len(comment_list)>0):
        data1.extend(comment_list) 
    data1.append(last_line)
        

    # now have the first version with continuation lines removed
    # process into f90 format - comment out common/save/

    data2=[]
    decl_info=[]
    use_list=[]
    
    for line in data1:

        # change comments c to !
        # comment out common, save
        # leave parameter for now
        # change declarations to public,type ::
        # extract array declaration information in every line
        # replace declarations by ':'
        if (len(line.strip())==0):
            # include blank lines from source
            data2.append('\n')
        elif ((line.lstrip().find('real')==0) or \
              (line.lstrip().find('integer')==0) or \
              (line.lstrip().find('logical')==0) or \
              (line.lstrip().find('character')==0) or \
              (line.lstrip().find('double precision')==0) or \
              (line.lstrip().find('complex')==0)):
            data2.append(process_line(line,decl_info))
        elif ((line.lstrip().find('parameter')==0)):
              data2.append(line.strip()+'\n')
        elif ((line.lstrip().find('include')==0)):
            include_name = ''
            line_data = line.lstrip().split("'")
            if (len(line_data) > 2):
                include_name = line_data[1]
            else:
                line_data = line.lstrip().split('"')
                if (len(line_data) > 2):
                    include_name = line_data[1]

            if (include_name != ''):
                # found an include name
                module_line = 'use mod_'+include_name.strip()+'\n'
                use_list.append(module_line)
            else:
                print "WARNING: Include statement found with no name available:"+line.strip()
            # comment out include
            data2.append('! '+line.strip()+'\n')
                
        elif ((line.lstrip().find('common')==0) or \
            (line.lstrip().find('save')==0)):
            data2.append('! '+line.strip()+'\n')
        elif (line.lstrip()[0]=='c'):
            # process comments last since some of the other lines can start with c
            data2.append('!'+line.strip()[1:]+'\n')
        else:
            data2.append(line.strip()+'\n')
            


    # construct fortran 90 output module
    #
    data3=[]
    allocate_data = []
    deallocate_data = []
    indent2 = '  '
    indent4 = '    '


    data3.append('module '+module_name.strip()+'\n')
    data3.append(indent2+'use debug_options\n')

    use_list = [indent2+s for s in use_list]
    data3.extend(use_list)

    data3.append(indent2+'implicit none\n')
    data3.append('\n')

    data2=format_data(data2)
    data2 = [indent2+s for s in data2]    
    data3.extend(data2)

    data3.append('\n')
    data3.append(indent2+'public :: allocate_'+module_name.strip()+',deallocate_'+module_name.strip()+'\n')
    data3.append('\n')
    data3.append('contains\n')
    data3.append('\n')
    data3.append(indent2+'subroutine allocate_'+module_name.strip()+'\n')
    data3.append(indent4+'use mod_params\n')
    data3.append(indent4+'use allocate_arrays\n')
    data3.append(indent4+'implicit none\n')
    data3.append(indent4+'integer :: ierr\n')
    data3.append('\n')
    data3.append(indent4+"call pr_trace('"+module_name.strip()+",'ALLOCATE')\n")
    data3.append('\n')
    
    allocate_data = format_data(allocate_dcl_data(decl_info))
    allocate_data = [indent4+s for s in allocate_data]
    data3.extend(allocate_data)
    
    data3.append('\n')
    data3.append(indent2+'end subroutine allocate_'+module_name.strip()+'\n')
    data3.append('\n')
    data3.append('\n')
    data3.append(indent2+'subroutine deallocate_'+module_name.strip()+'\n')
    data3.append(indent4+'implicit none\n')
    data3.append('\n')
    data3.append(indent4+"call pr_trace('"+module_name.strip()+",'DEALLOCATE')\n")
    data3.append('\n')
    
    deallocate_data = format_data(deallocate_dcl_data(decl_info))
    deallocate_data = [indent4+s for s in deallocate_data]
    data3.extend(deallocate_data)

    data3.append('\n')
    data3.append(indent2+'end subroutine deallocate_'+module_name.strip()+'\n')
    data3.append('\n')
    data3.append('end module '+module_name)

    
    #data3 = format_data(data3)
    

    fi.close()
        
    #outfilename = module_name+'.f90'
    fo=open(outfilename,"w")
    fo.writelines(data3)
    fo.close()


def allocate_dcl_data(decl_info):


    #
    # Need to identify the allocate_array fingerprint to use for different options
    #
    # 1D  - name, dim1, 'name', ierr)
    # 1DB - name, dima1, 'name', dima2, ierr)
    # 2D  - name, dim1, dim2, 'name', ierr)
    # 2DB - name, dima1, dima2, dimb1, dimb2, 'name',ierr
    # 3D  - name, dim1, dim2, dim3, 'name',ierr
    # 3DB - name, dima1, dima2, dimb1, dimb2, dimc1, dimc2, 'name', ierr
    # 4DB - name, dima1, dima2, dimb1, dimb2, dimc1, dimc2, dimd1, dimd2, 'name', ierr
    # 5DB - name, dima1, dima2, dimb1, dimb2, dimc1, dimc2, dimd1, dimd2, dime1, dime2,'name', ierr
    #
    
    allocate_data=[]
    for item in decl_info:
        dim_vals=[]
        name = item[0].strip()
        line = 'call allocate_array('+name.strip()
        ndims = len(item) -1
        for dims in item[1:]:
            pos = dims.strip().find(':')
            if (pos == -1):
                dim_vals.append(['1',dims.strip()])
            else:
                dim_vals.append([dims[:pos].strip(),dims[pos+1:]])
                
        if (len(dim_vals)==1):
            if (dim_vals[0][0].strip()=='1'):
                # 1D
                line += ','+dim_vals[0][1].strip()+",'"+name.strip()+"',ierr)\n"
            else:
                # 1DB
                line += ','+dim_vals[0][0].strip()+",'"+name.strip()+"',"+dim_vals[0][1].strip()+",ierr)\n"
        elif (len(dim_vals)==2):
            if (dim_vals[0][0].strip()=='1' and dim_vals[1][0].strip()=='1'):
                # 2D
                for dim in dim_vals:
                    line += ','+dim[1].strip()
                line += ",'"+name.strip()+"',ierr)\n"
            else:
                # 2DB
                for dim in dim_vals:
                    line += ','+dim[0].strip()+","+dim[1].strip()
                line += ",'"+name.strip()+"',ierr)\n"
        elif(len(dim_vals)==3):
            if (dim_vals[0][0].strip()=='1' and dim_vals[1][0].strip()=='1' and dim_vals[2][0].strip()=='1'):
                # 3D
                for dim in dim_vals:
                    line += ','+dim[1].strip()
                line += ",'"+name.strip()+"',ierr)\n"
            else:
                # 3DB
                for dim in dim_vals:
                    line += ','+dim[0].strip()+","+dim[1].strip()
                line += ",'"+name.strip()+"',ierr)\n"
        else:
            # 4DB, 5DB etc
            for dim in dim_vals:
                line += ','+dim[0].strip()+","+dim[1].strip()
            line += ",'"+name.strip()+"',ierr)\n"

        allocate_data.append(line)

    return allocate_data


   
def deallocate_dcl_data(decl_info):

    deallocate_data=[]
    for item in decl_info:
        name = item[0].strip()
        deallocate_data.append('if (allocated('+name+')) deallocate('+name+')\n')

    return deallocate_data




def breakline(line):
    # add fortran line continuation and formatting for long lines
    lines = []
    text_line = line.strip()
    indent = '     '
    if (len(text_line) < 80):
        lines = [text_line+'\n']
    else:
        # if the line is a long comment then every line broken also needs to be commented
        if (text_line.find('!')==0):
            comment_char='!'
        else:
            comment_char=''
        last_x = 0
        for i in range(0,len(text_line),80):
            if (i==0):
                insert_char = ''
            else:
                insert_char = comment_char+indent
            next_break = min(i+80,len(text_line))
            if (next_break==len(text_line)):
                lines.append(insert_char+text_line[i+last_x:]+'\n')
            else:
                # need to get more complicated about , ... do not want one inside ()
                x = text_line[next_break:].find(',')
                while text_line[next_break+x-1]==':':
                    x = text_line[next_break+x+1:].find(',')+x+1
                x=x+1
                lines.append(insert_char+text_line[i+last_x:next_break+x]+'&\n')
                last_x = x
    return lines                     
            

                
def format_data(data):
    data_out = []
    lines=[]
    for line in data:
        lines = breakline(line)
        data_out.extend(lines)
            
    return data_out
            
def process_line(line,decl_info):

    # Go through line finding array declarations and extracting information about each
    dcl_data=[]
    text_line = line.strip()
    
    start_type = 0
    end_type = text_line.find(' ')

    if (text_line.find('double')==0):
        end_type=(text_line(end_type+1).find(' '))+end_type

    new_line = 'public,'+line[0:end_type]+'::'

    ext_line = text_line[end_type+1:].split(',')

    in_array = False
    for item in ext_line:
        dc_start=item.find('(')
        if (in_array) :
            dc_end = item.find(')')
            if (dc_end!=-1):
                # found end of array
                dcl_item.append(item[0:dc_end].strip())
                dcl_data.append(dcl_item)
                new_line+=dcl_item[0]+'('
                for i in range(len(dcl_item)-2):
                    new_line+=':,'
                new_line+=':),'
                dcl_item=[]
                in_array=False
            else:
                dcl_item.append(item.strip())
                
        elif (dc_start > -1):
            # found start of array declaration
            dc_end = item.find(')')
            name = item[0:dc_start].strip()
            if (dc_end>-1):
                dim = item[dc_start+1:dc_end].strip()
                dcl_data.append([name,dim])
                new_line+=name+'(:),'
            else:
                in_array=True
                dim = item[dc_start+1:].strip()
                dcl_item=[name,dim]
        else:
            new_line+=item.strip()+','


    # remove trailing comma
    new_line=new_line[:len(new_line)-1]+'\n'
    # add extracted declaration data to collection
    decl_info.extend(dcl_data)
            
    return new_line
    
                
        

if __name__ == "__main__":
    main(sys.argv[1:])
