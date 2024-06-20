"""
Module documentation.

https://python4astronomers.github.io/files/asciifiles.html

"""

# Imports
import sys
#import os
from process_data import *

# Global variables

# Class declarations

# Function declarations

# ----------------------------------------------------------------------

def parse_buffer(buffer):

    print(buffer)

    print('splitting')    
    buffer = buffer.strip()

    columns = buffer.split()

    print(columns[0:8])
    print(columns[1])

    print(len(columns))

    file_path = ''.join(columns[8:len(columns)])
    print(file_path)

    file_tree = file_path.split("/")
    print(file_tree)

    file_name = file_tree[len(file_tree)-1]
    print(file_name)

    file_ext = file_name.split(".")
    file_ext = file_ext[len(file_ext)-1]
    print(file_ext)



    
    result = process_data(buffer)
    print(result)
    

# ----------------------------------------------------------------------

def load_index_5():

    print("Trying")
    
    f = open('./index_5.txt', 'r')

    print("Loaded")
    
#    data = f.read()   

    for buffer in f:
        print(repr(buffer))

    print("parsing")
        
    parse_buffer(buffer)



    f.close()        

# ----------------------------------------------------------------------

def load_index_5_old():

    print("Trying")
    
    f = open('./index_5.txt', 'r')

    print("Loaded")
    
#    data = f.read()   

    for line in f:
        print(repr(line))

    f.close()        

    print('closed')
    print(line)

    print('splitting')    
    line = line.strip()

    columns = line.split()

    print(columns[0:8])
    print(columns[1])

    print(len(columns))

    file_path = ''.join(columns[8:len(columns)])
    print(file_path)

    file_tree = file_path.split("/")
    print(file_tree)

    file_name = file_tree[len(file_tree)-1]
    print(file_name)

    file_ext = file_name.split(".")
    file_ext = file_ext[len(file_ext)-1]
    print(file_ext)

#    file_ext = file_tree[len(file_tree)
    
#    file_name_elements = file_name_tree[len(file_name.split(".")
#    print(file_name_elements)
    
#    print(data)




 
