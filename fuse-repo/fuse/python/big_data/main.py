"""
Module documentation.
"""

# Imports
import sys
#import os
from load_data import *

# Global variables

# Class declarations

# Function declarations

def main():
    args = sys.argv[1:]

    if not args:
        print('usage: [--flags options] [inputs] ')
        sys.exit(1)

    else:
        print("shit")

        load_index_5()
        
# Main body
if __name__ == '__main__':
    main()


