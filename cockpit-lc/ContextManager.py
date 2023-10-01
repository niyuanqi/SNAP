#################################################################
# Name:     ContextManager.py                                   #
# Author:   Yuan Qi Ni                                          #
# Date:     May, 24, 2016                                       #
# Function: Program allows temporary change of directories.     #
#################################################################

import os

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
