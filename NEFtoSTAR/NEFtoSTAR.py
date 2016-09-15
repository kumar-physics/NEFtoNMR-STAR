'''
Created on Sep 15, 2016

@author: kumaran
'''

import sys
sys.path.append('../PyNMRSTAR')
import bmrb

class NEFtoSTAR(object):
    '''
    classdocs
    '''


    def __init__(self, nefFile):
        '''
        Constructor
        '''
        print nefFile
        
if __name__=="__main__":
    nt=NEFtoSTAR('/home/kumaran/git/NEF/data/CCPN_1nk2_docr.nef')