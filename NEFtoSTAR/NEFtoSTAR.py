'''
Created on Sep 15, 2016

@author: kumaran
'''

import sys,csv,ntpath
sys.path.append('../PyNMRSTAR') #NMR-STAR and NEF-Parser added as a submodule and imported into this project. This is a separate git repository
import bmrb

class NEFtoSTAR(object):
    '''
    Translates NEF file to NMR-STAR file
    '''
    #mapping file btween NEF and NMR-STAR
    mapFile='NEF_NMRSTAR_equivalence.csv'
    #Standard atom names for 20 standard amino acids
    NMR_STAR_atom_names={'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG', 'H', 'HA', 'HB2', 'HB3', 'HG'],
                         'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2', 'H', 'HA', 'HB2', 'HB3', 'HD2'],
                         'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG', 'H', 'HA', 'HB2', 'HB3', 'HG'],
                         'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE21', 'HE22'],
                         'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE2', 'HE3', 'HZ1', 'HZ2', 'HZ3'],
                         'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'H', 'HA', 'HB', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HD11', 'HD12', 'HD13'],
                         'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3'],
                         'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'H', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23'],
                         'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
                         'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', 'H', 'HA', 'HB2', 'HB3', 'HD21', 'HD22'],
                         'GLY': ['N', 'CA', 'C', 'O', 'H', 'HA2', 'HA3'],
                         'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2'],
                         'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'H', 'HA', 'HB2', 'HB3', 'HG', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'],
                         'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE', 'HH11', 'HH12', 'HH21', 'HH22'],
                         'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HE1', 'HE3', 'HZ2', 'HZ3', 'HH2'],
                         'ALA': ['N', 'CA', 'C', 'O', 'CB', 'H1', 'HA', 'HB1', 'HB2', 'HB3', 'H2', 'H3'],
                         'VAL': ['N', 'CA', 'C', "O'", 'CB', 'CG1', 'CG2', 'H', 'HA', 'HB', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', "O''"],
                         'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE2'],
                         'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
                         'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE1', 'HE2', 'HE3']}
    
    def read_map_file(self):
        '''Reads the NEF_NMRSTAR_equivalence.csv file and create a mapping as a list'''
        with open(self.mapFile,'rb') as csvfile:
            spamreader = csv.reader(csvfile,delimiter=',')
            map_dat=[]
            for r in spamreader:
                if r.count('') != 10:
                    map_dat.append(r)
        self.map=map(list,zip(*map_dat))
    
    def translate(self):
        self.nefData=bmrb.Entry.from_file(self.nefFile)                     #parse the NEF file
        self.starData=bmrb.Entry.from_scratch(self.nefData.entry_id)        # creates an empty NMR-STAR data structure
        for saveframe in self.nefData:                                      # for every saveframe
            sf=bmrb.Saveframe.from_scratch(saveframe.name)                  # create an equivalent NMR-STAR saveframe 
            for tag in saveframe.tags:                                      # for every tag in the saveframe 
                nef_tag="%s.%s"%(saveframe.tag_prefix,tag[0])
                star_tag=self.map[1][self.map[0].index(nef_tag)]            # get the equivalent tag from the mapping 
                if star_tag!="":
                    if star_tag.split(".")[1]=="Sf_category":
                        ssf_category=self.map[1][self.map[0].index(tag[1])] # Only for saveframe category find the corresponding NMR-STAR category
                        sf.add_tag(star_tag,ssf_category)                   # add the tag value
                    else:
                        sf.add_tag(star_tag,tag[1])                         # for rest of the tags simply copy the tag value
            for loop in saveframe:                                          # for every loop in the saveframve
                ll=bmrb.Loop.from_scratch()                                 # create a NMR-STAR loop from the scratch 
                missing_col=[]                                              # to collect the information about missing columns in the loop
                for coln in loop.columns:                                   # for every columnn in the loop
                    nl_tag="%s.%s"%(loop.category,coln)                     # extract the column name tag
                    sl_tag=self.map[1][self.map[0].index(nl_tag)]           # find the equivalent NMR-STAR column name tag
                    if sl_tag!="":                                          # If there is no mapping add the column index to missing column 
                        ll.add_column(sl_tag)                               # otherwise add column to loop
                    else:
                        missing_col.append(loop.columns.index(coln))
                if len(missing_col)==0:                                     # If there is no missing column copy the data
                    for dat in loop.data:
                        ll.add_data(dat[:])
                sf.add_loop(ll)                                             # add the loop to saveframe
            self.starData.add_saveframe(sf)                                 # add the saveframe to data structure
        (file_path,file_name)=ntpath.split(self.nefFile)                        # write outout file
        if file_path=="": file_path="."
        out_file_name="%s.str"%(file_name.split(".nef")[0])
        outfile=file_path+"/"+out_file_name
        with open(outfile,'w') as strfile:
            strfile.write(str(self.starData))
    

    def __init__(self, nefFile):
        '''
        Constructor
        '''
        self.nefFile=nefFile
        self.read_map_file()
        
if __name__=="__main__":
    nt=NEFtoSTAR('/home/kumaran/git/NEF/data/CCPN_1nk2_docr.nef')
    nt.translate()