'''
Created on Sep 15, 2016

@author: kumaran
'''

import sys,csv,ntpath,re
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
                         'ALA': ['N', 'CA', 'C', 'O', 'CB', 'H', 'HA', 'HB1', 'HB2', 'HB3'],
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
        if ".nef" in self.inFile:                                          #lookup index(li) retrieve index(ri) 0-nef 1-star-auth-tag,2-star-tag 
            (li,ri)=(0,1)                                                   # for net to star li=0 and ri=1 and for star to nef li=1 and ri=0
        elif ".str" in self.inFile:
            (li,ri)=(1,0)
        else:
            print "Invalid input file"
            exit(1)
        self.inData=bmrb.Entry.from_file(self.inFile)                     #parse the NEF file
        self.outData=bmrb.Entry.from_scratch(self.inData.entry_id)        # creates an empty NMR-STAR data structure
        for saveframe in self.inData:                                      # for every saveframe
            sf=bmrb.Saveframe.from_scratch(saveframe.name)                  # create an equivalent NMR-STAR saveframe 
            for tag in saveframe.tags:                                      # for every tag in the saveframe 
                in_tag="%s.%s"%(saveframe.tag_prefix,tag[0])
                out_tag=self.map[ri][self.map[li].index(in_tag)]            # get the equivalent tag from the mapping 
                if out_tag!="":
                    if out_tag.split(".")[1]=="Sf_category" or out_tag.split(".")[1]=="sf_category":
                        out_category=self.map[ri][self.map[li].index(tag[1])] # Only for saveframe category find the corresponding NMR-STAR category
                        sf.add_tag(out_tag,out_category)                   # add the tag value
                    else:
                        sf.add_tag(out_tag,tag[1])                         # for rest of the tags simply copy the tag value
            for loop in saveframe:                                          # for every loop in the saveframve
                ll=bmrb.Loop.from_scratch()                                 # create a NMR-STAR loop from the scratch 
                missing_col=[]                                              # to collect the information about missing columns in the loop
                for coln in loop.columns:                                   # for every columnn in the loop
                    inl_tag="%s.%s"%(loop.category,coln)                     # extract the column name tag
                    outl_tag=self.map[ri][self.map[li].index(inl_tag)]           # find the equivalent NMR-STAR column name tag
                    if outl_tag!="":                                          # If there is no mapping add the column index to missing column 
                        ll.add_column(outl_tag)                               # otherwise add column to loop
                    else:
                        missing_col.append(loop.columns.index(coln))
                if len(missing_col)==0:                                     # If there is no missing column copy the data
                    for dat in loop.data:
                        ll.add_data(dat[:])
                if ll.data: sf.add_loop(ll)                                             # add the loop to saveframe
            self.outData.add_saveframe(sf)                                 # add the saveframe to data structure
        (file_path,file_name)=ntpath.split(self.inFile)                        # write output file
        if file_path=="": file_path="."
        if ".nef" in file_name:
            out_file_name="%s_.str"%(file_name.split(".nef")[0])
        else:
            out_file_name="%s_.nef"%(file_name.split(".str")[0])
        outfile=file_path+"/"+out_file_name
        with open(outfile,'w') as strfile:
            strfile.write(str(self.outData))
    
    def convert(self):
        if ".nef" in self.inFile:                                          #lookup index(li) retrieve index(ri) 0-nef 1-star-auth-tag,2-star-tag 
            (li,ri)=(0,1)
            self.nef2star=True                                                   # for net to star li=0 and ri=1 and for star to nef li=1 and ri=0
        elif ".str" in self.inFile:
            (li,ri)=(1,0)
            self.nef2star=False
        else:
            print "Invalid input file"
            exit(1)
        self.inData=bmrb.Entry.from_file(self.inFile)                     #parse the NEF file
        self.outData=bmrb.Entry.from_scratch(self.inData.entry_id)        # creates an empty NMR-STAR data structure
        for saveframe in self.inData:                                      # for every saveframe
            sf=bmrb.Saveframe.from_scratch(saveframe.name)                  # create an equivalent NMR-STAR saveframe 
            for tag in saveframe.tags:                                      # for every tag in the saveframe 
                map_auth_tag=False
                in_tag="%s.%s"%(saveframe.tag_prefix,tag[0])
                out_tag=self.map[ri][self.map[li].index(in_tag)]            # get the equivalent tag from the mapping 
                outs_tag=self.map[2][self.map[li].index(in_tag)]
                if out_tag!=outs_tag : map_auth_tag=True
                if out_tag!="":
                    if out_tag.split(".")[1]=="Sf_category" or out_tag.split(".")[1]=="sf_category":
                        out_category=self.map[ri][self.map[li].index(tag[1])] # Only for saveframe category find the corresponding NMR-STAR category
                        sf.add_tag(out_tag,out_category)                   # add the tag value
                    else:
                        sf.add_tag(out_tag,tag[1])                         # for rest of the tags simply copy the tag value
                        if map_auth_tag : sf.add_tag(outs_tag,tag[1])
            for loop in saveframe:                                          # for every loop in the saveframve
                ll=bmrb.Loop.from_scratch()                                 # create a NMR-STAR loop from the scratch 
                missing_col=[]
                star_col=[]                                              # to collect the information about missing columns in the loop
                for coln in loop.columns:
                    map_authl_tag=False                                   # for every columnn in the loop
                    inl_tag="%s.%s"%(loop.category,coln)                     # extract the column name tag
                    outl_tag=self.map[ri][self.map[li].index(inl_tag)]           # find the equivalent NMR-STAR column name tag
                    star_outl_tag=self.map[2][self.map[li].index(inl_tag)]
                    if outl_tag != star_outl_tag : map_authl_tag=True
                    if outl_tag!="":                                          # If there is no mapping add the column index to missing column 
                        ll.add_column(outl_tag)
                        if map_authl_tag:
                            ll.add_column(star_outl_tag)
                            star_col.append(loop.columns.index(coln))           # otherwise add column to loop
                    else:
                        missing_col.append(loop.columns.index(coln))
                if len(missing_col)==0:
                    al=[]
                    rl=[]
                    for k in star_col:
                        if "atom_name" in loop.columns[k]:
                            al.append(k)
                            rl.append([loop.columns.index(x) for x in loop.columns  if loop.columns[k].replace('atom_name','residue_type')==x][0] )   
                    print al,rl,star_col  ,loop.columns                             # If there is no missing column copy the data
                    for dat in loop.data:
                        if len(star_col)!=0:
                            dat2=dat[:]
                            if len(al)>0:
                                for (r,a) in zip(rl,al):
                                    dat2=dat[:]
                                    for k in star_col:
                                        if a==k and dat[:][r] in self.NMR_STAR_atom_names.keys():
                                            alist=self.get_atm_list(dat[:][r], dat[:][a])
                                            print alist
                                            if len(alist)>1:
                                                for at in alist:
                                                    dat2.insert(dat2.index(dat[:][k]),at)
                                            else:
                                                dat2.insert(dat2.index(dat[:][k]),dat[:][k])
                                        else:
                                            dat2.insert(dat2.index(dat[:][k]),dat[:][k])
                                        print dat2[:]
                                        ll.add_data(dat2[:])
                            else:
                                for k in star_col:
                                    dat2.insert(dat2.index(dat[:][k]),dat[:][k])
                        else:
                            print dat2[:]
                            ll.add_data(dat[:])
                if ll.data: sf.add_loop(ll)                                             # add the loop to saveframe
            self.outData.add_saveframe(sf)                                 # add the saveframe to data structure
        (file_path,file_name)=ntpath.split(self.inFile)                        # write output file
        if file_path=="": file_path="."
        if ".nef" in file_name:
            out_file_name="%s_.str"%(file_name.split(".nef")[0])
        else:
            out_file_name="%s_.nef"%(file_name.split(".str")[0])
        outfile=file_path+"/"+out_file_name
        with open(outfile,'w') as strfile:
            strfile.write(str(self.outData))
                        
    def get_atm_list(self,res,nefAtom):
        try:
            atms=self.NMR_STAR_atom_names[res]
            alist=[]
            try:
                refatm=re.findall(r'(\S+)([XY])([%*])|(\S+)([%*])|(\S+)([XY])',nefAtom)[0]  
                set=[refatm.index(i) for i in refatm if i!=""]
                if set==[0,1,2]:
                    pattern=re.compile(r'%s\d\d+'%(refatm[0]))
                    alist=[i for i in atms if re.search(pattern, i)]
                    if refatm[1]=="Y":
                        alist.reverse()
                elif set==[3,4]:
                    if refatm[4]=="%":
                        pattern=re.compile(r'%s\d+'%(refatm[3]))
                    elif refatm[4]=="*":
                        pattern=re.compile(r'%s\S+'%(refatm[3]))
                    else:
                        print "something wrong"
                    alist=[i for i in atms if re.search(pattern, i)]
                    
                elif set==[5,6]:
                    pattern=re.compile(r'%s\d+'%(refatm[5]))
                    alist=[i for i in atms if re.search(pattern, i)]
                    if len(alist)!=2:
                        alist=[]
                    elif refatm[6]=="Y":
                        alist.reverse()
                    elif refatm[6]=="X":
                        alist
                    else:
                        print "Something wrong"
                        
                else:
                    print "Wrong regular expression"
            except IndexError:
                #print nefAtom
                pass
            if len(alist)==0:
                if nefAtom in atms:
                    alist.append(nefAtom)
        except KeyError:
            print "Residue not found"
            alist=[]
        return alist                    
                    
                
    
    def test(self):
        for i in range(len(self.map[0])):
            if self.map[1][i]!=self.map[2][i]:
                print self.map[0][i],self.map[1][i],self.map[2][i]

    def __init__(self, inFile):
        '''
        Constructor
        '''
        self.inFile=inFile
        self.read_map_file()
        
if __name__=="__main__":
    nt=NEFtoSTAR('/home/kumaran/git/NEF/data/CCPN_1nk2_docr.nef')
    #nt.translate()
    nt.convert()