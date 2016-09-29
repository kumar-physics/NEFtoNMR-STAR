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
                         'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE1', 'HE2', 'HE3'],
                         'A': ["C3'", "O3'", "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4', "O2'", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", 'H8', 'H61', 'H62', 'H2', "HO5'"],
                         'C': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N1', 'C6', 'C5', 'C4', 'N4', 'N3', 'C2', 'O2', "O2'", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", 'H41', 'H42', 'H5', 'H6'],
                         'G': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4', "O2'", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", 'H8', 'H1', 'H21', 'H22'],
                         'U': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N1', 'C6', 'C5', 'C4', 'O4', 'N3', 'C2', 'O2', "O2'", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", 'H3', 'H5', 'H6', "HO3'"],
                         'DA': ["C3'", "O3'", "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4', "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", 'H8', 'H61', 'H62', 'H2', "HO5'"],
                         'DC': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N1', 'C6', 'C5', 'C4', 'N4', 'N3', 'C2', 'O2', "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", 'H41', 'H42', 'H5', 'H6'],
                         'DG': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4', "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", 'H8', 'H1', 'H21', 'H22'],
                         'DT': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N1', 'C6', 'C5', 'C4', 'O4', 'N3', 'C2', 'O2', 'C7', "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", 'H3', 'H71', 'H72', 'H73', 'H6', "HO3'"]}
    
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
        if ".nef" in self.inFile:
            self.nef2star=True
        elif ".str" in self.inFile:
            self.nef2star=False
        else:
            print "Unsupported input file or file extension may be wrong"
            exit(1)
        self.inData=bmrb.Entry.from_file(self.inFile)
        self.outData=bmrb.Entry.from_scratch(self.inData.entry_id)
        if self.nef2star:
            for saveframe in self.inData:
                sf=bmrb.Saveframe.from_scratch(saveframe.name)
                for tag in saveframe.tags:
                    in_tag="%s.%s"%(saveframe.tag_prefix,tag[0])
                    out_auth_tag=self.map[1][self.map[0].index(in_tag)]
                    out_star_tag=self.map[2][self.map[0].index(in_tag)]
                    if out_auth_tag!="" and out_star_tag!="":
                        if out_auth_tag==out_star_tag:
                            if out_star_tag.split(".")[1]=="Sf_category" or out_star_tag.split(".")[1]=="sf_category":
                                out_category=self.map[1][self.map[0].index(tag[1])]
                                sf.add_tag(out_star_tag,out_category)
                            else:
                                sf.add_tag(out_star_tag,tag[1])
                        else:
                            sf.add_tag(out_auth_tag,tag[1])
                            sf.add_tag(out_star_tag,tag[1])
                            print out_auth_tag,out_star_tag
                for loop in saveframe:
                    ll=bmrb.Loop.from_scratch()
                    missing_col=[]
                    auth_col=[]
                    for coln in loop.columns:
                        inl_tag="%s.%s"%(loop.category,coln)
                        outl_auth_tag=self.map[1][self.map[0].index(inl_tag)]
                        outl_star_tag=self.map[2][self.map[0].index(inl_tag)]
                        if outl_auth_tag!="" and outl_star_tag!="":
                            if outl_auth_tag!=outl_star_tag:
                                auth_col.append(loop.columns.index(coln))
                                ll.add_column(outl_auth_tag)
                                ll.add_column(outl_star_tag)
                            else:
                                ll.add_column(outl_star_tag)
                        else:
                            missing_col.append(loop.columns.index(coln))
                    atm_id=[i for i in range(len(ll.columns)) if "Atom_ID" in ll.columns[i]]
                    if sf.category=="assigned_chemical_shifts":
                        ll.add_column("_Atom_chem_shift.Ambiguity_code")
                    if sf.category=="general_distance_constraints":
                        ll.add_column("_Gen_dist_constraint.Member_logic_code")
                        const_id=1
                    
                    #print atm_id,sf.category
                    if len(missing_col)==0:
                        for dat in loop.data:
                            if len(auth_col)==0:
                                ll.add_data(dat[:])
                            else:
                                if sf.category=="assigned_chemical_shifts":
                                    atm_index=loop.columns.index("atom_name")
                                    res_index=loop.columns.index("residue_type")
                                    atm_list=self.get_atm_list(dat[:][res_index], dat[:][atm_index])
                                    if len(atm_list)==0:
                                        atm_list.append(dat[:][atm_index])
                                    for atm in atm_list:
                                        dat2=dat[:]
                                        for k in auth_col:
                                            if k != atm_index:
                                                dat2.insert(k+auth_col.index(k)+1,dat[:][k])
                                            else:
                                                dat2.insert(k+auth_col.index(k)+1,atm)
                                        if len(atm_list)==1:
                                            dat2.append('1')
                                        else:
                                            dat2.append('2')
                                        ll.add_data(dat2)
                                elif sf.category=="general_distance_constraints":
                                    atm_index_1=loop.columns.index("atom_name_1")
                                    res_index_1=loop.columns.index("residue_type_1")
                                    atm_index_2=loop.columns.index("atom_name_2")
                                    res_index_2=loop.columns.index("residue_type_2")
                                    atm_list_1=self.get_atm_list(dat[:][res_index_1], dat[:][atm_index_1])
                                    atm_list_2=self.get_atm_list(dat[:][res_index_2], dat[:][atm_index_2])
                                    if len(atm_list_1)==0:
                                        atm_list_1.append(dat[:][atm_index_1])
                                    if len(atm_list_2)==0:
                                        atm_list_2.append(dat[:][atm_index_2])
                                    if len(atm_list_1)==1 and len(atm_list_2)==1:
                                        dat2=dat[:]
                                        #print auth_col,ll.columns
                                        for k in auth_col:
                                            dat2.insert(k+auth_col.index(k)+1,dat[:][k])
                                        dat2[0]="%d"%(const_id)
                                        dat2.append('.')
                                        ll.add_data(dat2)
                                        const_id+=1
                                    else:
                                        for atm1 in atm_list_1:
                                            for atm2 in atm_list_2:
                                                dat2=dat[:]
                                                for k in auth_col:
                                                    if k==atm_index_1:
                                                        dat2.insert(k+auth_col.index(k)+1,atm1)
                                                    elif k==atm_index_2:
                                                        dat2.insert(k+auth_col.index(k)+1,atm2)
                                                    else:
                                                        dat2.insert(k+auth_col.index(k)+1,dat[:][k])
                                                dat2.append('OR')
                                                dat2[0]="%d"%(const_id)
                                                ll.add_data(dat2)
                                                const_id+=1
                                        
                                        
                                    
                                    
                                else:
                                    dat2=dat[:]
                                    #print auth_col,ll.columns
                                    for k in auth_col:
                                        dat2.insert(k+auth_col.index(k)+1,dat[:][k])
                                    ll.add_data(dat2[:])
                                    
                    if ll.data: sf.add_loop(ll)                      
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
                refatm=re.findall(r'(\S+)([XY])([%*])$|(\S+)([%*])$|(\S+)([XY]$)',nefAtom)[0]  
                set=[refatm.index(i) for i in refatm if i!=""]
                if set==[0,1,2]:
                    pattern=re.compile(r'%s\S\d+'%(refatm[0]))
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
                    pattern=re.compile(r'%s\S+'%(refatm[5]))
                    alist=[i for i in atms if re.search(pattern, i)]
                    if len(alist)!=2:
                        alist=[]
                    elif refatm[6]=="Y":
                        #alist.reverse()[]
                        alist=alist[-1:]
                    elif refatm[6]=="X":
                        alist=alist[:1]
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
            print "Residue not found",res,nefAtom
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
    nt=NEFtoSTAR('/home/kumaran/git/NEF/data_0_2/CCPN_H1GI_clean.nef')
    #nt.translate()
    nt.convert()
    