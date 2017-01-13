#!/usr/bin/env python

'''
Created on Sep 15, 2016

This script translates NEV v1.0 files to NMR-STAR v 3.2.0.1

This script depends on 

    * NMR-STAR parser provided along with this script 
    * NEF_NMRSTAR_equivalence.csv file

@author: kumaran
'''
import os
import sys,csv,ntpath,re,time,datetime,string
from Bio.MarkovModel import save
(script_path,script_name)=ntpath.split(os.path.realpath(__file__))
sys.path.append(script_path+'/../PyNMRSTAR') #NMR-STAR and NEF-Parser added as a submodule and imported into this project. This is a separate git repository
import bmrb












class NEFtoSTAR(object):
    '''
    Translates NEF file to NMR-STAR file
    '''
    #mapping file btween NEF and NMR-STAR
    
    mapFile=script_path+'/NEF_NMRSTAR_equivalence.csv'
    __version__=1.0
    def st(self,ts):
        return datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    
    #Standard atom names for 20 standard amino acids and 8 Nucleic acid 
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
    
    def __init__(self, inFile):
        self.inFile=inFile
        (self.inFilePath,self.inFileName)=ntpath.split(self.inFile)
        self.read_map_file()
        self.logFile=self.inFile.split("."+self.inFileName.split(".")[-1])[0]+"_.log"
        self.outFile=self.inFile.split("."+self.inFileName.split(".")[-1])[0]+"_.str"
        #print self.logFile
        self.convert()
        
    def read_map_file(self):
        '''Reads the NEF_NMRSTAR_equivalence.csv file and create a mapping as a list'''
        with open(self.mapFile,'rb') as csvfile:
            spamreader = csv.reader(csvfile,delimiter=',')
            map_dat=[]
            for r in spamreader:
                #print r
                if r[0][0]!='#':
                    map_dat.append(r)
        self.map=map(list,zip(*map_dat))
    
  
    
    def convert(self):
        self.logfile=open(self.logFile,'w')
        self.logfile.write("\n%s:%s\n"%(string.ljust("Input",25),self.inFile))
        self.logfile.write("%s:%s\n"%(string.ljust("Output",25),self.outFile))
        self.logfile.write("%s:%s\n"%(string.ljust("Log ",25),self.logFile))
        self.logfile.write("%s:%s\n"%(string.ljust("Translator Version",25),self.__version__))
        self.logfile.write("%s:%s\n"%(string.ljust("STAR parser Version",25),bmrb._VERSION))
        self.logfile.write("%s:%s\n\n"%(string.ljust("Date",25),self.st(time.time())))
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
                self.details={}
                if saveframe.get_tag("sf_category")[0] in self.map[0]:
                    sf=bmrb.Saveframe.from_scratch(saveframe.name)
                    for tag in saveframe.tags:
                        in_tag="%s.%s"%(saveframe.tag_prefix,tag[0])
                        try:
                            out_auth_tag=self.map[1][self.map[0].index(in_tag)]
                            out_star_tag=self.map[2][self.map[0].index(in_tag)]
                        except ValueError:
                            self.logfile.write("%s\tSoftware specific tag %s ignored\n"%(self.st(time.time()),in_tag))
                            self.details[tag[0]]=tag[1]
                            out_auth_tag=""
                            out_star_tag=""
                        if out_auth_tag!="" and out_star_tag!="":
                            if out_auth_tag==out_star_tag:
                                if out_star_tag.split(".")[1]=="Sf_category" or out_star_tag.split(".")[1]=="sf_category":
                                    out_category=self.map[1][self.map[0].index(tag[1])]
                                    sf.add_tag(out_star_tag,out_category)
                                else:
                                    if out_star_tag=="_Entry.NMR_STAR_version":
                                        sf.add_tag(out_star_tag,'3.2.0.1')
                                    else:
                                        sf.add_tag(out_star_tag,tag[1])
                            else:
                                sf.add_tag(out_auth_tag,tag[1])
                                sf.add_tag(out_star_tag,tag[1])
                                #print out_auth_tag,out_star_tag
                    for loop in saveframe:
                        ll=bmrb.Loop.from_scratch()
                        missing_col=[]
                        auth_col=[]
                        for coln in loop.columns:
                            inl_tag="%s.%s"%(loop.category,coln)
                            try:
                                outl_auth_tag=self.map[1][self.map[0].index(inl_tag)]
                                outl_star_tag=self.map[2][self.map[0].index(inl_tag)]
                            except ValueError:
                                self.logfile.write("%s\tSoftware specific tag %s ignored\n"%(self.st(time.time()),inl_tag))
                                self.details[inl_tag]=str(loop.get_tag(inl_tag))
                                outl_auth_tag=""
                                outl_star_tag=""
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
                        
                        for dat in loop.data:
                            if len(auth_col)==0:
                                if len(missing_col)!=0:
                                    dat3=dat[:]
                                    for m in sorted(missing_col,reverse=True):
                                        del dat3[m]
                                else:
                                    dat3=dat[:]
                                ll.add_data(dat3[:])
                            else:
                                if sf.category=="assigned_chemical_shifts":
                                    atm_index=loop.columns.index("atom_name")
                                    res_index=loop.columns.index("residue_name")
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
                                        if len(missing_col)!=0:
                                            dat3=dat2[:]
                                            for m in sorted(missing_col,reverse=True):
                                                del dat3[m]
                                        else:
                                            dat3=dat2[:]
                                        ll.add_data(dat3)
                                elif sf.category=="general_distance_constraints":
                                    atm_index_1=loop.columns.index("atom_name_1")
                                    res_index_1=loop.columns.index("residue_name_1")
                                    atm_index_2=loop.columns.index("atom_name_2")
                                    res_index_2=loop.columns.index("residue_name_2")
                                    atm_list_1=self.get_atm_list(dat[:][res_index_1], dat[:][atm_index_1])
                                    atm_list_2=self.get_atm_list(dat[:][res_index_2], dat[:][atm_index_2])
                                    if len(atm_list_1)==0:
                                        atm_list_1.append(dat[:][atm_index_1])
                                    if len(atm_list_2)==0:
                                        atm_list_2.append(dat[:][atm_index_2])
#                                     if len(atm_list_1)==1 and len(atm_list_2)==1:
#                                         dat2=dat[:]
#                                         #print auth_col,ll.columns
#                                         for k in auth_col:
#                                             dat2.insert(k+auth_col.index(k)+1,dat[:][k])
#                                         dat2[0]="%d"%(const_id)
#                                         dat2.append('.')
#                                         if len(missing_col)!=0:
#                                             dat3=dat2[:]
#                                             for m in sorted(missing_col,reverse=True):
#                                                 del dat3[m]
#                                         else:
#                                             dat3=dat2[:]
#                                         ll.add_data(dat3)
#                                         const_id+=1
#                                     else:
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
                                            if len(missing_col)!=0:
                                                dat3=dat2[:]
                                                for m in sorted(missing_col,reverse=True):
                                                    del dat3[m]
                                            else:
                                                dat3=dat2[:]
                                            ll.add_data(dat3)
                                            const_id+=1
                                        
                                        
                                    
                                    
                                else:
                                    dat2=dat[:]
                                    #print auth_col,ll.columns
                                    for k in auth_col:
                                        dat2.insert(k+auth_col.index(k)+1,dat[:][k])
                                    if len(missing_col)!=0:
                                        dat3=dat2[:]
                                        for m in sorted(missing_col,reverse=True):
                                            del dat3[m]
                                    else:
                                        dat3=dat2[:]
                                    ll.add_data(dat3[:])
                                       
                        if not self.is_empty(ll.data): sf.add_loop(ll)
                        #print sf
                        #print "%s.Details"%(sf.tag_prefix)
                        
                    if self.details: 
                        sf.add_tag("Details",str(self.details))
                        self.logfile.write("%s\tSoftware specific tags fond in Saveframe %s;Copied to %s.Details tag\n"%(self.st(time.time()),saveframe.name,sf.tag_prefix))
                                              
                    self.outData.add_saveframe(sf)
                else:
                    self.logfile.write("%s\tSoftware specific saveframe found; %s\n"%(self.st(time.time()),saveframe.name))
                    sf=bmrb.Saveframe.from_scratch(saveframe.name,tag_prefix="_Software_specific_NEF_info_list")
                    sf.add_tag("Sf_category","NEF_specific_software_saveframes")
                    sf.add_tag("Sf_framecode",saveframe.name)
                    sf.add_tag("nef_sf_category",saveframe.get_tag('sf_category')[0])
                    sf.add_tag("nef_tag_prefix",saveframe.tag_prefix)
                    ll=bmrb.Loop.from_scratch(category="_Software_specific_NEF_info")
                    ll.add_column("nef_tag")
                    ll.add_column("nef_data")
                    for d in saveframe.tags:
                        ll.add_data(d[:-1])
                    for loop in saveframe:
                        ll.add_data([loop.category,str(loop)])
                    if not self.is_empty(ll.data): sf.add_loop(ll)
                    self.outData.add_saveframe(sf)    
                    #print "Ignoring saveframe",s                                 # add the saveframe to data structure
#         (file_path,file_name)=ntpath.split(self.inFile)                        # write output file
#         if file_path=="": file_path="."
#         if ".nef" in file_name:
#             out_file_name="%s_.str"%(file_name.split(".nef")[0])
#             out_file_name2="%s_.nef"%(file_name.split(".nef")[0])
#         else:
#             out_file_name="%s_.nef"%(file_name.split(".str")[0])
#             out_file_name2="%s_.str"%(file_name.split(".str")[0])
#         outfile=file_path+"/"+out_file_name
#         outfile2=file_path+"/"+out_file_name2
        #print outfile2,outfile
        self.logfile.write("%s\tWritting output.........\n"%(self.st(time.time())))
        with open(self.outFile,'w') as strfile:
            strfile.write(str(self.outData))
        self.logfile.write("%s\tOutput written\n"%(self.st(time.time())))
        self.logfile.write("%s\tFinished sucessfully\n"%(self.st(time.time())))
        self.logfile.close()
#         with open(outfile2,'w') as strfile2:
#             strfile2.write(str(self.inData))
    
    def is_empty(self,input_list):
        """Recursively iterate through values in nested lists."""
        for item in input_list:
            if not isinstance(item, list) or not self.is_empty(item):
                return False
        return True
        
                                
                        
    
        
                    
                    
                    
                            
                        
                
           
              
    
                        
    def get_atm_list(self,res,nefAtom):
        try:
            atms=self.NMR_STAR_atom_names[res]
            alist=[]
            try:
                refatm=re.findall(r'(\S+)([xy])([%*])$|(\S+)([%*])$|(\S+)([xy]$)',nefAtom)[0]  
                set=[refatm.index(i) for i in refatm if i!=""]
                if set==[0,1,2]:
                    pattern=re.compile(r'%s\S\d+'%(refatm[0]))
                    alist=[i for i in atms if re.search(pattern, i)]
                    if refatm[1]=="y":
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
                    elif refatm[6]=="y":
                        #alist.reverse()[]
                        alist=alist[-1:]
                    elif refatm[6]=="x":
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
            self.logfile.write("%s\tResidue not found,%s,%s\n"%(self.st(time.time()),res,nefAtom))
            #print "Residue not found",res,nefAtom
            alist=[]
        return alist                    
                    
                
    
    def test(self):
        for i in range(len(self.map[0])):
            if self.map[1][i]!=self.map[2][i]:
                print self.map[0][i],self.map[1][i],self.map[2][i]

    
   
        
if __name__=="__main__":
    fname=sys.argv[1]
    #nt=NEFtoSTAR("/home/kumaran/git/NEF/specification/Commented_Example.nef")
    nt=NEFtoSTAR(fname)
    #neflist=["CCPN_1nk2_docr.nef", "CCPN_2mqq_docr.nef","CCPN_Commented_Example.nef","CCPN_Sec5Part3.nef","CCPN_2kko_docr.nef","CCPN_2mtv_docr.nef","CCPN_H1GI_clean.nef","CCPN_XplorNIH-simple.nef"]
    
#     for fname in neflist:
#         fname2="/home/kumaran/git/NEF/data_0_2/%s"%(fname)
#         print fname2
#         nt=NEFtoSTAR(fname2)
#     #nt.translate()
#         nt.convert()
#     