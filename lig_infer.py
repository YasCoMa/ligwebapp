
import os
import pandas as pd
import sqlite3
import json
import urllib.request

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from scipy.spatial import distance

import requests
from requests.adapters import HTTPAdapter 
from requests.packages.urllib3.util.retry import Retry
import urllib.parse
from contextlib import closing

from prepare_complete_network  import PrepareNetwork

class LigandInference:
    conn=None
    
    def __init__(self):
        try:
            self.conn = sqlite3.connect("bionetwork.db")
        except Error as e:
            print(e)
    
    def download_pdbs(self, l):
        r=True
        if(not os.path.isdir("structures_pdb")):
            os.system("mkdir structures_pdb")
        
        if(not os.path.isfile("structures_pdb/"+l+".pdb")):
            os.system("wget https://files.rcsb.org/download/"+l+".pdb -O structures_pdb/"+l+".pdb") 
            if(os.path.getsize("structures_pdb/"+l+".pdb")==0):
                r=False
        return r
        
    def _parse_line_pdb(self, line):
        l=line.replace("\n","")
                    
        atom=l[12:16].replace(" ","")
        residue=l[17:20].replace(" ","")
        chain=l[21].replace(" ","")
        pos=l[22:26].replace(" ","")
        x=float(l[30:38].replace(" ",""))
        y=float(l[38:46].replace(" ",""))
        z=float(l[46:54].replace(" ",""))
        coordinates = [x,y,z]
        
        return atom, residue, chain, pos, coordinates 
    
    def get_sequence_from_pdb(self, pdb):
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
     
        resmap={}
        g=open('structures_pdb/'+pdb+'.pdb','r')
        for line in g:
            l=line.replace('\n','')
            if(line.startswith('ATOM')):
                atom, residue, chain, pos, coordinates = self._parse_line_pdb(line)
                if(not chain in resmap.keys()):
                    resmap[chain]=set()
                    
                if(not residue+"-"+pos in resmap):
                    resmap[chain].add(residue+"-"+pos)
        g.close()
            
        result = {}
        for k in resmap.keys():
            seq=""
            for res in resmap[k]:
                if( res.split('-')[0] in d.keys() ):
                    seq+=d[res.split('-')[0]]
            result[k]=seq
        return result
    
    def read_proteins(self, protein_file):
        prots=[]
        f=open(protein_file, "r")
        for line in f:
            l=line.replace("\n","")
            prots.append(l)
        f.close()
        
        return prots
    
    def get_protein_seq(self, p):
        seq=''
        cur = self.conn.cursor()
        cur.execute("select sequence from proteins where id=?;", (p,) )

        rows = cur.fetchall()

        for row in rows:
            seq=row[0]
            
        return seq
        
    def get_family_annotation(self, p):
        pfams=[]
        cur = self.conn.cursor()
        cur.execute("select annotation from map_uniprot_pfam where protein=?;", (p,) )

        rows = cur.fetchall()

        for row in rows:
            pfams.append(row[0])
            
        return pfams
        
    def get_linked_pdbs(self, p):
        pdbs=set()
        cur = self.conn.cursor()
        cur.execute("select annotation from map_uniprot_pdbsum where protein=?;", (p,) )

        rows = cur.fetchall()

        for row in rows:
            pdbs.add(row[0])
            
        return pdbs
    
    def get_contacts_domain(self, info, seqprot, seqpdb):
        matrix = matlist.blosum62
        aln=pairwise2.align.globaldx(seqpdb, seqprot, matrix)
        mx=0
        for i in range(len(info)):
            if( info.iloc[0, 2]!=None and info.iloc[0, 3]!=None ):
                init=int( info.iloc[0, 2] )
                end=int( info.iloc[0, 3] )
                #als1=aln[0].seqA[init:end+1]
                #als2=aln[0].seqB[init:end+1]
                #print(aln[0])
                try:
                    als1=aln[0][0][init:end+1]
                    als2=aln[0][1][init:end+1]
                    
                    num=0
                    den=0
                    for j in range(len(als1)):
                        if(als1[j]==als2[j] and als1[j]!='-' and als2[j]!='-'):
                            num+=1
                            
                        if(als1[j]!='-'):
                            den+=1
                    if(den!=0):
                        ratio=num/den
                        if(ratio>mx):
                            mx=ratio
                except:
                    pass
        return mx
    
    def get_contacts_ligand(self, pdb, cutoff_ligand, info, fchain):
        ligands={}
        receptor={}
        
        garbage=["ALA","LYS","VAL","ARG","GLY","CYS","TYR", "PHE","LEU","ILE","MET", "SER","PRO","THR","HIS","ASN","GLN","GLU","ASP","TRP","NAD","NAI","NDC", "TAP","ADJ","NAJ","NDO","ZID","CAN","NAP","NDP", "CND","NAQ","NHD","DND","NAX","NHO","NAC","NDB","ODP","NAE","NBP", "PAD","NAH","NDA","SND", "FAD","FMN","6FA","FNS","FAA","MGD","FAB","RFL","FAE","FAS","FDA", "FMA", "AMP","GDP","UDP","ADP","GNP","UTP","ANP", "GTP","PSU", "ATP","2GP","CMP","TMP","CDP","TDP","CTP","TTP","GMP","UMP","SO4", "UNK", "UNX", "BMA","FUC","MAN","POP","BOG","GAL","MES","PYR","C8E","GLC","MPD", "SPM","CIT","GOL","MYR","TRS","CRY","HED","NAG","XYS","DTT","LDA","NGA", "EPE","LI1","PEG","F6P","MAL","PG4", "HEM","CD","NA", "MG", "ZN", "CL", "MN", "FE", "HOH", "PO4", "CA", "EDO"]
        
        init=0
        end=0
        if( info.iloc[0, 2]!=None and info.iloc[0, 3]!=None ):
            init=int( info.iloc[0, 2] )
            end=int( info.iloc[0, 3] )
            
        f=pdb+".pdb"
        g=open("structures_pdb/"+f,"r")
        for line in g:
            if(line.startswith("ATOM")):
                atom, residue, chain, pos, coordinates = self._parse_line_pdb(line)
                if(atom=='CA' and chain==fchain and int(pos)>=init and int(pos)<=end): # filter positions corresponding to the domain region
                    receptor[chain]=coordinates
                    
            if(line.startswith("HETATM")):
                atom, residue, chain, pos, coordinates = self._parse_line_pdb(line)
                if(not residue in garbage):
                    if (not residue in ligands.keys()):
                        ligands[residue]={}
                        
                    ligands[residue][chain]=coordinates
                
        g.close()  
        
        result={}
        for l in ligands.keys():
            result[l]={}
            for c in ligands[l]:
                passed=0
                if(c in receptor.keys()):
                    for coordp in receptor[c]:
                        summ=0
                        for coordl in ligands[l][c]:
                            try:
                                dist = distance.euclidean(coordp, coordl)
                                if(dist>cutoff_ligand):
                                    summ+=1
                            except:
                                pass
                        if(summ>0):
                            passed+=1
                    result[l][c]=[passed, summ]
                
        return result
    
    def get_info_ligand(self, lig):
        smi=''
        chemblid=''
        cur = self.conn.cursor()
        cur.execute("select smiles from ligands_info where ligand=?;", (lig,) )
        rows = cur.fetchall()
        if(len(rows)==0):
            link="https://data.rcsb.org/rest/v1/core/chemcomp/"+lig
            resp= urllib.request.urlopen(link)
            data=resp.read()
            encoding = resp.info().get_content_charset('utf-8')
            dat=json.loads(data.decode(encoding))
            if( 'pdbx_chem_comp_descriptor' in dat.keys()):
                if( len(dat['pdbx_chem_comp_descriptor'])>0 ):
                    for d in dat['pdbx_chem_comp_descriptor']:
                        if(d['type']=="SMILES"):
                            smi=d['descriptor']
            
            if( 'rcsb_chem_comp_related' in dat.keys() ):
                if( len(dat['rcsb_chem_comp_related'])>0 ):
                    for d in dat['rcsb_chem_comp_related']:
                        if(d['resource_name']=="ChEMBL"):
                            chemblid=d['resource_accession_code']
            
            #if(smi!='' and chemblid!=''):            
            sql = ''' insert into ligands_info values (?,?,?) '''
            cur = self.conn.cursor()
            cur.execute(sql, (lig, smi, chemblid) )
            self.conn.commit()
        else:
            for row in rows:
                smi=row[0]
                            
        return smi
    
    def get_zinc(self, smiles, tanimoto_zinc=60):
        table=[]
        try:
            smiles_2=urllib.parse.quote(smiles.encode("utf8"))
            link="https://zinc20.docking.org/substances/subsets/world.txt?ecfp4_fp-tanimoto-"+str(tanimoto_zinc)+"="+smiles_2+"&purchasability=for-sale&count=all"
            resp= urllib.request.urlopen(link)
            data=resp.read()
            encoding = resp.info().get_content_charset('utf-8')
            dat=data.decode(encoding)
            for line in dat.split("\n"):
                if(line!=""):
                    l=line.split("\t")
                    if(l[0]!=""):
                        table.append( (l[0], l[1]) )
        except:
            pass
        return table
        
    def get_zinc_old(self, smiles, tanimoto_zinc=60):
        '''Busca en zinc ligandos comerciables similares al smiles original con un minimo
        de 60% de similaridad por tanimoto'''
        retry_strategy = Retry(total=3, status_forcelist=[429, 500, 502, 503, 504], 
                            method_whitelist=["HEAD", "GET", "OPTIONS", "POST"])
        adapter = HTTPAdapter(max_retries=retry_strategy)
        http = requests.Session() 
        http.mount("https://", adapter) 
        http.mount("http://", adapter)
        smiles_2=urllib.parse.quote(smiles.encode("utf8"))
        #request_txt=f"https://zinc20.docking.org/substances.txt?ecfp4_fp-tanimoto-{tanimoto_zinc}={smiles_2}&purchasability=for-sale"
        request_txt=f"https://zinc20.docking.org/substances/subsets/world.txt?ecfp4_fp-tanimoto-{tanimoto_zinc}={smiles_2}&purchasability=for-sale&count=all"
        #request_txt="https://zinc20.docking.org/substances/subsets/world.txt?ecfp4_fp-tanimoto-60=c1ccccc1O&purchasability=for-sale&count=all"
        out=requests.get(request_txt, stream=True, allow_redirects=True)
        if out.status_code==200:
            with closing(out) as r, open('/tmp/zinc', 'w') as h:
                for x in r.iter_content(chunk_size=512*1024):
                    h.write(x.decode())
            with open('/tmp/zinc', 'r') as h:
                info=h.readlines()
            data=[]
            if info:
                if len(info)>1 and info!=None:
                    for zinc_smiles in info:
                        zinc_smiles=zinc_smiles.replace('\n','')
                        zinc_id=zinc_smiles.split('\t')[0]
                        smiles_id=zinc_smiles.split('\t')[1]
                        if zinc_id:
                            data.append((zinc_id, smiles_id))
            return data
        else:
            sys.stderr.write(f'Fail request to Zinc using SMILES: {request_txt} \n')
                                
    def tanimoto_calc(self, smi1, smi2):
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
        s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        return s
        
    def enrich_ligands_protein(self, pr, pfam_map, cutoff_ligand):
        seq=self.get_protein_seq(pr)
        
        pfams=self.get_family_annotation(pr)
        ligs=[]
        for pf in pfams:
            print('\tst1', pf)
            
            mapp=pfam_map[ pfam_map['PFAM_ACCESSION']==pf ][ ['PDB', 'CHAIN', 'PDB_START', 'PDB_END'] ]
            if(len(mapp)>0):
                linked_pdbs=self.get_linked_pdbs(pr)
                all_pdbs_from_family=set(mapp['PDB'].unique())
                filtered=all_pdbs_from_family-linked_pdbs
                print('number of proteins: ', len(filtered) )
                for p in filtered:
                    print('\tst2', p)
                    cur = self.conn.cursor()
                    cur.execute("select count(*) from inferred_pdb_pfam where protein=? and pdb=? and pfam=?;", (pr, p, pf,) )
                    count = cur.fetchone()[0]
                    if(count==0):
                        p=p.upper()
                        flag=self.download_pdbs(p)
                        if(flag):
                            seqpdb_chains=self.get_sequence_from_pdb(p)
                            chains=set(seqpdb_chains.keys())
                            #print(chains)
                            if(len(chains)>0):
                                info=mapp[mapp['PDB'].str.upper()==p]
                                chainsd=set(info['CHAIN'].unique())
                                if(len(chainsd)>0):
                                    chain=list(chains.intersection(chainsd))
                                    seqpdb=seqpdb_chains[chain[0]]
                                    info=info[ info['CHAIN']==chain[0] ]
                                    
                                    # Get conserved aas in the mapped domain region aligning structures sequence with seed protein sequence
                                    ratio_domain = self.get_contacts_domain(info, seq, seqpdb)
                                    #print(ratio_domain)
                                    sql = ''' insert into inferred_pdb_pfam values (?,?,?,?) '''
                                    cur = self.conn.cursor()
                                    cur.execute(sql, (pr, pf, p, ratio_domain,) )
                                    self.conn.commit()
                                    
                                    # Get ligands 
                                    # pr, pdb, ligand, chain, receptor_contacts, all_contacts
                                    result=self.get_contacts_ligand(p, cutoff_ligand, info, chain[0])
                                    #print(result)
                                    for lg in result.keys():
                                        if(not lg in ligs):
                                            ligs.append(lg)
                                            
                                        for c in result[lg]:
                                            if(c==chain[0]):
                                                sql = ''' insert into contact_ligands_inferred values (?,?,?,?,?,?) '''
                                                cur = self.conn.cursor()
                                                cur.execute(sql, (pr, p, lg, c, result[lg][c][0], result[lg][c][1],) )
                                                self.conn.commit()
                        
        
        smis={}
        for l in ligs:
            smis[l]=self.get_info_ligand(l)
            print('\tlig', l)
        
        """
        i=0
        for l1 in ligs:
            j=0
            for l2 in ligs:
                if(i<j):
                    cur.execute("select distinct ligand1, ligand2, tanimoto from ligands_alignment where ligand1=? and ligand2=?;", (l1, l2,) )
                    rows = cur.fetchall()
                    if( len(rows) == 0 ):
                        s1=smis[l1]
                        s2=smis[l2]
                        if(s1!='' and s2!=''):
                            try:
                                score=self.tanimoto_calc(s1, s2)
                                sql = ''' insert into ligands_alignment values (?,?,?) '''
                                cur = self.conn.cursor()
                                cur.execute(sql, (l1, l2, score) )
                                self.conn.commit()
                            except:
                                pass
                j+=1
            i+=1
        """
            
        for l in ligs:
            smi=smis[l]
            
            if(smi!=''):
                print('lig2', l)
                cur = self.conn.cursor()
                cur.execute("select count(*) from inferred_ligands_zinc where ligand_query=? ;", (l,) )
                count = cur.fetchone()[0]
                if(count==0):
                    dat=self.get_zinc(smi, 70)
                    print('\t\tzinc', smi)
                    if(dat!=None):
                        if(len(dat)>0):
                            for d in dat:
                                sql = ''' insert into inferred_ligands_zinc values (?,?,?) '''
                                cur = self.conn.cursor()
                                cur.execute(sql, (l, smi, d[1]) )
                                self.conn.commit()
                        
        
    def migrate_ligands_data(self):
        j=1
        df=pd.read_csv('aux_data/merged_ligand_descriptors.tsv', sep='\t')
        total=len(df['ligand'])
        for i in df.index:
            lig=df.loc[i, 'ligand']
            smi=df.loc[i, 'smiles']
            chemblid=df.loc[i, 'chembl']
            print(j, '/', total)
            
            cur = self.conn.cursor()
            cur.execute("select smiles from ligands_info where ligand=?;", (lig,) )
            rows = cur.fetchall()
            if(len(rows)==0):
                sql = ''' insert into ligands_info values (?,?,?) '''
                cur = self.conn.cursor()
                cur.execute(sql, (lig, smi, chemblid) )
                self.conn.commit()
            j+=1
    
    def run(self, proteins, cutoff_ligand, id_job, cutoff_num_contact, cutoff_perc_contact, email): # file with proteins uniprot ids
        pfam_map=pd.read_csv("pdb_pfam_mapping.txt", sep="\t")
        
        prots=self.read_proteins(proteins)
        
        for p in prots:
            print(p)
            self.enrich_ligands_protein(p, pfam_map, cutoff_ligand)
        
        # Draw graph
        a=PrepareNetwork()
        a.mount_all_graph(id_job, proteins, cutoff_num_contact, cutoff_perc_contact, email)
        
import sys   

a=LigandInference()
protein_file=sys.argv[1]
cutoff_num_contact=float(sys.argv[2])
cutoff_perc_contact=float(sys.argv[3])
id_job=sys.argv[4]
email=sys.argv[5]
a.run(protein_file, 10, id_job, cutoff_num_contact, cutoff_perc_contact, email)

#a.migrate_ligands_data()

#dat=a.get_zinc('COc1ccc2ccc(cc2c1OCC(N)=O)C(N)=N', 60)
#print(dat)

