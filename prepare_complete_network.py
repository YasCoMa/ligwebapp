import os
import sqlite3
import shutil

import networkx as nx
import pandas as pd

import smtplib
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email import encoders

class PrepareNetwork:
    conn=None
    
    def __init__(self):
        try:
            self.conn = sqlite3.connect("bionetwork.db")
        except Error as e:
            print(e)
            
    def load_annotations(self, properties, ids, g, proteins, id_job, eprops):
        out=open(id_job+"/protein_annotations.tsv","w")
        out.write("protein\tannotation_type\tvalues\n")
        cur = self.conn.cursor()
        for p in proteins:
            if(not p in ids):
                g.add_node(p)
                ids.add(p)
                properties['type'][p]='Protein'
                properties['title'][p]=p
                properties['label'][p]=p
            
            types=["Family", "Structure", "Drug"]
            info=['pfam', 'pdbsum', 'drugbank']
            j=0
            for i in info:
                cur.execute("select distinct annotation from map_uniprot_"+i+" where protein=?;", (p,) )
                rows = cur.fetchall()
                
                ans=set()
                for row in rows:
                    ans.add(row[0])
                    
                    if(not row[0] in ids):
                        g.add_node(row[0])
                        ids.add(row[0])
                        properties['type'][row[0]]=types[j]
                        properties['title'][row[0]]=row[0]
                        properties['label'][row[0]]=row[0]
                        
                    g.add_edge(p, row[0])
                    eprops[(p, row[0])]={ "identity": '-', "tanimoto": "-", "number_contacts": "-", "percentage_contacts": "-"  }
                j+=1
                
                if(len(rows)>0):
                    out.write("%s\t%s\t%s\n" %(p, i, (','.join(list(ans))) ) )
        
        out.close()
        
        return g, ids, properties, eprops
            
    def load_protein_alignment(self, properties, ids, g, proteins, eprops, id_job):
        out=open(id_job+"/proteins_alignment.tsv","w")
        out.write("protein_1\tprotein_2\tidentity\n")
        cur = self.conn.cursor()
        i=0
        for p1 in proteins:
            j=0
            for p2 in proteins:
                if(i<j):
                    cur.execute("select identity from proteins_alignment where idp1=? and idp2=?;", (p1, p2,) )
                    rows = cur.fetchall()

                    for row in rows:
                        out.write("%s\t%s\t%.6f\n" %(p1, p2, row[0]/100 ) )
                        
                        g.add_edge(p1, p2)
                        eprops[(p1, p2)]={ "identity": row[0]/100, "tanimoto": "-", "number_contacts": "-", "percentage_contacts": "-"  }
                j+=1
            i+=1
        out.close()
        return g, ids, properties, eprops
            
    def load_ligand_information_alignment(self, properties, ids, g, proteins, eprops, cutoff_num_contact, cutoff_perc_contact, id_job):
        reports={'direct_relation_protein_pdb_ligand.tsv': ['protein','ligand','pdb','chembl','smiles'], 'inferred_contacts_pdb_ligand.tsv': ['protein','ligand','pdb','count_protein_contacts_with_ligands','smiles','chembl'], 'inferred_conserved_contacts_pdb_pfam.tsv': ['protein','pdb','family', 'percentage_conserved_region'], 'inferred_ligands_zinc.tsv': ['ligand_id','smiles','zinc_id'], 'ligands_alignment.tsv': ['ligand_1','ligand_2','tanimoto']}
        for r in reports.keys():
            f=open(id_job+'/'+r, 'w')
            f.write('\t'.join(reports[r])+'\n')
            f.close()
        
        cur = self.conn.cursor()
        ligs=set()
        for p1 in proteins:
            cur.execute("select distinct p.ligand, p.pdb, chembl, smiles from map_uniprot_pdbsum as m, pdb_ligand as p, ligands_info as l where m.annotation=p.pdb and l.ligand=p.ligand and m.protein=? ;", (p1,) )
            rows = cur.fetchall()
            for row in rows:
                if(not row[0] in ids): # ligand
                    p=row[0]
                    ids.add(p)
                    if(not p in ligs):
                        ligs.add(p)
                    g.add_node(p)
                    properties['type'][p]='Ligand'
                    properties['title'][p]=row[3]
                    properties['label'][p]=row[2]
                    
                if(not row[1] in ids): # pdb
                    p=row[1]
                    ids.add(p)
                    g.add_node(p)
                    properties['type'][p]='Structure'
                    properties['title'][p]=p
                    properties['label'][p]=p
                    
                g.add_edge(row[0], row[1])
                eprops[(row[0], row[1])]={ "identity": '-', "tanimoto": "-", "number_contacts": "-", "percentage_contacts": "-"  }
                
                with open(id_job+'/direct_relation_protein_pdb_ligand.tsv','a') as out:
                    out.write( "%s\t%s\t%s\t%s\t%s\n" %(p1, row[0], row[1], row[2], row[3]) )
             
            cur.execute("select p.ligand, p.pdb, p.receptor_contacts, l.smiles, l.chembl from contact_ligands_inferred as p, ligands_info as l where p.protein=? and p.ligand=l.ligand;", (p1,) )
            rows = cur.fetchall()
            for row in rows:   
                if(not row[0] in ids): # ligand
                    p=row[0]
                    ids.add(p)
                    if(not p in ligs):
                        ligs.add(p)
                    g.add_node(p)
                    properties['type'][p]='Ligand'
                    properties['title'][p]=row[4]
                    properties['label'][p]=row[3]
            
                if(row[2]>=cutoff_num_contact):
                    g.add_edge(p1, row[0])
                    eprops[(p1, row[0])]= { "identity": "-", "tanimoto": "-", "number_contacts": row[2], "percentage_contacts": "-"  }
                
                with open(id_job+'/inferred_contacts_pdb_ligand.tsv','a') as out:
                    out.write( "%s\t%s\t%s\t%s\t%s\t%s\n" %(p1, row[0], row[1], row[2], row[3], row[4]) )
             
            cur.execute("select p.pdb, p.ratio, p.pfam from inferred_pdb_pfam as p where p.protein=?;", (p1,) )
            rows = cur.fetchall()
            gone=set()
            for row in rows:   
                if(not row[0] in ids): # pdb
                    p=row[0]
                    ids.add(p)
                    g.add_node(p)
                    properties['type'][p]='Structure'
                    properties['title'][p]=p
                    properties['label'][p]=p
                
                if(not row[1] in gone):
                    gone.add(row[0])
                    if( row[1] >= (cutoff_perc_contact/100) ):
                        g.add_edge(p1, row[0])
                        eprops[(p1, row[0])]= { "identity": "-", "tanimoto": "-", "number_contacts": "-", "percentage_contacts": row[1]  }
                
                with open(id_job+'/inferred_conserved_contacts_pdb_pfam.tsv','a') as out:
                    out.write( "%s\t%s\t%s\t%.6f\n" %(p1, row[0], row[2], row[1]*100 ) )
            
        i=0
        for l1 in ligs:
            cur.execute("select distinct ligand_query, smiles, zincid from inferred_ligands_zinc where ligand_query=?;", (l1,) )
            rows = cur.fetchall()
            for row in rows:
                if(not row[2] in ids): # pdb
                    p=row[0]
                    ids.add(row[2])
                    g.add_node(row[2])
                    properties['type'][row[2]]='ZINC Compound'
                    properties['title'][row[2]]=row[2]
                    properties['label'][row[2]]=row[2]
                g.add_edge(l1, row[2])
                eprops[(l1, row[2])]= { "identity": "-", "tanimoto": 0.7, "number_contacts": "-", "percentage_contacts": "-"  }
                
                with open(id_job+'/inferred_ligands_zinc.tsv','a') as out:
                    out.write( "%s\t%s\t%s\n" %(l1, row[0], row[1] ) )
                
            """
            j=0
            for l2 in ligs:
                cur.execute("select distinct ligand1, ligand2, tanimoto from ligands_alignment where ligand1=? and ligand2=?;", (l1, l2,) )
                rows = cur.fetchall()

                for row in rows:
                    g.add_edge(l1, l2)
                    eprops[(l1, l2)]= { "identity": "-", "tanimoto": row[2], "number_contacts": "-", "percentage_contacts": "-"  }
                
                    with open(id_job+'/ligands_alignment.tsv','a') as out:
                        out.write( "%s\t%s\t%.6f\n" %(l1, l2, row[2] ) )
                j+=1
                
            i+=1
            """
        
        return g, ids, properties, eprops
    
    def mount_all_graph(self, id_job, protein_file, cutoff_num_contact, cutoff_perc_contact, email):
        proteins=set()
        f=open(protein_file, "r")
        for line in f:
            p=line.replace('\n','')
            proteins.add(p)
        f.close()    
            
        id_job="data_users/"+id_job
        if(not os.path.isdir(id_job)):
            os.system("mkdir "+id_job)
        
        ids=set()
        properties={'type': {}, 'title': {}, 'label': {}}
        eprops={}
        g = nx.Graph()
        
        g, ids, properties, eprops = self.load_annotations(properties, ids, g, proteins, id_job, eprops)
        g, ids, properties, eprops = self.load_protein_alignment(properties, ids, g, proteins, eprops, id_job)
        g, ids, properties, eprops = self.load_ligand_information_alignment(properties, ids, g, proteins, eprops, cutoff_num_contact, cutoff_perc_contact, id_job)
        
        for t in properties.keys():
            nx.set_node_attributes(g, properties[t], name=t)
        nx.set_edge_attributes(g, eprops)
            
        df=pd.DataFrame()
        df["ids"]=list(properties["type"].keys())
        for p in properties.keys():
            df[p]=properties[p].values()
        df.to_csv(id_job+"/table_nodes.tsv")
        
        df=pd.DataFrame()
        keys=list(eprops.keys())
        c={"from": [], "to": []}
        if(len(keys)>0):
            for p in eprops[keys[0]].keys():
                c[p]=[]
        for k in eprops:
            c["from"].append(k[0])
            c["to"].append(k[1])
            if(len(keys)>0):
                for p in eprops[k].keys():
                    c[p]=eprops[k][p]
                
        for k in c.keys():
            df[k]=c[k]
        df.to_csv(id_job+"/table_edges.tsv")
        
        try:
            df = pd.DataFrame(dict(
                DEGREE_CENTRALITY      = nx.degree_centrality(g),
                EIGENVECTOR            = nx.eigenvector_centrality(g),
                #KATZ                   = nx.katz_centrality_numpy(g),
                CLOSENESS_CENTRALITY   = nx.closeness_centrality(g),
                BETWEENNESS_CENTRALITY = nx.betweenness_centrality(g),
                CLUSTCOEF              = nx.clustering(g),
            ))
            df.index += 1
            df.to_csv(id_job+'/report_nodes_main_centrality_metrics.tsv', sep="\t")
        except: 
            pass
            
        nx.write_graphml(g, id_job+"/complete_graph.graphml")
        
        shutil.make_archive(id_job, 'zip', id_job)
        self._send_success_email(email, id_job)
        
    def _send_success_email(self, dest, file) :
        name=file.split("/")[-1]
        
        fromaddr = 'ligqapp@outlook.fr'
        frompass = 'ligquba2022'
        
        toaddr = dest
        msg = MIMEMultipart()
        msg['From'] = fromaddr
        msg['To'] = toaddr
        
        msg['Subject'] = "Your Job in Ligand Inference App finished"
        body = "<p>Hello! <br /> <p>The complete generated network files related to your job are attached in this e-mail </p> <br /> Thanks for using our app."
        msg.attach(MIMEText(body, 'html'))
        
        part = MIMEBase("application", "octet-stream")
        part.set_payload(open(file + ".zip", "rb").read())
        encoders.encode_base64(part)
        part.add_header("Content-Disposition", "attachment; filename=\"%s.zip\"" % (name))
        msg.attach(part)
        
        #server = smtplib.SMTP('smtp.gmail.com', 587)
        server = smtplib.SMTP("smtp-mail.outlook.com", 587)
        server.ehlo()
        server.starttls()
        server.login(fromaddr, frompass)
        
        text = msg.as_string()
        server.sendmail(fromaddr, toaddr, text)
        server.quit()

        return None
        
