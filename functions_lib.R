#library(redland)
library(readr)
library(jsonlite)
library(dplyr)
library(RSQLite)
library(httr)
library(Biostrings)
library(ChemmineR)
#library(ChemmineOB)

#source("compare_contacts.R")

# initialize tables
init_db<-function(conn){
    temp <- dbGetQuery(conn, "select name from sqlite_schema where type ='table' and name not like 'sqlite_%';")
    info <- c('pfam', 'pdbsum', 'drugbank')
    for (i in info){
        if(! paste0("map_uniprot_", i) %in% temp$name ){
            dbExecute(conn, paste0("create table ", "map_uniprot_", i, " (id int auto_increment primary key, protein text, annotation text);") )
        }
    }
    
    tables=c("proteins_alignment", "ligands_alignment", "pdb_ligand", "proteins", "ligands_info", "inferred_pdb_pfam", "inferred_ligands_zinc", "contact_ligands_inferred", "ligands_broken")
    commands=c("create table proteins_alignment (idp1 text, idp2 text, identity decimal(10,5) );", 
               "create table ligands_alignment (ligand1 text, ligand2 text, tanimoto decimal(10,5) );",
               "create table pdb_ligand (pdb text, ligand text);",
               "create table proteins (id text, sequence text);",
               "create table ligands_info (ligand text, smiles text, chembl text);",
               "create table inferred_pdb_pfam (protein text, pfam text, pdb text, ratio decimal(10,5) );",
               "create table inferred_ligands_zinc (ligand_query text, smiles text, zincid text );",
               "create table contact_ligands_inferred (protein text, pdb text, ligand text, chain text, receptor_contacts decimal(10,5), all_contacts decimal(10,5));",
               "create table ligands_broken (ligand text);"
    )
    i=1
    for (t in tables){
        if(! t %in% temp$name){
            dbExecute(conn, commands[[i]])    
        }
        i=i+1
    }
}

# make sequence alignment
process_protein_sequence_aln <- function(conn, p1, p2 ){
    aas1=NULL
    aas2=NULL
    temp <- dbGetQuery(conn, paste0("select sequence from proteins where id='", p1,"';"))
    if(! is.na(temp$sequence[1])){
        aas1 = AAString(gsub(" ", "", temp$sequence[1]))
    }
    
    temp <- dbGetQuery(conn, paste0("select sequence from proteins where id='", p2,"';"))
    if(! is.na(temp$sequence[1])){
        aas2 = AAString(gsub(" ", "", temp$sequence[1]))
    }
    idd=0
    if( !is.null(aas1) && !is.null(aas2)){
        aln=pairwiseAlignment(aas1, aas2)
        idd=pid(aln)
        dfp=data.frame(idp1=c(p1), idp2=c(p2), identity=c(idd) )
        dbWriteTable(conn, "proteins_alignment", dfp, append = TRUE)
    }
    return(idd)
}

# Loading proteins from user's input file
process_protein_list <- function(nameFile, conn){
  proteins = read_lines(nameFile)
  folder="proteins/"
  for (p in proteins){
    if(!file.exists(paste0(folder, p, ".rdf"))){
        download.file(paste0('https://www.uniprot.org/uniprot/', p, ".rdf"), paste0(folder, p, ".rdf"), mode="wb")
      
        # feed proteins in database with sequence
        get_feed_info_protein(conn, p)
        
        # load mappings
        load_rdf_feed_annotation(conn, folder, p)
        
    }    
      
        # Get ligands inference
        #get_inference(conn, p)
    
  }
}

get_sequence <- function(p){
    aas = NULL
    folder="sequences/"
    if(!file.exists(paste0(folder, p, ".fasta"))){
        download.file(paste0('https://www.uniprot.org/uniprot/', p, ".fasta"), paste0(folder, p, ".fasta"), mode="wb")
    }
    if(file.info(paste0(folder, p, ".fasta"))$size!=0){
        aa=readAAStringSet(paste0(folder, p, ".fasta"))
        aas=aa[[1]]
    }
    return(aas)
}

get_feed_info_protein <- function(conn, p) {
    temp <- dbGetQuery(conn, paste0("select count(*) as cnt from proteins where id='", p,"';"))
    if(temp$cnt==0){
        obj=get_sequence(p)
        if(! is.null(obj)){
            seq=gsub(" ", "", toString( obj ))
            dfp=data.frame(id=c(p), sequence=c(seq) )
            dbWriteTable(conn, "proteins", dfp, append = TRUE)
        }
    }
}

load_rdf_feed_annotation <- function(conn, folder, p){
    info <- c('pfam', 'pdbsum', 'drugbank')
    annotations = read_lines( paste0(folder, p, ".rdf") )
    for (line in annotations){
        for (i in info){
            if( grepl( paste0('/',i,'/'), line, fixed = TRUE) && grepl( 'seeAlso', line, fixed = TRUE) ){
                ann=strsplit(line, '=')[[1]][[2]]
                ann= gsub(paste0("\"http://purl.uniprot.org/",i,"/"), "", ann )
                ann= gsub(paste0("\"/>"), "", ann )
                
                df=data.frame(protein=p, annotation=ann)
                dbWriteTable(conn, paste0("map_uniprot_", i), df, append = TRUE)
                temp <- dbGetQuery(conn, paste0("select distinct protein from map_uniprot_", i, ";"))
                
                if(i=='pdbsum'){
                    get_ligand_from_pdb2(conn, ann)
                }
            }   
        }
    }
}

load_rdf_feed_annotation2 <- function(conn, folder, p){

      # Load rdf file rdf file
  world <- new("World")
  storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
  model <- new("Model", world=world, storage, options="")
  parser <- new("Parser", world)
  parseFileIntoModel(parser, world, paste0(folder, p, ".rdf"), model)
  
  # Obtain annotation information from triples with resource being pfam, pdb or drugbank
  info <- c('pfam', 'pdbsum', 'drugbank')
  #info <- c('drugbank')
  for (i in info){
    query <- paste0("
      prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
      select ?protein ?info where {
        ?protein rdfs:seeAlso ?info .
        
        filter ( regex(str(?info), '/", i, "/')  ) .
      }
      ")
    query <- new("Query", world, query, base_uri=NULL, query_language="sparql", query_uri=NULL)
    queryResult <- executeQuery(query, model)
    results <- getResults(query, model, "json")
    if(! is.null(results) ){
      resultsj <- fromJSON(results) 
      resultsj <- resultsj$results$bindings %>% as.data.frame
      j=0
      proteins <- c()
      annotations <- c()
      for (ann in resultsj$info$value){
        a=strsplit(ann, '/')[[1]]
        ann=a[length(a)]
        proteins[j] <- p
        annotations[j] <- ann
        
        if(i=='pdbsum'){
            get_ligand_from_pdb2(conn, ann)
        }
        
        j=j+1
      }
      
      df=data.frame(protein=proteins, annotation=annotations)
      dbWriteTable(conn, paste0("map_uniprot_", i), df, append = TRUE)
      temp <- dbGetQuery(conn, paste0("select distinct protein from map_uniprot_", i, ";"))
      #print(temp)
    }
  }
  
}

# api pdb https://data.rcsb.org/redoc/index.html#tag/Entry-Service/operation/getEntryById

get_ligand_from_pdb2 <- function(conn, pdb){
    r2 <- GET( paste0("https://data.rcsb.org/rest/v1/core/entry/", pdb), accept_json(), add_headers(Accept = 'application/json') )
    nrl=c()
    i=1
    if(! is.null(content(r2)) ){
        ligands=content(r2)$pdbx_vrpt_summary$restypes_notchecked_for_bond_angle_geometry
        if(! is.null(ligands) ){
            if( length(ligands)>0 ){
                for ( l in ligands ){
                    if(!is.null(l) && ! l %in% nrl && !is.na(l)){
                        nrl[i]=l
                        i=i+1
                    }
                }
            }
        }
        
        ligands=content(r2)$rcsb_binding_affinity[[1]]$comp_id
        if(! is.null(ligands) ){
            if( length(ligands)>0 ){
                for ( l in ligands ){
                    if(!is.na(l) && !is.null(l) && ! l %in% nrl ){
                        nrl[i]=l
                        i=i+1
                    }
                }
            }
        }
    }
   # print(paste('ligs', nrl))
    
    if( length(nrl)>0 ){
        for (l in nrl){
            temp <- dbGetQuery(conn, paste0("select count(*) as cnt from pdb_ligand where pdb='", pdb,"' and ligand='", l,"';"))
            if(temp$cnt==0 && !is.na(l)){
                dfp=data.frame(pdb=c(pdb), ligand=c(l) )
                dbWriteTable(conn, "pdb_ligand", dfp, append = TRUE)
                
                #download_ligand(l)
                #if(file.exists( paste0("ligands/", l, "_ideal.sdf") ) ){
                get_info_ligand2(conn, l)
                #}
            }
        }
    }
}
#get_ligand_from_pdb2(conn, "7LYB")
#get_ligand_from_pdb2(conn, "1GIH")

#get_info_ligand2(conn, "1UP")
get_info_ligand2 <- function(conn, l){
    chemblid=''
    smi=''
    r2 <- GET( paste0("https://data.rcsb.org/rest/v1/core/chemcomp/", l), accept_json(), add_headers(Accept = 'application/json') )
    
    if(! is.null(content(r2)) ){
        if( length(content(r2)$pdbx_chem_comp_descriptor)>0 ){
            for (d in content(r2)$pdbx_chem_comp_descriptor){
                if(d$type=="SMILES"){
                    smi=d$descriptor    
                }
            }
        }
        if( length(content(r2)$rcsb_chem_comp_related)>0 ){
            for (d in content(r2)$rcsb_chem_comp_related){
                if(d$resource_name=="ChEMBL"){
                    chemblid=d$resource_accession_code    
                }
            }
        }
    }
    
    dfp=data.frame(ligand=c(l), smiles=c(smi), chembl=c(chemblid) )
    dbWriteTable(conn, "ligands_info", dfp, append = TRUE)
}

align_ligands2 <-function(conn, l1, l2){
    sim=0
    ap1=NULL
    ap2=NULL
    
    temp <- dbGetQuery(conn, paste0("select smiles from ligands_info where ligand='", l1,"';"))
    if(length(temp$smiles)!=0){
        sd=smiles2sdf(temp$smiles[1])
        tryCatch(  {
        ap1=sdf2ap(sd[[1]])
            },
            error= function(e){
                dfp=data.frame(ligand=c(l1) )
                dbWriteTable(conn, "ligands_broken", dfp, append = TRUE)
                print("Bad molecule formula")
            }
        )
    }
    
    temp <- dbGetQuery(conn, paste0("select smiles from ligands_info where ligand='", l2,"';"))
    if(length(temp$smiles)!=0){
        sd=smiles2sdf(temp$smiles[1])
        tryCatch(  {
        ap2=sdf2ap(sd[[1]])
            },
            error= function(e){
                dfp=data.frame(ligand=c(l2) )
                dbWriteTable(conn, "ligands_broken", dfp, append = TRUE)
                print("Bad molecule formula")
            }
        )
    }
    
    if( !is.null(ap1) && !is.null(ap2) ){ 
        sim=cmp.similarity(ap1, ap2)
        dfp=data.frame(ligand1=c(l1), ligand2=c(l2), tanimoto=c(sim) )
        dbWriteTable(conn, "ligands_alignment", dfp, append = TRUE)
            
    }
    return(sim)
}

library('visNetwork') 
get_network <- function(conn, nameFile, tanimoto_cutoff, idu, protein){
    proteins = read_lines(nameFile)
    if(protein!="All"){
        proteins=c(protein)
    }
    
    inde=1
    edges_to=c()
    edges_from=c()
    edges_titles=c()
    
    ind=1
    ids=c()
    titles=c() # label equals to title for now, after use it to show smiles in title and ligand id as label
    node_type=c()
    node_type_id=c()
    
    print("Get edges and nodes from protein annotations")
    # Get edges and nodes from protein annotations
    for (p in proteins){
        ids[ind]=p #paste0('p',ind)
        titles[ind]=p
        node_type[ind]="Protein"
        node_type_id[ind]=1
        ind=ind+1
        
        ts=c("Family", "Structure", "Drug")
        init=c('f', 's', 'd')
        j=1
        info <- c('pfam', 'pdbsum', 'drugbank')
        for (i in info){
            temp <- dbGetQuery(conn, paste0("select distinct annotation from map_uniprot_", i, " where protein='", p,"';"))
            for (t in temp$annotation){
                
                if(! t %in% titles){
                    ids[ind]=t #paste0(init[j],ind)
                    titles[ind]=t
                    node_type[ind]=ts[j]
                    node_type_id[ind]=j+1
                    ind=ind+1
                }
                edges_from[inde]=p #paste0('p',ind)
                edges_to[inde]=t #paste0(init[j],ind-1)
                edges_titles[inde]=''
                inde=inde+1
            }
            j=j+1
        }
    }
    
    print("Get edges from protein sequence alignment")
    # Get edges from protein sequence alignment
    i=0
    for (p1 in proteins) {
        k=0
        for (p2 in proteins){
            if(p1!=p2 && i<k){
                temp <- dbGetQuery(conn, paste0("select identity from proteins_alignment where idp1='", p1,"' and idp2='", p2,"';"))
                if( length(temp$identity) == 0){
                    idd=process_protein_sequence_aln(conn, p1, p2)
                }
                else{
                    idd=temp$identity[1]
                }
                
                edges_from[inde]=p1 #paste0('p',ind)
                edges_to[inde]=p2
                edges_titles[inde]=idd/100
                inde=inde+1
            }
            k=k+1
        }
        i=i+1
    }
    labels=titles
    
    print("Get edges from ligands alignment")
    # Get edges from ligands alignment
        # get list of non redundant ligands for all structures belonging to the seed proteins
    ligs=c()
    k=0
    for (p1 in proteins) {
        # map_uniprot_pdbsum, pdb_ligand, ligands_alignment
        temp <- dbGetQuery(conn, paste0("select distinct p.ligand, p.pdb, chembl, smiles from map_uniprot_pdbsum as m, pdb_ligand as p, ligands_info as l where m.annotation=p.pdb and l.ligand=p.ligand and m.protein='", p1,"' ;"))
        if(length(temp$ligand)>0){
            for (h in 1:length(temp$ligand)){
                if( ! is.na(temp$ligand[h]) ){
                    if( ! temp$ligand[h] %in% ligs){
                        val1 <- dbGetQuery(conn, paste0("select count(ligand) as cnt from ligands_broken where ligand='", temp$ligand[h],"' ;" ))
                        if(val1$cnt==0){
                            ligs[k]=temp$ligand[h] 
                        }
                    }
                    
                    if( ! temp$ligand[h] %in% ids){
                        ids[ind]=temp$ligand[h] #paste0(init[j],ind)
                        titles[ind]=temp$smiles[h]
                        labels[ind]=temp$chembl[h]
                        node_type[ind]="Ligand"
                        node_type_id[ind]=j+1
                        ind=ind+1 
                    }
                    k=k+1
                    
                    edges_from[inde]=temp$ligand[h] #paste0('p',ind)
                    edges_to[inde]=temp$pdb[h]
                    edges_titles[inde]=""
                    inde=inde+1
                }
            }
        }
    }
    
    
    if(FALSE){
        print("Get the tanimoto score of the ligands")
        # Get the tanimoto score of the ligands in database or calculate from loaded sdfs
        i=0
        for (l1 in ligs) {
            j=0
            for (l2 in ligs){ 
                if(!is.na(l1) && !is.na(l2)){
                    if(l1!=l2 && i<j){ 
                        idd=0
                        temp <- dbGetQuery(conn, paste0("select distinct ligand1, ligand2, tanimoto from ligands_alignment where ligand1='", l1,"' and ligand2='", l2, "';" ))
                        if( length(temp$tanimoto) == 0){
                            idd=align_ligands2(conn, l1, l2)
                        }
                        else{
                            idd=temp$tanimoto[1]
                        }
                        if(idd>=tanimoto_cutoff){
                            edges_from[inde]=temp$ligand1[1] #paste0('p',ind)
                            edges_to[inde]=temp$ligand2[1]
                            edges_titles[inde]=idd
                            inde=inde+1
                        }
                    }
                }
                j=j+1
            }
            i=i+1
        }
    }
    
    # get the inferred ligands and pdbs
    # for (p1 in proteins) {
    #     temp <- dbGetQuery(conn, paste0("select * from pdb_ligands_inferred as p, ligands_info as l where uniprot='", p1,"' and p.ligand=l.ligand;"))
    #     temp<- na.omit(temp)
    #     if(length(temp$ligand)>0){
    #         for (h in 1:nrow(temp)){
    #             if( ! temp$ligand[h] %in% ids){
    #                 ids[ind]=temp$ligand[h] #paste0(init[j],ind)
    #                 titles[ind]=temp$smiles[h]
    #                 labels[ind]=temp$chembl[h]
    #                 node_type[ind]="Ligand"
    #                 node_type_id[ind]=j+1
    #                 ind=ind+1 
    #             }
    #             
    #             if( ! temp$pdb[h] %in% ids){
    #                 ids[ind]=temp$pdb[h] #paste0(init[j],ind)
    #                 titles[ind]=temp$pdb[h]
    #                 labels[ind]=temp$pdb[h]
    #                 node_type[ind]="Structure"
    #                 node_type_id[ind]=j+1
    #                 ind=ind+1 
    #             }
    #             
    #             if(temp$number_contact[h]>=num_contact){
    #                 edges_from[inde]=p1 #paste0('p',ind)
    #                 edges_to[inde]=temp$ligand[h]
    #                 edges_titles[inde]=temp$number_contact[h]
    #                 inde=inde+1
    #             }
    #             
    #             if(temp$percent_contact[h]>=perc_contact){
    #                 edges_from[inde]=p1 #paste0('p',ind)
    #                 edges_to[inde]=temp$pdb[h]
    #                 edges_titles[inde]=temp$percent_contact[h]
    #                 inde=inde+1
    #             }
    #         }
    #     }
    # }
    
    print("Preparing table nodes")
    types=c("Protein", "Family", "Structure", "Drug","Ligand")
    colors=c("cadetblue", "darkolivegreen", "lightcoral", "lavenderblush", "khaki")
    # color options: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    # filter expr[expr$cell_type == "hesc", ] or expr[expr$cell_type %in% c("hesc", "bj fibroblast"), ]
    
    # filter_nodes=c()
    # for(t in types){
    #     i=1
    #     j=1
    #     for(d in ids){
    #         if(node_type[i]==t && j<=100){
    #             filter_nodes[j] = d
    #             j=j+1
    #         }
    #         i=i+1
    #     }
    # }
    
    nodes <- data.frame(id=ids, title=titles, label=labels, type=node_type, type_id=node_type_id)
    nodes$shape <- "dot"  
    nodes$shadow <- TRUE
    nodes$borderWidth <- 1
    nodes$color.background <- colors[node_type_id]
    nodes$color.border <- "black"
    
    #nodes$color.highlight.background <- "lighblue"
    #nodes$color.highlight.border <- "darkblue"
    
    print("Preparing table edges")
    edges=data.frame(from=edges_from, to=edges_to)
    edges$width=1
    edges$value=edges_titles
    edges$title=edges_titles
    edges$label=edges_titles
    edges$color="black"

    #write.csv(nodes, paste0(idu,"/nodes.csv") )
    #write.csv(edges, paste0(idu,"/edges.csv") )
    
    #nodes=nodes[nodes$id %in% filter_nodes, ]
    #edges=edges[ which(edges$from %in% filter_nodes & edges$to %in% filter_nodes) , ]
    
    return( list(nodes, edges) )
}

render_network<-function(nodes, edges){
    visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE,  selectedBy = "type")
    
    #return (list(nodes, edges))
}

#l<-"JY"
#download_ligand(l)
download_ligand <- function(l){
    if(!file.exists( paste0("ligands/", l, "_ideal.sdf") ) ){
        r <- GET( paste0("https://wwwdev.ebi.ac.uk/pdbe/download/api/pdb/compound/sdf?combined=false&conformer=ideal&id=", l), accept_json(), add_headers(Accept = 'application/json') )
        if(! is.null(content(r)) ){
            print(content(r))
            if( ! grepl("503 Service", content(r), fixed=TRUE) ){
                url<-content(r)$url
                download.file(url, "temp.tar.gz")
                untar("temp.tar.gz", files=paste0("temp/", l, "_ideal.sdf"))
                file.copy( paste0("temp/", l, "_ideal.sdf"), 'ligands')
                unlink('temp', recursive=TRUE)
                file.remove('temp.tar.gz')
            }
        }
    }
}

get_info_ligand <- function(conn, l){
    if(file.exists( paste0("ligands/", l, "_ideal.sdf") ) ){
        sdfset <- read.SDFset( paste0("ligands/", l, "_ideal.sdf") )
        smiles <- sdf2smiles(sdfset)
        smi=smiles[[1]]
        
        chemblid=''
        r2 <- GET( paste0("https://www.ebi.ac.uk/pdbe/api/pdb/compound/mappings/", l), accept_json(), add_headers(Accept = 'application/json') )
        if( is.null(content(r2)) ){
            chemblid=content(r2)[[l]][[1]]$chembl_id
        }
        dfp=data.frame(ligand=c(l), smiles=c(smi), chembl=c(chemblid) )
        dbWriteTable(conn, "ligands_info", dfp, append = TRUE)
    }
}

align_ligands <-function(conn, l1, l2){
    sim=0
    if( file.exists( paste0("ligands/", l1, "_ideal.sdf") ) &&   file.exists( paste0("ligands/", l2, "_ideal.sdf") ) ){
        sdfset <- read.SDFset( paste0("ligands/", l1, "_ideal.sdf") )
        ap1 <- sdf2ap(sdfset[[1]])
        
        sdfset <- read.SDFset( paste0("ligands/", l2, "_ideal.sdf") )
        app2 <- sdf2ap(sdfset[[1]])
        
        sim=cmp.similarity(ap1, ap2)
        dfp=data.frame(ligand1=c(l1), ligand2=c(l2), tanimoto=c(sim) )
        dbWriteTable(conn, "ligands_alignment", dfp, append = TRUE)
    }
    return(sim)
}
