library(bio3d)
library(seqinr)
library(Biostrings)

# '7JZV', 'A', 'PF12820', 10, 'MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRRTPSTEKKVDLNADPLCERKEWNKQKLPCSENPRDTEDVPWITLNSSIQKVNEWFSRSDELLGSDDSHDGESESNAKVADVLDVLNEVDEYSGSSEKIDLLASDPHEALICKSERVHSKSVESNIEDKIFGKTYRKKASLPNLSHVTENLIIGAFVTEPQIIQERPLTNKLKRKRRPTSGLHPEDFIKKADLAVQKTPEMINQGTNQTEQNGQVMNITNSGHENKTKGDSIQNEKNPNPIESLEKESAFKTKAEPISSSISNMELELNIHNSKAPKKNRLRRKSSTRHIHALELVVSRNLSPPNCTELQIDSCSSSEEIKKKKYNQMPVRHSRNLQLMEGKEPATGAKKSNKPNEQTSKRHDSDTFPELKLTNAPGSFTKCSNTSELKEFVNPSLPREEKEEKLETVKVSNNAEDPKDLMLSGERVLQTERSVESSSISLVPGTDYGTQESISLLEVSTLGKAKTEPNKCVSQCAAFENPKGLIHGCSKDNRNDTEGFKYPLGHEVNHSRETSIEMEESELDAQYLQNTFKVSKRQSFAPFSNPGNAEEECATFSAHSGSLKKQSPKVTFECEQKEENQGKNESNIKPVQTVNITAGFPVVGQKDKPVDNAKCSIKGGSRFCLSSQFRGNETGLITPNKHGLLQNPYRIPPLFPIKSFVKTKCKKNLLEENFEEHSMSPEREMGNENIPSTVSTISRNNIRENVFKEASSSNINEVGSSTNEVGSSINEIGSSDENIQAELGRNRGPKLNAMLRLGVLQPEVYKQSLPGSNCKHPEIKKQEYEEVVQTVNTDFSPYLISDNLEQPMGSSHASQVCSETPDDLLDDGEIKEDTSFAENDIKESSAVFSKSVQKGELSRSPSPFTHTHLAQGYRRGAKKLESSEENLSSEDEELPCFQHLLFGKVNNIPSQSTRHSTVATECLSKNTEENLLSLKNSLNDCSNQVILAKASQEHHLSEETKCSASLFSSQCSELEDLTANTNTQDPFLIGSSKQMRHQSESQGVGLSDKELVSDDEERGTGLEENNQEEQSMDSNLGEAASGCESETSVSEDCSGLSSQSDILTTQQRDTMQHNLIKLQQEMAELEAVLEQHGSQPSNSYPSIISDSSALEDLRNPEQSTSEKAVLTSQKSSEYPISQNPEGLSADKFEVSADSSTSKNKEPGVERSSPSKCPSLDDRWYMHSCSGSLQNRNYPSQEELIKVVDVEEQQLEESGPHDLTETSYLPRQDLEGTPYLESGISLFSDDPESDPSEDRAPESARVGNIPSSTSALKVPQLKVAESAQSPAAAHTTDTAGYNAMEESVSREKPELTASTERVNKRMSMVVSGLTPEEFMLVYKFARKHHITLTNLITEETTHVVMKTDAEFVCERTLKYFLGIAGGKWVVSYFWVTQSIKERKMLNEHDFEVRGDVVNGRNHQGPKRARESQDRKIFRGLEICCYGPFTNMPTDQLEWMVQLCGASVVKELSSFTLGTGVHPIVVVQPDAWTEDNGFHAIGQMCEAPVVTREWVLDSVALYQCQELDTYLIPQIPHSHY'
pfam_map <- read.table("pdb_pfam_mapping.txt", skip=1, sep="\t", header = 1)
garbage <- c("ALA","LYS","VAL","ARG","GLY","CYS","TYR", "PHE","LEU","ILE","MET",
             "SER","PRO","THR","HIS","ASN","GLN","GLU","ASP","TRP","NAD","NAI","NDC",
             "TAP","ADJ","NAJ","NDO","ZID","CAN","NAP","NDP",
             "CND","NAQ","NHD","DND","NAX","NHO","NAC","NDB","ODP","NAE","NBP",
             "PAD","NAH","NDA","SND", "FAD","FMN","6FA","FNS","FAA","MGD","FAB","RFL","FAE","FAS","FDA",
             "FMA", "AMP","GDP","UDP","ADP","GNP","UTP","ANP", "GTP","PSU",
             "ATP","2GP","CMP","TMP","CDP","TDP","CTP","TTP","GMP","UMP","SO4", "UNK",
             "UNX", "BMA","FUC","MAN","POP","BOG","GAL","MES","PYR","C8E","GLC","MPD",
             "SPM","CIT","GOL","MYR","TRS","CRY","HED","NAG","XYS","DTT","LDA","NGA",
             "EPE","LI1","PEG","F6P","MAL","PG4", "HEM","CD","NA", "MG", "ZN", "CL",
             "MN", "FE", "HOH", "PO4", "CA", "EDO")

contacts_with_domain <- function(pdb_name, ch, pfam_domain, contact_cutoff, seq_prot){
  
  # pdb_name = pdb entry code
  # ch = chain
  # pfam_domain = pfam entry
  # contact cutoff = cutoff in Angstroms for a residue with an atom within 
  # the range to be considered a contact
  # seq_prot = protein sequence
  
  result <- data.frame(pdb=NA, ligand=NA, percent_conserved=NA,
                       pfam=NA,domain_contacts=NA)
  
  pdb <- tryCatch(read.pdb(pdb_name), error=function(e) NULL)
  
  if(is.null(pdb) != TRUE){
    ind_pdb <- atom.select(pdb, chain=ch, type="ATOM", elety="CA")
    sequence <- paste(pdbseq(pdb, ind_pdb), collapse = "")
    #print( paste(pdb_name, pfam_domain, ch))
    subset <- pfam_map[which(pfam_map$PDB == tolower(pdb_name) & pfam_map$PFAM_ACCESSION == pfam_domain & 
                               pfam_map$CHAIN == ch),]
    
    ### trim and create distance matrix
    
    pdb <- trim.pdb(pdb, atom.select(pdb, chain=ch))
    dmat <- dm(pdb, mask.lower=FALSE, grp=FALSE, all.atom=TRUE)
    ligand_contacts <- data.frame(dmat[,])
    
    pdb_atom <- pdb$atom
    
    pdb_het <- pdb_atom[which(pdb_atom$type == "HETATM"),]
   
    if(nrow(pdb_het) != 0){
      pdb_het <- pdb_het[-which(pdb_het$resid %in% garbage),]
      ligands <- unique(sort(pdb_het$resid))
      
      if(length(ligands) != 0){
        
        for(lig in ligands){
          ### store data
          pdb_het_lig <- pdb_het[which(pdb_het$resid == lig),]
          
          tmp <- which(pdb_atom$resid == lig)
          
          cont <- which(ligand_contacts[,tmp] <= contact_cutoff)
          in_cont <- unique(sort(pdb_atom$resno[cont]))
          ### remove ligand from index list
          in_cont <- in_cont[-which(in_cont == unique(sort(
            pdb_atom$resno[which(pdb_atom$resid == lig)])))]
          
          
          ### align PDB and ref-sequence
          pdb_seq <- paste(pdbseq(pdb), collapse = "")
          
          aas1 <- AAString(pdb_seq)
          aas2 <- AAString(seq_prot)
          
          aln<-pairwiseAlignment(aas1, aas2)
          s1=pattern(aln)
          s2=subject(aln)
          
          m1 <- strsplit(as.character(s1), split="")[[1]][subset[1,3]:subset[1,4]]
          m2 <- strsplit(as.character(s2), split="")[[1]][subset[1,3]:subset[1,4]]
          
          percent_con <- length(which(m1 == m2 & m2 != '-'))/length(m1[which(m1 != "-")])
          
          tmp <- data.frame(pdb=pdb_name, ligand=lig, percent_conserved=percent_con,
                           pfam=pfam_domain,domain_contacts=length(in_cont))
          result <- rbind(result, tmp) 
          #return(tmp)
        }
      }
      else{
        print("No ligands")
      } 
    } 
    else{
      print("no HETATM")
    }
    
    result <- na.omit(result)
    if(nrow(result) != 0){
      return(result)
    }
    
  }
  else{
    return("pdb not found")
  }
}

get_pdb_info <- function(pfams){
    
    ### subj from format_querry_uniprot
    pfam_vec_final <- c()
    
    pfam_vec <- pfams
    pdb_vec <- pfam_map$PDB[which(pfam_map$PFAM_ACCESSION %in% pfam_vec)]
    chain_vec <- pfam_map$CHAIN[which(pfam_map$PFAM_ACCESSION %in% pfam_vec)]
    
    result <- list(cid, pfam_vec, pdb_vec, chain_vec)
    return(result) 
}

#pfams<- c("PF02782")

get_inference <- function(conn, cid){
    temp=dbGetQuery(conn, paste0("select annotation from map_uniprot_pfam where protein='",cid,"';"))
    pfams=c()
    i=1
    for (t in temp$annotation){
        pfams[i]=t
        i=i+1
    }
    
    temp=dbGetQuery(conn, paste0("select sequence from proteins where id='",cid,"';"))
    if(length(pfams)>0 && nrow(temp)>0){
        seq=temp$sequence[1]
        
        pdb_stuff <- get_pdb_info(pfams)
        #print(pdb_stuff[[3]])
        gone=c()
        i=1
        
        final_res <- data.frame(pdb=NA, ligand=NA, percent_conserved=NA, pfam=NA, domain_contacts=NA)
        for(current in 1:length(pdb_stuff[[3]])){
            
            pdb_name <- toupper(pdb_stuff[[3]][current])
            if(! pdb_name %in% gone){
                gone[i]=pdb_name
                i=i+1
            
                ch <- toString(pdb_stuff[[4]][current])
                pfam_domain <- pdb_stuff[[2]]
                
                result <- contacts_with_domain(pdb_name, ch, pfam_domain, 10, seq )
                final_res<-rbind(final_res, result)
            }
        }
        #print(final_res)
        
        ligs=c()
        j=1
        for (i in 1:nrow(final_res)){
            if(!final_res$ligand[i] %in% ligs){
                if( !is.na(final_res$ligand[[i]]) ){
                    ligs[j]=final_res$ligand[[i]]
                    get_info_ligand2(conn, ligs[j])
                    j=j+1
                }
            }
            
            dfp=data.frame(uniprot=c(cid), pfam=c(final_res$pfam[i]), pdb=c(final_res$pdb[i]), ligand=c(final_res$ligand[i]), percent_contact=c(final_res$percent_conserved[i]), number_contact=c(final_res$domain_contacts, percent_contact=c(final_res$percent_conserved[i])[i]) )
            dbWriteTable(conn, "pdb_ligands_inferred", dfp, append = TRUE)
        }
    }
    #return(final_res)
}

#res <- get_inference(conn, "A0A7I0DYN5")
