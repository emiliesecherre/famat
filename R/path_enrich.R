#' @export path_enrich
#############################FUNCTIONS###########################
##perform pathway enrichment analysis on user's gene list
gene_path<-function(gene, db){
    gostres<-gprofiler2::gost(gene, "hsapiens", sources=db, significant=FALSE)
    res<-gostres$result

    ##only keep pathway ids
    res[,9]<-stringr::str_sub(res[, 9], 6)
    res=res[,c(11,9)];names(res)=c("name","id")
    return(res)
}

##perform pathway enrichment analysis on user's metabolites list
meta_path<-function(meta, db){
    ##convert kegg ids in pubchem ids
    meta_pubchem<-mapper[which(mapper[, 1] %in% meta), 2]

    ##perform pathway enrichment analysis on user's metabolites list
    pss<-MPINet::getPSS(meta_pubchem, plot=FALSE)
    if(db == "KEGG"){
        medium<-MPINet::identifypathway(meta_pubchem, pss, "KEGG")
    }
    else if(db == "REAC"){
        medium<-MPINet::identifypathway(meta_pubchem, pss, "Reactome")
    }
    else if (db == "WP"){
        medium<-MPINet::identifypathway(meta_pubchem, pss, "Wikipathways")
    }

    ##extract informations about pathways from MPINet results
    resm<-vector()
    for (i in seq_len(length(medium))) {
        path<-medium[[i]]
        if(db == "KEGG"){
            resm[i]<-path[[2]][[1]]
        }
        else if (db == "REAC" || db == "WP"){
            resm[i]<-path[[1]][[1]]
        }
    }
    resm=as.data.frame(resm);names(resm)="name"
    return(resm)
}

##put pathway enrichment analysis results together and add pathways ids
path_enrich<-function(source, metabo, genes){
    ##genes and metabolites lists
    genes_list<-genes[, 1]
    meta_list<-metabo[, 1]
    ##perform pathway enrichment analysis on user's gene list
    resgene<-gene_path(genes_list, source)
    ##perform pathway enrichment analysis on user's metabolites list
    resmeta<-meta_path(meta_list, source)
    ##add pathways ids
    resmeta<-cbind(resmeta, id=rep(NA, nrow(resmeta)))
    if(source == "KEGG"){
        prev_db<-KEGGREST::keggList("pathway", "hsa")
        idskegg=stringr::str_sub(names(prev_db), 9, nchar(names(prev_db)))
        idskegg=paste("hsa:",idskegg,sep="")
        keggdb<-as.list(idskegg)
        nameskegg<-unlist(stringr::str_split(unname(prev_db),
                                            " - Homo sapiens (.)human(.)"))
        names(keggdb)<-nameskegg[!(nameskegg %in% "")]

        resgene$id<-paste("hsa:", resgene$id, sep="")
        for (n in seq_len(nrow(resmeta))){
            id<-keggdb[[resmeta[n,1]]]
            if(!is.null(id)){resmeta[n, 2]<-id}
        }
    }
    else if(source == "REAC"){
        readb<-unlist(as.list(reactome.db::reactomePATHID2NAME))
        readb<-readb[which(stringr::str_sub(readb, 1, 14)%in%"Homo sapiens: ")]
        pathids=names(readb);readb<-stringr::str_sub(readb, 15, nchar(readb))
        names(readb)=pathids
        for (n in seq_len(nrow(resmeta))){
            id<-names(readb[which(readb %in% resmeta[n, 1])])
            if (length(id)>0){resmeta[n, 2]<-id}
        }
    }
    else if (source == "WP"){
        wpdb<-rWikiPathways::listPathways("Homo sapiens")
        wpdb<-wpdb[, c(3, 1)]
        for (r in seq_len(nrow(resmeta))){
            pathname=stringr::str_to_upper(resmeta[r, 1])
            id<-wpdb[which(pathname == stringr::str_to_upper(wpdb$name)), 2]
            if(length(id)>0){resmeta[r, 2]<-id[1]}
        }
        resgene$id<-paste("WP", resgene$id, sep="")
    }
    return(list(resmeta, resgene, genes, metabo))
}

