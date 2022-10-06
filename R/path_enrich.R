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
    ##extract informations about pathways from MPINet results
    #pss<-getPSS(meta_pubchem, plot=FALSE)#removed MPINet::
    if(db == "KEGG"){
        medium<-identifypathway(meta_pubchem, "KEGG")#removed MPINet::
        resm<-vapply(medium, function(x){x[[2]][[1]]}, character(1))
    }else if(db == "REAC"){
        medium<-identifypathway(meta_pubchem, "Reactome")#removed MPINet::
        resm<-vapply(medium, function(x){x[[1]][[1]]}, character(1))
    }else if (db == "WP"){
        medium<-identifypathway(meta_pubchem, "Wikipathways")#removed MPINet::
        resm<-vapply(medium, function(x){x[[1]][[1]]}, character(1))
    }

    resm=as.data.frame(resm);names(resm)="name"
    return(resm)
}

##put pathway enrichment analysis results together and add pathways ids
path_enrich<-function(source, metabo, genes){
    ##perform pathway enrichment analysis on user's gene list
    resgene<-gene_path(genes, source)
    ##perform pathway enrichment analysis on user's metabolites list
    resmeta<-meta_path(metabo, source)
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
        resmeta[,2]=apply(resmeta, 1, function(x){
            id=keggdb[[x[1]]]
            if(!is.null(id)){id}
            else{NA}
        })
    }
    else if(source == "REAC"){
        readb<-unlist(as.list(reactome.db::reactomePATHID2NAME))
        readb<-readb[which(stringr::str_sub(readb, 1, 14)%in%"Homo sapiens: ")]
        pathids=names(readb);readb<-stringr::str_sub(readb, 15, nchar(readb))
        names(readb)=pathids
        resmeta[,2]=apply(resmeta, 1, function(x){
            id=names(readb[which(readb %in% x[1])])
            if(length(id)>0){id}
            else{NA}
        })
    }
    else if (source == "WP"){
        wpdb<-rWikiPathways::listPathways("Homo sapiens")
        wpdb<-wpdb[, c(3, 1)]
        ids=apply(resmeta, 1, function(x){
            pathname=stringr::str_to_upper(x[1])
            id<-wpdb[which(pathname == stringr::str_to_upper(wpdb$name)), 2]
            if(length(id)>0){id[1]}
            else{NA}
        })
        resmeta[,2]=unname(unlist(ids))
        resgene$id<-paste("WP", resgene$id, sep="")
    }
    return(list(resmeta, resgene, genes, metabo, resmeta))
}
