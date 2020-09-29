#' @export interactions
##############################FONCTIONS ################################
##remove duplicated values in a vector
rm_vector<-function(vec){
    dp<-which(duplicated(vec))
    if (length(dp)>0){
        vec<-vec[-dp]
    }
    return(vec)
}

##remove duplicated values in a dataframe
rm_df<-function(df, col){
    dp<-which(duplicated(df[, col]))
    if (length(dp)>0){
        df<-df[-dp, ]
    }
    return(df)
}

##gather all pathways obtained in pathway enrichment analysis
get_pathways<-function(resgk, resgr, resgw, resmk, resmr, resmw){
    pathway_tot<-data.frame(name=character(), id=character())
    pathway_tot<-rbind(pathway_tot,resmk,resmr,resmw,resgk,resgr,resgw)
    pathway_tot<-rm_df(pathway_tot,c(seq_len(ncol(pathway_tot))))
    pathway_tot<-pathway_tot[!is.na(pathway_tot[,2]),]
    pathway_tot<-pathway_tot[!(pathway_tot[,2] %in% c("0000000","hsa:00000")),]
    return(pathway_tot)
}

##function of data_size, fills the dataframe
elem_size<-function(size, i, j, user_elem, path_elem){
    size[i, j+2]<-length(user_elem)
    if (length(user_elem)>0){
        size[i, j+3]<-paste(user_elem, collapse=" # ")
        size[i, j]<-length(intersect(user_elem, path_elem))
        if (length(intersect(user_elem, path_elem))>0){
            elements=intersect(user_elem, path_elem)
            size[i, j+1]<-paste(elements, collapse=" # ")
        }
        else{size[i, j+1]<-NA}
    }
    else{
        size[i, j+3]<-NA
        size[i, j]<-0
    }
    return(size)
}

#dataframe showing pathways and which genes/metabolites they have
data_size<-function(inter, meta_list, gene_list, mapkegg, pathk, pathr, pathp){
    size<-data.frame(path=character(), nb_gene_query=integer(),
                        gene_que=character(), nb_gene_tot=integer(),
                        genes=character(), nb_meta_query=integer(),
                        meta_que=character(), nb_meta_tot=integer(),
                        meta=character())

    pathtot_ids<-c(pathk$id, pathr$id, pathp$id)
    for (k in seq_len(length(pathtot_ids))){
        path_inter<-inter[which(inter$path == pathtot_ids[k]),]
        elements<-rm_vector(c(path_inter[, 1], path_inter[, 3]))

        meta<-elements[which(elements %in% mapkegg$name)]
        gene<-elements[which(!(elements %in% mapkegg$name))]

        size[k, 1]<-pathtot_ids[k]

        size<-elem_size(size, k, 2, gene, gene_list)
        size<-elem_size(size, k, 6, meta, meta_list)
    }
    return(size)
}

#add direct interactions types to direct interactions dataframe
interactions_type<-function(inter, meta, genes){
    inter<-cbind(inter, type=rep(NA, nrow(inter)))
    for (d in seq_len(nrow(inter))){
        tag<-""
        if(inter[d, 1] %in% meta && inter[d, 3] %in% meta){
            tag<-"m/m"
        }
        else if((inter[d, 1] %in% meta && inter[d, 3] %in% genes) ||
                (inter[d, 1] %in% genes && inter[d, 3] %in% meta)){
            tag<-"g/m"
        }
        else if (inter[d, 1] %in% genes && inter[d, 3] %in% genes){
            tag<-"g/g"
        }
        inter[d, 5]<-tag
    }
    return(inter)
}

##for every pathway, determine how many direct interactions an element from
##user's list have with others elements in the pathway.
centrality<-function(interac, list_elem){
    central<-list()
    inter_allpath<-rm_vector(interac[, 4])
    for(p in seq_len(length(inter_allpath))){
        path_inter<-interac[interac[, 4] == inter_allpath[p], ]
        elements<-rm_vector(c(path_inter[path_inter[, 1] %in% list_elem, 1],
                            path_inter[path_inter[, 3] %in% list_elem, 3]))
        if (length(elements)>0){
            for ( e in seq_len(length(elements))){
                path_elem<-rbind(path_inter[path_inter[, 1] == elements[e], ],
                            path_inter[path_inter[, 3]==elements[e], ])
                central[[inter_allpath[p]]][[elements[e]]]<-nrow(path_elem)
            }
        }
    }
    return(central)
}

##regroup interactions without taking pathways into account
filter_inter<-function(inter){
    tagged<-cbind(inter[, c(seq_len(4))], tag=rep(NA, nrow(inter)), inter[, 5])
    for (d in seq_len(nrow(inter))){
        tagged[d, 5]<-paste(tagged[d, 1]," / ", tagged[d, 3], sep="")
    }
    tagged<-rm_df(tagged, c(seq_len(5)))

    elements<-tagged[,c(1, 3)]
    filtered<-elements[!duplicated(lapply(as.data.frame(t(elements)), sort)),]
    no_path<-data.frame(from=character(), link=character(), to=character(),
                        path=character(), tag=character(), type=character())
    for(t in seq_len(nrow(filtered))){
        filtered_rows<-tagged[intersect(
            rm_vector(c(which(tagged[, 1] %in% filtered[t, 1]),
                        which(tagged[, 3] %in% filtered[t, 1]))),
            rm_vector(c(which(tagged[, 1] %in% filtered[t, 2]),
                        which(tagged[, 3] %in% filtered[t, 2])))), ]
        no_path[t, 1]<-filtered[t, 1]; no_path[t, 3]<-filtered[t, 2]
        no_path[t, 2]<-paste(rm_vector(filtered_rows[, 2]), collapse=", ")
        no_path[t, 4]<-paste(rm_vector(filtered_rows[, 4]), collapse=", ")
        no_path[t, 5]<-paste(c(filtered[t, 1], " / ",
                                filtered[t, 2]), collapse="")
        no_path[t, 6]<-filtered_rows[1, 6]
    }
    no_path<-rm_df(no_path, c(seq_len(4)))
    return(list(tagged, no_path))
}

##Input : path_enrich results with KEGG, Reactome and Wikipathways.
interactions<-function(listk, listr, listw){
    resmetak<-listk[[1]]; resgenek<-listk[[2]]; resmetar<-listr[[1]]
    resgener<-listr[[2]]; resmetaw<-listw[[1]]; resgenew<-listw[[2]]
    genes<-listk[[3]]; meta<-listk[[4]]
    pathtot<-get_pathways(resgenek, resgener, resgenew,
                            resmetak, resmetar, resmetaw)
    pathtotr<-pathtot[which(stringr::str_sub(pathtot$id, 1, 1)=="R"), ]
    pathtotp<-pathtot[which(stringr::str_sub(pathtot$id, 1, 1)=="W"), ]
    pathtotk<-pathtot[which(stringr::str_sub(pathtot$id, 1, 1)=="h"), ]
    ##mappings for metabolites(chebi,kegg,name) and genes(gene symbol,name)
    keggchebi<-KEGGREST::keggConv("chebi", "compound")
    keggchebi<-rm_df(data.frame(kegg=names(keggchebi), chebi=keggchebi),
                    c(seq_len(2)))
    keggname<-KEGGREST::keggList("compound")
    keggname<-data.frame(kegg=names(keggname), name=keggname)
    for (t in seq_len(nrow(keggname))){
        keggname[t, 2]<-stringr::str_split(keggname[t, 2], ";")[[1]][1]
    }
    keggchebiname<-merge(keggchebi, keggname, by="kegg")
    keggchebiname$chebi<-stringr::str_to_upper(keggchebiname$chebi)
    keggchebiname<-keggchebiname[which(!is.na(keggchebiname[,1])), ]
    keggchebiname<-keggchebiname[which(!is.na(keggchebiname[,2])), ]
    keggchebiname<-keggchebiname[which(!is.na(keggchebiname[,3])), ]
    ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
    namegeneid<-biomaRt::getBM(
        attributes=c("hgnc_symbol", "entrezgene_description"), mart=ensembl)
    interac<-interac[which(interac[, 4] %in% pathtot[, 2]), ]
    ##Data informations
    meta<-rm_vector(keggname[which(keggname$kegg
                                    %in% paste("cpd:", meta, sep="")), 2])
    size<-data_size(interac, meta, genes, keggchebiname,
                    pathtotk, pathtotr, pathtotp)
    list_elem<-c(meta, genes)
    interac<-rm_df(rbind(interac[which(interac[, 1] %in% list_elem), ],
                        interac[which(interac[, 3] %in% list_elem), ]))
    central<-centrality(interac, list_elem)
    interac<-interac[intersect(which(interac[, 1] %in% list_elem),
                                which(interac[, 3] %in% list_elem)), ]
    interac<-interactions_type(interac, meta, genes)
    list_filter<-filter_inter(interac)
    tagged<-list_filter[[1]]; no_path<-list_filter[[2]]

    return(list(size, pathtot, tagged, namegeneid, keggchebiname, central,
                                                no_path, genes, meta, ensembl))
}
