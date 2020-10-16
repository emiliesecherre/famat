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
elem_size<-function(user_elem, path_elem){
    sumary=vector()
    if (length(user_elem)>0){
        sumary=c(sumary,length(intersect(user_elem, path_elem)))
        if (length(intersect(user_elem, path_elem))>0){
            elements=intersect(user_elem, path_elem)
            sumary=c(sumary,paste(elements, collapse=" # "))
        }
        else{sumary=c(sumary,NA)}
        sumary=c(sumary,length(user_elem))
        sumary=c(sumary,paste(user_elem, collapse=" # "))
    }else{
        sumary=c(sumary,0,NA,length(user_elem),NA)
    }
    return(sumary)
}

#dataframe showing pathways and which genes/metabolites they have
data_size<-function(inter, meta_list, gene_list, mapkegg, pathk, pathr, pathp){
    pathtot_ids<-c(pathk$id, pathr$id, pathp$id)
    pre_size=lapply(pathtot_ids, function(x){
        path_inter<-inter[which(inter$path == x),]
        elements<-rm_vector(c(path_inter[, 1], path_inter[, 3]))
        meta<-elements[which(elements %in% mapkegg$name)]
        gene<-elements[which(!(elements %in% mapkegg$name))]

        c(x,elem_size(gene, gene_list),elem_size(meta, meta_list))
    })
    size=data.frame(matrix(unlist(pre_size), nrow=length(pre_size), byrow=TRUE))
    names(size)=c("path", "nb_gene_query","gene_que", "nb_gene_tot","genes",
                    "nb_meta_query","meta_que", "nb_meta_tot","meta")
    size[,c(2,4,6,8)]=apply(size[,c(2,4,6,8)],2,as.integer)
    return(size)
}

#add direct interactions types to direct interactions dataframe
interactions_type<-function(inter, meta, genes){
    tags=apply(inter,1,function(x){
        tag<-""
        if(x[1] %in% meta && x[3] %in% meta){
            tag<-"m/m"
        }
        else if((x[1] %in% meta && x[3] %in% genes) ||
                (x[1] %in% genes && x[3] %in% meta)){
            tag<-"g/m"
        }
        else if (x[1] %in% genes && x[3] %in% genes){
            tag<-"g/g"
        }
        tag
    })
    inter<-cbind(inter, type=unname(tags))
    return(inter)
}

##for every pathway, determine how many direct interactions an element from
##user's list have with others elements in the pathway.
centrality_calc<-function(interac, list_elem){
    inter_allpath<-rm_vector(interac[, 4])
    central=vapply(inter_allpath,function(x){
        path_inter<-interac[interac[, 4] == x, ]
        elements<-rm_vector(c(path_inter[path_inter[, 1] %in% list_elem, 1],
                                path_inter[path_inter[, 3] %in% list_elem, 3]))
        if (length(elements)>0){
            pre_central=vapply(elements,function(y){
                path_elem<-rbind(path_inter[path_inter[, 1] == y, ],
                                    path_inter[path_inter[, 3]==y, ])
                list(nrow(path_elem))
            }, list(1))
            list(pre_central)
        }
    }, list(1))
    return(central)
}

##regroup interactions without taking pathways into account
filter_inter<-function(inter){
    tags=apply(inter,1,function(x){
        paste(x[1]," / ", x[3], sep="")
    })
    tagged<-cbind(inter[, c(seq_len(4))], tag=unname(tags), inter[, 5])
    tagged<-rm_df(tagged, c(seq_len(5)))

    elements<-tagged[,c(1, 3)]
    filtered<-elements[!duplicated(lapply(as.data.frame(t(elements)), sort)),]
    no_path=apply(filtered,1,function(x){
        filtered_rows<-tagged[intersect(
            rm_vector(c(which(tagged[, 1] %in% x[1]),
                        which(tagged[, 3] %in% x[1]))),
            rm_vector(c(which(tagged[, 1] %in% x[2]),
                        which(tagged[, 3] %in% x[2])))), ]
        c(x[1],paste(rm_vector(filtered_rows[, 2]), collapse=", "),x[2],
            paste(rm_vector(filtered_rows[, 4]), collapse=", "),
            paste(c(x[1], " / ", x[2]), collapse=""),filtered_rows[1, 6])
    })
    no_path=as.data.frame(t(no_path))
    names(no_path)=names(tagged)=c("from", "link", "to", "path", "tag", "type")
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
    name=apply(keggname,1,function(x){stringr::str_split(x[2], ";")[[1]][1]})
    keggname[,2]=unname(name)
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
    central<-centrality_calc(interac, list_elem)
    interac<-interac[intersect(which(interac[, 1] %in% list_elem),
                                which(interac[, 3] %in% list_elem)), ]
    interac<-interactions_type(interac, meta, genes)
    list_filter<-filter_inter(interac)
    tagged<-list_filter[[1]]; no_path<-list_filter[[2]]

    return(list(size, pathtot, tagged, namegeneid, keggchebiname, central,
                                                no_path, genes, meta, ensembl))
}
