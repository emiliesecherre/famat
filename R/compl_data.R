#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr desc
#' @export compl_data
###############################FUNCTIONS##########################
##build walks using 2-2 relations between items
##for every parent-child relation, if a parent of the parent term is found,
##put this parent before the parent term. Same is done for children.
build_walk<-function(df){
    walks<-apply(df, 1, paste, collapse=">")
    walks<-as.data.frame(unname(walks)); names(walks)<-"hier"
    save_walks<-walks[-seq_len(nrow(walks)),]; stop<-0
    while(stop == 0){
        save_walks_child<-apply(walks,1,function(s){
            walk<-as.data.frame(unname(s))
            #parents
            names(walk)=names(walks)
            parents<-df[which(df[,2] %in%
                                stringr::str_split(walk[1,1],">")[[1]][1]),1]
            if(length(parents)>0){
                pre_save_walks<-vapply(parents, function(c){
                    paste(c,walk[1,1],sep=">")
                }, character(1))
                unname(pre_save_walks)
            }
            else{walk}
        })
        save_walks_parent<-apply(walks,1,function(s){
            walk<-as.data.frame(unname(s))
            #child
            child<-stringr::str_split(walk[1,1],">")[[1]]
            children<-df[which(df[,1] %in% child[length(child)]),2]
            if(length(children)>0){
                pre_save_walks<-vapply(children, function(c){
                    paste(walk[1,1],c,sep=">")
                }, character(1))
                unname(pre_save_walks)
            }
            else{walk}
        })
        save_walks<-c(unname(unlist(save_walks_parent)),
                        unname(unlist(save_walks_child)))
        save_walks<-as.data.frame(rm_vector(save_walks))
        names(save_walks)="hier"
        if(length(save_walks[!(save_walks[,1] %in% walks[,1]),1]) == 0){
            stop<-1
        }
        walks<-save_walks
        save_walks<-walks[-seq_len(nrow(walks)),]
    }
    walks<-apply(walks, 1, function(w){
        if(length(walks[stringr::str_detect(walks[,1], w),1])<=1){w}
    })
    walks<-as.data.frame(unname(unlist(walks)))
    names(walks)<-c("hier")
    return(walks)
}

#gather go terms hierarchies, their index, go terms and genes
list_go<-function(hgo){
    root_index<-which(!(stringr::str_sub(hgo[, 1], 1, 1) == "_"))

    golist<-vapply(root_index, function(r){
        if (r!=root_index[length(root_index)]){
            hiera_index<-c(r:(root_index[which(root_index %in% r)+1]-1))
        }
        else{
            hiera_index<-c(r:nrow(hgo))
        }
        genes_hiera<-vapply(hiera_index, function(h){
            list(stringr::str_split(hgo[h, 3], ", ")[[1]])
        },list(1))
        genes_hiera<-unname(unlist(genes_hiera))
        go_terms<-vapply(hiera_index, function(h){
            gomed<-stringr::str_split(hgo[h, 1], "__")[[1]]
            gomed<-gomed[gomed!=""]
            list(stringr::str_split(gomed, " / ")[[1]])
        },list(1))
        go_terms<-unname(unlist(go_terms))
        genes_hiera<-genes_hiera[!is.na(genes_hiera)]
        genes_hiera<-rm_vector(genes_hiera[genes_hiera!=""])
        go_terms<-go_terms[!is.na(go_terms)]
        go_terms<-rm_vector(go_terms[go_terms!=""])
        list(list(goterm=go_terms, index=hiera_index, gene=genes_hiera))
    }, list(1))
    return(golist)
}

##for every walk ending by a same leaf, if one or more walk contain a node
##required in hierarchies but not present in leaves, remove other walks
sort_by_leaves<-function(df, notin_leaf){
    df[]<-lapply(df, as.character)
    pre_walk<-build_walk(df)
    pre_walk<-rm_vector(pre_walk$hier)
    walk<-vapply(pre_walk, function(p){
        leaf<-stringr::str_sub(p, nchar(p)-9, nchar(p))
        leaf_walks<-pre_walk[stringr::str_sub(pre_walk, nchar(pre_walk)-9,
                                                nchar(pre_walk)) == leaf]
        if(length(leaf_walks)>1){
            match<-vapply(notin_leaf, function(n){
                list(leaf_walks[stringr::str_detect(leaf_walks,n)])
            }, list(1))
            match<-unname(unlist(match))
            if(length(match) == 0){#keep the smaller walk
                match<-leaf_walks[stringr::str_count(leaf_walks, ">") == min(
                    stringr::str_count(leaf_walks, ">"))][1]
            }
        }
        else{match<-leaf_walks}
        list(match)
    }, list(1))
    walk<-unname(unlist(walk))
    walk<-rm_vector(walk);walk<-walk[!is.na(walk)]

    max_len<-max(stringr::str_count(walk,">"))+1
    walk<-as.data.frame(walk)
    walk<-tidyr::separate(walk, 1, as.character(c(seq_len(max_len))),
                                sep=">", extra="drop", fill="right")
    return(walk)
}

##if there is only one different node between two walks, merge walks
merge_walks<-function(walk){
    merged<-walk[-seq_len(nrow(walk)),]; b<-1
    while (nrow(walk)>=b){
        go_terms<-walk[b, ];go_terms<-go_terms[!is.na(go_terms)]
        leaf<-go_terms[length(go_terms)]; index<-which(go_terms == leaf)
        studied<-remove<-walk[which(walk[, index] == leaf), ]
        if(nrow(studied)>1){
            stop<-0 ; bye<-0
            while (bye == 0 && nrow(studied)>1){
                leaf_len<-nrow(studied); a<-1
                while(a<=nrow(studied) && nrow(studied)>1){#first walk
                    i<-1
                    while(i<=nrow(studied) && nrow(studied)>1){#second
                        stop<-1
                        dif<-vapply(seq_len(length(studied[a, ])), function(v){
                            if(!is.na(studied[a, v]) &&
                                !is.na(studied[i, v])){
                                term_one<-stringr::str_split(
                                    studied[a, v]," / ")[[1]]
                                term_two<-stringr::str_split(
                                    studied[i, v]," / ")[[1]]
                                #difference between terms
                                if(length(intersect(term_one, term_two)) == 0){
                                    list(v)
                                }
                                else{list(NULL)}
                            }
                            else{list(NULL)}
                        }, list(1))
                        dif<-unname(unlist(dif))
                        if(length(dif) == 1){
                            studied[a, dif]<-paste(studied[a, dif]," / ",
                                                    studied[i, dif], sep="")
                            studied<-studied[-i, ]; stop<-0}
                        else if (length(dif) == 0 && a!=i){
                            studied<-studied[-i, ];stop<-0
                        }
                        if(stop == 1){i<-i+1}
                        else{i<-1}
                    }
                    if(stop == 1){a<-a+1}
                }
                if(nrow(studied) == leaf_len){bye<-1}
            }}
        merged<-rbind(merged, studied);walk<-dplyr::setdiff(walk, remove)}
    return(merged)
}

##sort GO terms walks from the hierarchy
sort_go_steps<-function(links, res_enrich){
    leaves<-as.vector(unique(links[which(!(links$to %in% links$from)), 2]))
    notin_leaf<-dplyr::setdiff(res_enrich[, 1], leaves)

    walk<-sort_by_leaves(links, notin_leaf)
    #walks included into others ?
    toremove<-apply(walk,1,function(w){
        walk_one<-unname(w)
        walk_one<-walk_one[!is.na(walk_one)]
        comp<-apply(walk,1,function(x){
            walk_two<-unname(x); walk_two<-walk_two[!is.na(walk_two)]
            common<-intersect(walk_one, walk_two)
            if (length(common) == min(length(walk_one), length(walk_two)) &&
                ((length(common) != length(walk_one)) ||
                (length(common) != length(walk_two)))){
                if(length(walk_one)>length(walk_two)){
                    nas<-ncol(walk)-length(walk_two)
                    c(walk_two,rep(NA,nas))
                }
                else if(length(walk_two)>length(walk_one)){
                    nas<-ncol(walk)-length(walk_one)
                    c(walk_one,rep(NA,nas))
                }
            }
        })
        if(!is.null(comp)){
            unname(unlist(comp))
        }
    })
    if(!is.null(toremove)){
        toremove<-as.data.frame(matrix(unname(unlist(toremove)),
                                    ncol=ncol(walk), byrow=TRUE))
        toremove<-rm_df(toremove, seq_len(ncol(toremove)))
        names(toremove)=names(walk)
        walk<-dplyr::setdiff(walk,toremove)
    }
    walk<-merge_walks(walk)
    return(walk)
}

##build a hierarchy in tree view form
##for ordered walks beginning with the same root, this function goes from the
##last walk to the second. For a node, if in the walk above the walk studied
##there is the same node at the same position in the hierarchy, remove this
##node. In the dataframe, the position is the column. So, nodes in the df
##will look like a tree.
tree_view<-function(walks){
    row.names(walks)=as.character(seq_len(nrow(walks)))
    num<-as.integer(row.names(walks))
    arr_walks<-dplyr::arrange(walks,desc(num))
    ref<-arr_walks[2:nrow(arr_walks),]
    pre_walks<-vapply(sort(as.integer(row.names(ref))),function(r){
        temp<-arr_walks[r,]
        refe<-ref[r,]
        replace<-vapply(seq_len(ncol(refe)),function(y){
            if(!is.na(temp[y]) && !is.na(refe[y])
                && temp[y] == refe[y]){
                list(NA)
            }
            else{list(temp[y])}
        }, list(1))
        unname(unlist(replace))
    }, character(ncol(ref)))
    pre_walks<-as.data.frame(t(pre_walks))
    names(pre_walks)=names(arr_walks)
    pre_walks<-rbind(pre_walks,arr_walks[nrow(arr_walks),])
    pre_walks<-dplyr::arrange(pre_walks,desc(num))
    row.names(pre_walks)=seq_len(nrow(pre_walks))

    #one node per line
    treeview<-apply(pre_walks,1,function(x){
        pre_tree<-vapply(seq_len(length(x)), function(y){
            if(!is.na(x[y])){list(c(rep(NA,y-1),x[y],rep(NA,length(x)-y)))}
            else{list(NULL)}
        }, list(1))
        unlist(pre_tree)
    })
    treeview<-unname(unlist(treeview))
    treeview<-as.data.frame(matrix(treeview, ncol=ncol(walks), byrow=TRUE))
    return(treeview)
}

#add Go terms names and genes from user's list to GO hierarchies
info_go<-function(go_tab, go_genelist, res_enrich){
    map_names<-as.list(GO.db::GOTERM)

    go_tab<-apply(go_tab, 1, function(g){
        go_terms<-stringr::str_split(g[1], " / ")[[1]]
        go_terms<-c(stringr::str_split(go_terms[1], "__")[[1]],
                    go_terms[2:length(go_terms)])
        go_terms<-go_terms[go_terms!=""]
        go_terms<-rm_vector(go_terms[!is.na(go_terms)])

        go_names<-vapply(go_terms, function(t){
            if (is.null(map_names[[t]])){
                list(res_enrich[res_enrich[, 1] == t, 2])
            }
            else{list(map_names[[t]]@Term)}
        }, list(1))
        go_names<-unname(unlist(go_names))

        go_genes<-vapply(go_terms, function(t){
            pre_go_genes<-vector()
            pre_go_genes<-go_genelist[go_genelist$go_id == t, 1]
            pre_go_genes<-rm_vector(pre_go_genes[!is.na(pre_go_genes)])
            list(paste(pre_go_genes, collapse=", "))
        }, list(1))
        go_genes<-unname(unlist(go_genes))

        go_names<-go_names[!is.na(go_names)]
        go_genes<-go_genes[go_genes!=""]

        c(g[1],paste(go_names,collapse=" / " ), paste(go_genes,collapse=" / "))
    })
    go_tab<-as.data.frame(matrix(unlist(unname(go_tab)), ncol=3, byrow=TRUE))
    return(go_tab)
}

#with enriched GO terms, build GO terms hierarchies
hieraGO<-function(type, res_enrich, go_genelist){
    #parent terms of enriched GO terms
    if(type == "MF"){ancestors <- as.list(GO.db::GOMFPARENTS)}
    else if (type == "BP"){ancestors <- as.list(GO.db::GOBPPARENTS)}
    ancestors <- ancestors[!is.na(ancestors)]
    ancestors<-ancestors[!(ancestors %in% "all")]
    if(type == "MF"){ancestors[["GO:0003674"]]<-vector()}
    else if (type == "BP"){ancestors[["GO:0008150"]]<-vector()}
    onto<-ontologyIndex::ontology_index(parents=ancestors)
    parents<-ontologyIndex::get_ancestors(onto, res_enrich[, 1])
    candidate<-c(res_enrich[, 1], parents) #Go terms + ancestors
    children<-names(ancestors)
    go_links<-vapply(children, function(child){
        parent<-unname(ancestors[[child]])
        toadd<-vapply(parent, function(s){
            if (s %in% candidate && child %in% candidate){list(c(s,child))}
            else{list(NULL)}
        }, list(1))
        list(toadd)
    }, list(1))
    go_links<-as.data.frame(matrix(unname(unlist(go_links)),ncol=2,byrow=TRUE))
    names(go_links)=c("from","to")
    walks<-sort_go_steps(go_links, res_enrich)
    #build walks
    pre_walks<-apply(walks, 1, function(wk){
        wk<-wk[!is.na(wk)]
        paste(wk, collapse=">")
    })
    pre_walks<-unname(pre_walks)
    pre_walks<-sort(pre_walks) #useful to build tree view
    pre_walks<-as.data.frame(pre_walks)
    walks<-tidyr::separate(pre_walks, 1, as.character(c(seq_len(ncol(walks)))),
                            sep=">", extra="drop", fill="right")
    walks<-walks[,-1] #remove root term of BP/MF hierarchies
    treeview<-tree_view(walks)

    #reduce treeview to a single column, with "__" instead of NAs
    names(treeview)<-as.character(c(seq_len(ncol(treeview))))

    go_tab<-apply(treeview, 1, function(t){
        index<-which(!is.na(t))
        space<-paste(rep("__", index-1), collapse="")
        node<-t[index]
        paste(space, node, sep="")
    })
    go_tab<-as.data.frame(go_tab)
    go_tab<-info_go(go_tab, go_genelist, res_enrich)
    names(go_tab)<-c("goterm", "go_name", "genes")
    return(go_tab)
}

##find which elements are found in the same pathways, and put them together
##find which pathways contain the same elements also
##if element=T, also return user's elements which aren't in pathways
cluster_items<-function(items, element){
    items[is.na(items)]<-"N";items<-items[-1, ]
    elem_names<-colnames(items)

    dist_data<-vapply(elem_names, function(e){#distances between elements
        index<-which(elem_names %in% e)
        first_item<-items[, index]
        if(index != length(elem_names)){
            sub_elem<-elem_names[(index+1):length(elem_names)]
            diff<-vapply(sub_elem, function(s){
                sec_index<-which(elem_names %in% s)
                sec_item<-items[, sec_index]
                dist<-length(which(first_item != sec_item))
                list(c(e,s,dist))
            }, list(1))
            list(diff)
        }
        else{list(NULL)}
    }, list(1))
    dist_data<-unname(unlist(dist_data))
    dist_data<-as.data.frame(matrix(dist_data, ncol=3, byrow=TRUE))
    names(dist_data)=c("first","second","dist")
    notin_path<-vapply(elem_names, function(e){
        index<-which(elem_names %in% e)
        first_item<-items[, index]
        nb_path<-length(first_item[first_item %in% "X"])
        if(element == TRUE && nb_path == 0){list(e)}
        else{list(NULL)}
    }, list(1))
    notin_path<-unname(unlist(notin_path))

    ##order elements by distances
    dist_data<-dist_data[order(dist_data$dist),]
    cluster<-apply(dist_data, 1, function(d){
        c(d[1], d[2])
    })
    cluster<-rm_vector(as.character(cluster))
    if (element == TRUE){return(list(cluster, notin_path))}
    else{return(cluster)}
}

#filter entire pathways hierarchy to build a hierarchy concerning our pathways
sort_hiera<-function(pathways){
    kegg_path<-pathways[stringr::str_sub(pathways, 1, 3) == "hsa"]
    kegg_path<-paste(stringr::str_sub(kegg_path,1,3),
                        stringr::str_sub(kegg_path,5),sep="")
    path_walks_k<-vapply(kegg_path, function(x){
        list(kegg_hiera[stringr::str_detect(kegg_hiera[,1],x),1])
    }, list(1))
    path_walks_k<-as.data.frame(sort(unlist(path_walks_k)))
    if(ncol(path_walks_k) == 0){path_walks_k<-data.frame(walks=character())}

    wp_path<-pathways[stringr::str_sub(pathways, 1, 2) == "WP"]
    path_walks_w<-vapply(wp_path, function(x){
        list(wp_hiera[stringr::str_detect(wp_hiera[,1],x),1])
    }, list(1))
    path_walks_w<-as.data.frame(sort(unlist(path_walks_w)))
    if(ncol(path_walks_w) == 0){path_walks_w<-data.frame(walks=character())}

    first_walks_r<-rea_hiera[2:nrow(rea_hiera),1]
    path_walks_r<-vapply(first_walks_r, function(x){
        rea_walks<-stringr::str_split(x, ">")[[1]]
        if(length(pathways[pathways %in% rea_walks])>0){
            rea_walks<-rm_vector(rea_walks[c(1, which(rea_walks%in%pathways))])
            list(paste(rea_walks, collapse=">"))
        }
        else{list(NULL)}
    }, list(1))
    path_walks_r<-rm_vector(unname(unlist(path_walks_r)))
    path_walks_r<-path_walks_r[stringr::str_detect(path_walks_r, ">")]

    final_walks_r<-vapply(path_walks_r, function(x){
        dupl<-which(stringr::str_detect(path_walks_r, x))
        dupl<-dupl[-which(dupl == which(path_walks_r == x))]
        if(length(dupl) == 0){list(x)}
        else{list(NULL)}
    }, list(1))
    final_walks_r<-rm_vector(unname(unlist(final_walks_r)))
    final_walks_r<-as.data.frame(sort(final_walks_r))

    names(final_walks_r)<-names(path_walks_w)<-names(path_walks_k)<-"walks"
    path_walks<-rbind(final_walks_r, path_walks_k,path_walks_w)
    max<-max(stringr::str_count(path_walks[,1],">"))+1
    return(list(path_walks, max))
}

#add informations about pathway hierarchies to the final heatmap
#informations are : ratio of user's elements / total elements,
#names and ids of pathways in hierarchies
hiera_info<-function(pathways, size, sorted_pathways,
                                        treeview, no_path, list_elem){
    pathidtoname <- as.list(reactome.db::reactomePATHID2NAME)
    heatmap<-apply(treeview, 1, function(t){
        index<-which(!is.na(t))
        space<-paste(rep("__", index-1), collapse="")
        node<-t[index]
        if(stringr::str_sub(node, 1, 3) == "hsa"){
            node<-paste("hsa:", stringr::str_sub(node, 4, nchar(node)), sep="")
        }
        name<-pathways[pathways[,2] == node, 1]
        htmp<-vector()
        if(length(name)>0){
            htmp<-c(htmp,paste(space,
                            pathways[pathways[,2] == node, 1], sep=""))
        }
        else if(stringr::str_sub(node, 1, 5) == "R-HSA"){
            htmp<-c(htmp,paste(space, stringr::str_sub(pathidtoname[[node]],
                                    15, nchar(pathidtoname[[node]])), sep=""))
        }
        else{htmp<-c(htmp,paste(space, node, sep=""))}
        htmp<-c(htmp,node,
                paste("'", size[size$path == node, 2], "/",
                        size[size$path == node, 4], sep=""),
                paste("'", size[size$path == node, 6], "/",
                        size[size$path == node, 8], sep=""),NA)
        htmp
    })
    heatmap<-as.data.frame(t(heatmap))
    heatmap<-rbind(rep(NA, ncol(heatmap)), heatmap)
    colnames(heatmap)<-c("path_name", "path_id", "meta_ratio", "gene_ratio",
                            "blank")
    heatmap[which(heatmap[, 3] == "'/"), 3]<-"'0/0"
    heatmap[which(heatmap[, 4] == "'/"), 4]<-"'0/0";tags<-no_path$tag
    elements<-c(list_elem, tags); elements_df<-t(as.data.frame(elements))
    elements_df<-as.data.frame(elements_df); colnames(elements_df)<-elements
    toadd<-as.data.frame(matrix(seq_len(ncol(elements_df)),
                    ncol=ncol(elements_df), nrow=nrow(heatmap)-1,byrow=TRUE))
    names(toadd)=names(elements_df)
    elements_df<-rbind(elements_df, toadd)
    heatmap<-cbind(heatmap, elements_df)
    return(heatmap)
}

##cluster hierarchies using their elements
##to do so, all roots of hierarchies are used, and elements from all
##the hierarchy pathways are added to the root
cluster_hiera<-function(heatmap, size, tagged, no_path){
    help_hm<-heatmap<-heatmap[-1, ]; root<-help_hm[-(seq_len(nrow(help_hm))), ]
    index<-c(which(stringr::str_sub(help_hm[,1], 1, 1)!="_"), nrow(help_hm)+1)
    root<-vapply(index, function(i){
        if(!(i == index[length(index)])){
            hiera<-help_hm[c(i:(index[which(index == i)+1]-1)), ]
            hiera<-rbind(rep(NA, ncol(hiera)), hiera); hiera[1, 1]<-hiera[2, 2]
            hiera<-apply(hiera, 2, function(h){
                if(is.na(h[1]) ){
                    if("X" %in% h){"X"}
                    else{NA}
                }
                else{h[1]}
            })
            list(unname(hiera))
        }
        else{list(NULL)}
    }, list(1))
    root[length(root)]=NULL
    root<-data.frame(matrix(unlist(root), nrow=length(root), byrow=TRUE))
    names(root)=names(help_hm); root_data<-root[seq_len(nrow(root)), ]
    root_data<-t(root_data);colnames(root_data)=root_data[1,]
    clusters<-cluster_items(root_data, FALSE)
    heatmap<-help_hm[-(seq_len(nrow(heatmap))), ]
    heatmap<-vapply(clusters, function(c){
        root_id<-which(help_hm[, 2] == c)
        next_root<-index[which(index == root_id)+1]
        if(!is.na(next_root)){
            row<-rbind(help_hm[root_id:(next_root-1),],rep(NA, ncol(heatmap)))
        }
        list(t(row))
    }, list(1))
    heatmap<-as.data.frame(matrix(unname(unlist(heatmap)),
                                    ncol=ncol(help_hm), byrow=TRUE))
    names(heatmap)=names(help_hm)
    root_ids<-c(which(!stringr::str_detect(heatmap[,1],"__")), nrow(heatmap)+1)
    hierapath<-vapply(root_ids, function(r){
        next_root<-root_ids[which(root_ids == r)+1]
        if(!is.na(next_root)){
            element<-root[which(root[,1] %in% heatmap[r,2]),]
            element<-t(element)[,1]; element<-element[-1]
            element<-names(element[!is.na(element)])
            list(list(name=heatmap[r:(next_root-2),2],
                        index=c(r:(next_root-2)),elem=element))
        }
        else{list(NULL)}
    }, list(1))
    hierapath[length(hierapath)]=NULL
    return(list(heatmap, hierapath))
}

##cluster hierarchies and elements
cluster_htmp<-function(heatmap, tags, size, tagged, no_path){
    elem_data<-heatmap[,c(6:(ncol(heatmap)-length(tags)))]
    elem_data<-as.matrix(elem_data);tags<-sort(tags)
    listele<-cluster_items(elem_data, TRUE)
    cluster_elem<-save_cluster_elem<-listele[[1]];notin_path<-listele[[2]]
    cluster_elem<-c(cluster_elem, tags) #add interactions
    cluster_elem<-cluster_elem[!(cluster_elem %in% notin_path)]
    heatmap<-heatmap[,c("path_name", "path_id", "meta_ratio", "gene_ratio",
                        "blank", cluster_elem)]
    listhiera<-cluster_hiera(heatmap, size, tagged, no_path)
    heatmap<-listhiera[[1]]; hierapath<-listhiera[[2]]
    row.names(heatmap)<-c(seq_len(nrow(heatmap)))
    hierapath<-lapply(hierapath,function(x){
        index<-x[["index"]]
        list(name=x[["name"]],index=c(index,index[length(index)]+1),
                elem=x[["elem"]])
    })
    heatmap<-apply(heatmap,1,function(x){#pathway with direct interaction ?
        if(x[2] %in% tagged$path){
            c(x[1]<-paste(x[1],"(*)", sep=" "),x[2:length(x)])
        }
        else{x}
    })
    heatmap<-as.data.frame(t(heatmap))
    names(heatmap)<-c("path_name", "path_id", "meta_ratio", "gene_ratio",
                        "blank", cluster_elem)
    return(list(heatmap, notin_path, hierapath, save_cluster_elem))
}

##build heatmap of hierarchies of pathways and elements included in them
final_tab<-function(build_hm, pathways, size, sorted_path, no_path,
                    list_elem, tagged){
    heatmap<-hiera_info(pathways, size, sorted_path, build_hm,
                        no_path, list_elem)
    sub_htmp<-heatmap[2:nrow(heatmap),]; tags<-no_path$tag #direct interactions
    pre_htmp<-apply(sub_htmp,1,function(h){
        path<-h[2]
        pre_elem<-c(size[size$path %in% path, 3], size[size$path %in% path, 7])
        pre_elem<-pre_elem[!is.na(pre_elem)]
        if(length(pre_elem)>0){
            element<-vapply(pre_elem, function(e){
                list(stringr::str_split(e," # ")[[1]])}, list(1))
            element<-rm_vector(unname(unlist(as.list(element))))
            htmp_elem<-h[6:(length(h)-length(tags))]
            first_vec<-vapply(htmp_elem,function(t){
                if(names(htmp_elem[which(htmp_elem %in% t)]) %in% element){
                    list("X")
                }
                else{list(NA)}
            }, list(1))
            first_vec<-unname(unlist(first_vec))
            htmp_tags<-h[(length(h)-length(tags)+1):length(h)]
            sec_vec<-vapply(htmp_tags, function(m){
                path_inter<-tagged[tagged$path == path,]
                path_inter<-path_inter[path_inter$tag ==
                                    names(htmp_tags[which(htmp_tags %in% m)]),]
                if(nrow(path_inter)>0){list("X")}
                else{list(NA)}
            }, list(1))
            sec_vec<-unname(unlist(sec_vec))
            htmp_row<-c(h[seq_len(5)],first_vec,sec_vec)
        }
        else{htmp_row<-c(h[seq_len(5)],rep(NA, ncol(sub_htmp)-5))}
    })
    pre_htmp<-as.data.frame(t(pre_htmp));names(pre_htmp)=names(heatmap)
    heatmap<-rbind(heatmap[1,],pre_htmp)

    ##CLUSTERS
    listhtmp<-cluster_htmp(heatmap, tags, size, tagged, no_path)
    heatmap<-listhtmp[[1]]; notin_path<-listhtmp[[2]]; hierapath<-listhtmp[[3]]
    save_cluster_elem<-listhtmp[[4]]
    return(list(heatmap, notin_path, hierapath, save_cluster_elem))
}

##perform go term enrichment analysis
enr_go<-function(genes){
    ##entrez genes ids for user's genes
    genes_entrez<-gprofiler2::gconvert(genes, organism="hsapiens",
                                        target='ENTREZGENE_ACC')$target
    allResBP<-clusterProfiler::enrichGO(genes_entrez, keyType="ENTREZID",
                            'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)@result
    allResBP<-allResBP[c(seq_len(20)), ]
    allResMF<-clusterProfiler::enrichGO(genes_entrez, keyType="ENTREZID",
                            'org.Hs.eg.db', ont="MF", pvalueCutoff=0.01)@result
    allResMF<-allResMF[c(seq_len(20)), ]
    topgo<-rbind(allResBP[, c(1, 2)], allResMF[, c(1, 2)])

    go_genelist<-go_gene[go_gene$hgnc_symbol %in% genes, ]
    go_genelist<-go_genelist[go_genelist$go_id %in% topgo[,1], ]
    return(list(go_genelist, allResBP, allResMF))
}

#build dataframes showing informations about elements and their interactions
infos_elem<-function(genes, notin_path, meta, keggchebiname, no_path,
                        go_genelist){
    #genes informations
    pre_genetab<-vapply(genes,function(g){
        gene_name<-namegeneid[namegeneid$hgnc_symbol == g, 2]
        if(length(gene_name)>0){c(g, gene_name)}
        else{c(g, "" )}
    }, character(2))
    pre_genetab<-as.data.frame(t(pre_genetab))
    genetab<-pre_genetab[which(!(pre_genetab[,1] %in% notin_path)),]
    gene_notin<-pre_genetab[which(pre_genetab[,1] %in% notin_path),]
    row.names(genetab)=seq_len(nrow(genetab))
    row.names(gene_notin)=seq_len(nrow(gene_notin))
    names(genetab)<-names(gene_notin)<-c("gene_symbol","name")
    genetab<-genetab[order(genetab[,1]), ]
    gene_notin<-gene_notin[order(gene_notin[,1]), ]

    #metabolites informations
    metatab<-vapply(meta, function(m){
        c(m, paste(keggchebiname[keggchebiname[,3] == m, 2],
                    collapse=", "))
    }, character(2))
    metatab<-as.data.frame(t(metatab))
    row.names(metatab)=seq_len(nrow(metatab))
    names(metatab)<-c("name", "chebi_id")
    metatab<-metatab[order(metatab$name), ]

    #interactions informations
    intetab<-apply(no_path, 1, function(p){
        if(p[1] %in% go_genelist$hgnc_symbol ||
            p[3] %in% go_genelist$hgnc_symbol){
            goterm<-1
        }
        else{goterm<-0}
        c(p[5], p[1], p[2], p[3], goterm, p[4], p[6])
    })
    intetab<-as.data.frame(t(intetab))
    row.names(intetab)=seq_len(nrow(intetab))
    names(intetab)=c("tag", "first_item", "link", "sec_item",
                        "go", "path", "type")
    return(list(genetab, metatab, intetab, gene_notin))
}

#reactome pathways types
type_reactome<-function(sorted_path){
    rea_path<-sorted_path[stringr::str_sub(sorted_path, 1, 3) == "R-H"]
    mapnameid <- as.list(reactome.db::reactomePATHID2NAME) #id-name mapping
    rea_types<-vapply(rea_path,function(r){
        walk<-rea_hiera[stringr::str_detect(rea_hiera[,1],r),1]
        roots<-vapply(walk, function(w){
            chose_root<-stringr::str_split(w,">")[[1]]
            root<-chose_root[1]
            if(chose_root[1] == "R-HSA-1430728" && length(chose_root)>1){
                root<-chose_root[2]
            }
            list(root)
        }, list(1))
        roots<-unname(unlist(roots))
        names<-vapply(roots,function(o){
            if(!is.null(mapnameid[[o]])){
                list(stringr::str_sub(mapnameid[[o]],15))
            }
            else{list("unknown")}
        }, list(1))
        names<-rm_vector(unname(unlist(names)))
        if(length(names)>1 && length(names[names %in% "unknown"])>0){
            names<-names[!(names %in% "unknown")]
        }
        toadd<-vapply(names, function(n){
            c(r,n)
        }, character(2))
        list(as.character(toadd))
    }, list(1))
    rea_types<-unname(unlist(rea_types))
    rea_types<-as.data.frame(matrix(rea_types, ncol=2, byrow=TRUE))
    rea_types<-rm_df(rea_types,c(seq_len(ncol(rea_types))))
    names(rea_types)=c("id","root")
    row.names(rea_types)=seq_len(nrow(rea_types))
    return(rea_types)
}

#kegg pathways types
type_kegg<-function(sorted_path){
    kegg_path<-sorted_path[stringr::str_sub(sorted_path, 1, 3) == "hsa"]
    kegg_types<-vapply(kegg_path, function(k){
        path<-paste(stringr::str_sub(k, 1, 3),
                    stringr::str_sub(k, 5), sep="")
        hiera<-kegg_hiera[stringr::str_detect(kegg_hiera[, 1], path), ]
        if(length(hiera)>0){
            pre_types<-vapply(hiera,function(h){
                root<-stringr::str_split(h, ">")[[1]][1]
                if(root %in% c("Metabolism", "Genetic Information Processing",
                                    "Cellular Processes", "Organismal Systems",
                                    "Human Diseases", "Drug Development")){
                    root<-stringr::str_split(h, ">")[[1]][2]
                }
                c(k,root)
            }, character(2))
            as.character(pre_types)
        }
        else if (path == "hsa01100"){
            c(k, "Global")
        }
        else{
            c(k, "unknown")
        }
    }, character(2))
    kegg_types<-unname(unlist(kegg_types))
    kegg_types<-as.data.frame(matrix(kegg_types, ncol=2, byrow=TRUE))
    names(kegg_types)=c("id","root")
    return(kegg_types)
}

#wikipathways pathways types
type_wp<-function(sorted_path){
    wp_path<-sorted_path[stringr::str_sub(sorted_path, 1, 2) == "WP"]
    wp_types<-vapply(wp_path, function(w){
        hiera<-wp_hiera[stringr::str_detect(wp_hiera[, 1], w), ]
        if(length(hiera)>0){
            pre_types<-vapply(hiera,function(h){
                walk<-stringr::str_split(h, ">")[[1]]
                root<-walk[1]
                if(root%in%c("classic metabolic pathway", "regulatory pathway")
                    && length(walk)>2){
                    root<-walk[2]
                }
                c(w,root)
            }, character(2))
            list(as.character(pre_types))
        }
        else{
            list(c(w, "unknown"))
        }
    }, list(1))
    wp_types<-unname(unlist(wp_types))
    if(is.null(wp_types)){
        wp_types<-data.frame(id=character(), root=character())
    }
    else{
        wp_types<-as.data.frame(matrix(wp_types, ncol=2, byrow=TRUE))
        names(wp_types)=c("id","root")
    }
    return(wp_types)
}

#pathways types=roots of pathways hierarchy
type_path<-function(sorted_path, hierapath){
    kegg_type<-type_kegg(sorted_path)#kegg types
    rea_types<-type_reactome(sorted_path)#Reactome types
    wp_types<-type_wp(sorted_path)#wikipathways types

    if(nrow(kegg_type)>0){
        kegg_type<-apply(kegg_type,1,function(k){
            if(length(which(k[2] %in% types_links[,1]))>0){
                c(k[1], types_links[which(types_links[, 1] %in% k[2]), 2])
            }
            else{k}
        })
        kegg_type<-as.data.frame(t(kegg_type))
    }

    if(nrow(wp_types)>0){
        wp_types<-apply(wp_types,1,function(w){
            if(length(which(w[2] %in% types_links[,3]))>0){
                c(w[1], types_links[which(types_links[, 3] %in% w[2]), 2])
            }
            else{w}
        })
        wp_types<-as.data.frame(t(wp_types))
    }

    names(rea_types)=names(kegg_type)=names(wp_types)=c("id","root")
    types<-rbind(rea_types, kegg_type, wp_types)

    ##list of concerned hierarchies by pathways types
    path_types<-unique(types$root)
    hieratypes<-vapply(path_types, function(p){
        type_path<-types[types[, 2] %in% p, 1]#concerned pathways
        index<-lapply(hierapath, function(h){
            if(length(intersect(type_path, h[["name"]]))>0){h[["index"]]}
        })
        name<-lapply(hierapath, function(h){
            if(length(intersect(type_path, h[["name"]]))>0){h[["name"]]}
        })
        list(list(index=unlist(index), name=unname(unlist(name))))
    }, list(1))
    return(list(types, hieratypes))
}

#gene types=uniprot keywords from BP and MF
type_elem<-function(genetab){
    pre_genetype<-uni_terms[uni_terms[,2] %in% genetab[,1],]
    genetype<-vapply(unique(pre_genetype[,1]), function(u){
        list(pre_genetype[which(pre_genetype[,1] %in% u), 2])
    }, list(1))

    genes_with_kw<-rm_vector(unname(unlist(genetype)))
    if(length(genetab[,1][!(genetab[,1] %in% genes_with_kw)])>0){
        genetype[["unknown"]]<-genetab[,1][!(genetab[,1] %in% genes_with_kw)]
    }
    return(genetype)
}

#build centrality dataframe
build_centrality<-function(htmap, central, size){
    centrality<-htmap[-(seq_len(nrow(htmap))), c(6:ncol(htmap))]
    centrality<-apply(htmap, 1, function(s){ #centrality
        col<-vapply(seq_len(ncol(centrality)), function(c){
            element<-colnames(centrality[c])
            numerator<-central[[s[2]]][[element]]
            denominator=(size[size[, 1] %in% s[2], 4] +
                                size[size[, 1] %in% s[2],8])
            centralratio<-round(numerator/denominator*100, digits=1)
            if(s[element] == 1){
                if(length(centralratio)>0){list(centralratio+10)}
                else{list(NA)}
            }
            else{list(NA)}
        }, list(1))
    })
    centrality<-as.data.frame(matrix(unlist(centrality),
                                    ncol=ncol(htmap[-(seq_len(nrow(htmap))),
                                            c(6:ncol(htmap))]), byrow=TRUE))
    names(centrality)=names(htmap[-(seq_len(nrow(htmap))),
                                    c(6:ncol(htmap))])
    return(centrality)
}

#build sub dataframe
build_sub<-function(htmap,central, size, tagged){
    sub<-htmap[-(seq_len(nrow(htmap))), c(6:ncol(htmap))]
    sub<-apply(htmap, 1, function(s){ #sub
        col<-vapply(seq_len(ncol(sub)), function(c){
            path<-stringr::str_split(s[1], "__")[[1]]
            element<-colnames(sub[c])
            numerator<-central[[s[2]]][[element]]
            denominator=(size[size[, 1] %in% s[2], 4] +
                                size[size[, 1] %in% s[2], 8])
            centralratio<-round(numerator/denominator*100, digits=1)
            if(s[element] == 1){
                elem_inter<-c(tagged[which(tagged[, 4] %in% s[2]), 1],
                                tagged[which(tagged[, 4] %in% s[2]), 3])
                if(element %in% elem_inter){
                    with<-tagged[which(tagged[, 4] %in% s[2]), ]
                    with<-paste(rm_vector(c(with[which(with[,1]%in%element),3],
                                    with[which(with[,3] %in% element), 1])),
                                collapse=", ")
                    list(paste("x : ",element,"\ny : ",path[length(path)],
                            "\nId : ", s[2], "\nGene ratio : ",
                            s[3], "\nMeta ratio : ", s[4],
                            "\nCentralite : ", paste(
                                central[[s[2]]][[element]], " (",
                                centralratio,"%)",sep=""),"\nInteract with : ",
                            with, sep=""))
                }
                else{
                    list(paste("x : ", element, "\ny : ", path[length(path)],
                            "\nId : ", s[2], "\nGene ratio : ",
                            s[3],"\nMeta ratio : ", s[4],
                            "\nCentralite : ", paste(
                                central[[s[2]]][[element]], " (",
                                centralratio, "%)", sep=""), sep=""))
                }
            }
            else{list(NA)}
        }, list(1))
    })
    sub<-as.data.frame(matrix(unlist(sub),
                                ncol=ncol(htmap[-(seq_len(nrow(htmap))),
                                            c(6:ncol(htmap))]), byrow=TRUE))
    names(sub)=names(htmap[-(seq_len(nrow(htmap))),c(6:ncol(htmap))])
    return(sub)
}

#build inter_values dataframe
build_inter_values<-function(htmap, tagged, gene_list, meta_list){
    inter_values<-htmap[-(seq_len(nrow(htmap))), c(6:ncol(htmap))]
    inter_values<-apply(htmap, 1, function(s){ #inter_values
        col<-vapply(seq_len(ncol(inter_values)), function(c){
            element<-colnames(inter_values[c])
            if(s[element] == 1){
                elem_inter<-c(tagged[which(tagged[, 4] %in% s[2]), 1],
                                tagged[which(tagged[, 4] %in% s[2]), 3])
                if(element %in% elem_inter){
                    if(element %in% gene_list){list(3)}
                    else if(element %in% meta_list){list(2)}
                }
                else{list(1)}
            }
            else{list(NA)}
        }, list(1))
    })
    inter_values<-as.data.frame(matrix(unlist(inter_values),
                                    ncol=ncol(htmap[-(seq_len(nrow(htmap))),
                                            c(6:ncol(htmap))]), byrow=TRUE))
    names(inter_values)=names(htmap[-(seq_len(nrow(htmap))),c(6:ncol(htmap))])
    return(inter_values)
}

##Dataframes of values to color the heatmap. Values are centrality, up/down
##quantitative values and direct interactions.
##Also, subtitles of the heatmap are built.
values_shiny<-function(htmap, central, size, tagged, gene_list, meta_list){
    centrality<-build_centrality(htmap, central, size)
    sub<-build_sub(htmap,central, size, tagged)
    inter_values<-build_inter_values(htmap, tagged, gene_list, meta_list)

    sub[is.na(sub)]<-""
    sub<-cbind(htmap[, c(seq_len(5))], sub); sub<-as.matrix(sub)
    centrality<-cbind(htmap[, c(seq_len(5))], centrality)
    centrality[,c(seq_len(5))]<-100; centrality[is.na(centrality)]<-0
    centrality<-as.matrix(centrality)
    centrality<-mapply(centrality, FUN=as.integer)
    centrality<-matrix(data=centrality, ncol=ncol(sub), nrow=nrow(sub))
    inter_values<-cbind(htmap[, c(seq_len(5))], inter_values)
    inter_values[,c(seq_len(5))]<-100; inter_values[is.na(inter_values)]<-0
    inter_values<-as.matrix(inter_values)
    inter_values<-mapply(inter_values, FUN=as.integer)
    inter_values<-matrix(data=inter_values, ncol=ncol(sub), nrow=nrow(sub))
    return(list(centrality, inter_values, sub))
}

#filter pathways regarding user's element ratio and direct interactions
filter_path<-function(tagged,size){
    path_inter<-as.vector(tagged[,4])
    sorted_path<-apply(size,1,function(x){ #sort pathways obtained
        path_elem<-as.integer(x[4])+as.integer(x[8])
        if (path_elem>0){
            num<-as.integer(x[2])+as.integer(x[6])
            if(num/path_elem>=0.2){x[1]}
        }
    })
    sorted_path<-unname(unlist(sorted_path))
    sorted_path<-rm_vector(c(sorted_path[!is.na(sorted_path)],path_inter))
    return(sorted_path)
}

compl_data<-function(listparam){
    size<-listparam[[1]]; pathways<-listparam[[2]]; tagged<-listparam[[3]];
    keggchebiname<-listparam[[4]];
    central<-listparam[[5]]; no_path<-listparam[[6]];
    gene_list<-rm_vector(listparam[[7]]);
    meta_list<-rm_vector(listparam[[8]]);list_elem<-c(gene_list, meta_list)

    sorted_path<-filter_path(tagged,size)
    listpath<-sort_hiera(sorted_path)
    path_walks<-listpath[[1]]; max<-listpath[[2]]
    path_walks<-tidyr::separate(path_walks, 1, as.character(c(seq_len(max))),
                                sep=">", extra="drop", fill="right")
    treeview<-tree_view(path_walks);names(treeview)<-c(seq_len(ncol(treeview)))
    listtab<-final_tab(treeview, pathways, size, sorted_path, no_path,
                        list_elem, tagged)
    heatmap<-listtab[[1]]; notin_path<-listtab[[2]]; hierapath<-listtab[[3]]
    save_cluster_elem<-listtab[[4]]
    listgo<-enr_go(gene_list) #go terms enrichment
    go_genelist<-listgo[[1]]; allResBP<-listgo[[2]]; allResMF<-listgo[[3]]
    listelm<-infos_elem(gene_list, notin_path, meta_list, keggchebiname,
                        no_path, go_genelist)
    genetab<-listelm[[1]]; metatab<-listelm[[2]]; intetab<-listelm[[3]]
    gene_notin<-listelm[[4]];heatmap[heatmap == "X"]<-1
    heatmap[is.na(heatmap)]<-0; heatmap[heatmap[,1] == 0, 1]<-""
    listype<-type_path(sorted_path, hierapath)
    types<-listype[[1]]; hierabrite<-listype[[2]]
    genetype<-type_elem(genetab)
    listparam[[7]]<-rm_vector(listparam[[7]])
    listparam[[8]]<-rm_vector(listparam[[8]])
    listval<-values_shiny(heatmap, central, size, tagged, gene_list, meta_list)
    centrality<-listval[[1]]; inter_values<-listval[[2]]; sub<-listval[[3]]
    intetab<-cbind(intetab, cat=rep(NA, nrow(intetab)),
                    catnumb=rep(NA, nrow(intetab)))
    intetab<-apply(intetab, 1, function(i){
        path_cat<-stringr::str_split(i[6], ", ")[[1]]
        cate<-paste(rm_vector(types[which(types$id
                                          %in% path_cat),2]), collapse=", ")
        c(i[seq_len(7)], cate, stringr::str_count(cate, ", ")+1)
    })
    intetab<-as.data.frame(t(intetab))
    names(intetab)<-c("tag", "first_item", "link", "sec_item", "go", "path",
                        "type", "cat", "catnumb")
    gobp_tab<-hieraGO("BP", allResBP, go_genelist)
    gomf_tab<-hieraGO("MF", allResMF, go_genelist)
    gomflist<-list_go(gomf_tab);gobplist<-list_go(gobp_tab)
    return(list(heatmap, meta_list, allResBP, go_genelist, allResMF, types,
                genetype, metatab, genetab, intetab, gomf_tab, gobp_tab,
                gene_list, gomflist, gobplist, hierabrite,
                hierapath, save_cluster_elem, centrality, inter_values,
                gene_notin, sub))
}
