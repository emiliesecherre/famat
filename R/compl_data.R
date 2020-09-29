#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export compl_data
###############################FUNCTIONS##########################
##build walks using 2-2 relations between items
##for every parent-child relation, if a parent of the parent term is found,
##put this parent before the parent term. Same is done for children.
build_walk<-function(df){
    walks<-data.frame(hier=character())
    for (n in seq_len(nrow(df))){walks[n, 1]<-paste(df[n, ], collapse=">")}
    i<-1
    while(i<=nrow(walks)){
        start<-1
        while(start == 1){ #as long as a walk is modified
            start<-0
            studied<-stringr::str_split(walks[i, 1], ">")[[1]]#walk studied
            parent<-df[df[, 2] == studied[1], 1]
            if(length(parent)>1){
                walks[i, 1]<-paste(c(parent[1], studied), collapse=">")
                for (m in 2:length(parent)){#duplicate walk if multiple parents
                    walks<-rbind(walks,
                                paste(c(parent[m], studied), collapse=">"))
                }
                start<-1
            }
            else if (length(parent) == 1){
                walks[i, 1]<-paste(c(parent[1], studied), collapse=">")
                start<-1
            }

            child<-df[df[, 1] == studied[length(studied)], 2]
            if(length(child)>1){
                walks[i, 1]<-paste(c(studied, child[1]), collapse=">")
                for (m in 2:length(child)){#duplicate walk if multiple children
                    walks<-rbind(walks,
                                paste(c(studied, child[m]), collapse=">"))
                }
                start<-1
            }
            else if (length(child) == 1){
                walks[i, 1]<-paste(c(studied, child[1]), collapse=">")
                start<-1
            }
        }
        check<-as.data.frame(rm_df(walks, c(1)))
        if(nrow(check) == nrow(walks)){
            i<-i+1
        }
        else{ #duplicated walks : go back to be sure to treat all walks
            nbd<-which(duplicated(walks[, 1]));nbd<-nbd[nbd<i]
            walks<-as.data.frame(rm_df(walks, c(1)))
            i<-i-length(nbd)
        }
    }
    names(walks)<-c("hier")
    return(walks)
}

#gather go terms hierarchies, their index, go terms and genes
listgo<-function(hgo){
    golist<-list()
    root_index<-which(!(stringr::str_sub(hgo[, 1], 1, 1)=="_"))
    for (g in seq_len(length(root_index))){
        if (root_index[g]!=root_index[length(root_index)]){
        hiera_index<-c(root_index[g]:(root_index[g+1]-1))
    }
    else{
        hiera_index<-c(root_index[g]:nrow(hgo))
    }
    genes_hiera<-go_terms<-vector()
    for (h in seq_len(length(hiera_index))){
        genes_hiera<-c(genes_hiera,
                        stringr::str_split(hgo[hiera_index[h], 3], ", ")[[1]])
        gomed<-stringr::str_split(hgo[hiera_index[h], 1], "__")[[1]]
        gomed<-gomed[gomed!=""]
        gomed<-stringr::str_split(gomed, " / ")[[1]]
        go_terms<-c(go_terms, gomed)
    }
    genes_hiera<-genes_hiera[!is.na(genes_hiera)]
    genes_hiera<-rm_vector(genes_hiera[genes_hiera!=""])
    go_terms<-go_terms[!is.na(go_terms)]
    go_terms<-rm_vector(go_terms[go_terms!=""])
    golist[[g]]<-list(goterm=go_terms, index=hiera_index, gene=genes_hiera)
    }
    return(golist)
}

##for every walk ending by a same leaf, if one or more walk contain a node
##required in hierarchies but not present in leaves, remove other walks
sort_by_leaves<-function(df, notin_leaf){
    df[]<-lapply(df, as.character)
    pre_walk<-build_walk(df)
    pre_walk<-rm_vector(pre_walk$hier)
    walk<-vector()
    b<-1
    while (b<=length(pre_walk)){
        len<-length(pre_walk)#to know if walks are removed
        leaf<-stringr::str_sub(pre_walk[b],
                                nchar(pre_walk[b])-9, nchar(pre_walk[b]))
        leaf_walks<-pre_walk[stringr::str_sub(pre_walk, nchar(pre_walk)-9,
                                                nchar(pre_walk)) == leaf]
        match<-vector()
        if(length(leaf_walks)>1){
            for (n in seq_len(length(notin_leaf))){
                match<-c(match, leaf_walks[stringr::str_detect(leaf_walks,
                                                            notin_leaf[n])])
            }
            if(length(match) == 0){#keep the smaller walk
                match<-leaf_walks[stringr::str_count(leaf_walks, ">") == min(
                    stringr::str_count(leaf_walks, ">"))][1]
            }
        }
        else{match=leaf_walks}
        walk<-c(walk, match)
        pre_walk<-dplyr::setdiff(pre_walk, leaf_walks)#remove analyzed walks
        if(length(pre_walk) == len){b<-b+1}
    }
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
                        stop<-1 ; dif<-vector()
                        for (v in seq_len(length(studied[a, ]))){
                            if(!is.na(studied[a, v]) &&
                                !is.na(studied[i, v])){
                                term_one<-stringr::str_split(studied[a, v],
                                                                " / ")[[1]]
                                term_two<-stringr::str_split(studied[i, v],
                                                                " / ")[[1]]
                                if(length(intersect(term_one, term_two)) == 0){
                                    dif<-c(dif, v)} #difference between terms
                            }
                        }
                        if(length(dif) == 1){
                            studied[a, dif]<-paste(studied[a, dif]," / ",
                                                    studied[i, dif], sep="")
                            studied<-studied[-i, ]; stop<-0}
                        else if (length(dif) == 0 && a!=i){
                            for (v in seq_len(length(studied[a, ]))){
                                if(!is.na(studied[a, v]) &&
                                    !is.na(studied[i, v])){
                                    goterms<-rm_vector(c(term_one, term_two))
                                    studied[a, dif]<-paste(goterms,
                                                            collapse=" / ")
                                }
                            }
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
    colnames(links)<-c("from", "to")
    leaves<-as.vector(unique(links[which(!(links$to %in% links$from)), 2]))
    notin_leaf<-dplyr::setdiff(res_enrich[, 1], leaves)

    walk<-sort_by_leaves(links, notin_leaf)
    #walks included into others ?
    p<-1
    while(p<=nrow(walk)){
        walk_one<-walk[p, ]
        walk_one<-walk_one[!is.na(walk_one)]
        h<-1
        while (h<=nrow(walk)){
            walk_two<-walk[h, ]; walk_two<-walk_two[!is.na(walk_two)]
            common<-intersect(walk_one, walk_two)
            if (length(common) == min(length(walk_one), length(walk_two))){
                if(length(walk_one)>length(walk_two)){
                    walk<-walk[-h, ]
                }
                else if(length(walk_two)>length(walk_one)){
                    walk<-walk[-p, ]
                }
                else{h<-h+1}
            }
            else{h<-h+1}
        }
        final_walk<-walk[p, ]; final_walk<-final_walk[!is.na(final_walk)]
        if(length(intersect(final_walk, walk_one)) == length(walk_one)){p<-p+1}
    }
    walk<-merge_walks(walk)
    return(walk)
}

##build a hierarchy in tree view form
##for ordered walks beginning with the same root, this function goes from the
##last walk to the second. For a node, if in the walk above the walk studied
##there is the same node at the same position in the hierarchy, remove this
##node. IN the dataframe, the position is the column. So, nodes in the df
##will look like a tree.
tree_view<-function(walks){
    pre_walks<-treeview<-walks[-(seq_len(nrow(walks))), ]
    go_terms_roots<-rm_vector(walks[, 1])

    for (v in seq_len(length(go_terms_roots))){
        roots_walks<-walks[walks[, 1] == go_terms_roots[v], ]

        if(nrow(roots_walks)>1){
            r<-nrow(roots_walks)
            while (r >= 2){
                c<-1
                stop<-0
                while(stop == 0 && c<=ncol(roots_walks)){
                    if(!is.na(roots_walks[r, c]) && !is.na(roots_walks[r-1, c])
                        && roots_walks[r, c] == roots_walks[r-1, c]){
                        roots_walks[r, c]<-NA
                    }
                    else{
                        stop<-1
                    }
                    c<-c+1
                }
                r<-r-1
            }
        }
        pre_walks<-rbind(pre_walks, roots_walks)
    }

    #add a space between each hierarchy
    for (f in seq_len(nrow(pre_walks))){
        for (c in seq_len(ncol(pre_walks))){
            if(!is.na(pre_walks[f, c])){
                treeview<-rbind(treeview, rep(NA, ncol(pre_walks)))
                treeview[nrow(treeview), c]<-pre_walks[f, c]
            }
        }
    }
    return(treeview)
}

#add Go terms names and genes from user's list to GO hierarchies
info_go<-function(go_tab, go_genelist, res_enrich){
    map_names<-as.list(GO.db::GOTERM)
    go_tab<-cbind(go_tab,name=rep(NA,nrow(go_tab)), genes=rep(NA,nrow(go_tab)))

    for (g in seq_len(nrow(go_tab))){
        go_terms<-stringr::str_split(go_tab[g, 1], " / ")[[1]]
        go_terms<-c(stringr::str_split(go_terms[1], "__")[[1]],
                                        go_terms[2:length(go_terms)])
        go_terms<-go_terms[go_terms!=""]
        go_terms<-rm_vector(go_terms[!is.na(go_terms)])

        go_names<-vector()
        pre_go_genes<-vector()
        go_genes<-vector()
        for (l in seq_len(length(go_terms))){
            if (is.null(map_names[[go_terms[l]]])){
                go_names<-c(go_names,
                            res_enrich[res_enrich[, 1] == go_terms[l], 2])
            }
            else{go_names<-c(go_names, map_names[[go_terms[l]]]@Term)}
            pre_go_genes<-go_genelist[go_genelist$go_id == go_terms[l], 1]
            pre_go_genes<-rm_vector(pre_go_genes[!is.na(pre_go_genes)])
            go_genes<-c(go_genes, paste(pre_go_genes, collapse=", "))
        }
        go_names<-go_names[!is.na(go_names)]
        go_tab[g, 2]<-paste(go_names, collapse=" / " )
        go_genes<-go_genes[go_genes!=""]
        go_tab[g, 3]<-paste(go_genes, collapse=" / ")
    }
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
    go_links<-data.frame(parent=character(), child=character())
    for (x in seq_len(length(ancestors))){
        child<-names(ancestors[x]); parent<-ancestors[[x]]
        for (p in seq_len(length(parent))){
            nr<-nrow(go_links)+1
            if (unname(parent[p]) %in% candidate && child %in% candidate){
                go_links[nr, 1]<-unname(parent[p])
                go_links[nr, 2]<-child
            }
        }
    }
    walks<-sort_go_steps(go_links, res_enrich)
    #build walks
    pre_walks<-vector()
    for (b in seq_len(nrow(walks))){
        wk<-walks[b, ];wk<-wk[!is.na(wk)]
        pre_walks<-c(pre_walks, paste(wk, collapse=">"))
    }
    pre_walks<-sort(pre_walks) #useful to build tree view
    pre_walks<-as.data.frame(pre_walks)
    walks<-tidyr::separate(pre_walks, 1, as.character(c(seq_len(ncol(walks)))),
                            sep=">", extra="drop", fill="right")
    walks<-walks[,-1] #remove root term of BP/MF hierarchies
    treeview<-tree_view(walks)

    #reduce treeview to a single column, with "__" instead of NAs
    names(treeview)<-as.character(c(seq_len(ncol(treeview))))
    go_tab<-data.frame(name=character())
    for (n in seq_len(nrow(treeview))){
        index<-which(!is.na(treeview[n,]))
        space<-paste(rep("__", index-1), collapse="")
        go_tab<-rbind(go_tab, rep(NA, ncol(go_tab)))
        node<-treeview[n, !is.na(treeview[n, ])]
        go_tab[n, 1]<-paste(space, node, sep="")
    }
    go_tab<-info_go(go_tab, go_genelist, res_enrich)

    return(go_tab)
}

##find which elements are found in the same pathways, and put them together
##find which pathways contain the same elements also
##if element=T, also return user's elements which aren't in pathways
cluster_items<-function(items, element){
    items[is.na(items)]<-"N";items<-items[-1, ]
    elem_names<-colnames(items)
    dist_data<-data.frame(one=character(), two=character(), dist=integer())
    nb_path<-notin_path<-vector() #elements not in pathways
    for (e in seq_len(length(elem_names))){ #distances between elements
        first_item<-items[, e]
        nb_path<-c(nb_path, length(first_item[first_item %in% "X"]))
        if(element == TRUE && nb_path[e] == 0){
            notin_path<-c(notin_path, elem_names[e])
        }
        if(e != length(elem_names)){
            for(f in (e+1):length(elem_names)){
                dist<-0
                sec_item<-items[,f]
                for(i in seq_len(length(sec_item))){
                    if(first_item[i] != sec_item[i]){
                        dist<-dist+1
                    }
                }
                nr<-nrow(dist_data)+1
                dist_data[nr, 1]<-elem_names[e]
                dist_data[nr, 2]<-elem_names[f];dist_data[nr, 3]<-dist
            }
        }
    }

    ##order elements by distances
    dist_data<-dist_data[order(dist_data$dist),]
    cluster<-vector()
    for (d in seq_len(nrow(dist_data))){
        cluster<-c(cluster, dist_data[d,1], dist_data[d,2])
    }
    cluster<-rm_vector(cluster)
    if (element == TRUE){return(list(cluster, notin_path))}
    else{return(cluster)}
}

#filter entire pathways hierarchy to build a hierarchy concerning our pathways
sort_hiera<-function(pathways){
    kegg_path<-pathways[stringr::str_sub(pathways, 1, 3) == "hsa"]
    kegg_path<-paste(stringr::str_sub(kegg_path,1,3),
                    stringr::str_sub(kegg_path,5),sep="")
    path_walks_k<-vector()
    for(k in seq_len(length(kegg_path))){
        toadd<-kegg_hiera[stringr::str_detect(kegg_hiera[,1],kegg_path[k]),1]
        path_walks_k<-c(path_walks_k,toadd)
    }
    path_walks_k<-as.data.frame(sort(path_walks_k))

    wp_path<-pathways[stringr::str_sub(pathways, 1, 2) == "WP"]
    path_walks_w<-vector()
    for (w in seq_len(length(wp_path))){
        toadd<-wp_hiera[stringr::str_detect(wp_hiera[,1],wp_path[w]),1]
        path_walks_w<-c(path_walks_w,toadd)
    }
    path_walks_w<-as.data.frame(sort(path_walks_w))

    path_walks_r<-rea_hiera[,1];t<-1
    while (t<=length(path_walks_r)){
        rea_walks<-stringr::str_split(path_walks_r[t], ">")[[1]]
        if(length(rea_walks %in% pathways)>0){
            rea_walks<-rm_vector(rea_walks[c(1, which(rea_walks%in%pathways))])
            path_walks_r[t]<-paste(rea_walks, collapse=">")
            t<-t+1
        }
        else{path_walks_r<-path_walks_r[-t]}
    }
    path_walks_r<-path_walks_r[stringr::str_detect(path_walks_r, ">")]
    path_walks_r<-rm_vector(path_walks_r);u<-1
    while(u<=length(path_walks_r)){
        dupl<-which(stringr::str_detect(path_walks_r, path_walks_r[u]))
        dupl<-dupl[-which(dupl == u)]
        if(length(dupl)>0){path_walks_r<-path_walks_r[-u]}
        else{u<-u+1}
    }
    path_walks_r<-as.data.frame(sort(path_walks_r))

    names(path_walks_r)<-names(path_walks_w)<-names(path_walks_k)<-"walks"
    path_walks<-rbind(path_walks_r, path_walks_k,path_walks_w)
    max=max(stringr::str_count(path_walks[,1],">"))+1
    return(list(path_walks, max))
}

#add informations about pathway hierarchies to the final heatmap
#informations are : ratio of user's elements / total elements,
#names and ids of pathways in hierarchies
hiera_info<-function(heatmap, pathways, size, sorted_pathways,
                                                treeview, no_path, list_elem){
    pathidtoname <- as.list(reactome.db::reactomePATHID2NAME)
    for (n in seq_len(nrow(treeview))){
        index<-which(!is.na(treeview[n,]))
        space<-paste(rep("__", index-1), collapse="")
        heatmap<-rbind(heatmap, rep(NA, ncol(heatmap)))
        node<-treeview[n, !is.na(treeview[n, ])]
        if(stringr::str_sub(node, 1, 3) == "hsa"){
            node<-paste("hsa:", stringr::str_sub(node, 4, nchar(node)), sep="")
        }
        name<-pathways[pathways[,2] == node, 1]
        if(length(name)>0){
            heatmap[n, 1]<-paste(space,
                                    pathways[pathways[,2] == node, 1], sep="")
        }
        else if(stringr::str_sub(node, 1, 5) == "R-HSA"){
            heatmap[n, 1]<-paste(space, stringr::str_sub(pathidtoname[[node]],
                                    15, nchar(pathidtoname[[node]])), sep="")
        }
        else{heatmap[n, 1]<-paste(space, node, sep="")}
        heatmap[n, 2]<-node
        heatmap[n, 3]<-paste("'", size[size$path == node, 2], "/",
                                size[size$path == node, 4], sep="")
        heatmap[n, 4]<-paste("'", size[size$path == node, 6], "/",
                                size[size$path == node, 8], sep="")
    }
    heatmap<-rbind(rep(NA, ncol(heatmap)), heatmap); hiera_nodes<-heatmap[, 2]
    hiera_nodes<-rm_vector(hiera_nodes[!is.na(hiera_nodes)])
    for (m in seq_len(length(sorted_pathways))){
        if(!(sorted_pathways[m] %in% hiera_nodes)){ #pathways not in hierarchy
            path_entry<-size[size$path == sorted_pathways[m], ]
            heatmap<-rbind(heatmap,
                        c(pathways[pathways[,2] == sorted_pathways[m], 1],
                        sorted_pathways[m], paste("'", path_entry[, 2], "/",
                        path_entry[, 4], sep=""), paste("'", path_entry[, 6],
                        "/", path_entry[, 8], sep=""), NA))
        }
    }
    colnames(heatmap)<-c("un", "deux", "trois", "quatre", "cinq")
    heatmap[which(heatmap[, 3] == "'/"), 3]<-"'0/0"
    heatmap[which(heatmap[, 4] == "'/"), 4]<-"'0/0";tags<-no_path$tag
    elements<-c(list_elem, tags); elements_df<-t(as.data.frame(elements))
    elements_df<-as.data.frame(elements_df); colnames(elements_df)<-elements
    for (d in seq_len((nrow(heatmap)-1))){
        elements_df<-rbind(elements_df, rep(NA, ncol(elements_df)))
    }
    heatmap<-cbind(heatmap, elements_df)
    return(heatmap)
}

#help building hierapath list
hierapath_build<-function(to_build, heatmap, size, tagged, no_path,
                                                            elements, index){
    to_build<-rbind(to_build, heatmap[index ,])
    pre_elements<-as.vector(t(size[which(size[, 1] == heatmap[index, 2]),
                            c(3, 7)]))
    pre_elements<-pre_elements[!is.na(pre_elements)]
    if(length(pre_elements)>0){
        for(e in seq_len(length(pre_elements))){
            elements<-c(elements,
                        stringr::str_split(pre_elements[e], " # ")[[1]])
        }
    }
    interactions<-tagged[which(tagged$path == heatmap[index, 2]), 5]
    interactions<-interactions[interactions %in% no_path$tag]
    elements<-c(elements, interactions)

    return(list(to_build, elements))
}

##cluster hierarchies using their elements
##to do so, all roots of hierarchies are used, and elements from all
##the hierarchy pathways are added to the root
cluster_hiera<-function(heatmap, size, tagged, no_path){
    help_hm<-heatmap<-heatmap[-1, ]
    root<-help_hm[-(seq_len(nrow(help_hm))), ]
    index<-which(stringr::str_sub(help_hm$un, 1, 1)!="_")
    for (t in seq_len(length(index))){
        if(index[t] == index[length(index)] && index[t]!=nrow(help_hm)){
            hiera<-help_hm[c(index[t]:nrow(help_hm)), ]}
        else if(index[t]+1 == index[t+1] || index[t] == nrow(help_hm)){
            hiera<-help_hm[index[t], ]}
        else{hiera<-help_hm[c(index[t]:(index[t+1]-1)), ]}
        hiera<-rbind(rep(NA, ncol(hiera)), hiera)
        hiera[1, 1]<-hiera[2, 2]
        for (c in 2:ncol(hiera)){ #add elements to roots
            if("X" %in% hiera[, c]){hiera[1, c]<-"X"}
            else{hiera[1, c]<-NA}
        }
        root<-rbind(root, hiera[1, ])
    }
    root_data<-root[seq_len(nrow(root)), ]
    root_data<-t(root_data);colnames(root_data)=root_data[1,]
    clusters<-cluster_items(root_data, FALSE)
    heatmap<-help_hm[-(seq_len(nrow(heatmap))), ]
    hierapath<-list(); a<-1 #hierarchies list
    for (r in seq_len(length(clusters))){
        root_id<-which(help_hm[, 2] == clusters[r])
        next_root<-index[which(index==root_id)+1]
        if(is.na(next_root)){next_root=nrow(help_hm)}
        if (root_id == nrow(help_hm) || next_root == root_id+1){#one path
            element<-vector()
            listbuid<-hierapath_build(heatmap, help_hm, size, tagged, no_path,
                                    element, root_id)
            heatmap<-listbuid[[1]]; element<-listbuid[[2]]
            hierapath[[a]]<-list(name=help_hm[root_id, 2],index=nrow(heatmap),
                                elem=element)
            a<-a+1}
        else{ #hierarchy
            nb_path<-next_root-root_id-1;path_id<-path_index<-element<-vector()
            for (d in 0:nb_path){
                path_id<-c(path_id, help_hm[root_id+d, 2])
                listbuid<-hierapath_build(heatmap, help_hm, size, tagged,
                                        no_path, element, root_id+d)
                heatmap<-listbuid[[1]]; element<-listbuid[[2]]
                path_index<-c(path_index, nrow(heatmap))
            }
            hierapath[[a]]<-list(name=path_id, index=path_index,
                                elem=rm_vector(element))
            a<-a+1}
        heatmap<-rbind(heatmap, rep(NA, ncol(heatmap)))}
    return(list(heatmap, hierapath))
}

##build heatmap of hierarchies of pathways and elements included in them
final_tab<-function(build_hm, pathways, size, sorted_path, no_path,
                                                            list_elem, tagged){
    heatmap<-data.frame(name=character(), id=character(), ratiog=character(),
                        ratiom=character(), void=character())
    heatmap<-hiera_info(heatmap, pathways, size, sorted_path, build_hm,
                        no_path, list_elem)
    tags<-no_path$tag #direct interactions
    for (p in 2:nrow(heatmap)){ #if an element is in a pathway : "X"
        path<-heatmap[p, 2]
        element<-vector()
        pre_elem<-c(size[size$path == path, 3], size[size$path == path, 7])
        pre_elem<-pre_elem[!is.na(pre_elem)]
        if(length(pre_elem)>0){
            for (f in seq_len(length(pre_elem))){
                element<-c(element, stringr::str_split(pre_elem[f]," # ")[[1]])
                element<-rm_vector(element)
            }
            for (e in 6:(ncol(heatmap)-length(tags))){ #all except interactions
                if(heatmap[1, e] %in% element){
                    heatmap[p, e]<-"X" #if the element is in the pathway
                }
            }
            for (i in (ncol(heatmap)-length(tags)+1):ncol(heatmap)){
                path_inter<-tagged[tagged$path == path,]
                path_inter<-path_inter[path_inter$tag == heatmap[1, i], ]
                if(nrow(path_inter)>0){
                    heatmap[p, i]<-"X" #if the interaction is in the pathway
                }
    }}}
    ##CLUSTERS
    elem_data<-heatmap[,c(6:(ncol(heatmap)-length(tags)))]
    elem_data<-as.matrix(elem_data);tags<-sort(tags)
    listele<-cluster_items(elem_data, TRUE)
    cluster_elem<-save_cluster_elem<-listele[[1]];notin_path<-listele[[2]]
    cluster_elem<-c(cluster_elem, tags) #add interactions
    cluster_elem<-cluster_elem[!(cluster_elem %in% notin_path)]
    heatmap<-heatmap[,c("un", "deux", "trois", "quatre", "cinq", cluster_elem)]
    listhiera<-cluster_hiera(heatmap, size, tagged, no_path);
    heatmap<-listhiera[[1]]; hierapath<-listhiera[[2]]
    row.names(heatmap)<-c(seq_len(nrow(heatmap)))
    for (h in seq_len(length(hierapath))){
        index<-hierapath[[h]][["index"]]
        hierapath[[h]][["index"]]<-c(index, index[length(index)]+1)
    }
    for (p in seq_len(nrow(heatmap))){ #if the pathway has a direct interaction
        if(heatmap[p, 2] %in% tagged$path){heatmap[p, 1]<-paste(heatmap[p, 1],
                                                            "(*)", sep=" ")}
    }
    return(list(heatmap, notin_path, hierapath, save_cluster_elem))
}

##perform go term enrichment analysis
enr_go<-function(genes,ensembl){
    ##gene-go terms mapping
    go_gene<-biomaRt::getBM(attributes=c("hgnc_symbol", 'go_id'), mart=ensembl)
    go_gene<-go_gene[!(go_gene$go_id == ""),]
    go_gene<-go_gene[!(go_gene$hgnc_symbol == ""),]

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
                                                    go_genelist, namegeneid){
    #genes informations
    genetab<-data.frame(id=character(), name=character())#for genes
    gene_notin<-data.frame(id=character(), name=character())#for notin_path
    for (g in seq_len(length(genes))){ #build genetab and gene_notin
        gene_name<-namegeneid[namegeneid$hgnc_symbol == genes[g], 2]
        if (genes[g] %in% notin_path){
            if(length(gene_name)>0){
                gene_notin<-rbind(gene_notin, c(genes[g], gene_name))
            }
            else{gene_notin<-rbind(gene_notin, c(genes[g], "" ))}
        }
        else{
            if(length(gene_name)>0){
                genetab<-rbind(genetab, c(genes[g], gene_name))
            }
            else{genetab<-rbind(genetab, c(genes[g], "" ))}
        }
    }
    names(genetab)<-names(gene_notin)<-c("id","name")
    genetab<-genetab[order(genetab$id), ]
    gene_notin<-gene_notin[order(gene_notin$id), ]

    #metabolites informations
    metatab<-data.frame(id=character(), name=character())
    for (m in seq_len(length(meta))){#build metatab
        chebi=paste(keggchebiname[keggchebiname[,3] == meta[m], 2],
                    collapse=", ")
        metatab<-rbind(metatab, c(meta[m], chebi))
    }
    names(metatab)<-c("id", "name")
    metatab<-metatab[order(metatab$name), ]

    #interactions informations
    intetab<-data.frame(tag=character(), idu=character(), link=character(),
            idd=character(), go=integer(), path=character(), type=character())
    no_path<-no_path[order(no_path$tag), ]
    for (d in seq_len(nrow(no_path))){
        if(no_path[d, 1] %in% go_genelist$hgnc_symbol ||
            no_path[d, 3] %in% go_genelist$hgnc_symbol){
            goterm<-1
        }
        else{goterm<-0}
        intetab<-rbind(intetab, c(no_path[d, 5],no_path[d, 1], no_path[d, 2],
                        no_path[d, 3], goterm, no_path[d, 4], no_path[d, 6]))
    }
    names(intetab)=c("tag", "idu", "link", "idd", "go", "path", "type")
    return(list(genetab, metatab, intetab, gene_notin))
}

#reactome pathways types
type_reactome<-function(sorted_path){
    rea_path<-sorted_path[stringr::str_sub(sorted_path, 1, 3) == "R-H"]
    mapnameid <- as.list(reactome.db::reactomePATHID2NAME) #id-name mapping
    rea_types<-data.frame(id<-character(), root<-character())

    for(f in seq_len(length(rea_path))){
        walk<-rea_hiera[stringr::str_detect(rea_hiera[,1],rea_path[f]),1]
        roots<-vector()
        for(w in seq_len(length(walk))){
            chose_root<-stringr::str_split(walk[w],">")[[1]]
            root<-chose_root[1]
            if(chose_root[1] == "R-HSA-1430728" && length(chose_root)>1){
                root<-chose_root[2]
            }
            roots<-rm_vector(c(roots, root))
        }
        names<-vector()
        for(r in seq_len(length(root))){
            if(!is.null(mapnameid[[root[r]]])){
                names<-c(names, stringr::str_sub(mapnameid[[root[r]]],15))
                names<-rm_vector(names)
            }
            else{names<-c(names, "unknown")}
        }
        if(length(names)>1 && length(names[names %in% "unknown"])>0){
            names<-names[!(names %in% "unknown")]
        }
        for(n in seq_len(length(names))){
            nr<-nrow(rea_types)+1
            rea_types[nr,1]<-rea_path[f]
            rea_types[nr,2]<-names[n]
        }
    }
    return(rea_types)
}

#kegg pathways types
type_kegg<-function(sorted_path){
    kegg_path<-sorted_path[stringr::str_sub(sorted_path, 1, 3) == "hsa"]
    kegg_types<-data.frame(id=character(), root=character())
    for (k in seq_len(length(kegg_path))){
        path<-paste(stringr::str_sub(kegg_path[k], 1, 3),
                    stringr::str_sub(kegg_path[k], 5), sep="")
        hiera<-kegg_hiera[stringr::str_detect(kegg_hiera[, 1], path), ]
        if(length(hiera)>0){
            for(h in seq_len(length(hiera))){
                nr<-nrow(kegg_types)+1
                root<-stringr::str_split(hiera[h], ">")[[1]][1]
                if(root %in% c("Metabolism", "Genetic Information Processing",
                                "Cellular Processes", "Organismal Systems",
                                "Human Diseases", "Drug Development")){
                    root<-stringr::str_split(hiera[h], ">")[[1]][2]
                }
                kegg_types[nr, 1]<-kegg_path[k]
                kegg_types[nr, 2]<-root
            }
        }
        else if (path == "hsa01100"){
            nr<-nrow(kegg_types)+1
            kegg_types[nr, 1]<-kegg_path[k]
            kegg_types[nr, 2]<-"Global"
        }
        else{
            nr<-nrow(kegg_types)+1
            kegg_types[nr, 1]<-kegg_path[k]
            kegg_types[nr, 2]<-"unknown"
        }
    }
    return(kegg_types)
}

#wikipathways pathways types
type_wp<-function(sorted_path){
    wp_path<-sorted_path[stringr::str_sub(sorted_path, 1, 2) == "WP"]
    wp_types<-data.frame(id=character(), root=character())
    for (w in seq_len(length(wp_path))){
        hiera<-wp_hiera[stringr::str_detect(wp_hiera[, 1], wp_path[w]), ]
        if(length(hiera)>0){
            for(h in seq_len(length(hiera))){
                nr<-nrow(wp_types)+1
                walk<-stringr::str_split(hiera[h], ">")[[1]]
                root<-walk[1]
                if(root%in%c("classic metabolic pathway", "regulatory pathway")
                                                            && length(walk)>2){
                    root<-walk[2]
                }
                wp_types[nr, 1]<-wp_path[w]
                wp_types[nr, 2]<-root
            }
        }
        else{
            nr<-nrow(wp_types)+1
            wp_types[nr, 1]<-wp_path[w]
            wp_types[nr, 2]<-"unknown"
        }
    }
    return(wp_types)
}

#pathways types = roots of pathways hierarchy
type_path<-function(sorted_path, hierapath){
    kegg_type<-type_kegg(sorted_path)#kegg types
    rea_types<-type_reactome(sorted_path)#Reactome types
    wp_types<-type_wp(sorted_path)#wikipathways types

    for (c in seq_len(nrow(types_links))){#convert kegg types to reactome ones
        kegg_type[which(kegg_type[,2]%in%types_links[c,1]),2]<-types_links[c,2]
        wp_types[which(wp_types[,2]%in%types_links[c,3]),2]<-types_links[c,2]
    }
    names(rea_types)=names(kegg_type)=names(wp_types)=c("id","root")
    types<-rbind(rea_types, kegg_type, wp_types)

    ##list of concerned hierarchies by pathways types
    path_types<-unique(types$root)
    hieratypes<-list()
    for (u in seq_len(length(path_types))){
        type_path<-types[types[, 2] %in% path_types[u], 1]#concerned pathways
        index<-name<-vector()
        for (h in seq_len(length(hierapath))){
            if(length(intersect(type_path, hierapath[[h]][["name"]]))>0){
                index<-rm_vector(c(index, hierapath[[h]][["index"]]))
                name<-rm_vector(c(name, hierapath[[h]][["name"]]))
            }
        }
        hieratypes[[path_types[u]]]<-list(index=index, name=name)
    }
    return(list(types, hieratypes))
}

##build genetype list using uniprot keywords from BP and MF
##if a gene has only BP keywords, collect these keywords and see if others
##genes are concerned by them
build_genetype=function(genetab, map_gene_uni, typemf, typebp){
    genetype<-list()
    for (g in seq_len(nrow(genetab))){#if genes are only concerned by BP
        if(genetab[g, 1] %in% map_gene_uni$hgnc_symbol){
            gene_kw<-stringr::str_split(UniprotR::GetProteinAnnontate(
                    map_gene_uni[which(map_gene_uni$hgnc_symbol %in%
                                    genetab[g, 1]), 2] , "keywords"), ";")[[1]]
            if(length(gene_kw[gene_kw %in% typemf]) == 0){
                for(u in seq_len(length(gene_kw))){
                    if((gene_kw[u] %in% typebp) == TRUE){
                        genetype[[gene_kw[u]]]<-rm_vector(c(
                            genetype[[gene_kw[u]]], genetab[g, 1]))
                    }
                }
            }
        }
    }
    genes_bp<-names(genetype) #BP keywords from previous loop
    for (g in seq_len(nrow(genetab))){ #genes concerned by MF terms
        if(genetab[g, 1] %in% map_gene_uni$hgnc_symbol){
            gene_kw<-stringr::str_split(UniprotR::GetProteinAnnontate(
                    map_gene_uni[which(map_gene_uni$hgnc_symbol %in%
                                    genetab[g, 1]), 2], "keywords"), ";")[[1]]
            for(u in seq_len(length(gene_kw))){
                if((gene_kw[u] %in% typemf) == TRUE ||
                    (gene_kw[u] %in% genes_bp) == TRUE){
                    genetype[[gene_kw[u]]]<-rm_vector(c(genetype[[gene_kw[u]]],
                                                        genetab[g, 1]))
                }
            }
        }
    }
    genes_with_kw<-rm_vector(unname(unlist(genetype)))
    if(length(genetab$id[!(genetab$id %in% genes_with_kw)])>0){
        genetype[["unknown"]]<-genetab$id[!(genetab$id %in% genes_with_kw)]
    }
    return(genetype)
}

#gene types = uniprot keywords from BP and MF
type_elem<-function(genetab,ensembl){
    uniprotkw<-readLines('https://www.uniprot.org/docs/keywlist', warn=FALSE)
    mf_keywords<-uniprotkw[which(stringr::str_sub(uniprotkw, 1, 24)
                                %in% "HI   Molecular function:")]
    typemf<-vector()
    for (t in seq_len(length(mf_keywords))){
        parse_data<-rm_vector(stringr::str_split(mf_keywords[t],
                                        "HI   Molecular function: ")[[1]][2])
        parse_data<-stringr::str_sub(parse_data, 1, nchar(parse_data)-1)
        typemf<-rm_vector(c(typemf, stringr::str_split(parse_data, "; ")[[1]]))
    }

    bp_keywords<-uniprotkw[which(stringr::str_sub(uniprotkw, 1, 24)
                                %in% "HI   Biological process:")]
    typebp<-vector()
    for (t in seq_len(length(bp_keywords))){
        parse_data<-rm_vector(stringr::str_split(bp_keywords[t],
                                        "HI   Biological process: ")[[1]][2])
        parse_data<-stringr::str_sub(parse_data, 1, nchar(parse_data)-1)
        typebp<-rm_vector(c(typebp, stringr::str_split(parse_data, "; ")[[1]]))
    }
    map_gene_uni<-biomaRt::getBM(attributes=c("hgnc_symbol",'uniprotswissprot'),
                        filters="hgnc_symbol", values=genetab$id, mart=ensembl)
    map_gene_uni<-rm_df(map_gene_uni[map_gene_uni[, 2]!="", ])
    genetype<-build_genetype(genetab, map_gene_uni, typemf, typebp)
    return(genetype)
}

#help building dataframes for values_shiny
build_values_shiny=function(htmap,central, size, tagged, gene_list, meta_list){
    sub<-centrality<-htmap[-(seq_len(nrow(htmap))), c(6:ncol(htmap))]
    inter_values<-htmap[-(seq_len(nrow(htmap))), c(6:ncol(htmap))]
    for (l in seq_len(nrow(htmap))){
        for (c in seq_len(ncol(sub))){
            path<-stringr::str_split(htmap[l, 1], "__")[[1]]
            element<-colnames(sub[c])
            numerator=central[[htmap[l, 2]]][[element]]
            denominator=(size[size[, 1] %in% htmap[l, 2], 4]+
                            size[size[, 1] %in% htmap[l, 2], 8])
            centralratio<-round(numerator/denominator*100, digits=1)
            if(htmap[l, element] == 1){
                if(length(centralratio)>0){centrality[l, c]<-centralratio+10}
                sub[l, c]<-paste("x : ", element, "\ny : ", path[length(path)],
                                "\nId : ", htmap[l, 2], "\nGene ratio : ",
                                htmap[l, 3],"\nMeta ratio : ", htmap[l, 4],
                                "\nCentralite : ", paste(
                                    central[[htmap[l, 2]]][[element]], " (",
                                    centralratio, "%)", sep=""), sep="")
                elem_inter<-c(tagged[which(tagged[, 4] %in% htmap[l, 2]), 1],
                                tagged[which(tagged[, 4] %in% htmap[l, 2]), 3])
                if(element %in% elem_inter){
                    with<-tagged[which(tagged[, 4] %in% htmap[l, 2]), ]
                    with<-paste(rm_vector(c(with[which(with[,1]%in%element),3],
                            with[which(with[,3]%in%element),1])),collapse=", ")
                    sub[l,c]<-paste("x : ",element,"\ny : ",path[length(path)],
                                "\nId : ", htmap[l, 2], "\nGene ratio : ",
                                htmap[l, 3], "\nMeta ratio : ", htmap[l, 4],
                                "\nCentralite : ", paste(
                                    central[[htmap[l, 2]]][[element]], " (",
                                centralratio,"%)",sep=""),"\nInteract with : ",
                                with, sep="")
                    if(element %in% gene_list){inter_values[l, c]<-3}
                    else if(element %in% meta_list){inter_values[l, c]<-2}
                }
                else{inter_values[l, c]<-1}
            }
        }
    }
    return(list(sub,centrality,inter_values))
}

##Dataframes of values to color the heatmap. Values are centrality, up/down
##quantitative values and direct interactions.
##Also, subtitles of the heatmap are built.
values_shiny<-function(htmap, central, size, tagged, gene_list, meta_list){
    listtab=build_values_shiny(htmap,central,size,tagged, gene_list, meta_list)
    sub=listtab[[1]]; centrality=listtab[[2]]; inter_values=listtab[[3]]

    sub[is.na(sub)]<-""; sub<-rbind(sub, rep("", ncol(sub)))
    sub<-cbind(htmap[, c(seq_len(5))], sub); sub<-as.matrix(sub)
    centrality<-rbind(centrality, rep(0, ncol(centrality)))
    centrality<-cbind(htmap[, c(seq_len(5))], centrality)
    centrality[,c(seq_len(5))]<-100; centrality[is.na(centrality)]<-0
    centrality<-as.matrix(centrality)
    centrality<-mapply(centrality, FUN=as.integer)
    centrality<-matrix(data=centrality, ncol=ncol(sub), nrow=nrow(sub))
    inter_values<-rbind(inter_values, rep(0, ncol(inter_values)))
    inter_values<-cbind(htmap[, c(seq_len(5))], inter_values)
    inter_values[,c(seq_len(5))]<-100; inter_values[is.na(inter_values)]<-0
    inter_values<-as.matrix(inter_values)
    inter_values<-mapply(inter_values, FUN=as.integer)
    inter_values<-matrix(data=inter_values, ncol=ncol(sub), nrow=nrow(sub))
    return(list(centrality, inter_values, sub))
}

#filter pathways regarding user's element ratio and direct interactions
filter_path=function(tagged,size){
    path_inter<-as.vector(tagged[,4])
    sorted_path<-vector()
    for (f in seq_len(nrow(size))){ #sort pathways obtained
        path_elem<-size[f, 4]+size[f, 8]
        if (path_elem>0){
            num<-size[f, 2]+size[f, 6]
            if(num/path_elem>=0.2){sorted_path<-c(sorted_path, size[f, 1])}
        }
    }
    sorted_path<-rm_vector(c(sorted_path, path_inter))
    return(sorted_path)
}

compl_data<-function(listparam){
    size<-listparam[[1]]; pathways<-listparam[[2]]; tagged<-listparam[[3]];
    namegeneid<-listparam[[4]]; keggchebiname<-listparam[[5]];
    central<-listparam[[6]]; no_path<-listparam[[7]]; ensembl=listparam[[10]];
    gene_list<-rm_vector(listparam[[8]]);
    meta_list<-rm_vector(listparam[[9]]);list_elem<-c(gene_list, meta_list)

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
    listgo<-enr_go(gene_list,ensembl) #go terms enrichment
    go_genelist<-listgo[[1]]; allResBP<-listgo[[2]]; allResMF<-listgo[[3]]
    listelm<-infos_elem(gene_list, notin_path, meta_list, keggchebiname,
                                            no_path, go_genelist, namegeneid)
    genetab<-listelm[[1]]; metatab<-listelm[[2]]; intetab<-listelm[[3]]
    gene_notin<-listelm[[4]];heatmap[heatmap == "X"]<-1
    heatmap[is.na(heatmap)]<-0; heatmap[heatmap[,1] == 0, 1]<-""
    listype<-type_path(sorted_path, hierapath)
    types<-listype[[1]]; hierabrite<-listype[[2]]
    genetype<-type_elem(genetab,ensembl)
    listparam[[8]]<-rm_vector(listparam[[8]])
    listparam[[9]]<-rm_vector(listparam[[9]])
    listval<-values_shiny(heatmap, central, size, tagged, gene_list, meta_list)
    centrality<-listval[[1]]; inter_values<-listval[[2]]; sub<-listval[[3]]
    intetab<-cbind(intetab, cat=rep(NA, nrow(intetab)))
    for (d in seq_len(nrow(intetab))){#interactions categories
        path_cat<-stringr::str_split(intetab[d, 6], ", ")[[1]]
        intetab[d,8]<-paste(rm_vector(types[which(types$id %in% path_cat), 2]),
                            collapse=", ")
    }
    gobp_tab<-hieraGO("BP", allResBP, go_genelist)
    gomf_tab<-hieraGO("MF", allResMF, go_genelist)
    gomflist<-listgo(gomf_tab);gobplist<-listgo(gobp_tab)
    return(list(heatmap, meta_list, allResBP, go_genelist, allResMF, types,
        genetype, metatab, genetab, intetab, gomf_tab, gobp_tab, gene_list,
        gomflist, gobplist, namegeneid,hierabrite, hierapath,save_cluster_elem,
        centrality, inter_values, gene_notin, sub))
}
