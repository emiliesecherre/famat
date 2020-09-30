#' @importFrom magrittr %>%
#' @export rshiny
#####UI#####
ui<-function(types, genetype, gotermsgene){
    ui <- shinydashboard::dashboardPage(
        skin = "blue",
        shinydashboard::dashboardHeader(title="Famat"),
        shinydashboard::dashboardSidebar(
            shinydashboard::sidebarMenu(
                shinydashboard::menuItem("Elements", tabName = "elements",
                                            icon=shiny::icon("check-square")),
                shinydashboard::menuItem("Pathways", tabName = "pathways",
                                            icon = shiny::icon("align-left")),
                shinydashboard::menuItem("GO Molecular Function",
                                        tabName = "gomf",
                                        icon = shiny::icon("project-diagram")),
                shinydashboard::menuItem("GO Biological Process",
                                        tabName = "gobp",
                                        icon = shiny::icon("project-diagram")),
                shinydashboard::menuItem("History", tabName = "hist",
                                                icon = shiny::icon("history")),
                shinydashboard::menuItem("Elements not in pathways",
                                        tabName = "ncments",
                                        icon = shiny::icon("window-close")),
                shiny::fluidRow(
                    shiny::column(3, align="center", offset=0,
                                        shiny::actionButton("reset", "Reset"))
                ),
                shiny::fluidRow(
                    shiny::column(12,offset=0,shiny::verbatimTextOutput('aff'))
                )
            )
        ),
        shinydashboard::dashboardBody(
            shiny::tags$head(shiny::tags$style(shiny::HTML('
                /* logo */
                .skin-blue .main-header .logo {
                    background-color:#48c9b0 ;
                }
                /* logo when hovered */
                .skin-blue .main-header .logo:hover {
                    background-color:#48c9b0 ;
                }

                /* navbar (rest of the header) */
                .skin-blue .main-header .navbar {
                    background-color:#48c9b0 ;
                }
                /* body */
                .content-wrapper, .right-side {
                    background-color: #fdfefe;
                }
            '))),
            shinydashboard::tabItems(
                shinydashboard::tabItem(tabName = "elements",
                    shiny::div(
                        style="left:240px; right:0px; bottom=0px; top=450px;
                                position:fixed; cursor:inherit; z-index: 2;",
                        shinyBS::bsCollapse(id = "collapseEx",
                            shinyBS::bsCollapsePanel("Filters",
                                shiny::fluidRow(
                                    shiny::column(8, align="center", offset=0,
                                        shiny::radioButtons("mode", "Mode :",
                                        choices=c("a AND b","a OR b","a NOT b"),
                                        selected="a AND b", inline=TRUE)),
                                    shiny::column(4, align="center", offset=0,
                                        shiny::actionButton("elemfilters",
                                        "Apply filters"))
                                ),
                                shiny::fluidRow(
                                    shiny::verbatimTextOutput('walk')
                                ), style = "info"
                            )
                        )
                    ),
                    shiny::fluidRow(
                        shiny::verbatimTextOutput("secblank")
                    ),
                    shiny::fluidRow(
                        shinydashboard::box(title="Genes",
                                DT::dataTableOutput('x2'),width=5),
                        shinydashboard::box(title="Metabolites",
                                DT::dataTableOutput('x3'),width=3),
                        shinydashboard::box(title="Interactions",
                                DT::dataTableOutput('x4'), width=4)
                    )
                ),
                shinydashboard::tabItem(tabName = "pathways",
                    shiny::div(
                        style="left:240px; right:0px; top=10px; position:fixed;
                        cursor:inherit; z-index: 2;",
                        shinyBS::bsCollapse(id = "collapseExample",
                            shinyBS::bsCollapsePanel("Filters",
                                shiny::fluidRow(
                                    shiny::column(4, align="center", offset=0,
                                    shiny::selectInput("pathtype",
                                        "Categories : ",
                                        choices = c("all", unique(types$root)),
                                        selected="all", multiple=FALSE)),
                                    shiny::column(4, align="center", offset=0,
                                        shiny::selectInput("typegene",
                                            "Genes types : ",choices = c("all",
                                            names(genetype)), selected="all",
                                            multiple=TRUE)),
                                    shiny::column(4, align="center", offset=0,
                                        shiny::selectInput("gotype",
                                            "Enriched GO terms : ",
                                            choices = c("all", gotermsgene),
                                            selected="all", multiple=TRUE))
                                ),
                                shiny::fluidRow(
                                    shiny::column(4, align="center", offset=0,
                                        shiny::selectInput("intetype",
                                            "Interactions types : ",
                                            choices = c("all", "g/g",
                                                "g/m", "m/m"), selected="all",
                                            multiple=TRUE)),
                                    shiny::column(4, align="center", offset=0,
                                        shiny::actionButton("pathfilters",
                                            "Apply filters")),
                                    shiny::column(4, align="center", offset=0,
                                        shiny::radioButtons("view", "View :",
                                            choices=c("default", "centrality",
                                                "interactions"),
                                            selected="default", inline=FALSE))
                                ), style = "info"
                            )
                        )
                    ),
                    shiny::fluidRow(
                        shiny::verbatimTextOutput("firstblank")
                    ),
                    shiny::fluidRow(
                        shinydashboard::box(title="Pathways",
                            shiny::div(style = 'overflow-x: scroll',
                                plotly::plotlyOutput("x1", height = "100%")),
                            width="100%", height = "100%")
                    )
                ),
                shinydashboard::tabItem(tabName = "gomf",
                    shiny::fluidRow(
                        shinydashboard::box(title="GO Molecular Function",
                                            DT::dataTableOutput('x5'))
                    )
                ),
                shinydashboard::tabItem(tabName = "gobp",
                    shiny::fluidRow(
                        shinydashboard::box(title="GO Biological Process",
                                            DT::dataTableOutput('x6'))
                    )
                ),
                shinydashboard::tabItem(tabName = "hist",
                    shiny::fluidRow(
                        shinydashboard::box(title="History",
                                            DT::dataTableOutput('x7'))
                    )
                ),
                shinydashboard::tabItem(tabName = "ncments",
                    shiny::fluidRow(
                        shinydashboard::box(title="Elements not in pathways",
                                            DT::dataTableOutput('x8'))
                    )
                )
            )
        )
    )
    return(ui)
}

rshiny=function(listdata){
    heatmap<-listdata[[1]];meta_list<-listdata[[2]];allResBP<-listdata[[3]];
    go_genelist<-listdata[[4]];allResMF<-listdata[[5]];types<-listdata[[6]];
    genetype<-listdata[[7]];metatab<-listdata[[8]];genetab<-listdata[[9]];
    intetab<-listdata[[10]];gomf_tab<-listdata[[11]];gobp_tab<-listdata[[12]];
    genes<-listdata[[13]];gomflist<-listdata[[14]];gobplist<-listdata[[15]];
    namegeneid<-listdata[[16]];hierabrite<-listdata[[17]];
    hierapath<-listdata[[18]];save_cluster_elem<-listdata[[19]]
    centrality<-listdata[[20]];inter_values<-listdata[[21]]
    gene_notin<-listdata[[22]];sub<-listdata[[23]]

    gotermsgene<-c(allResBP[allResBP[, 1] %in% go_genelist$go_id, 2],
                    allResMF[allResMF[, 1] %in% go_genelist$go_id, 2])
    server<-function(input, output, session) {
        v<-shiny::reactiveValues(selecgo=gotermsgene, info_bubble="",
                                selec_genetype=names(genetype),rows=vector(),
                                selec_pathtype=unique(types$root),
                                selec_goterm=vector(), a_not_b=vector(),
                                selec_inter=c("g/g","g/m","m/m"),
                                column=vector(), select_intertype=vector(),
                                heatmap_shiny=heatmap, metatab_shiny=metatab,
                                genetab_shiny=genetab, intetab_shiny=intetab,
                                x1_selected_rows=seq_len(nrow(heatmap)),
                                x2_selected_rows=seq_len(nrow(genetab)),
                                x3_selected_rows=seq_len(nrow(metatab)),
                                x4_selected_rows=seq_len(nrow(intetab)),
                                walk=vector(),mode="a AND b",history=list(),
                                gomf=gomf_tab, gobp=gobp_tab,
                                elements=vector(),suspended = TRUE,
                                histo_tab=data.frame(histo=character()))

        #reset button action
        shiny::observeEvent(input$reset, {
            v$x1_selected_rows <- NULL
            v$x2_selected_rows <- NULL
            v$x3_selected_rows <- NULL
            v$x4_selected_rows <- NULL
            v$selecgo<-gotermsgene
            v$selec_genetype<-names(genetype)
            v$selec_pathtype<-unique(types$root)
            v$selec_inter<-c("g/g", "g/m", "m/m")
            v$heatmap_shiny<-heatmap
            v$column<-v$selec_goterm<-vector()
            v$select_intertype<-v$rows<-vector()
            v$metatab_shiny<-metatab
            v$genetab_shiny<-genetab
            v$intetab_shiny<-intetab
            v$gomf<-gomf_tab
            v$gobp<-gobp_tab
            v$elements<-vector()
            v$history<-list()
            v$histo_tab<-data.frame(histo=character())
            v$mode<-"a AND b"
            v$info_bubble<-""
            v$walk<-vector()
            v$a_not_b<-vector()

            shiny::updateSelectInput(session, "typegene", "Genes types : ",
                        choices = c("all", v$selec_genetype), selected="all")
            shiny::updateSelectInput(session, "pathtype", "Categories : ",
                        choices = c("all", v$selec_pathtype), selected="all")
            shiny::updateSelectInput(session,"gotype", "Enriched GO terms : ",
                        choices = c("all", v$selecgo), selected="all")
            shiny::updateSelectInput(session, "intetype",
                "Interactions types : ", choices = c("all", v$selec_inter),
                selected="all")
        })

        #click on heatmap cell : pop up with pathway informations
        observer <- shiny::observeEvent(plotly::event_data("plotly_click",
                                        source = "x1"), suspended = TRUE, {
            v$x1_selected_rows <- plotly::event_data("plotly_click",
                                                        source = "x1")[["y"]]
            v$x1_selected_rows<-rev(c(seq_len(nrow(
                                        v$heatmap_shiny))))[v$x1_selected_rows]
            v$x1_selected_rows<-v$heatmap_shiny[v$x1_selected_rows, 1]
            if(v$x1_selected_rows!=""){
                sd <- which(heatmap$un %in% v$x1_selected_rows)
                popgene<-genetab[genetab[, 1] %in%rm_vector(
                            colnames(heatmap[sd, ])[which(heatmap[sd, ] == 1,
                                                arr.ind = TRUE)[, "col"]]), ]
                popmeta<-metatab[metatab[, 1] %in%rm_vector(
                            colnames(heatmap[sd, ])[which(heatmap[sd, ] == 1,
                                                arr.ind = TRUE)[, "col"]]), ]
                popinte<-intetab[intetab[, 1] %in%rm_vector(
                            colnames(heatmap[sd, ])[which(heatmap[sd, ] == 1,
                                                arr.ind = TRUE)[, "col"]]), ]

                if(stringr::str_sub(heatmap[sd, 2], 1, 1) == "R"){
                    href<-paste("https://reactome.org/PathwayBrowser/#/",
                                                    heatmap[sd, 2], sep="")
                    url<-a("Visualize pathway", href=href)
                }
                else if(stringr::str_sub(heatmap[sd, 2], 1, 1) == "h"){
                    link<-
            "https://www.genome.jp/kegg-bin/show_pathway?orgs=hsa eco&mapno="
                    href<-paste(link, stringr::str_sub(heatmap[sd, 2], 5,
                                            nchar(heatmap[sd, 2])), sep="")
                    url<-a("Visualize pathway", href=href)
                }
                else if(stringr::str_sub(heatmap[sd, 2], 1, 1) == "W"){
                    link<-"https://www.wikipathways.org/index.php/Pathway:"
                    href<-paste(link, heatmap[sd, 2], sep="")
                    url<-a("Visualize pathway", href=href)
                }

                shiny::showModal(shiny::modalDialog(
                    shiny::fluidRow(
                        shiny::h3(heatmap[sd, 1]),
                        shiny::h3(shiny::tagList(paste(heatmap[sd, 2], " : ",
                                                                sep=""), url)),
                        shinydashboard::box(title="Genes", DT::renderDataTable(
                            DT::datatable(popgene, rownames = FALSE,
                                selection ='none', options = list(pageLength =
                                                    nrow(popgene), dom = 't'))
                        ), width="100%"),
                        shinydashboard::box(title="Metabolites",
                                                        DT::renderDataTable(
                            DT::datatable(popmeta, rownames = FALSE,
                                selection ='none', options = list(pageLength =
                                                    nrow(popmeta), dom = 't'))
                        ), width="100%"),
                        shinydashboard::box(title="Interactions",
                                                        DT::renderDataTable(
                            DT::datatable(popinte, rownames = FALSE,
                                selection ='none', options = list(pageLength =
                                                    nrow(popinte), dom = 't'))
                        ), width="100%")
                    )
                ))
            }
        })

        #click on go terms dataframes : pop up with informations
        #about the go terms in th hierarchy node selected
        shiny::observeEvent(input$x5_rows_selected | input$x6_rows_selected, {
            goterms<-vector()
            if(length(v$gomf[input$x5_rows_selected, 1])>0){
                goterms<-stringr::str_split(v$gomf[input$x5_rows_selected, 1],
                                                                    " / ")[[1]]
            }
            else if(length(v$gobp[input$x6_rows_selected, 1])>0){
                goterms<-stringr::str_split(v$gobp[input$x6_rows_selected, 1],
                                                                    " / ")[[1]]
            }
            if(length(goterms)>0){
                if(length(goterms)>1){
                    goterms<-c(goterms[2:length(goterms)], stringr::str_split(
                                                        goterms[1], "__")[[1]])
                }
                else{goterms<-stringr::str_split(goterms, "__")[[1]]}
                goterms<-goterms[!(goterms %in% "")]
            }

            #genes related to go terms + pop up
            geneterms<-go_genelist[go_genelist[, 2] %in% goterms, 1]
            if(length(geneterms)>0){
                tabgeneterms<-namegeneid[namegeneid[, 1] %in% geneterms, ]
                shiny::showModal(shiny::modalDialog(
                    shinydashboard::box(title="Genes", DT::renderDataTable(
                        DT::datatable(tabgeneterms, rownames = FALSE,
                            selection = 'none', options = list(pageLength =
                                            nrow(tabgeneterms), dom = 't'))
                    ), width="100%")
                ))
            }
        }, ignoreNULL=TRUE)

        #click on history dataframe
        shiny::observeEvent(input$x7_rows_selected, {
            r<-input$x7_rows_selected
            v$x2_selected_rows <- NULL
            v$x3_selected_rows <- NULL
            v$x4_selected_rows <- NULL
            v$column<-v$walk<-v$a_not_b<-v$rows<-v$elements<-vector()
            v$heatmap_shiny<-v$history[[r]][["path"]]
            v$gomf<-v$history[[r]][["mf"]]
            v$gobp<-v$history[[r]][["bp"]]
            v$genetab_shiny<-v$history[[r]][["tg"]]
            v$metatab_shiny<-v$history[[r]][["tm"]]
            v$intetab_shiny<-v$history[[r]][["ti"]]
            v$selecgo<-v$history[[r]][["sgo"]]
            v$selec_genetype<-v$history[[r]][["segene"]]
            v$selec_pathtype<-v$history[[r]][["spath"]]
            v$selec_inter<-v$history[[r]][["sinter"]]

            shiny::updateSelectInput(session, "pathtype", "Categories : ",
                        choices=c("all", v$selec_pathtype), selected="all")
            shiny::updateSelectInput(session, "typegene", "Genes types : ",
                        choices=c("all", v$selec_genetype), selected="all")
            shiny::updateSelectInput(session,"gotype","Enriched GO terms : ",
                                choices = c("all", v$selecgo), selected="all")
            shiny::updateSelectInput(session, "intetype",
                    "Interactions types : ", choices = c("all", v$selec_inter),
                                                                selected="all")
        })

        #click on filter pathway button : filter the heatmap rows and columns
        shiny::observeEvent(input$pathfilters, {
            #elements
            coltemp<-colnames(v$heatmap_shiny)[colnames(v$heatmap_shiny)
                                                                    %in% genes]
            v$column<-rm_vector(c(v$column, v$selec_goterm))
            if(length(v$column) == 0){
                v$column<-colnames(v$heatmap_shiny[, 6:ncol(v$heatmap_shiny)])
            }
            colselec<-v$column
            colnotselec<-coltemp[!(coltemp %in% v$column)]
            v$column<-c(v$column[v$column %in% coltemp],
                        colnames(v$heatmap_shiny)[!(colnames(v$heatmap_shiny)
                                                                %in% genes)])

            #interactions
            noninte<-intetab[intetab[, 2] %in% colnotselec, 4]
            noninte<-c(noninte, intetab[intetab[, 4] %in% colnotselec, 4])
            v$column<-v$column[!(v$column %in% noninte)]

            #pathways
            v$rows<-rm_vector(c(v$rows, v$select_intertype))
            if(length(v$rows) == 0){
                v$rows<-row.names(v$heatmap_shiny)
            }

            v$heatmap_shiny<-v$heatmap_shiny[row.names(v$heatmap_shiny)[
                row.names(v$heatmap_shiny) %in% as.character(v$rows)],
                colnames(v$heatmap_shiny)[colnames(v$heatmap_shiny)
                                                            %in% v$column]]

            #update elements dataframes with elements
            #in the final heatmap pathways
            v$x2_selected_rows <- NULL
            v$x3_selected_rows <- NULL
            v$x4_selected_rows <- NULL
            v$elements<-vector()

            final_elements<-colnames(v$heatmap_shiny[,6:ncol(v$heatmap_shiny)])
            selected_genes<-final_elements[final_elements %in% genes]
            selected_meta<-final_elements[final_elements %in% meta_list]
            selected_inter<-final_elements[final_elements %in% intetab$tag]
            v$genetab_shiny<-genetab[genetab[, 1] %in% selected_genes, ]
            v$metatab_shiny<-metatab[metatab[, 1] %in% selected_meta, ]
            v$intetab_shiny<-intetab[intetab[, 1] %in% selected_inter, ]

            #histo
            if(length(v$history) == 10){
                v$history<-v$history[-1]
                for(d in seq_len(nrow(v$histo_tab)-1)){
                    v$histo_tab[d, 1]<-v$histo_tab[d+1, 1]
                }
                v$histo_tab<-v$histo_tab[-10, ]
            }

            v$history[[length(v$history)+1]]<-list(path=v$heatmap_shiny,
                mf=v$gomf, bp=v$gobp, tg=v$genetab_shiny, tm=v$metatab_shiny,
                ti=v$intetab_shiny, sgo=v$selecgo, spath=v$selec_pathtype,
                sinter=v$selec_inter, segene=v$selec_genetype)
            v$histo_tab[length(v$history), 1]<-paste("Updated : ",
                                paste(input$pathtype, collapse=" "), sep="")

            v$rows<-v$column<-vector()

            #update available choices on filters
            hm_elements<-names(v$heatmap_shiny)
            v$selec_pathtype<-v$selec_genetype<-vector()
            v$selecgo<-v$selec_inter<-vector()

            pathselected<-v$heatmap_shiny[row.names(v$heatmap_shiny), 2]
            for(h in seq_len(length(hierabrite))){ #path
                if(length(hierabrite[[h]][["name"]][hierabrite[[h]][["name"]]
                                                        %in% pathselected])>0){
                    v$selec_pathtype<-rm_vector(c(v$selec_pathtype,
                                                        names(hierabrite)[h]))
                }
            }
            gomfid<-vector()
            gobpid<-vector()
            #types d'interactions, types genes et go
            for(i in seq_len(length(hierapath))){
                if(length(hierapath[[i]][["name"]][hierapath[[i]][["name"]]
                                                        %in% pathselected])>0){
                    hm_genes<-hm_elements[hm_elements %in%
                            hierapath[[i]][["elem"]][hierapath[[i]][["elem"]]
                                                                %in% genes]]
                    hm_inter<-hm_elements[hm_elements %in%
                            hierapath[[i]][["elem"]][hierapath[[i]][["elem"]]
                                                            %in% intetab[, 1]]]

                    if(length(hm_genes)>0){
                        for(g in seq_len(length(genetype))){
                            if(length(genetype[[g]][genetype[[g]] %in%
                                                                hm_genes])>0){
                                v$selec_genetype<-rm_vector(c(v$selec_genetype,
                                                        names(genetype)[g]))
                            }
                        }
                        if(length(hm_genes[hm_genes %in%
                                                go_genelist$hgnc_symbol] )>0){
                            goenr<-go_genelist[go_genelist$hgnc_symbol %in%
                                                    hm_genes[hm_genes %in%
                                                    go_genelist$hgnc_symbol],2]
                            if(length(goenr[goenr %in% allResBP[, 1]])>0){
                                goenrbp<-goenr[goenr %in% allResBP[, 1]]
                                for(b in seq_len(length(gobplist))){
                                    if (length(gobplist[[b]][["goterm"]]
                                                [gobplist[[b]][["goterm"]]
                                                            %in% goenrbp])>0){
                                        gobpid<-c(gobpid,
                                                    gobplist[[b]][["index"]])
                                    }
                                }
                                v$selecgo<-rm_vector(c(v$selecgo,
                                        allResBP[allResBP[,1] %in% goenrbp,2]))
                            }
                            else if(length(goenr[goenr %in% allResMF[, 1]])>0){
                                goenrmf<-goenr[goenr %in% allResMF[, 1]]
                                for(m in seq_len(length(gomflist))){
                                    if (length(gomflist[[m]][["goterm"]][
                                                    gomflist[[m]][["goterm"]]
                                                            %in% goenrmf])>0){
                                        gomfid<-c(gomfid,
                                                    gomflist[[m]][["index"]] )
                                    }
                                }
                                v$selecgo<-rm_vector(c(v$selecgo,
                                    allResMF[allResMF[, 1] %in% goenrmf, 2]))
                            }
                        }
                    }
                    if (length(hm_inter)>0){
                        v$selec_inter<-rm_vector(c(v$selec_inter,
                                    intetab[intetab[, 1] %in% hm_inter, 6]))
                    }
                }
            }
            gomfid<-gomfid[gomfid %in% as.integer(row.names(v$gomf))]
            gobpid<-gobpid[gobpid %in% as.integer(row.names(v$gobp))]
            v$gomf<-v$gomf[as.character(gomfid), ]
            v$gobp<-v$gobp[as.character(gobpid), ]

            shiny::updateSelectInput(session, "gotype", "Enriched GO terms : ",
                                choices = c("all", v$selecgo), selected="all")
            shiny::updateSelectInput(session, "typegene", "Genes types : ",
                        choices = c("all", v$selec_genetype), selected="all")
            shiny::updateSelectInput(session, "pathtype", "Categories : ",
                        choices = c("all", v$selec_pathtype), selected="all")
            shiny::updateSelectInput(session, "intetype",
                    "Interactions types : ", choices = c("all", v$selec_inter),
                                                                selected="all")
        })

        #pathways categories selected
        shiny::observeEvent(input$pathtype, {
            root<-input$pathtype
            if(!("all" %in% root)){
                v$rows<-hierabrite[[root[1]]][["index"]]
            }
            else{
                v$rows<-vector()
            }
        }, ignoreNULL=FALSE)

        #interactions categories selected
        shiny::observeEvent(input$intetype, {
            selected<-input$intetype
            if(length(selected) == 0){
                shiny::updateSelectInput(session, "intetype",
                    "Interactions types : ", choices = c("all", v$selec_inter),
                                                                selected="all")
                selected<-"all"
            }
            if(!("all" %in% selected)){
                v$select_intertype<-vector()
                selected_inter<-rm_vector(intetab[intetab[,7] %in% selected,4])
                for (h in seq_len(length(hierapath))){
                    if(length(hierapath[[h]][["elem"]][hierapath[[h]][["elem"]]
                                                    %in% selected_inter])>0){
                        v$select_intertype<-c(v$select_intertype,
                                                    hierapath[[h]][["index"]])
                    }
                }
            }
            else{
                v$select_intertype<-vector()
            }
        }, ignoreNULL=FALSE)

        #genes categories selected
        shiny::observeEvent(input$typegene, {
            selected<-input$typegene
            if(length(selected) == 0){
                shiny::updateSelectInput(session, "typegene", "Genes types : ",
                        choices = c("all", v$selec_genetype), selected="all")
                selected="all"
            }
            if(!("all" %in% selected)){
                v$column<-vector()
                for (s in seq_len(length(selected))){
                    v$column<-c(v$column, genetype[[selected[s]]])
                }
                v$column<-rm_vector(v$column)
            }
            else{
                v$column<-vector()
            }
        }, ignoreNULL=FALSE)

        #go terms selected
        shiny::observeEvent(input$gotype, {
            selected<-input$gotype
            if(length(selected) == 0){
                shiny::updateSelectInput(session, "gotype",
                    "Enriched GO terms : ", choices = c("all", v$selecgo),
                                                            selected="all")
                selected="all"
            }
            if(!("all" %in% selected)){
                goidterms<-c(allResBP[allResBP[, 2] %in% selected, 1],
                                    allResMF[allResMF[, 2] %in% selected, 1])
                v$selec_goterm<-rm_vector(go_genelist[go_genelist$go_id %in%
                                                                goidterms, 1])
            }
            else{
                v$selec_goterm<-vector()
            }
        }, ignoreNULL=FALSE)

        #find pathways containing elements selected or not
        shiny::observeEvent(input$elemfilters, {
            if(length(v$walk)>0){
                #read walk to determine pathways
                prev_elem_path<-vector()
                a<-1
                while(a<=length(v$walk)){
                    if(a == 1){
                        prev_elem_path<-rm_vector(as.vector(v$heatmap_shiny[
                            v$heatmap_shiny[, v$walk[a]] %in% c(1), "deux"]))
                        if(length(v$walk) == 1){
                            elem_path<-prev_elem_path
                        }
                        a<-a+1
                    }

                    else if(length(v$walk)!=1){
                        elem_path<-rm_vector(as.vector(v$heatmap_shiny[
                            v$heatmap_shiny[, v$walk[a+1]] %in% c(1), "deux"]))
                        if(v$walk[a] == "&"){
                            elem_path<-intersect(elem_path, prev_elem_path)
                        }
                        else if(v$walk[a] == "|"){
                            elem_path<-rm_vector(c(elem_path, prev_elem_path))
                        }
                        else{
                            elem_path<-dplyr::setdiff(prev_elem_path,elem_path)
                        }
                        prev_elem_path<-elem_path
                        a<-a+2
                    }
                }

                #use pathways to find which hierarchies to show,
                #and which elements
                index<-vector()
                hm_elements<-vector()
                for (b in seq_len(length(hierapath))){
                    if(length(which(names(table(elem_path %in%
                                    hierapath[[b]][["name"]])) == TRUE))>0){
                        hiera_elements<-hierapath[[b]][["elem"]]
                        if(length(hiera_elements[hiera_elements %in%
                                                        v$a_not_b]) == 0){
                            hm_elements<-c(hm_elements,
                                                    hierapath[[b]][["elem"]])
                            index<-c(index, hierapath[[b]][["index"]])
                        }
                    }
                }
                index<-sort(index)
                hm_elements<-rm_vector(hm_elements)

                if(length(which(colnames(v$heatmap_shiny)%in%hm_elements))>0){
                    v$heatmap_shiny<-v$heatmap_shiny[c(as.character(index)),
                            c(1, 2, 3, 4, 5, which(colnames(v$heatmap_shiny)
                                                            %in% hm_elements))]
                }
                else{ #si aucun match entre les colonnes
                    v$heatmap_shiny<-v$heatmap_shiny[-c(seq_len(
                                                    nrow(v$heatmap_shiny))), ]
                }

                v$rows<-v$column<-vector()

                #GOTERMS
                genes_in_walk<-v$walk[v$walk %in% genes]
                walk_only_genes<-vector()
                if (length(genes_in_walk)>0){
                    wi<-which(v$walk %in% genes_in_walk)
                    walk_only_genes<-c(walk_only_genes, v$walk[wi[1]])
                    if (length(wi)>1){
                        for (w in 2:length(wi)){
                            walk_only_genes<-c(walk_only_genes,
                                            v$walk[wi[w]-1], v$walk[wi[w]])
                        }
                    }

                    g<-1
                    while (g<=length(walk_only_genes)){
                        if(g == 1){
                            prev_goterms_walk<-rm_vector(as.vector(go_genelist[
                                                    go_genelist$hgnc_symbol ==
                                                        walk_only_genes[g],2]))
                            prev_walk_gomf<-prev_goterms_walk[prev_goterms_walk
                                                            %in% allResMF[, 1]]
                            prev_walk_gobp<-prev_goterms_walk[prev_goterms_walk
                                                            %in% allResBP[, 1]]
                            g<-g+1
                        }

                        else if(length(walk_only_genes)!=1){
                            goterms_walk<-rm_vector(as.vector(go_genelist[
                                                    go_genelist$hgnc_symbol ==
                                                    walk_only_genes[g+1], 2]))
                            walk_gomf<-goterms_walk[goterms_walk %in%
                                                                allResMF[, 1]]
                            walk_gobp<-goterms_walk[goterms_walk %in%
                                                                allResBP[, 1]]
                            if(walk_only_genes[g] == "&"){
                                walk_gomf<-intersect(walk_gomf, prev_walk_gomf)
                                walk_gobp<-intersect(walk_gobp, prev_walk_gobp)
                            }
                            else if(walk_only_genes[g] == "|"){
                                walk_gomf<-rm_vector(c(walk_gomf,
                                                            prev_walk_gomf))
                                walk_gobp<-rm_vector(c(walk_gobp,
                                                            prev_walk_gobp))
                            }
                            else{
                                walk_gomf<-dplyr::setdiff(prev_walk_gomf,
                                                                    walk_gomf)
                                walk_gobp<-dplyr::setdiff(prev_walk_gobp,
                                                                    walk_gobp)
                            }
                            prev_walk_gomf<-walk_gomf
                            prev_walk_gobp<-walk_gobp
                            g<-g+2
                        }
                    }

                    hiera_gomf<-vector()
                    for (m in seq_len(length(gomflist))){
                        if(length(which(names(table(prev_walk_gomf %in%
                                    gomflist[[m]][["goterm"]])) == TRUE))>0){
                            hierago_gene<-gomflist[[m]][["gene"]]
                            if(length(hierago_gene[hierago_gene %in%
                                                            v$a_not_b]) == 0){
                                hiera_gomf<-c(hiera_gomf,
                                                    gomflist[[m]][["index"]])
                            }
                        }
                    }
                    hiera_gomf<-sort(hiera_gomf)
                    hiera_gomf<-hiera_gomf[hiera_gomf %in%
                                                as.integer(row.names(v$gomf))]
                    v$gomf<-v$gomf[as.character(hiera_gomf), ]

                    hiera_gobp<-vector()
                    for (b in seq_len(length(gobplist))){
                        if(length(which(names(table(prev_walk_gobp %in%
                                    gobplist[[b]][["goterm"]])) == TRUE))>0){
                            hierago_gene<-gobplist[[b]][["gene"]]
                            if(length(hierago_gene[hierago_gene %in%
                                                            v$a_not_b]) == 0){
                                hiera_gobp<-c(hiera_gobp,
                                                    gobplist[[b]][["index"]])
                            }
                        }
                    }
                    hiera_gobp<-sort(hiera_gobp)
                    hiera_gobp<-hiera_gobp[hiera_gobp %in%
                                                as.integer(row.names(v$gobp))]
                    v$gobp<-v$gobp[as.character(hiera_gobp), ]
                }
            }
            else if (length(v$a_not_b)>0 && length(v$walk) == 0){
                #pathways NOT
                elem_path<-rm_vector(v$heatmap_shiny[!(v$heatmap_shiny[, 2]
                                                                %in% c(0)), 2])

                index<-vector()
                hm_elements<-vector()
                for (b in seq_len(length(hierapath))){
                    if(length(which(names(table(elem_path %in%
                                    hierapath[[b]][["name"]])) == TRUE))>0){
                        hiera_elements<-hierapath[[b]][["elem"]]
                        if(length(hiera_elements[hiera_elements %in%
                                                        v$a_not_b]) == 0){
                            hm_elements<-c(hm_elements,
                                                    hierapath[[b]][["elem"]])
                            index<-c(index, hierapath[[b]][["index"]])
                        }
                    }
                }
                index<-sort(index)
                hm_elements<-rm_vector(hm_elements)

                v$heatmap_shiny<-v$heatmap_shiny[c(as.character(index)),
                            c(1, 2, 3, 4, 5, which(colnames(v$heatmap_shiny)
                                                            %in% hm_elements))]

                #GO NOT
                prev_walk_gomf<-allResMF[, 1]
                prev_walk_gobp<-allResBP[, 1]

                hiera_gomf<-vector()
                for (m in seq_len(length(gomflist))){
                    if(length(which(names(table(prev_walk_gomf %in%
                                    gomflist[[m]][["goterm"]])) == TRUE))>0){
                        hierago_gene<-gomflist[[m]][["gene"]]
                        if(length(hierago_gene[hierago_gene %in%
                                                            v$a_not_b]) == 0){
                            hiera_gomf<-c(hiera_gomf, gomflist[[m]][["index"]])
                        }
                    }
                }
                hiera_gomf<-sort(hiera_gomf)
                hiera_gomf<-hiera_gomf[hiera_gomf %in%
                                                as.integer(row.names(v$gomf))]
                v$gomf<-v$gomf[as.character(hiera_gomf), ]

                hiera_gobp<-vector()
                for (b in seq_len(length(gobplist))){
                    if(length(which(names(table(prev_walk_gobp %in%
                                    gobplist[[b]][["goterm"]])) == TRUE))>0){
                        hierago_gene<-gobplist[[b]][["gene"]]
                        if(length(hierago_gene[hierago_gene %in% v$a_not_b])
                                                                        == 0){
                            hiera_gobp<-c(hiera_gobp, gobplist[[b]][["index"]])
                        }
                    }
                }
                hiera_gobp<-sort(hiera_gobp)
                hiera_gobp<-hiera_gobp[hiera_gobp %in%
                                                as.integer(row.names(v$gobp))]
                v$gobp<-v$gobp[as.character(hiera_gobp), ]
            }

            #elements des hiera
            v$x2_selected_rows <- NULL
            v$x3_selected_rows <- NULL
            v$x4_selected_rows <- NULL
            v$elements<-vector()

            final_elements<-colnames(v$heatmap_shiny[,6:ncol(v$heatmap_shiny)])
            selected_genes<-final_elements[final_elements %in% genes]
            selected_meta<-final_elements[final_elements %in% meta_list]
            selected_inter<-final_elements[final_elements %in% intetab$tag]
            v$genetab_shiny<-genetab[genetab[, 1] %in% selected_genes, ]
            v$metatab_shiny<-metatab[metatab[, 1] %in% selected_meta, ]
            v$intetab_shiny<-intetab[intetab[, 1] %in% selected_inter, ]

            #histo
            if(length(v$history) == 10){
                v$history<-v$history[-1]
                for(d in seq_len(nrow(v$histo_tab)-1)){
                    v$histo_tab[d, 1]<-v$histo_tab[d+1, 1]
                }
                v$histo_tab<-v$histo_tab[-10, ]
            }

            v$history[[length(v$history)+1]]<-list(path=v$heatmap_shiny,
                mf=v$gomf, bp=v$gobp, tg=v$genetab_shiny, tm=v$metatab_shiny,
                ti=v$intetab_shiny, sgo=v$selecgo, spath=v$selec_pathtype,
                                sinter=v$selec_inter, segene=v$selec_genetype)
            v$histo_tab[length(v$history), 1]<-paste("Pathways focus : ",
                                        paste(v$walk, collapse=" "), sep="")
            v$walk<-v$a_not_b<-vector()

            #update available choices on filters
            hm_elements<-names(v$heatmap_shiny)
            v$selec_pathtype<-v$selec_genetype<-vector()
            v$selecgo<-v$selec_inter<-vector()

            pathselected<-v$heatmap_shiny[row.names(v$heatmap_shiny), 2]
            for(h in seq_len(length(hierabrite))){ #path
                if(length(hierabrite[[h]][["name"]][hierabrite[[h]][["name"]]
                                                        %in% pathselected])>0){
                    v$selec_pathtype<-rm_vector(c(v$selec_pathtype,
                                                        names(hierabrite)[h]))
                }
            }
            #types d'interactions, types genes et go
            for(i in seq_len(length(hierapath))){
                if(length(hierapath[[i]][["name"]][hierapath[[i]][["name"]]
                                                        %in% pathselected])>0){
                    hm_genes<-hm_elements[hm_elements %in%
                            hierapath[[i]][["elem"]][hierapath[[i]][["elem"]]
                                                                %in% genes]]
                    hm_inter<-hm_elements[hm_elements %in%
                            hierapath[[i]][["elem"]][hierapath[[i]][["elem"]]
                                                        %in% intetab[, 1]]]

                    if(length(hm_genes)>0){
                        for(g in seq_len(length(genetype))){
                            if(length(genetype[[g]][genetype[[g]] %in%
                                                                hm_genes])>0){
                                v$selec_genetype<-rm_vector(c(v$selec_genetype,
                                                        names(genetype)[g]))
                            }
                        }
                        if(length(hm_genes[hm_genes %in%
                                                go_genelist$hgnc_symbol])>0){
                            goenr<-go_genelist[go_genelist$hgnc_symbol %in%
                                                hm_genes[hm_genes %in%
                                                go_genelist$hgnc_symbol], 2]
                            if(length(goenr %in% allResBP[, 1])>0){
                                v$selecgo<-rm_vector(c(v$selecgo,
                                                allResBP[allResBP[, 1] %in%
                                                                goenr, 2]))
                            }
                            else if(length(goenr %in% allResMF[, 1])>0){
                                v$selecgo<-rm_vector(c(v$selecgo,
                                                allResMF[allResMF[, 1] %in%
                                                                goenr, 2]))
                            }
                        }
                    }
                    if (length(hm_inter)>0){
                        v$selec_inter<-rm_vector(c(v$selec_inter,
                                                    intetab[intetab[, 1] %in%
                                                                hm_inter, 6]))
                    }
                }
            }

            shiny::updateSelectInput(session, "gotype", "Enriched GO terms : ",
                                choices = c("all", v$selecgo), selected="all")
            shiny::updateSelectInput(session, "typegene", "Genes types : ",
                            choices=c("all", v$selec_genetype), selected="all")
            shiny::updateSelectInput(session, "pathtype", "Categories : ",
                            choices=c("all", v$selec_pathtype), selected="all")
            shiny::updateSelectInput(session, "intetype",
                        "Interactions types : ",choices=c("all",v$selec_inter),
                                                                selected="all")
        })

        #selected elements : build a walk
        shiny::observeEvent(input$x2_rows_selected | input$x3_rows_selected |
                                                    input$x4_rows_selected, {
            v$x2_selected_rows <- input$x2_rows_selected
            v$x3_selected_rows <- input$x3_rows_selected
            v$x4_selected_rows <- input$x4_rows_selected

            element<-c(v$genetab_shiny[v$x2_selected_rows, 1],
                        v$metatab_shiny[v$x3_selected_rows, 1],
                        v$intetab_shiny[v$x4_selected_rows, 1])

            add_elem<-dplyr::setdiff(element, v$elements)
            remove_elem<-dplyr::setdiff(v$elements, element)

            #traitement remove_elem sur walk
            if(length(remove_elem)>0 && length(v$a_not_b[v$a_not_b %in%
                                                        remove_elem])>0){
                v$a_not_b<-v$a_not_b[!(v$a_not_b %in% remove_elem)]
            }
            else if(length(remove_elem)>0 && length(v$walk) == 1){
                v$walk<-vector()
            }
            else if(length(remove_elem)>0 && which(v$walk %in% remove_elem)
                                                                        == 1){
                v$walk<-v$walk[-c(1, 2)]
            }
            else if(length(remove_elem)>0){
                i<-which(v$walk %in% remove_elem)
                v$walk<-v$walk[-c(i, i-1)]
            }

            if(length(element)>0){
                if(length(add_elem)>0){

                    if(length(v$walk) == 0){
                        if(input$mode == "a NOT b"){
                            v$a_not_b<-c(v$a_not_b, add_elem)
                        }
                        else{
                            v$walk<-c(v$walk, add_elem)
                        }
                    }
                    else{
                        if(input$mode == "a AND b"){
                            op<-"&"
                            v$walk<-c(v$walk, op, add_elem)
                        }
                        else if(input$mode == "a OR b"){
                            op<-"|"
                            v$walk<-c(v$walk, op, add_elem)
                        }
                        else{
                            v$a_not_b<-c(v$a_not_b, add_elem)
                        }
                    }
                }
            }
            v$elements<-element
        })

        #info bubble --> genes
        shiny::observeEvent(input$geneindex, {
            gi<-as.vector(unname(t(v$genetab_shiny[input$geneindex+1, ])))
            v$info_bubble<-paste(gi[1], "\nName : ", gi[2], sep="")
        })

        #info bubble --> meta
        shiny::observeEvent(input$metaindex, {
            mi<-as.vector(unname(t(v$metatab_shiny[input$metaindex+1, ])))
            v$info_bubble<-paste(mi[1], "\nName : ", mi[2], sep="")
        })

        #info bubble --> interactions
        shiny::observeEvent(input$inteindex, {
            ii<-as.vector(unname(t(v$intetab_shiny[input$inteindex+1, ])))
            v$info_bubble<-paste(ii[1], "\nCompo 1 : ", ii[2], "\nCompo 2 : ",
                                    ii[4], "\nGO term : ", ii[5], "\nPath : ",
                                    ii[6], "\nType : ", ii[7], "\nLien : ",
                                    paste(stringr::str_split(ii[3], ", ")[[1]],
                                                        collapse="\n"), sep="")
        })


        ##RENDER##
        output$x1<- plotly::renderPlotly({
            input$reset
            input$view

            y<-rev(v$heatmap_shiny$un)
            if(length(y)>0){
                tempor<-data.matrix(v$heatmap_shiny[,
                                    which(colnames(v$heatmap_shiny) %in%
                                                    save_cluster_elem)])
                if(input$view == "default"){
                    data<-tempor
                    gamme<-grDevices::colorRampPalette(c("ghostwhite","blue4"))
                }
                else if(input$view == "centrality"){
                    data<-centrality[as.integer(row.names(v$heatmap_shiny)),
                                    which(colnames(sub) %in% colnames(tempor))]
                    gamme<-grDevices::colorRampPalette(c("ghostwhite", "blue",
                                        "green", "gold", "darkorange" , "red"))
                }
                else if (input$view == "interactions"){
                    data<-inter_values[as.integer(row.names(v$heatmap_shiny)),
                                    which(colnames(sub) %in% colnames(tempor))]
                    gamme<-grDevices::colorRampPalette(c("ghostwhite", "blue",
                                                                "green","red"))
                }
                subtitles<-sub[as.integer(row.names(v$heatmap_shiny)),
                                which(colnames(sub) %in% colnames(tempor))]

                p <- plotly::plot_ly(
                    x=colnames(v$heatmap_shiny[, which(colnames(v$heatmap_shiny)
                    %in% save_cluster_elem)]), y = seq_along(y),
                    z = apply(data, 2, rev), type = "heatmap", source = "x1",
                    height=(30*nrow(v$heatmap_shiny)+100),
                    width=(6*ncol(v$heatmap_shiny)+1200), colors=gamme(100),
                    showscale=FALSE, hoverinfo='text',
                    text= apply(subtitles, 2, rev)
                ) %>% plotly::layout(yaxis = list(side = "right",
                                    tickvals = seq_along(y), ticktext = y),
                                    margin =list(l=0, r=1200, b=100, t=50))
            }
            else{ #no heatmap printed
                p <- plotly::plotly_empty(type = "scatter", mode = 'lines')
            }
            ## resume observer only if suspended
            if(v$suspended) {
                observer$resume()
                v$suspended <- FALSE
            }
            return(p)
        })
        output$x2 <- DT::renderDataTable({
            input$x7_rows_selected
            input$reset

            DT::datatable(v$genetab_shiny, rownames = FALSE,
                selection=list(mode ='multiple', selected = v$x2_selected_rows,
                target ='row'), options=list(pageLength=nrow(v$genetab_shiny),
                                                dom = 't', rowCallback=DT::JS(
                'function(row, data) {
                    $(row).mouseenter(function(){
                        var gene_index = $(this)[0]._DT_RowIndex
                        /* console.log(hover_index); */
                        Shiny.onInputChange("geneindex", gene_index);
                    });
                }'
            )))
        })
        output$x3 <- DT::renderDataTable({
            input$x7_rows_selected
            input$reset

            DT::datatable(v$metatab_shiny, rownames = FALSE,
                selection=list(mode ='multiple', selected=v$x3_selected_rows,
                target ='row'), options=list(pageLength=nrow(v$metatab_shiny),
                                                dom = 't', rowCallback=DT::JS(
                'function(row, data) {
                    $(row).mouseenter(function(){
                        var meta_index = $(this)[0]._DT_RowIndex
                        /* console.log(hover_index); */
                        Shiny.onInputChange("metaindex", meta_index);
                    });
                }'
            )))
        })
        output$x4 <- DT::renderDataTable({
            input$x7_rows_selected
            input$reset
            DT::datatable(v$intetab_shiny[, c(1, 5, 8, 7)], rownames = FALSE,
                selection =list(mode ='multiple', selected=v$x4_selected_rows,
                target ='row'), options=list(pageLength=nrow(v$intetab_shiny),
                                                            rowCallback=DT::JS(
                'function(row, data) {
                    $(row).mouseenter(function(){
                        var inte_index = $(this)[0]._DT_RowIndex
                        /* console.log(hover_index); */
                        Shiny.onInputChange("inteindex", inte_index);
                    });
                }'
            )))
        })
        output$x5 <- DT::renderDataTable({
            input$x5_rows_selected
            input$reset
            DT::datatable(v$gomf, rownames = FALSE, selection="single",
                    options = list(pageLength = nrow(v$gomf), dom = 't'))
        })
        output$x6 <- DT::renderDataTable({
            input$x6_rows_selected
            input$reset
            DT::datatable(v$gobp, rownames = FALSE, selection="single",
                    options = list(pageLength = nrow(v$gobp), dom = 't'))
        })
        output$x7 <- DT::renderDataTable({
            input$x7_rows_selected
            input$reset
            DT::datatable(v$histo_tab, rownames = FALSE, selection="single",
                    options = list(pageLength = nrow(v$histo_tab), dom = 't'))
        })
        output$x8 <- DT::renderDataTable({
            DT::datatable(gene_notin, rownames = FALSE, selection="none",
                    options = list(pageLength = nrow(gene_notin), dom = 't'))
        })
        output$aff<-shiny::renderText({
            v$info_bubble
        })
        output$walk<-shiny::renderText({
            input$x2_rows_selected
            input$x3_rows_selected
            paste("walk :", paste(v$walk, collapse=" "), "\nexcluded :",
                                    paste(v$a_not_b, collapse=" "), sep=" ")
        })
        output$firstblank<-shiny::renderText({
            paste(rep("\n", 2), collapse="")
        })
        output$secblank<-shiny::renderText({
            paste(rep("\n", 2), collapse="")
        })
    }

    shiny::shinyApp(ui = ui(types, genetype, gotermsgene), server = server)
}

