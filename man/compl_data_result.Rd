\name{compl_data_result}
\alias{compl_data_result}
\docType{data}
\title{
    Output of \code{compl_data} function
}
\description{
    Enriched and structured informations about genes, metabolites, their
    interactions, pathways and enriched GO terms. All these informations
    are then used by the Shiny interface (\code{rshiny} function). GO
    terms enrichment analysis is performed using clusterProfiler.
}
\usage{data("compl_data_result")}
\format{
    List of 25 items.
    \describe{
        \item{\code{heatmap}}{
            Dataframe heatmap-like, with in abscissa elements of pathways
            ("X" is written if an element is present in a pathway), and with
            in ordinate hierarchies of pathways
        }
        \item{\code{meta_list}}{
            User's metabolites given in \code{path_enrich} function
        }
        \item{\code{allResBP}}{
            Results of Go BP terms enrichment analysis performed by
            clusterProfileR (20 best)
        }
        \item{\code{go_genelist}}{
            Dataframe containing enriched GO terms per genes of user's list
        }
        \item{\code{allResMF}}{
            Results of Go MF terms enrichment analysis performed by
            clusterProfileR (20 best)
        }
        \item{\code{types}}{
            Dataframe containing categories of pathways per pathways
        }
        \item{\code{genetype}}{
            List of genes categories, containing genes of user's list
        }
        \item{\code{metatab}}{
            Dataframe of metabolites and informations about them, as
            names and chebi ids
        }
        \item{\code{genetab}}{
            Dataframe of genes and informations about them, as gene
            symbols and names
        }
        \item{\code{intetab}}{
            Dataframe of direct interactions between elements of user's list,
            and informations about them as elements in the interaction, how
            they interact, and which pathways are concerned
        }
        \item{\code{gomf_tab}}{
            Dataframe of Go MF terms hierarchies containing our enriched GO
            terms, plus description of these GO terms and genes of user's
            list concerned by enriched GO terms
        }
        \item{\code{gobp_tab}}{
            Dataframe of Go BP terms hierarchies containing our enriched GO
            terms, plus description of these GO terms and genes of user's
            list concerned by enriched GO terms
        }
        \item{\code{gene_list}}{
            User's genes given in \code{path_enrich} function
        }
        \item{\code{gomflist}}{
            List containing GO MF terms hierarchies, with indices in the joliMF
            dataframe and genes concerned by the hierarchie
        }
        \item{\code{gobplist}}{
            List containing GO BP terms hierarchies, with indices in the joliMF
            dataframe and genes concerned by the hierarchy
        }
        \item{\code{namegeneid}}{
            Dataframe containing all human genes symbols and names
        }
        \item{\code{hierabrite}}{
            List of pathways categories, containing pathways concerned by a
            category and their indices in "trait" dataframe
        }
        \item{\code{hierapath}}{
            List of hierarchies of pathways, containing indices of pathways
            in "trait" dataframe and elements contained in the hierarchy.
        }
        \item{\code{save_cluster_elem}}{Vector of clustered elements}
        \item{\code{centrality}}{
            Matrix trait-like with values of centralities (number of direct
            interactions between an element of user's list and other elements
            of the pathway) instead of "X"
        }
        \item{\code{inter_values}}{
            Matrix trait-like with values representing direct interactions
            (3/2/1, respectively for genes/metabolites implicated in a
            direct interaction, and for elements not implicated in a direct
            interaction) instead of "X"
        }
        \item{\code{gene_notin}}{
            Dataframe of genes which aren't in pathways and informations about
            them, as gene symbols and names
        }
        \item{\code{sub}}{
            Matrix trait-like with pathway and element informations instead
            of "X"
        }
        \item{\code{quanti}}{
            Matrix trait-like with user's quantitative values instead
            of "X"
        }
        \item{\code{updown}}{
            Matrix trait-like with user's "UP"/"DOWN" values instead
            of "X"
        }
    }
}
\source{
    \code{compl_data} function
}
\references{
    Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He. clusterProfiler:
    an R package for comparing biological themes among gene clusters. OMICS:
    A Journal of Integrative Biology 2012, 16(5):284-287
}
\keyword{datasets}
