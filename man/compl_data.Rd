\name{compl_data}
\alias{compl_data}
\title{Data preparation for Shiny interface}
\description{
    Complete and prepare data obtained with \code{interactions} function, to
    use it in the Shiny interface. GO terms enrichment analysis is performed
    using clusterProfiler.
}
\usage{
    compl_data(listparam)
}
\arguments{
    \item{listparam}{
        Output from \code{interactions} function
    }
}
\value{
    A list containing :
    \item{heatmap}{
        Dataframe heatmap-like, with in abscissa elements of pathways ("X" is
        written if an element is present in a pathway), and with in ordinate
        hierarchies of pathways
    }
    \item{meta_list}{User's metabolites given in \code{path_enrich} function}
    \item{allResBP}{
        Results of Go BP terms enrichment analysis performed by clusterProfileR
        (20 best)
    }
    \item{go_genelist}{
        Dataframe containing enriched GO terms per genes of user's list
    }
    \item{allResMF}{
        Results of Go MF terms enrichment analysis performed by clusterProfileR
        (20 best)
    }
    \item{types}{Dataframe containing categories of pathways per pathways}
    \item{genetype}{List of genes categories, containing genes of user's list}
    \item{metatab}{
        Dataframe of metabolites and informations about them, as names and
        chebi ids
    }
    \item{genetab}{
        Dataframe of genes and informations about them, as gene symbols
        and names
    }
    \item{intetab}{
        Dataframe of direct interactions between elements of user's list, and
        informations about them as elements in the interaction, how they
        interact, and which pathways are concerned
    }
    \item{gomf_tab}{
        Dataframe of Go MF terms hierarchies containing our enriched GO terms,
        plus description of these GO terms and genes of user's list concerned
        by enriched GO terms
    }
    \item{gobp_tab}{
        Dataframe of Go BP terms hierarchies containing our enriched GO terms,
        plus description of these GO terms and genes of user's list concerned
        by enriched GO terms
    }
    \item{gene_list}{User's genes given in \code{path_enrich} function}
    \item{gomflist}{
        List containing GO MF terms hierarchies, with indices in the joliMF
        dataframe and genes concerned by the hierarchie
    }
    \item{gobplist}{
        List containing GO BP terms hierarchies, with indices in the joliMF
        dataframe and genes concerned by the hierarchy
    }
    \item{namegeneid}{Dataframe containing all human genes symbols and names}
    \item{hierabrite}{
        List of pathways categories, containing pathways concerned by a
        category and their indices in "trait" dataframe
    }
    \item{hierapath}{
        List of hierarchies of pathways, containing indices of pathways in
        "trait" dataframe and elements contained in the hierarchy.
    }
    \item{save_cluster_elem}{Vector of clustered elements}
    \item{centrality}{
        Matrix trait-like with values of centralities (number of direct
        interactions between an element of user's list and other elements of
        the pathway) instead of "X"
    }
    \item{inter_values}{
        Matrix trait-like with values representing direct interactions (3/2/1,
        respectively for genes/metabolites implicated in a direct interaction,
        and for elements not implicated in a direct interaction) instead of "X"
    }
    \item{gene_notin}{
        Dataframe of genes which aren't in pathways and informations about them,
        as gene symbols and names
    }
    \item{sub}{
        Matrix trait-like with pathway and element informations instead of "X"
    }
    This list is used by \code{rshiny} function.
}
\references{
    Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He. clusterProfiler:
    an R package for comparing biological themes among gene clusters. OMICS:
    A Journal of Integrative Biology 2012, 16(5):284-287
}
\author{
    Emilie Secherre
    \email{emisecherre@gmail.com}
}
\seealso{
    \code{\link{interactions}}
    \code{\link{rshiny}}
}
\examples{
    ## load example data
    data(interactions_result)

    compl_data_result=compl_data(interactions_result)
}
