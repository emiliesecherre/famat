\name{interactions_result}
\alias{interactions_result}
\docType{data}
\title{
    Output of \code{interactions} function
}
\description{
    List containing informations about interactions between genes and
    metabolites, centrality and pathways. Direct interactions extraction
    was performed using BioPax, KGML and GPML files parsed with PaxtoolsR,
    graphite and author's parsers.
}
\usage{data("interactions_result")}
\format{
    List of 10 items.
    \describe{
        \item{\code{size}}{
            Dataframe containing pathways, genes and metabolites in pathways
            (from the list or not), and number of elements in pathways
        }
        \item{\code{pathtot}}{
            Dataframe containing pathways names and ids from pathway enrichment
            analysis on Reactome, Kegg and Wikipathways pathways
        }
        \item{\code{tagged}}{
            Dataframe containing direct interactions between elements from the
            user's list per pathways
        }
        \item{\code{namegeneid}}{
            Dataframe containing all human genes symbols and names
        }
        \item{\code{keggchebiname}}{
            Dataframe containing all human metabolites ids (kegg and chebi)
            and names
        }
        \item{\code{central}}{
            List of pathways, each pathway containing the number of direct
            interactions between an element of user's list and other elements
            in the pathway
        }
        \item{\code{no_path}}{
            Dataframe containing direct interactions between elements from the
            user's list, but not per pathways
        }
        \item{\code{genes}}{User's genes given in \code{path_enrich} function}
        \item{\code{meta}}{
            User's metabolites given in \code{path_enrich} function
        }
        \item{\code{ensembl}}{
            Result of useMart("ensembl", dataset="hsapiens_gene_ensembl")
            function of BiomaRt
        }
    }
}
\source{
    \code{interactions} function.
}
\references{
    Luna, A., Babur, O., Aksoy, A. B, Demir, E., Sander, C. (2015). “PaxtoolsR:
    Pathway Analysis in R Using Pathway Commons.” Bioinformatics.

    Sales G, Calura E, Cavalieri D, Romualdi C (2012). “graphite - a
    Bioconductor package to convert pathway topology to gene network.” BMC
    Bioinformatics. https://bmcbioinformatics.biomedcentral.com/articles/10 .
    1186/1471-2105-13-20.
}
\keyword{datasets}
