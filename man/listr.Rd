\name{listr}
\alias{listr}
\docType{data}
\title{
    Pathway enrichment analysis results for Reactome pathways.
}
\description{
    Results of pathways enrichment analysis on the list of genes and
    metabolites, using Reactome pathways. Pathways enrichment analysis
    is performed using MPINet for metabolites and gprofiler2 for genes.
}
\usage{data("listr")}
\format{
    A list of 4 items.
    \describe{
        \item{resmeta}{
            Pathway enrichment analysis results for metabolites.
            A data frame with 278 rows and 2 variables.
            \describe{
                \item{name}{
                    Name of pathways resulting in metabolites pathway
                    enrichment analysis.
                }
                \item{id}{
                    Identifiers of pathways resulting in metabolites
                    pathway enrichment analysis.
                }
            }
        }
        \item{resgene}{
            Pathway enrichment analysis results for genes.
            A data frame with 27 rows and 2 variables.
            \describe{
                \item{name}{
                    Name of pathways resulting in genes pathway
                    enrichment analysis.
                }
                \item{id}{
                    Identifiers of pathways resulting in genes
                    pathway enrichment analysis.
                }
            }
        }
        \item{gened}{Vector containing user's genes (ACAA1, SLC6A12)}
        \item{metad}{Vector containing user's metabolites (C00002, C00719)}
    }
}
\references{
    Yanjun Xu, Chunquan Li and Xia Li (2013). MPINet: The package can implement
    the network-based metabolite pathway identification of pathways.. R package
    version 1.0. https://CRAN.R-project.org/package=MPINet

    Liis Kolberg and Uku Raudvere (2020). gprofiler2: Interface to the
    'g:Profiler' Toolset. R package version 0.2.0.
    https://CRAN.R-project.org/package=gprofiler2
}
\source{
    \code{path_enrich} function.
}
\keyword{datasets}
