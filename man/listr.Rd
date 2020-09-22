\name{listr}
\alias{listr}
\docType{data}
\title{
    Output of \code{path_enrich} function for REACTOME pathways
}
\description{
    Results of pathways enrichment analysis on the list of genes and
    metabolites, using REACTOME pathways. Pathways enrichment analysis
    is performed using MPINet for metabolites and gprofiler2 for genes.
}
\usage{data("listr")}
\format{
    A list of 4 items.
    \describe{
        \item{\code{resmeta}}{
            Results of metabolites pathway enrichment analysis
        }
        \item{\code{resgene}}{Results of genes pathway enrichment analysis}
        \item{\code{gened}}{
            Dataframe containing genes and values related to them
        }
        \item{\code{metad}}{
            Dataframe containing metabolites and values related to them
        }
    }
}
\source{
    \code{path_enrich} function.
}
\references{
    Yanjun Xu, Chunquan Li and Xia Li (2013). MPINet: The package can implement
    the network-based metabolite pathway identification of pathways.. R package
    version 1.0. https://CRAN.R-project.org/package=MPINet

    Liis Kolberg and Uku Raudvere (2020). gprofiler2: Interface to the
    'g:Profiler' Toolset. R package version 0.2.0.
    https://CRAN.R-project.org/package=gprofiler2
}
\keyword{datasets}
