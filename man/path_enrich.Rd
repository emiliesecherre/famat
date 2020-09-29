\name{path_enrich}
\alias{path_enrich}
\title{Pathway enrichment analysis}
\description{
    Perform a pathway enrichment analysis using a list of genes and a list of
    metabolites. Pathways enrichment analysis is performed using MPINet for
    metabolites and gprofiler2 for genes.
}
\usage{
    path_enrich(source, metabo, genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{source}{
        Pathways database used, either Kegg ("KEGG"), Reactome ("REAC") or
        Wikipathays ("WP")
    }
    \item{metabo}{
        Dataframe with three columns : the first column contain the list of
        metabolites, the second some quantitative data about the metabolites,
        the last one words "DOWN" or "UP" depending on the metabolites
        concentration behavior in a certain condition. Last two columns can
        contain only/some NAs. All metabolites ids are KEGG Compound ids.
    }
    \item{genes}{
        Dataframe with three columns : the first column contain the list of
        genes, the second some quantitative data about the genes, the last
        one words "DOWN" or "UP" depending on the genes expression behavior
        in a certain condition. Last two columns can contain only/some NAs.
        All genes ids are gene symbol.
    }
}
\value{
    A list containing :
    \item{resmeta}{Results of metabolites pathway enrichment analysis}
    \item{resgene}{Results of genes pathway enrichment analysis}
    \item{genes}{Vector containing genes}
    \item{metabo}{Vector containing metabolites}
    This list is used by \code{interactions} function.
}
\references{
    Yanjun Xu, Chunquan Li and Xia Li (2013). MPINet: The package can implement
    the network-based metabolite pathway identification of pathways.. R package
    version 1.0. https://CRAN.R-project.org/package=MPINet

    Liis Kolberg and Uku Raudvere (2020). gprofiler2: Interface to the
    'g:Profiler' Toolset. R package version 0.2.0.
    https://CRAN.R-project.org/package=gprofiler2
}
\author{
    Emilie Secherre
    \email{emisecherre@gmail.com}
}
\seealso{
    \code{\link{interactions}}
}
\examples{
    ## load example data
    data(genes)
    data(meta)

    ## perform pathway enrichment analysis on Reactome pathways
    listr=path_enrich("REAC", meta, genes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
