\name{interactions}
\alias{interactions}
\title{Interactions between genes and metabolites}
\description{
    Gather informations about direct interactions between genes and metabolites
    inside a pathway, and about pathways themselves. These informations are
    direct interactions between these two elements and number of relations
    between an element from the list provided by the user and other elements of
    the pathway (centrality). Direct interactions extraction was performed
    using BioPax, KGML and GPML files parsed with PaxtoolsR, graphite and
    author's parsers.
}
\usage{
    interactions(listk, listr, listw)
}
\arguments{
    \item{listk}{
        Output from \code{path_enrich} function, with "KEGG" argument.
    }
    \item{listr}{
        Output from \code{path_enrich} function, with "REAC" argument.
    }
    \item{listw}{
        Output from \code{path_enrich} function, with "WP" argument.
    }
}
\value{
    A list containing :
    \item{size}{
        Dataframe containing pathways, genes and metabolites in pathways (from
        the list or not), and number of elements in pathways
    }
    \item{pathtot}{
        Dataframe containing pathways names and ids from pathway enrichment
        analysis on Reactome, Kegg and Wikipathways pathways
    }
    \item{tagged}{
        Dataframe containing direct interactions between elements from the
        user's list per pathways
    }
    \item{keggchebiname}{
        Dataframe containing all human metabolites ids (kegg and chebi) and
        names
    }
    \item{central}{
        List of pathways, each pathway containing the number of direct
        interactions between an element of user's list and other elements
        in the pathway
    }
    \item{no_path}{
        Dataframe containing direct interactions between elements from the
        user's list, but not per pathways
    }
    \item{genes}{User's genes given in \code{path_enrich} function}
    \item{meta}{User's metabolites given in \code{path_enrich} function}
    This list is used by \code{compl_data} function.
}
\references{
    Luna, A., Babur, O., Aksoy, A. B, Demir, E., Sander, C. (2015).“PaxtoolsR:
    Pathway Analysis in R Using Pathway Commons.” Bioinformatics.

    Sales G, Calura E, Cavalieri D, Romualdi C (2012). “graphite - a
    Bioconductor package to convert pathway topology to gene network.” BMC
    Bioinformatics. https://bmcbioinformatics.biomedcentral.com/articles/10 .
    1186/1471-2105-13-20.
}
\author{
    Emilie Secherre
    \email{emisecherre@gmail.com}
}
\seealso{
    \code{\link{path_enrich}}
    \code{\link{compl_data}}
}
\examples{
    ## load example data
    data(listk)
    data(listr)
    data(listw)

    interactions_result=interactions(listk,listr,listw)
}
