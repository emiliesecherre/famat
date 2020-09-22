\name{genes}
\alias{genes}
\docType{data}
\title{
    List of genes.
}
\description{
    Dataframe contaning a list of genes, quantitative data about them, and
    mention "DOWN" or "UP" depending on gene expression in a certain condition.
}
\usage{data("genes")}
\format{
    A data frame with 2 observations on the following 3 variables.
    \describe{
        \item{\code{gene_symbol}}{List of genes}
        \item{\code{data}}{Quantitative data}
        \item{\code{up_down}}{"UP" or "DOWN"}
    }
}
\source{
    Sample of data from a study on chickens, under heat-stress condition.
}
\keyword{datasets}
