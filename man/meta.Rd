\name{meta}
\alias{meta}
\docType{data}
\title{
    List of metabolites.
}
\description{
    Dataframe contaning a list of metabolites, quantitative data about them,
    and mention "DOWN" or "UP" depending on metabolite levels in a certain
    condition.
}
\usage{data("meta")}
\format{
    A data frame with 2 observations on the following 3 variables.
    \describe{
        \item{\code{meta_id}}{List of metabolites}
        \item{\code{data}}{Quantitative data}
        \item{\code{up_down}}{"UP" or "DOWN"}
    }
}
\source{
    Sample of data from a study on chickens, under heat-stress condition.
}
\keyword{datasets}
