\name{rshiny}
\alias{rshiny}
\title{Shiny interface}
\description{
    Visualize and filter all functionnal informations gathered by famat using
    a Shiny interface.
}
\usage{
rshiny(listdata)
}
\arguments{
    \item{listdata}{
        Output from \code{compl_data} function
    }
}
\value{
    Shiny interface
}
\references{
    Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson
    (2020). shiny: Web Application Framework for R. R package version 1.5.0.
    https://CRAN.R-project.org/package=shiny
}
\author{
    Emilie Secherre
    \email{emisecherre@gmail.com}
}
\seealso{
    \code{\link{compl_data}}
}
\examples{
    ## load example data
    data(compl_data_result)

    \dontrun{rshiny(compl_data_result)}
}
