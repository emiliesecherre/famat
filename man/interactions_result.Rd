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
    List of 8 items.
    \describe{
        \item{size}{
            Description on which elements (from user's list or not) are
            contained in pathways from pathway enrichment analysis.
            A data frame with 286 rows and 9 variables.
            \describe{
                \item{path}{
                    Pathways obtained throught pathways enrichment analysis on
                    KEgg, Reactome and Wikipathways pathways.
                }
                \item{nb_gene_query}{
                    Number of user's genes contained in the pathway (0--2).
                }
                \item{gene_que}{
                    User's genes contained in the pathway
                    (ACAA1, SLC6A12, ACAA1 # SLC6A12).
                }
                \item{nb_gene_tot}{
                    Total number of genes contained in the pathway
                    (0--2075).
                }
                \item{genes}{
                    All the genes contained in the pathway.
                }
                \item{nb_meta_query}{
                    Number of user's metabolites contained in the pathway
                    (0--2).
                }
                \item{meta_que}{
                    User's metabolites contained in the pathway
                    (Betaine, ATP, Betaine # ATP, ATP # Betaine).
                }
                \item{nb_meta_tot}{
                    Total number of metabolites contained in the pathway
                    (0--915).
                }
                \item{meta}{
                    All the metabolites contained in the pathway.
                }
            }
        }
        \item{pathtot}{
            All results of pathways enrichment analysis performed on Kegg,
            Reactome and Wikipathways pathways.
            A data frame with 286 rows and 2 variables.
            \describe{
                \item{name}{
                    Name of pathways resulting in genes pathway
                    enrichment analysis performed on Kegg, Reactome
                    and Wikipathways.
                }
                \item{id}{
                    Identifiers of pathways resulting in genes pathway
                    enrichment analysis performed on Kegg, Reactome
                    and Wikipathways.
                }
            }
        }
        \item{tagged}{
            Description of all direct interactions between user's elements in
            pathways.
            A data frame with 11 rows and 6 variables.
            \describe{
                \item{from}{
                    First element of the direct interaction (Betaine, SLC6A12)
                }
                \item{link}{
                    Description of how the two elements interact
                    (Control(In: ACTIVATION of BiochemicalReaction),
                    controls-transport-of-chemical).
                }
                \item{to}{
                    Second element of the direct interaction (Betaine, SLC6A12)
                }
                \item{path}{
                    Pathway containing the direct interaction (R-HSA-112310,
                    R-HSA-112315, R-HSA-112316, R-HSA-382551, R-HSA-425366,
                    R-HSA-425393, R-HSA-425407, R-HSA-888590, R-HSA-352230,
                    R-HSA-442660, R-HSA-888593)
                }
                \item{tag}{
                    Summary of elements concerned by the interaction
                    (Betaine / SLC6A12, SLC6A12 / Betaine)
                }
                \item{type}{
                    Interaction type, can be gene/gene, metabolite/metabolite,
                    or gene/metabolite (g/m)
                }
            }
        }
        \item{keggchebiname}{
            Dataframe containing all human metabolites ids (kegg and chebi)
            and names.
            A data frame with 16075 rows and 3 variables.
            \describe{
                \item{kegg}{
                    Kegg_compound identifiers of all human metabolites.
                }
                \item{chebi}{
                    Chebi identifiers of all human metabolites.
                }
                \item{name}{
                    Names of all human metabolites.
                }
            }
        }
        \item{central}{
            List of pathways, each pathway containing the number of direct
            interactions between an element of user's list and other elements
            in the pathway.
            A list of 138 items.
        }
        \item{no_path}{
            Dataframe containing direct interactions between elements from the
            user's list, but not per pathways.
            A data frame with 1 rows and 6 variables.
            \describe{
                \item{from}{
                    First element of the direct interaction (Betaine, SLC6A12)
                }
                \item{link}{
                    Description of how the two elements interact
                    (Control(In: ACTIVATION of BiochemicalReaction),
                    controls-transport-of-chemical).
                }
                \item{to}{
                    Second element of the direct interaction (Betaine, SLC6A12)
                }
                \item{path}{
                    Pathways containing the direct interaction ("R-HSA-112310,
                    R-HSA-112315, R-HSA-112316, R-HSA-382551, R-HSA-425366,
                    R-HSA-425393, R-HSA-425407, R-HSA-888590, R-HSA-352230,
                    R-HSA-442660, R-HSA-888593")
                }
                \item{tag}{
                    Summary of elements concerned by the interaction
                    (Betaine / SLC6A12)
                }
                \item{type}{
                    Interaction type, can be gene/gene, metabolite/metabolite,
                    or gene/metabolite (g/m)
                }
            }
        }
        \item{genes}{Vector containing user's genes (ACAA1, SLC6A12)}
        \item{meta}{Vector containing user's metabolites (C00002, C00719)}
    }
}
\source{
    \code{interactions} function.
}
\references{
    Luna, A., Babur, O., Aksoy, A. B, Demir, E., Sander, C. (2015).“PaxtoolsR:
    Pathway Analysis in R Using Pathway Commons.” Bioinformatics.

    Sales G, Calura E, Cavalieri D, Romualdi C (2012). “graphite - a
    Bioconductor package to convert pathway topology to gene network.” BMC
    Bioinformatics. https://bmcbioinformatics.biomedcentral.com/articles/10 .
    1186/1471-2105-13-20.
}
\keyword{datasets}
