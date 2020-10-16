\name{compl_data_result}
\alias{compl_data_result}
\docType{data}
\title{
    Output of \code{compl_data} function
}
\description{
    Enriched and structured informations about genes, metabolites, their
    interactions, pathways and enriched GO terms.
}
\usage{data("compl_data_result")}
\format{
    List of 23 items.
    \describe{
        \item{heatmap}{
            Dataframe heatmap-like, with in abscissa elements of pathways
            ("X" is written if an element is present in a pathway), and with
            in ordinate hierarchies of pathways.
            A data frame with 16 rows and 9 variables.
            \describe{
                \item{path_name}{
                    Hierarchies of pathways obtained by pathways enrichment
                    analysis. Pathways are given by their name.
                }
                \item{path_id}{
                    Identifiers of pathways in hierarchies.
                }
                \item{meta_ratio}{
                    Metabolites ratio, so the part of user's metabolites in
                    the total number of metabolites in the pathway.
                }
                \item{gene_ratio}{
                    Genes ratio, so the part of user's genes in
                    the total number of genes in the pathway.
                }
                \item{blank}{
                    Just an empty column to separate heatmap data from pathways
                    data (NA).
                }
                \item{SLC6A12}{
                    The value "1" shows that the element SLC6A12 is in the
                    pathway on the row, "0" shows it is not in this pathway
                    (0--1).
                }
                \item{Betaine}{
                    The value "1" shows that the element Betaine is in the
                    pathway on the row, "0" shows it is not in this pathway
                    (0--1).
                }
                \item{ATP}{
                    The value "1" shows that the element ATP is in the
                    pathway on the row, "0" shows it is not in this pathway
                    (0--1).
                }
                \item{Betaine / SLC6A12}{
                    The value "1" shows that the interaction Betaine / SLC6A12
                    is in the pathway on the row, "0" shows it is not in this
                    pathway (0--1).
                }
            }
        }
        \item{meta_list}{
            Vector containing user's metabolites (C00002, C00719)
        }
        \item{allResBP}{
            Results of Go BP terms enrichment analysis performed by
            clusterProfileR (20 best).
            A data frame with 20 rows and 9 variables.
            \describe{
                \item{ID}{
                    Identifiers of enriched GO terms.
                }
                \item{Description}{
                    Names of enriched GO terms.
                }
                \item{GeneRatio}{
                    Number of user's genes concerned by the enriched GO term,
                    by total number of user's genes.
                }
                \item{BgRatio}{
                    Number of genes concerned by the enriched GO term, by the
                    total number of annotated genes in the database.
                }
                \item{pvalue}{
                    Pvalue of the go term enrichment analysis
                    (0.001392161--0.011536208).
                }
                \item{p.adjust}{
                    Adjusted pvalue of the go term enrichment analysis
                    (0.01713291--0.02670275).
                }
                \item{qvalue}{
                    Qvalue of the go term enrichment analysis
                    (0.001803464--0.002810816).
                }
                \item{geneID}{
                    Entrez Gene identifiers of all genes concerned by the
                    enriched GO term (30, 6539).
                }
                \item{Count}{
                    Number of genes concerned by the enriched GO term described
                    by an Entrez Gene identifier(1).
                }
            }
        }
        \item{go_genelist}{
            Dataframe containing enriched GO terms per genes of user's list.
            A data frame with 13 rows and 2 variables.
            \describe{
                \item{hgnc_symbol}{
                    Gene symbol of the gene concerned by an enriched GO term
                    (SLC6A12, ACAA1).
                }
                \item{go_id}{
                    Identifier of the GO term concerning the gene (GO:0005328,
                    GO:0015293, GO:0003333, GO:0015812, GO:0015171, GO:0042165,
                    GO:0008028, GO:0006635, GO:0006625, GO:0033540, GO:0036109,
                    GO:0008206, GO:0000038)
                }
            }
        }
        \item{allResMF}{
            Results of Go MF terms enrichment analysis performed by
            clusterProfileR (20 best)
            A data frame with 20 rows and 9 variables.
            \describe{
                \item{ID}{
                    Identifiers of enriched GO terms.
                }
                \item{Description}{
                    Names of enriched GO terms.
                }
                \item{GeneRatio}{
                    Number of user's genes concerned by the enriched GO term,
                    by total number of user's genes.
                }
                \item{BgRatio}{
                    Number of genes concerned by the enriched GO term, by the
                    total number of annotated genes in the database.
                }
                \item{pvalue}{
                    Pvalue of the go term enrichment analysis
                    (0.001392161--0.011536208).
                }
                \item{p.adjust}{
                    Adjusted pvalue of the go term enrichment analysis
                    (0.01713291--0.02670275).
                }
                \item{qvalue}{
                    Qvalue of the go term enrichment analysis
                    (0.001803464--0.002810816).
                }
                \item{geneID}{
                    Entrez Gene identifiers of all genes concerned by the
                    enriched GO term (30, 6539).
                }
                \item{Count}{
                    Number of genes concerned by the enriched GO term described
                    by an Entrez Gene identifier(1).
                }
            }
        }
        \item{types}{
            Dataframe containing categories of pathways per pathways.
            A data frame with 12 rows and 2 variables.
            \describe{
                \item{id}{
                    Identifier of pathways from pathways enrichment analysis.
                }
                \item{root}{
                    Name of the pathway category concerning a pathway.
                }
            }
        }
        \item{genetype}{
            List of genes categories, containing genes of user's list.
            A list of 3 items.
        }
        \item{metatab}{
            Dataframe of metabolites and informations about them, as
            names and chebi ids.
            A data frame with 2 rows and 2 variables.
            \describe{
                \item{id}{
                    Name of user's metabolites (ATP, Betaine).
                }
                \item{name}{
                    Chebi identifier of user's metabolites
                    (CHEBI:15422, CHEBI:17750).
                }
            }
        }
        \item{genetab}{
            Dataframe of genes and informations about them, as gene
            symbols and names.
            A data frame with 1 rows and 2 variables.
            \describe{
                \item{id}{
                    Gene symbols of user's genes contained in pathways
                    (SLC6A12).
                }
                \item{name}{
                    Name of user's genes contained in pathways
                    (solute carrier family 6 member 12).
                }
            }
        }
        \item{intetab}{
            Dataframe of direct interactions between elements of user's list,
            and informations about them as elements in the interaction, how
            they interact, and which pathways are concerned
            A data frame with 1 row and 8 variables.
            \describe{
                \item{tag}{
                    Summary of elements concerned by the interaction
                    (Betaine / SLC6A12)
                }
                \item{first_item}{
                    First element of the direct interaction (Betaine)
                }
                \item{link}{
                    Description of how the two elements interact
                    (Control(In: ACTIVATION of BiochemicalReaction),
                    controls-transport-of-chemical).
                }
                \item{sec_item}{
                    Second element of the direct interaction (SLC6A12)
                }
                \item{go}{
                    Value "1" means that a gene of the interaction is concerned
                    by an enriched GO term, "0" means no element is concerned
                    by an enriched GO term (1).
                }
                \item{path}{
                    Pathways containing the direct interaction ("R-HSA-112310,
                    R-HSA-112315, R-HSA-112316, R-HSA-382551, R-HSA-425366,
                    R-HSA-425393, R-HSA-425407, R-HSA-888590, R-HSA-352230,
                    R-HSA-442660, R-HSA-888593")
                }
                \item{type}{
                    Interaction type, can be gene/gene, metabolite/metabolite,
                    or gene/metabolite (g/m)
                }
                \item{cat}{
                    Categories of pathways containing the direct interaction
                    (Neuronal System, Transport of small molecules)
                }
            }

        }
        \item{gomf_tab}{
            Dataframe of Go MF terms hierarchies containing our enriched GO
            terms, plus description of these GO terms and genes of user's
            list concerned by enriched GO terms.
            A data frame with 93 rows and 3 variables.
            \describe{
                \item{goterm}{
                    Hierarchies of enriched GO terms.
                }
                \item{go_name}{
                    Names of GO terms.
                }
                \item{genes}{
                    Genes concerned by GO terms.
                }
            }
        }
        \item{gobp_tab}{
            Dataframe of Go BP terms hierarchies containing our enriched GO
            terms, plus description of these GO terms and genes of user's
            list concerned by enriched GO terms.
            A data frame with 107 rows and 3 variables.
            \describe{
                \item{goterm}{
                    Hierarchies of enriched GO terms.
                }
                \item{go_name}{
                    Names of GO terms.
                }
                \item{genes}{
                    Genes concerned by GO terms.
                }
            }
        }
        \item{gene_list}{Vector containing user's genes (ACAA1, SLC6A12)}
        \item{gomflist}{
            List containing GO MF terms hierarchies, with indices in the joliMF
            dataframe and genes concerned by the hierarchie.
            A list of 3 items.
        }
        \item{gobplist}{
            List containing GO BP terms hierarchies, with indices in the joliMF
            dataframe and genes concerned by the hierarchy.
            A list of 5 items.
        }
        \item{namegeneid}{
            Dataframe containing all human genes symbols and names.
            A data frame with 39946 rows and 2 variables.
            \describe{
                \item{hgnc_symbol}{
                    Gene symbols of all human genes.
                }
                \item{entregene_description}{
                    Names of all human genes.
                }
            }
        }
        \item{hierabrite}{
            List of pathways categories, containing pathways concerned by a
            category and their indices in "trait" dataframe.
            A list of 3 items.
        }
        \item{hierapath}{
            List of hierarchies of pathways, containing indices of pathways
            in "trait" dataframe and elements contained in the hierarchy.
            A list of 3 items.
        }
        \item{save_cluster_elem}{Vector of clustered elements}
        \item{centrality}{
            Matrix heatmap-like with values of centralities (number of direct
            interactions between an element of user's list and other elements
            of the pathway) instead of "X". Other cells contain the value "0"
            (0--65).
            An integer matrix with 16 rows and 9 columns.
        }
        \item{inter_values}{
            Matrix heatmap-like with values representing direct interactions
            (3/2/1, respectively for genes/metabolites implicated in a
            direct interaction, and for elements not implicated in a direct
            interaction) instead of "X". Other cells contain the value "0"
            (0--3).
            An integer matrix with 16 rows and 9 columns.
        }
        \item{gene_notin}{
            Dataframe of genes which aren't in pathways and informations about
            them, as gene symbols and names.
            A data frame with 1 row and 2 variables.
            \describe{
                \item{id}{
                    Gene symbols of genes (ACAA1).
                }
                \item{name}{
                    Names of genes (acetyl-CoA acyltransferase 1).
                }
            }
        }
        \item{sub}{
            Matrix heatmap-like with pathway and element informations instead
            of "X". Cells with no informations contain only "".
            A character matrix with 16 rows and 9 columns.
        }
    }
}
\source{
    \code{compl_data} function
}
