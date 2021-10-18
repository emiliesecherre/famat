genes=c("ACAA1","SLC6A12")
meta=c("C00002","C00719")
listk=path_enrich("KEGG",meta,genes)
listr=path_enrich("REAC",meta,genes)
listw=path_enrich("WP",meta,genes)
etape_inter=interactions(listk, listr, listw)
etape_donnees=compl_data(etape_inter)
heatmap=etape_donnees[[1]];go_genelist=etape_donnees[[4]];types=etape_donnees[[6]]
genetype=etape_donnees[[7]];gomf_tab=etape_donnees[[11]]
gomflist=etape_donnees[[14]]; hierabrite=etape_donnees[[16]]
hierapath=etape_donnees[[17]];inter_values=etape_donnees[[20]];gene_notin=etape_donnees[[21]]
sub=etape_donnees[[22]]

test_that("informations are good", {
  expect_equal(nrow(heatmap),16)
#  expect_equal(nrow(go_genelist),13)
  expect_equal(nrow(types),12)
  expect_equal(length(genetype),3)
#  expect_equal(nrow(gomf_tab),85)
#  expect_equal(length(gomflist),3)
  expect_equal(length(hierabrite),3)
  expect_equal(length(hierapath),3)
  expect_equal(nrow(inter_values),16)
  expect_equal(nrow(gene_notin),1)
  expect_equal(nrow(sub),16)
})
