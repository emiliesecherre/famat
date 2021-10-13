genes=c("ACAA1","SLC6A12")
meta=c("C00002","C00719")
listk=path_enrich("KEGG",meta,genes)
listr=path_enrich("REAC",meta,genes)
listw=path_enrich("WP",meta,genes)
etape_inter=interactions(listk, listr, listw)
size=etape_inter[[1]];tagged=etape_inter[[3]];central=etape_inter[[5]];no_path=etape_inter[[6]]

test_that("interactions are well extracted", {
#  expect_equal(nrow(size),282)
  expect_equal(nrow(tagged),11)
#  expect_equal(length(central),135)
  expect_equal(nrow(no_path),1)
})

