genes=c("ACAA1","SLC6A12")
meta=c("C00002","C00719")
listk=path_enrich("KEGG",meta,genes)
resmeta=listk[[1]];resgene=listk[[2]]

test_that("pathways are found", {
  #expect_equal(nrow(resmeta),7)
  expect_equal(nrow(resgene),11)
})
