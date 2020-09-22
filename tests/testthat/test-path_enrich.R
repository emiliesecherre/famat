setwd(stringr::str_sub(system.file("NAMESPACE",package="famat"),1,32))
genes=data.frame(c("ACAA1","SLC6A12"),c(1.740000e-09,1.136095e-03),c("DOWN","DOWN"))
meta=data.frame(c("C00002","C00719"),c(0.8929425,3.0269316),c("DOWN","UP"))
listk=path_enrich("KEGG",meta,genes)
resmeta=listk[[1]];resgene=listk[[2]]

test_that("pathways are found", {
  expect_equal(nrow(resmeta),7)
  expect_equal(nrow(resgene),11)
})
