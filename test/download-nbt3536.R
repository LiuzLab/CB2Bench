shrna.link <- "http://www.nature.com/nbt/journal/v34/n6/extref/nbt.3536-S2.xlsx"
crispr.link <- "http://www.nature.com/nbt/journal/v34/n6/extref/nbt.3536-S3.xlsx"
crispri.link <- "http://www.nature.com/nbt/journal/v34/n6/extref/nbt.3536-S4.xlsx"

fetch_file <- function(link) {
  require(readxl)
  tmp <- tempfile()
  download.file(link, tmp)
  read_xlsx(tmp)
}
df.shrna <- fetch_file(shrna.link)
df.crispr <-fetch_file(crispr.link)
df.crispri <-fetch_file(crispri.link)

shRNA.RT112 <- select(df.shrna, 2, 3, 1, 4, 5:7, 8:10)
shRNA.UMUC3 <- select(df.shrna, 2, 3, 1, 4, 17:22)

crispr.RT112 <- select(df.crispr, 2, 3, 1, 4, 5:7, 8:10)
crispr.UMUC3 <- select(df.crispr, 2, 3, 1, 4, 17:22)

crispri.RT112 <- select(df.crispri, 2, 3, 1, 4, 5:7, 8:10)


dataset <- list("shRNA.RT112" = shRNA.RT112,
                "shRNA.UMUC3" = shRNA.UMUC3,
                "CRISPR.RT112"= crispr.RT112,
                "CRIPSR.UMUC3" = crispr.UMUC3,
                "CRISPRi.RT112" = crispri.RT112)

for(d in names(dataset)) {
  colnames(dataset[[d]]) <- c("X", "gene", "sgRNA", "class", "B1", "B2", "B3", "A1", "A2", "A3")
  dataset[[d]]$class[dataset[[d]]$class==0] <- "inactive"
  dataset[[d]]$class[dataset[[d]]$class==1] <- "decreasing"

}
save(dataset, file="inst/extdata/nature-biotech.Rdata")
