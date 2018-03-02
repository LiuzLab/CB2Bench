load("inst/extdata/nature-biotech.Rdata")

for(d in names(dataset)) {
  file_path <- sprintf("inst/extdata/%s.csv", d)
  print(file_path)
  dataset[[d]] %>% select(-X, -class) %>% write_csv(file_path)
}
#
# for(d in names(dataset)) {
#   file_path <- sprintf("inst/extdata/%s_essentiality.csv", d)
#   print(file_path)
#   dataset[[d]]%>% write_csv(file_path)
# }
