library(tidyverse)

cal_bg <- function(lib,sample,nick,size){
  f = sprintf("%s/%s_bwa_sort.bmat",lib,sample)
  bmat_f <-  read_tsv(f)
  bmat_f_bg <- filter(bmat_f,
                      chr_index > 15 & chr_index < max(chr_index)-15+1) %>% filter(
                        chr_index < nick-size | chr_index > nick+size)%>% mutate(
                          mut_ratio = mut_num *100 / (`A`+`G`+`C`+`T`)) %>% filter(
                            mut_ratio <= 10)
  bg = sum(bmat_f_bg$mut_ratio) / nrow(bmat_f_bg)
  
  return(list("BG"=bg,"BG_file"=bmat_f_bg))
}

Res <- cal_bg(lib="FANCF-off-PE3-rep3",sample = "F-AS-mis3-10",nick = 164,size = 25)
Res$BG
Res_f <- test$BG_file
write_csv(Res_f,file = "Test.csv")
