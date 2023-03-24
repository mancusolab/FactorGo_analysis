#### Map genes to HUGO approved gene names

library(tidyverse)
setwd("./data/PanUKBB/3_ldprune/sig_pval5e-08_nct2_maf1e-02_1kg/varinfo")

# Reference (all approved gene symbol)
ref <- read_tsv("./data/GENE_MAP/hgnc_complete_set_20221219.txt") %>%
  mutate(chr = as.integer(sub("q.*|p.*", "", location)),
         Approved_symbol_lower = tolower(symbol))

ref_approved <- ref %>% # filter(Status == "Approved") %>%
  mutate(Approved_symbol_lower = tolower(name)) # 42571 genes

# query for 11583 genes

query <- read_tsv("allvar.pruned.closesttss") %>%
  distinct(Gene_name, chr) %>%
  mutate(Gene_name_lower = tolower(Gene_name),
         name_1kg = Gene_name)

# 10911 found
query_found <- query %>%
  inner_join(ref, by = c("Gene_name_lower" = "Approved_symbol_lower", "chr")) %>%
  select(name_1kg, chr:name, omim_id) %>%
  mutate(notfound=0) %>%
  select(-symbol)

# 672 not found
query_unfound <- query %>%
  anti_join(ref, by = c("Gene_name_lower" = "Approved_symbol_lower", "chr")) %>%
  mutate(old_name = NA,
         alias_name = NA)

# look for if symbol is previously used
for (idx in seq_along(query_unfound$Gene_name)){
  tmp <- ref %>%
    filter(grepl(query_unfound$Gene_name[idx], prev_symbol, ignore.case = TRUE))
  for (each_row in 1:nrow(tmp)){
    # remove white space and parse by |
    gene_list <- tolower(unlist(strsplit(gsub("\\s", "", tmp$prev_symbol[each_row]), split = "\\|")))
    # print(gene_list)

    # check exact match
    if (query_unfound$Gene_name_lower[idx] %in% gene_list){
      query_unfound$old_name[idx] <- tmp$symbol[each_row]
    }
  }
}

# look for if symbol is alias name
for (idx in seq_along(query_unfound$Gene_name)){
  tmp <- ref %>%
    filter(grepl(query_unfound$Gene_name[idx], alias_symbol))
  for (each_row in 1:nrow(tmp)){
    gene_list <- tolower(unlist(strsplit(gsub("\\s", "", tmp$alias_symbol[each_row]), split = "\\|")))
    # print(gene_list)
    if (query_unfound$Gene_name_lower[idx] %in% gene_list){
      query_unfound$alias_name[idx] <- tmp$symbol[each_row]
    }
  }
}

# 19 unfound
query_unfound2 <- query_unfound %>%
  mutate(Gene_name = ifelse(!is.na(old_name), old_name,
                            ifelse(!is.na(alias_name), alias_name, Gene_name)),
         notfound = ifelse((is.na(old_name) & is.na(alias_name)), 1, 0))

allfound <- query_unfound2 %>% left_join(ref, by = c("Gene_name" = "symbol", "chr")) %>%
  select(colnames(query_found)) %>%
  bind_rows(query_found) %>%
  select(-Gene_name_lower)

# go to https://www.genenames.org/
# search for these 19 names
# found 5
# 14 either pseudogene; withdrawn;long-noncoding RNA; RNA
allfound <- allfound %>% mutate(
  Gene_name = case_when(
    Gene_name == "SF3B14" ~ "SF3B6",
    Gene_name == "SGK223" ~ "PRAG1",
    Gene_name == "HDGFRP3" ~ "HDGFL3",
    Gene_name == "SGK494" ~ "RSKR",
    Gene_name == "HDGFRP2" ~ "HDGFL2",
    TRUE ~ Gene_name
  ),
  notfound = ifelse(name_1kg %in% c("SF3B14", "SGK223", "HDGFRP3", "SGK494", "HDGFRP2"), 0, notfound)
)

sum(allfound$notfound) # 14

# back to old file

query <- read_tsv("allvar.pruned.closesttss") %>%
  rename(name_1kg = Gene_name)

out <- query %>%
  left_join(allfound, by=c("name_1kg", "chr"))

all.equal(out$SNP,query$SNP)

out %>% write_tsv("./allvar.pruned.closesttss.hugo")
