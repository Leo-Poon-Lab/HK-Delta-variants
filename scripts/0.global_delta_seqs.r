library(tidyverse)
library(Biostrings)
library(lubridate)

# filter 617 sequences from all data
# untar("../data/metadata_tsv_2021_06_28.tar.xz", files = "metadata.tsv", exdir = "../data/")
# data <- read_tsv("../data/metadata.tsv", col_types = cols(.default = "c"))

# data_delta <- data %>% filter(`Pango lineage` %in% c("B.1.617.1", "B.1.617.2", "B.1.617.3"))
# table(data_delta$`Pango lineage`)
# table(data_delta$`Location`)

# seqs <- readDNAStringSet("../data/mmsa_2021-06-27/2021-06-27_masked.fa")

# sum(names(seqs) %in% data_delta$`Accession ID`)
# seqs_delta <- seqs[names(seqs) %in% data_delta$`Accession ID`]
# seqs_delta <- seqs_delta[order(as.numeric(gsub("\\D", "", names(seqs_delta))))]

# data_delta <- data_delta %>% filter(`Accession ID` %in% names(seqs_delta))
# data_delta <- data_delta %>% arrange(order(as.numeric(gsub("\\D", "", `Accession ID`))))

# writexl::write_xlsx(data_delta, "../results/data_delta.xlsx")
# writeXStringSet(seqs_delta, "../results/seqs_delta.fasta")

# sapply(c("617.1", "617.2", "617.3"), function(lineage_i){
# 	print(lineage_i)
# 	data_tmp <- data_delta %>% filter(grepl(lineage_i, `Pango lineage`, fixed = T))
# 	seqs_tmp <- seqs_delta[names(seqs_delta) %in% data_tmp$`Accession ID`]
# 	writexl::write_xlsx(data_tmp, paste0("../results/data_delta_", lineage_i, ".xlsx"))
# 	writeXStringSet(seqs_tmp, paste0("../results/seqs_delta_", lineage_i, ".fasta"))
# })

# global 617 subset
ref_seq <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
seqs_delta_glo_617_1 <- readDNAStringSet("../results/seqs_delta_617.1.fasta")
seqs_delta_glo_617_2 <- readDNAStringSet("../results/seqs_delta_617.2.fasta")
seqs_delta_glo_617_3 <- readDNAStringSet("../results/seqs_delta_617.3.fasta")

seqs_delta_glo <- c(seqs_delta_glo_617_1, seqs_delta_glo_617_2, seqs_delta_glo_617_3)

unique(seqs_delta_glo)

metadata <- readxl::read_excel("../results/data_delta.xlsx")
data_snps <- read_csv("../results/delta_hk.snpeff.csv")
data_snps <- data_snps %>% filter(effect == "missense_variant")

# global tree
metadata$year <- year(ymd(metadata$`Collection date`))
metadata$week <- week(ymd(metadata$`Collection date`))
metadata$country <- sapply(metadata$Location, function(x){
	strsplit(x, " / ", fixed = T)[[1]][2]
})
table(metadata$country)
set.seed(2021)
metadata <- metadata %>% filter(!grepl("Hong Kong", Location)) %>% group_by(country, year, week) %>% mutate(idx = sample(n())) %>% ungroup()
(metadata_sub <- metadata %>% filter(idx<=3))

sort(names(table(metadata$country)))

seq_global <- seqs_delta_glo[names(seqs_delta_glo) %in% metadata_sub$`Accession ID`] 
seq_hk_delta_ref <- readDNAStringSet("../results/seqs_hk_delta_ref.fasta")
seq_out <- c(seq_hk_delta_ref, seq_global)
seq_out <- subseq(seq_out, 1, width(seq_global)[1])
writeXStringSet(seq_out, "../results/global_hk_617_subset.fasta")
write_csv(metadata_sub, "../results/metadata_global_617_subset.csv")


# drop outliers
seq_hk_delta_ref <- readDNAStringSet("../results/seqs_hk_delta_ref.fasta")
metadata_sub <- read_csv("../results/metadata_global_617_subset.csv")
df_date_hk <- read_csv("../results/seqs_hk_delta_ref_date.csv")
seq_out <- readDNAStringSet("../results/global_hk_617_subset.fasta")

metadata_sub <- metadata_sub %>% filter(!`Accession ID` %in% paste0("EPI_ISL_", c("2661593", "2661578", "2362678", "2434976", "2484895", "2661566", "2688555", "2558719", "1969250", "2648241")))
seq_out <- seq_out[names(seq_out) %in% metadata_sub$`Accession ID`]
seq_out <- c(seq_hk_delta_ref, seq_out)
seq_out <- subseq(seq_out, 1, 29891)
writeXStringSet(seq_out, "../results/global_hk_617_subset.fasta")
write_csv(metadata_sub, "../results/metadata_global_617_subset.csv")

df_date_glo <- metadata_sub %>% select(`Accession ID`, `Collection date`)
names(df_date_glo) <- c("name", "date")
df_date_hk$date <- as.character(df_date_hk$date)

df_date_all <- bind_rows(df_date_hk, df_date_glo)
write_csv(df_date_all, "../results/seqs_all_delta_ref_date.csv")
