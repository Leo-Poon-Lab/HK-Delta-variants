library(tidyverse)
library(Biostrings)

seqs_hk <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome_caseid.fasta")
ref_seq <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
data_meta <- read_csv("../../2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 20000)

system("pangolin ../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome_caseid.fasta --outfile ../results/lineage.csv")

df_pan <- read_csv("../results/lineage.csv")
df_pan_delta <- df_pan %>% filter(lineage %in% c("B.1.617.1", "B.1.617.2", "B.1.617.3"))
df_pan_delta <- df_pan_delta %>% arrange(taxon)
df_pan_delta <- df_pan_delta %>% filter(!grepl("testindex", taxon)) 
df_pan_delta <- df_pan_delta %>% filter(!grepl("Nepal", taxon)) 
sort(df_pan_delta$taxon)
df_pan_delta$case_id <- sapply(strsplit(df_pan_delta$taxon, "_"), function(x){x[length(x)]})
df_pan_delta <- df_pan_delta %>% group_by(case_id) %>% mutate(n = n()) %>% ungroup()
df_pan_delta %>% filter(n>1) %>% .$taxon
df_pan_delta <- df_pan_delta %>% filter(!(n>1 & grepl("iseq", taxon)))
df_pan_delta <- df_pan_delta %>% group_by(case_id) %>% mutate(n = n()) %>% ungroup()
df_pan_delta$repeated_samples <- df_pan_delta$n > 1

table(df_pan_delta$lineage)

seqs_hk_delta <- seqs_hk[names(seqs_hk) %in% df_pan_delta$taxon]
case_id_delta <- sapply(strsplit(names(seqs_hk_delta), "_"), function(x){x[length(x)]})
data_meta_delta <- data_meta %>% filter(case_id %in% case_id_delta)
sort(data_meta_delta$case_id)
case_id_delta %in% data_meta_delta$case_id

case_id_delta[!case_id_delta %in% data_meta_delta$case_id]
writexl::write_xlsx(data_meta_delta, "../results/metadata_hk_delta.xlsx")
writeXStringSet(seqs_hk_delta, "../results/seqs_hk_delta.fasta")

source("./translate_mod_v3.r")

sapply(c("617.1", "617.2", "617.3"), function(lineage_i){
	print(lineage_i)
	df_pan_tmp <- df_pan_delta %>% filter(grepl(lineage_i, lineage, fixed = T))
	seqs_tmp <- seqs_hk[names(seqs_hk) %in% df_pan_tmp$taxon]
	case_id_tmp <- sapply(strsplit(names(seqs_tmp), "_"), function(x){x[length(x)]})
	data_tmp <- data_meta %>% filter(case_id %in% case_id_tmp)

	seq_s <- get_nt(seqs_tmp, "s")
	seq_s_aa <- get_aa(seqs_tmp, "s")

	writexl::write_xlsx(data_tmp, paste0("../results/metadata_hk_delta_", lineage_i, ".xlsx"))
	writeXStringSet(seqs_tmp, paste0("../results/seqs_hk_delta_", lineage_i, ".fasta"))
	writeXStringSet(seq_s, paste0("../results/seqs_hk_delta_", lineage_i, "_spike_nt.fasta"))
	writeXStringSet(seq_s_aa, paste0("../results/seqs_hk_delta_", lineage_i, "_spike_aa.fasta"))
})


# focus on Delta 1 and 2.
ref_seq <- readDNAStringSet("./reference.fasta")
seqs_hk_617 <- readDNAStringSet("../results/seqs_hk_delta.fasta")
case_id_t <- sapply(strsplit(names(seqs_hk_617), "_"), function(x){x[length(x)]})
seqs_hk_617 <- seqs_hk_617[!duplicated(case_id_t)]
case_id_t <- sapply(strsplit(names(seqs_hk_617), "_"), function(x){x[length(x)]})
seqs_hk_617 <- seqs_hk_617[case_id_t %in% data_meta_delta$case_id]
case_id_t <- sapply(strsplit(names(seqs_hk_617), "_"), function(x){x[length(x)]})

df_Date <- data_meta %>% filter(case_id %in% case_id_t) %>% select(case_id, `Report date`)
df_Date <- left_join(tibble(name = names(seqs_hk_617), case_id = as.numeric(case_id_t)), df_Date)
names(df_Date)[3] <- "date"

df_Date %>% select(name, date) %>% bind_rows(tibble(name = "MN908947_3", date = lubridate::ymd("2019-12-26"))) %>% write_csv("../results/seqs_hk_delta_ref_date.csv")

seqs_hk_617_ref <- c(ref_seq, seqs_hk_617) 
cov_check <- apply(alphabetFrequency(seqs_hk_617_ref), 1, function(x){
	sum(x[1:4])
})
cov_check <- round(cov_check/29903*100, 2)
quantile(cov_check)
sort(cov_check*29903/100)

writeXStringSet(seqs_hk_617_ref, "../results/seqs_hk_delta_ref.fasta")
