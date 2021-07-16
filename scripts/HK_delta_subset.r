library(tidyverse)
library(Biostrings)

seqs_hk <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome_caseid.fasta")
ref_seq <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
data_meta <- read_csv("../data/cleaned_metadata.csv", guess_max = 20000)

system("pangolin ../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome_caseid.fasta --outfile ../results/lineage.csv")

df_pan <- read_csv("../results/lineage.csv")
df_pan_delta <- df_pan %>% filter(lineage %in% c("B.1.617.1", "B.1.617.2", "B.1.617.3"))
sort(df_pan_delta$taxon)
seq_id_iseq <- sort(df_pan_delta$taxon)[117:123]
df_pan_delta <- bind_rows(df_pan_delta %>% filter(!grepl("iseq", taxon)), df_pan_delta %>% filter(taxon %in% seq_id_iseq))

table(df_pan_delta$lineage)

seqs_hk_delta <- seqs_hk[names(seqs_hk) %in% df_pan_delta$taxon]
case_id_delta <- sapply(strsplit(names(seqs_hk_delta), "_"), function(x){x[length(x)]})
data_meta_delta <- data_meta %>% filter(case_id %in% case_id_delta)
sort(data_meta_delta$case_id)

case_id_delta[!case_id_delta %in% data_meta_delta$case_id]
writexl::write_xlsx(data_meta_delta, "../results/metadata_hk_delta.xlsx")
writeXStringSet(seqs_hk_delta, "../results/seqs_hk_delta.fasta")

source("./translate_mod_v2.r")
orf <- read_csv("../../2020-09-01_COVID_NGS_pipeline/scripts/ORF_SCoV2.csv")
ref_seq_s <- subseq(ref_seq, orf$start[orf$sequence=="S"], orf$stop[orf$sequence=="S"])
ref_seq_s_aa <- translateGappedAln(ref_seq_s)
writeXStringSet(ref_seq_s_aa, "../results/MN908947_3_spike_aa.fasta")

sapply(c("617.1", "617.2", "617.3"), function(lineage_i){
	print(lineage_i)
	df_pan_tmp <- df_pan_delta %>% filter(grepl(lineage_i, lineage, fixed = T))
	seqs_tmp <- seqs_hk[names(seqs_hk) %in% df_pan_tmp$taxon]
	case_id_tmp <- sapply(strsplit(names(seqs_tmp), "_"), function(x){x[length(x)]})
	data_tmp <- data_meta %>% filter(case_id %in% case_id_tmp)

	seq_s <- subseq(seqs_tmp, orf$start[orf$sequence=="S"], orf$stop[orf$sequence=="S"])
	seq_s <- chartr("-", "N", seq_s)
	seq_s_aa <- translateGappedAln(seq_s)

	writexl::write_xlsx(data_tmp, paste0("../results/metadata_hk_delta_", lineage_i, ".xlsx"))
	writeXStringSet(seqs_tmp, paste0("../results/seqs_hk_delta_", lineage_i, ".fasta"))
	writeXStringSet(seq_s, paste0("../results/seqs_hk_delta_", lineage_i, "_spike_nt.fasta"))
	writeXStringSet(seq_s_aa, paste0("../results/seqs_hk_delta_", lineage_i, "_spike_aa_tmp.fasta"))

	file_in <- paste0("../results/seqs_hk_delta_", lineage_i, "_spike_aa_tmp.fasta")
	file_out <- paste0("../results/seqs_hk_delta_", lineage_i, "_spike_aa.fasta")
	system(paste0("mafft --6merpair --thread -10 --keeplength --addfragments ", file_in, " ../results/MN908947_3_spike_aa.fasta > ", file_out)) #MAFFT
	file.remove(file_in)
})


