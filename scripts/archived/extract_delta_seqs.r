library(tidyverse)
library(Biostrings)

untar("../data/metadata_tsv_2021_06_28.tar.xz", files = "metadata.tsv", exdir = "../data/")
data <- read_tsv("../data/metadata.tsv", col_types = cols(.default = "c"))

data_delta <- data %>% filter(`Pango lineage` %in% c("B.1.617.1", "B.1.617.2", "B.1.617.3"))
table(data_delta$`Pango lineage`)
table(data_delta$`Location`)

seqs <- readDNAStringSet("../data/mmsa_2021-06-27/2021-06-27_masked.fa")

sum(names(seqs) %in% data_delta$`Accession ID`)
seqs_delta <- seqs[names(seqs) %in% data_delta$`Accession ID`]
seqs_delta <- seqs_delta[order(as.numeric(gsub("\\D", "", names(seqs_delta))))]

data_delta <- data_delta %>% filter(`Accession ID` %in% names(seqs_delta))
data_delta <- data_delta %>% arrange(order(as.numeric(gsub("\\D", "", `Accession ID`))))

writexl::write_xlsx(data_delta, "../results/data_delta.xlsx")
writeXStringSet(seqs_delta, "../results/seqs_delta.fasta")

sapply(c("617.1", "617.2", "617.3"), function(lineage_i){
	print(lineage_i)
	data_tmp <- data_delta %>% filter(grepl(lineage_i, `Pango lineage`, fixed = T))
	seqs_tmp <- seqs_delta[names(seqs_delta) %in% data_tmp$`Accession ID`]
	writexl::write_xlsx(data_tmp, paste0("../results/data_delta_", lineage_i, ".xlsx"))
	writeXStringSet(seqs_tmp, paste0("../results/seqs_delta_", lineage_i, ".fasta"))
})
