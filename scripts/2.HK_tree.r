library(tidyverse)
library(Biostrings)
library(lubridate)
library(ggtree)
library(ggrepel)
library(colorspace)
library(ggnewscale)

data_meta <- read_csv("../../2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 20000)
df_date <- read_tsv("../results/dates.tsv")
seqs_hk_617_1 <- readDNAStringSet("../results/seqs_hk_delta_617.1.fasta")
seqs_hk_617_2 <- readDNAStringSet("../results/seqs_hk_delta_617.2.fasta")
date_meta_glo <- read_csv("../results/metadata_global_617_subset.csv")
load("./colors.rdata")

# system(paste0('~/softwares/iqtree-2.1.3-MacOSX/bin/iqtree2 -T 10 -m JC --redo -s ../results/seqs_hk_delta_ref.fasta -o "MN908947_3"'))
# ancestral sequences
# system("treetime --aln ../results/seqs_hk_delta_ref.fasta --tree ../results/seqs_hk_delta_ref.fasta.treefile --dates ../results/seqs_hk_delta_ref_date.csv --outdir ../results/")
# system("treetime ancestral --aln ../results/seqs_hk_delta_ref.fasta --tree ../results/timetree.nexus --outdir ../results/")

# tree <- ape::read.tree("../results/seqs_hk_delta_ref.fasta.treefile")
tree <- treeio::read.beast("../results/annotated_tree.nexus")
treetext <- readLines("../results/timetree.nexus")
treetext <- treetext[which.max(sapply(treetext, nchar))]
date_t <- strsplit(treetext, "date=", fixed = T)[[1]]
date_t <- sapply(date_t, function(x){
	as.numeric(strsplit(x, "]", fixed = T)[[1]][1])
}, USE.NAMES = F)

mrsd <- max(date_t, na.rm = T)
(p <- ggtree(tree, mrsd=decimal2Date(mrsd)) + theme_tree2())

clean_p <- function(p){
	p$data$case_id <- sapply(p$data$label, function(x){
		if(is.na(x)){return(x)}
		if(x == "MN908947_3"){return(x)}
		tmp <- strsplit(x, "_")[[1]]
		tmp[length(tmp)]	
	})

	p$data$is_617_1 <- p$data$label %in% names(seqs_hk_617_1)
	p$data$is_617_2 <- p$data$label %in% names(seqs_hk_617_2)
	p$data$lineage <- NA
	# p$data$lineage[p$data$is_617_1] <- "B.1.617.1"
	# p$data$lineage[p$data$is_617_2] <- "B.1.617.2"

	p$data$date <- sapply(p$data$case_id, function(x){
		if(is.na(x)){return(x)}
		if(x == "MN908947_3"){return("2019-12-26")}
		data_meta %>% filter(case_id == x) %>% .$`Report date` %>% as.character()
	})

	p$data$metaCluster <- sapply(p$data$case_id, function(x){
		if(is.na(x)){return(x)}
		if(x == "MN908947_3"){return("None")}
		tmp <- data_meta %>% filter(case_id == x) %>% .$metaCluster
		if(length(tmp)==0){return(NA)}else if(is.na(tmp)){return("Sporadic")}else{return(tmp)}
	})
	l_t <- unique(p$data$metaCluster[!is.na(p$data$metaCluster)])
	p$data$Cluster <- factor(p$data$metaCluster, levels = l_t, labels = c(l_t[1:2], paste0("Cluster_", 1:(length(l_t)-2))))

	p$data$country <- sapply(p$data$case_id, function(x){
		if(is.na(x)){return(x)}
		if(x == "MN908947_3"){return("China")}
		tmp <- data_meta %>% filter(case_id == x) %>% .$`Country of importation` %>% as.character()
		if(length(tmp)==0){return(NA)}else{return(tmp)}
	})
	p$data$country[is.na(p$data$country)] <- "Local"

	p$data$whp_id <- sapply(p$data$label, function(x){
		if(is.na(x)){return(x)}
		if(x == "MN908947_3"){return(x)}
		strsplit(x, "_")[[1]][1]
	})
	p$data$label_new <- paste0(p$data$country, "_", p$data$date, "_", p$data$Cluster)
	return(p)
}
p <- clean_p(p)
table(p$data$country)

p$data$label_new[p$data$label_new=="China_2019-12-26_None"] <- "MN908947_3"
p$data$country[p$data$label_new=="MN908947_3"] <- NA

# annotate with snpEFF
source("./helper/prepare_vcf.r")
# prepare_vcf(p, "../results/delta_hk")
# system("java -jar ~/softwares/snpEff/snpEff.jar download -c ~/softwares/snpEff/snpEff.config -v MN908947.3")
# system("java -jar ~/softwares/snpEff/snpEff.jar MN908947.3 -fastaProt ../results/delta_hk.snpeff.vcf.faa -csvStats ../results/delta_hk.snpeff.vcf.stats ../results/delta_hk.vcf > ../results/delta_hk.snpeff.vcf")
# system("Rscript ./helper/Parse_SnpEff.r ../results/delta_hk.snpeff.vcf ../results/delta_hk.snpeff.csv")

df_ann <- read_csv("../results/delta_hk.snpeff.csv")
df_ann$annot <- NA
df_ann$annot[df_ann$in_gene] <- df_ann$mut_aa[df_ann$in_gene]
df_ann$annot <- gsub("p.", "", df_ann$annot, fixed = T)
df_ann$mut <- gsub("MUT=", "", df_ann$mut, fixed = T)
df_ann$annot[!df_ann$in_gene] <- paste0("nt", df_ann$mut[!df_ann$in_gene])
df_ann_s <- df_ann %>% filter(gene == "S") # only consider spike gene
df_ann_ns <- df_ann %>% filter(!gene == "S") # consider other genes

p$data$mutations_others <- sapply(p$data$mutations, function(x){
	if(length(x)==0){return(NA)}
	if(all(x=="")){return(NA)}
	tmp <- sapply(x, function(y){
		mut <- df_ann_ns$annot[df_ann_ns$mut == y]
		check <- df_ann_ns$effect[df_ann_ns$mut == y]
		if(length(check)==0){return(NA)}
		if(check == "missense_variant"){return(mut)} else {return(NA)} # only consider missense_variant
	})
	tmp <- tmp[!grepl("-", tmp)] # ignore deletion
	tmp <- tmp[!is.na(tmp)]
	if(length(tmp)==0){return(NA)}
	paste0(tmp, collapse = " ")
})

p$data$mutations_spike <- sapply(p$data$mutations, function(x){
	if(length(x)==0){return(NA)}
	if(all(x=="")){return(NA)}
	tmp <- sapply(x, function(y){
		mut <- df_ann_s$annot[df_ann_s$mut == y]
		check <- df_ann_s$effect[df_ann_s$mut == y]
		if(length(check)==0){return(NA)}
		if(check == "missense_variant"){return(mut)} else {return(NA)} # only consider missense_variant
	})
	tmp <- tmp[!grepl("-", tmp)] # ignore deletion
	tmp <- tmp[!is.na(tmp)]
	if(length(tmp)==0){return(NA)}
	paste0(tmp, collapse = " ")
})

nodes <- unique(tree@phylo$edge[,1])
check <- sapply(nodes, function(x) {
	tmp <- phangorn::Children(tree@phylo, x)
	all(tmp <= length(tree@phylo$tip.label))
})
end_nodes <- nodes[check]
color_others <- "#a3a3a3"
names(color_others) <- "Others"

p1 <- p + 
	geom_tiplab(aes(color = country, label = label_new), size = 2, offset = 0.01)+
	geom_tippoint(aes(color = country), size = 1)+
	# scale_color_discrete_sequential(name = "Country",palette = "Batlow", na.translate=FALSE, guide = guide_legend(override.aes=aes(size = 3)))+
	scale_color_manual(name = "Country", values = c(colors_countries_t, color_others), na.value = "#a3a3a3", guide = guide_legend(override.aes=aes(size = 4, alpha = 1, shape = 19)))+
	# geom_nodelab(aes(label = mutations_new), size = 1, nudge_y = 0.2, color = "#a45e67")+
	# scale_color_manual(values = colors_t)+
	xlim(2020,2022.2)+
	# geom_treescale(0.0001, 20, offset = 1, fontsize = 2)+
	theme(legend.key.width=unit(0.1,"cm"))

# system(paste0('~/softwares/iqtree-2.1.3-MacOSX/bin/iqtree2 -T 16 -m JC --redo -s ../results/global_hk_617_subset.fasta -o "MN908947_3"'))
# dir.create("../results/global_tree")
# system("treetime --aln ../results/global_hk_617_subset.fasta --tree ../results/global_hk_617_subset.fasta.treefile --dates ../results/seqs_all_delta_ref_date.csv --outdir ../results/global_tree/")
# system("treetime ancestral --aln ../results/global_hk_617_subset.fasta --tree ../results/global_tree/timetree.nexus --outdir ../results/global_tree/")

tree2_raw <- treeio::read.beast("../results/global_tree/annotated_tree.nexus")
p2 <- ggtree(tree2_raw)

tree2_raw <- treeio::read.beast("../results/global_tree/annotated_tree.nexus")
treetext <- readLines("../results/global_tree/timetree.nexus")
treetext <- treetext[which.max(sapply(treetext, nchar))]
date_t <- strsplit(treetext, "date=", fixed = T)[[1]]
date_t <- sapply(date_t, function(x){
	as.numeric(strsplit(x, "]", fixed = T)[[1]][1])
}, USE.NAMES = F)

mrsd <- max(date_t, na.rm = T)
p2 <- ggtree(tree2_raw, mrsd=decimal2Date(mrsd)) + theme_tree2()

# prepare_vcf(p2, "../results/delta_global_617_sub")
# system("java -jar ~/softwares/snpEff/snpEff.jar MN908947.3 -fastaProt ../results/delta_glo.snpeff.vcf.faa -csvStats ../results/delta_glo.snpeff.vcf.stats ../results/delta_global_617_sub.vcf > ../results/delta_global_617_sub.snpeff.vcf")
# system("Rscript ./helper/Parse_SnpEff.r ../results/delta_global_617_sub.snpeff.vcf ../results/delta_global_617_sub.snpeff.csv")

tree2 <- tree2_raw@phylo
tip_edge_length <- setNames(tree2$edge.length[sapply(1:length(tree2$tip.label), function(x,y) which (y==x),y=tree2$edge[,2])], tree2$tip.label)
# tip_to_rmv <- names(tip_edge_length[tip_edge_length>0.0005])
# tree2 <- ape::drop.tip(tree2, tip_to_rmv) 
# tree2_raw@phylo <- tree2

df_ann <- read_csv("../results/delta_global_617_sub.snpeff.csv")
df_ann$annot <- NA
df_ann$annot[df_ann$in_gene] <- df_ann$mut_aa[df_ann$in_gene]
df_ann$annot <- gsub("p.", "", df_ann$annot, fixed = T)
df_ann$mut <- gsub("MUT=", "", df_ann$mut, fixed = T)
df_ann$annot[!df_ann$in_gene] <- paste0("nt", df_ann$mut[!df_ann$in_gene])
df_ann_s <- df_ann %>% filter(gene == "S") # only consider spike gene
df_ann_ns <- df_ann %>% filter(!gene == "S") # consider other genes

p2$data$mutations_others <- sapply(p2$data$mutations, function(x){
	if(length(x)==0){return(NA)}
	if(all(x=="")){return(NA)}
	tmp <- sapply(x, function(y){
		mut <- df_ann_ns$annot[df_ann_ns$mut == y]
		check <- df_ann_ns$effect[df_ann_ns$mut == y]
		if(length(check)==0){return(NA)}
		# if(check %in% c("missense_variant", "synonymous_variant")){return(mut)} else {return(NA)} 
		if(check == "missense_variant"){return(mut)} else {return(NA)} # only consider missense_variant
	})
	tmp <- tmp[!grepl("-", tmp)] # ignore deletion
	tmp <- tmp[!is.na(tmp)]
	if(length(tmp)==0){return(NA)}
	paste0(tmp, collapse = " ")
})
p2$data$mutations_spike <- sapply(p2$data$mutations, function(x){
	if(length(x)==0){return(NA)}
	if(all(x=="")){return(NA)}
	tmp <- sapply(x, function(y){
		mut <- df_ann_s$annot[df_ann_s$mut == y]
		check <- df_ann_s$effect[df_ann_s$mut == y]
		if(length(check)==0){return(NA)}
		# if(check %in% c("missense_variant", "synonymous_variant")){return(mut)} else {return(NA)} 
		if(check == "missense_variant"){return(mut)} else {return(NA)} # only consider missense_variant
	})
	tmp <- tmp[!grepl("-", tmp)] # ignore deletion
	tmp <- tmp[!is.na(tmp)]
	if(length(tmp)==0){return(NA)}
	paste0(tmp, collapse = " ")
})

p2$data <- left_join(p2$data, p1$data %>% select(label, metaCluster:label_new, starts_with("mutations_")))

# p1$data %>% filter(grepl("3209V", mutations_others))
# p1$data %>% filter(grepl("2046L", mutations_others))
# p1 <- flip(p1, 156, 178)

p1$data %>% filter(grepl("385K", mutations_others))
p1$data %>% filter(grepl("222V", mutations_spike))
p1 <- flip(p1, 157, 167)

p1$data %>% filter(grepl("95I", mutations_spike))
phangorn::Siblings(tree@phylo, 183)
p1 <- flip(p1, 192, 183)

d1 <- p1$data
d2 <- p2$data

d2$x_ori <- d2$x
d2$y_ori <- d2$y
# scale_t <- 1.5/(max(d2$x, na.rm = T)-min(d2$x, na.rm = T))
# d2$x <- d2$x*scale_t
x_adj <- 0.3
d2$x <- (max(d1$x, na.rm = T) + x_adj) + (d2$x_ori - min(d2$x_ori, na.rm = T))# move to right 
d2$y <- d2$y_ori*(max(d1$y, na.rm = T)/max(d2$y_ori, na.rm = T)) # rescale y

metadata_sub <- read_csv("../results/metadata_global_617_subset.csv")
metadata_sub$label <- metadata_sub$`Accession ID`

d2 <- left_join(d2, metadata_sub %>% select(label, `Pango lineage`))
d2$lineage <- d2$`Pango lineage`
d2 <- d2 %>% select(-`Pango lineage`)

d2$country_glo <- sapply(d2$label, function(x) {
	tmp <- date_meta_glo %>% filter(`Accession ID` == x) %>% .$Location
	if(length(tmp)==0){return(NA)}
	strsplit(tmp, " / ", fixed = T)[[1]][2]
})
d2$country_glo[!d2$country_glo %in% names(colors_countries_t)] <- "Others"
d2$country[is.na(d2$country)] <- d2$country_glo[is.na(d2$country)]

breaks_hk <- seq(Date2decimal("2021-01-01"), Date2decimal("2021-07-06"), 14/365*2)
labels_hk <- decimal2Date(breaks_hk)
labels_hk <- paste0(year(labels_hk), "-", month(labels_hk, label = T), "-", day(labels_hk))

breaks_glo_t <- seq(Date2decimal("2021-01-01")-14/365*2, Date2decimal("2021-07-06"), 14/365*2)
breaks_glo <- (max(d1$x, na.rm = T) + x_adj) + (breaks_glo_t - min(d2$x_ori, na.rm = T))
labels_glo <- decimal2Date(breaks_glo_t)
labels_glo <- paste0(year(labels_glo), "-", month(labels_glo, label = T), "-", day(labels_glo))

breaks_t <- c(breaks_hk, breaks_glo)
labels_t <- c(labels_hk, labels_glo)

pp <- p1 + 
	# xlim(2020,2023) +
	# guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))+
	scale_x_continuous(limits = c(2020.5,2024.2), breaks = breaks_t, labels = labels_t)+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))+
	NULL

dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))

# id_tips_617_1 <- d2$label[which(d2$lineage == "B.1.617.1")]
# id_tips_617_2 <- d2$label[which(d2$lineage == "B.1.617.2")]
# id_tips_617_3 <- d2$label[which(d2$lineage == "B.1.617.3")]
# node_617_1 <- phangorn::mrca.phylo(tree2, id_tips_617_1)
# node_617_2 <- phangorn::mrca.phylo(tree2, id_tips_617_2)
# node_617_3 <- phangorn::mrca.phylo(tree2, id_tips_617_3)

nodes2 <- unique(tree2_raw@phylo$edge[,1])
# dist_node <- RRphylo::distNodes(tree2_raw@phylo, 1, 1)
# check <- sapply(nodes2, function(x) {
# 	# tmp <- phangorn::Descendants(tree2_raw@phylo, x, "tips")
# 	dist_node[names(dist_node[,2])==x,1]>15
# })
check <- sapply(nodes2, function(x) {
	tmp <- phangorn::Children(tree2_raw@phylo, x)
	all(tmp <= length(tree2_raw@phylo$tip.label))
})
end_nodes2 <- nodes2[check]

d1$nudge_x = -0.1
d2$nudge_x = 0.1
dd_label <- bind_rows(d1 %>% filter(!node %in% end_nodes), d2 %>% filter(!node %in% end_nodes2)) %>%
	filter(!is.na(label)) %>%
	filter(!isTip) %>%
	# filter(x<breaks_hk[6] | (x>2022 & x<breaks_glo[1])) %>% 
	select(x, y, nudge_x, starts_with("mutations_")) %>% 
	pivot_longer(starts_with("mutations_")) %>% 
	# filter(nchar(value)<40) %>% 
	unique()

# https://www.medrxiv.org/content/10.1101/2021.07.13.21260273v1
label_india_paper <- c("A1306S", "P2046L", "P2287S", "V2930L", "T3255I", "T3446A", "G5063S", "P5401L", "A6319V", "G215C", "P309L", "A3209V", "V3718A", "G5063S", "P5401L", "L116F", "A3209V", "V3718A", "T3750I", "G5063S", "P5401L", "A222V", "P309L", "D2980N", "F3138S", "K77T")

dd_label$in_ind_paper <- sapply(dd_label$value, function(x) {
	tmp <- strsplit(x, " ", fixed = T)[[1]]
	check <- any(tmp %in% label_india_paper)
	if(!check){return(NA)}
	return(paste(tmp[tmp %in% label_india_paper], collapse = " "))
})
dd_label <- dd_label %>% filter(!is.na(in_ind_paper))
dd_label$value <- dd_label$in_ind_paper

p_co <- pp + geom_line(aes(x, y, group=label, color = country), data= dd %>% filter(!grepl("^EP", label)) %>% filter(isTip), alpha = 0.5, show.legend = F, size = 0.6)+
	geom_tiplab(aes(color = country, label = label_new), size = 2, offset = 0.01, data = d2)+
	geom_tippoint(aes(color = country), size = 0.5, data = d2, na.rm = T) +
	# new_scale_color() +
	geom_tree(data=d2, aes(color = country), size = 0.2, show.legend = F) +
	# scale_color_manual(name = "Lineage", values = c("#7fc97f", "#80b1d3", "#beaed4"), guide = guide_legend(override.aes=aes(size = 4, alpha = 1, shape = 95)), na.value = "black", limits = paste0("B.1.617.", 1:3))+
	new_scale_color() +
	geom_text_repel(aes(label = value, color = name), size = 1.5, nudge_y = -0.3, bg.color = "white", bg.r = 0.1, min.segment.length = 0.003, segment.size= 0.2, data = dd_label)+
	scale_color_manual(name = "Missense variant", limits = c("mutations_spike", "mutations_others"), labels = c("Spike", "Other genes"), values = c("#b2182b", "#2166ac"), guide = guide_legend(override.aes=aes(size = 4, alpha = 1, shape = 15)))+
	# guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1, shape = 15)))+
	# theme(legend.key.width=unit(0.1,"cm"))+
	NULL

data_t <- d2 %>% filter(grepl("^EPI", label))
# data_t %>% filter(lineage == "B.1.617.1") %>% .$y
data_t$lineage[data_t$lineage == "B.1.617.1" & data_t$y > 48] <- NA
# data_t %>% filter(lineage == "B.1.617.2") %>% .$y
data_t$lineage[data_t$lineage == "B.1.617.2" & data_t$y < 20] <- NA

# data_t %>% filter(lineage == "B.1.617.3") %>% .$y
data_t$lineage[data_t$lineage == "B.1.617.3" & data_t$y < 20] <- NA

label_lineage <- function(plot, data, lineage_t, label_t, nudge_x, nudege_x_text, x_pos, color, alpha, size_text) {
	if(is.na(x_pos)){
		x_s <- sort(data %>% filter(lineage == lineage_t) %>% .$x)
		x_pos <- min(x_s) + nudge_x
	}
	y_s <- sort(data %>% filter(lineage == lineage_t) %>% .$y)
	plot + 
	geom_segment(x = x_pos, xend = x_pos, y = min(y_s), yend = max(y_s), color = color, alpha = alpha, size = 0.3) + 
	geom_text(x = x_pos + nudege_x_text, y = (max(y_s) + min(y_s))/2, color = color, label = label_t, angle = 270, size = size_text)
}

x_tail <- 2024.25
p_co_label <- label_lineage(plot = p_co, data = data_t, lineage_t = "B.1.617.1", label_t = "Kappa (B.1.617.1)", x_pos = x_tail-0.05, color = "#a3a3a3", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
p_co_label <- label_lineage(plot = p_co_label, data = data_t, lineage_t = "B.1.617.2", label_t = "Delta (B.1.617.2)", x_pos = x_tail+0.06, color = "#a3a3a3", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
p_co_label <- label_lineage(plot = p_co_label, data = data_t, lineage_t = "B.1.617.3", label_t = "B.1.617.3", x_pos = x_tail-0.15, color = "#a3a3a3", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)

label_clade <- function(plot, data, y1, y2, col_name, label_t, nudge_x, nudege_x_text, x_pos, color, alpha, size_text) {
	# y1_idx <- grepl(y_max, data[[col_name]])
	# y1 <- max(data$y[y1_idx], na.rm = T)
	# y2_idx <-  grepl(y_min, data[[col_name]])
	# y2 <- min(data$y[y2_idx], na.rm = T)

	y_s <- sort(data$y[data$y<=y1 & data$y>=y2])
	plot + 
	geom_segment(x = x_pos, xend = x_pos, y = min(y_s), yend = max(y_s), color = color, alpha = alpha, size = 0.3) + 
	geom_text(x = x_pos + nudege_x_text, y = (max(y_s) + min(y_s))/2, color = color, label = label_t, angle = 270, size = size_text, alpha = alpha, family = "Helvetica")
}

x_tail <- 2022.3
# p1$data %>% filter(isTip) %>% filter(grepl("India_2021-03-27_Sporadic", label_new))
# p1$data %>% filter(isTip) %>% filter(grepl("India_2021-04-17_Cluster_14", label_new))

p_co_label_1 <- label_clade(plot = p_co_label, data = d1, y1=43, y2=2, label_t = "Kappa (B.1.617.1)", x_pos = x_tail-0.32, color = "#a3a3a3", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
p_co_label_1 <- label_clade(plot = p_co_label_1, data = d1, y1=73, y2=61, label_t = "Delta II", x_pos = x_tail, color = "black", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
p_co_label_1 <- label_clade(plot = p_co_label_1, data = d1, y1=60, y2=45, label_t = "Delta III", x_pos = x_tail, color = "black", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
p_co_label_1 <- label_clade(plot = p_co_label_1, data = d1, y1=115, y2=74, label_t = "Delta I", x_pos = x_tail+0.2, color = "black", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
p_co_label_1 <- label_clade(plot = p_co_label_1, data = d1, y1=44, y2=44, label_t = "Delta IV", x_pos = x_tail-0.30, color = "black", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)

x_tail <- 2024.25
y1_t <- d2 %>% filter(isTip) %>% filter(grepl("2415482", label)) %>% .$y
y2_t <- d2 %>% filter(isTip) %>% filter(grepl("2284893", label)) %>% .$y
p_co_label_1 <- label_clade(plot = p_co_label_1, data = d1, y1=115, y2=y2_t, label_t = "Delta I", x_pos = x_tail-0.05, color = "black", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
y1_t <- d2 %>% filter(isTip) %>% filter(grepl("2641979", label)) %>% .$y
y2_t <- d2 %>% filter(isTip) %>% filter(grepl("2577027", label)) %>% .$y
p_co_label_1 <- label_clade(plot = p_co_label_1, data = d1, y1=y1_t, y2=y2_t, label_t = "Delta III", x_pos = x_tail-0.05, color = "black", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
y1_t <- d2 %>% filter(isTip) %>% filter(grepl("2577027", label)) %>% .$y
y2_t <- data_t %>% filter(lineage == "B.1.617.3") %>% .$y %>% max()
p_co_label_1 <- label_clade(plot = p_co_label_1, data = d1, y1=y1_t, y2=y2_t, label_t = "Delta IV", x_pos = x_tail-0.05, color = "black", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)
y1_t <- d2 %>% filter(isTip) %>% filter(grepl("2652216", label)) %>% .$y
y2_t <- d2 %>% filter(isTip) %>% filter(grepl("2155363", label)) %>% .$y
p_co_label_1 <- label_clade(plot = p_co_label_1, data = d1, y1=y1_t, y2=y2_t, label_t = "Delta II", x_pos = x_tail-0.05, color = "black", alpha = 0.8, nudege_x_text = 0.05, size_text = 2.5)

# annotate local cases
p_co_label_local <- p_co_label_1 +
	geom_label_repel(label = "Local case A", 
		color = "#a3a3a3", size = 3, nudge_x = -0.7, nudge_y = 5, bg.color = "white", bg.r = 0.1, min.segment.length = 0.003, segment.size= 0.5,
		# arrow = arrow(length = unit(0.01, "npc")),
		data = dd %>% filter(country == "Local" & isTip & x<2023 & grepl("Cluster_32", label_new)))+
	geom_label_repel(label = "Local case C", 
		color = "#a3a3a3", size = 3, nudge_x = -0.7, nudge_y = 5, bg.color = "white", bg.r = 0.1, min.segment.length = 0.003, segment.size= 0.5,
		# arrow = arrow(length = unit(0.01, "npc")),
		data = dd %>% filter(country == "Local" & isTip & x<2023 & grepl("Cluster_27", label_new)))

ggsave("../results/global_Imported_617_cotree.pdf", width = 10, height = 10*sqrt(2), plot = p_co_label_local)

# comapre collection date 
d2 %>% filter(!isTip) %>% filter(grepl("L116F", mutations_others))
tip_b1 <- tree2$tip.label[phangorn::Descendants(tree2, 1981)[[1]]]
date_b1 <- metadata_sub %>% filter(`Accession ID` %in% tip_b1)

d2 %>% filter(!isTip) %>% filter(grepl("3750I", mutations_others))
tip_b2 <- tree2$tip.label[phangorn::Descendants(tree2, 1784)[[1]]]
date_b2 <- metadata_sub %>% filter(`Accession ID` %in% tip_b2)

d2 %>% filter(!isTip) %>% filter(grepl("3646A", mutations_others))
tip_a <- tree2$tip.label[phangorn::Descendants(tree2, 2229)[[1]]]
date_a <- metadata_sub %>% filter(`Accession ID` %in% tip_a)

date_b <- bind_rows(date_b1, date_b2)

decimal2Date(median(Date2decimal(date_a$`Collection date`), na.rm = T))
decimal2Date(median(Date2decimal(date_b$`Collection date`), na.rm = T))
decimal2Date(median(Date2decimal(date_b1$`Collection date`), na.rm = T))
decimal2Date(median(Date2decimal(date_b2$`Collection date`), na.rm = T))


wilcox.test(Date2decimal(date_a$`Collection date`), Date2decimal(date_b$`Collection date`))
wilcox.test(Date2decimal(date_b1$`Collection date`), Date2decimal(date_b2$`Collection date`))
wilcox.test(Date2decimal(date_b1$`Collection date`), Date2decimal(date_b2$`Collection date`))
wilcox.test(Date2decimal(date_b1$`Collection date`), Date2decimal(date_a$`Collection date`))
