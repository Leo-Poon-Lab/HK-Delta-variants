library(tidyverse)
library(Biostrings)
library(ggplot2)
library(scico)
library(colorspace)
library(ggtree)
library(ggrepel)
library(ggnewscale)
library(phangorn)
library(patchwork)
library(lubridate)
library(writexl)
library(RColorBrewer)

data_meta <- read_csv("../../2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 20000)
data_delta <- readxl::read_excel("../results/metadata_hk_delta.xlsx")
seqs_hk_617 <- readDNAStringSet("../results/seqs_hk_delta.fasta")
df_lin <- read_csv("../results/lineage.csv")

tail(data_meta)
data_meta %>% filter(`Report date`>= ymd("2021-03-27")) %>% .$Classification %>% table()
data_meta %>% filter(`Report date`>= ymd("2021-03-27")) %>% summarise(n =n(), n_imported = sum(Classification == "Imported"), p_imported = n_imported/n*100)
# HK adapted an elimination strategy to control COVID-19, and 83.8% (N=433) of all RT-PCR positive COVID-19 cases were imported cases 

df_lin <- df_lin %>% filter(lineage != "None")
df_lin$case_id <- sapply(df_lin$taxon, function(x){
	tmp <- strsplit(x, "_")[[1]]
	tmp <- tmp[length(tmp)]
})
df_lin$case_id <- as.numeric(df_lin$case_id)
df_lin <- df_lin %>% filter(!is.na(case_id))
(df_lin_dup <- df_lin %>% select(case_id, lineage) %>% unique() %>% filter(case_id %in% case_id[duplicated(case_id)])) # WHP4602_11751 potential problematic sample
df_lin <- df_lin %>% filter(taxon != "WHP4602_11751")
(df_lin_dup <- df_lin %>% select(case_id, lineage) %>% unique() %>% filter(case_id %in% case_id[duplicated(case_id)]))
df_lin <- df_lin %>% filter(!(case_id %in% df_lin_dup$case_id & grepl("iseq", taxon)))

data_meta <- left_join(data_meta, df_lin %>% select(case_id, lineage) %>% unique())

data_meta_imported <- data_meta %>% filter(Classification == "Imported")
table(data_meta$Classification)
table(data_meta_imported$Classification)

data_meta_imported$sequenced <- data_meta_imported$case_id %in% df_lin$case_id
data_meta_imported$lineage_raw <- data_meta_imported$lineage
data_meta_imported$lineage[!data_meta_imported$sequenced] <- "Not sequenced"

data_meta_imported %>% filter(`Report date`>= ymd("2021-03-27") & `Report date`<= ymd("2021-07-16")) %>% mutate(Month = month(`Report date`, label = T, abbr = F), Lineage = lineage_raw) %>% filter(!is.na(Lineage)) %>% group_by(Month, Lineage) %>% summarize(N=n()) %>% arrange(Month, desc(N)) %>% write_xlsx("../results/lineage_month.xlsx")

# we sequenced 48.96% (N=212) of all imported cases with Illumina NovaSeq/iSeq platform (Appendix Methods), and 53.27% (N=114) of them are Delta or Kappa variants.
(num_seqed <- data_meta_imported %>% filter(`Report date`>= ymd("2021-03-27") & `Report date`<= ymd("2021-07-16")) %>% .$sequenced %>% table())
round(num_seqed[names(num_seqed)=="TRUE"]/sum(num_seqed)*100,2) # proportion of sequenced cases after 2021-03-27
num_lin <- data_meta_imported %>% filter(`Report date`>= ymd("2021-03-27") & `Report date`<= ymd("2021-07-16")) %>% filter(sequenced) %>% .$lineage %>% table()
num_lin/sum(num_lin)
sum(num_lin[grep("B.1.617", names(num_lin))])
round(sum(num_lin[grep("B.1.617", names(num_lin))])/sum(num_lin)*100, 2)

# lineages_minor <- names(table(data_meta_imported$lineage)[table(data_meta_imported$lineage)<15])
lineages_minor <- names(table(data_meta_imported$lineage)) 
lineages_minor <- lineages_minor[!lineages_minor %in% c("Not sequenced", "B.1.617.1", "B.1.617.2")]
data_meta_imported$lineage[data_meta_imported$lineage %in% lineages_minor] <- "Other lineages"

tmp <- table(data_meta_imported %>% filter(`Report date`>= ymd("2021-03-27")) %>% .$lineage)

colors <- scico(length(tmp)-2, palette = "lapaz", direction = -1)
colors <- c(colors, "#c7c7c7", "#636363")
names(colors) <- names(tmp)
tmp <- tmp[order(tmp, decreasing = T)]
tmp <- tmp[c(3:length(tmp), 2, 1)]
labels <- paste0(names(tmp), " (N=", tmp, ")")
# labels[length(tmp)-1] <- "Other lineages (N < 15)"
data_meta_imported$lineage <- factor(data_meta_imported$lineage, levels = names(tmp), labels = labels)
colors <- colors[match(names(tmp), names(colors))]
names(colors) <- labels

color_lin <- c("#b2df8a", "#a6cee3", "#fb9a99")
names(color_lin) <- c("B.1.617.1", "B.1.617.2", "B.1.617.3")
colors[grepl("617.1", names(colors))] <- color_lin[["B.1.617.1"]]
colors[grepl("617.2", names(colors))] <- color_lin[["B.1.617.2"]]
names(colors) <- gsub("B.1.617.1", "Kappa", names(colors))
names(colors) <- gsub("B.1.617.2", "Delta", names(colors))

data_meta_imported$lineage_voc <- as.character(data_meta_imported$lineage)
data_meta_imported$lineage_voc <- gsub("B.1.617.1", "Kappa", data_meta_imported$lineage_voc)
data_meta_imported$lineage_voc <- gsub("B.1.617.2", "Delta", data_meta_imported$lineage_voc)

# incidence of imported cases in 2021 colored by lineages
p_1 <- data_meta_imported %>% filter(`Report date` >= "2021-03-27") %>% ggplot()+
	geom_histogram(aes(x = `Report date`, fill = lineage_voc), binwidth = 1)+
	scale_x_date(date_breaks = "1 week", date_labels = "%b-%d", expand = c(0, 0), limits = c(ymd("2021-03-26"), NA), guide = guide_axis(n.dodge = 2))+
	theme_minimal()+
	# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	scale_fill_manual(name = "Lineage", values = colors)+
	# scale_color_manual(values = colors)+
	ylab("No. of imported cases")+
	NULL
ggsave("../results/Imported_2021_lineage.pdf", width = 10, height = 6, plot = p_1)


# Inbound cases by countries
df_plot <- data_meta_imported %>% filter(`Report date` >= "2021-01-01" & sequenced) %>% filter(grepl("617", lineage)
) %>% filter(Classification=="Imported")
df_plot$id <- as.character(seq_len(nrow(df_plot)))
countries_t <- names(sort(table(df_plot$`Country of importation`), decreasing = T))

# colors_countries_t <- c("#76c763",
# "#7c4db9",
# "#bead5a",
# "#c05b91",
# "#c3573b",
# "#8eafbd",
# "#4d403f")
colors_countries_t <- c("#b65a49","#738488","#99c161","#8f50ab")
colors_countries_t <- c(colors_countries_t, "#e41a1c", "#377eb8", "#ff7f00")
colors_countries_t <- c(colors_countries_t, "#b3e2cd", "#fdcdac", "#fabed4", "#fff2ae")
names(colors_countries_t) <- countries_t
labels_colors_countries_t <- sapply(names(colors_countries_t), function(x) {
	cnt <- sum(df_plot$`Country of importation` == x)
	paste0(x, " (N=", cnt, ")")
})

save(colors_countries_t, countries_t, file="./colors.rdata")

# annotate control measures
# https://www.news.gov.hk/eng/2021/04/20210415/20210415_172942_093.html
# https://www.news.gov.hk/chi/2021/04/20210401/20210401_182214_863.html
# https://www.info.gov.hk/gia/general/202104/18/P2021041800768.htm?fontSize=1
#  https://travelbans.org/asia/hong-kong/
date_start <- ymd(c("2021-04-20", "2021-05-01", "2021-06-25", "2021-07-01", "2021-05-26"))
events <- c("India", "Nepal", "Indonesia", "The UK ban", "The UK resume")
color_country <- colors_countries_t[c(1,2,3,4,4)]
ploicy <- c(rep("Ban", 4), "Resume")
df_sus <- tibble(`Report date`=date_start, events=events, color_country=color_country, `Flight policy`=ploicy, lineage_voc="Delta (N=70)")
df_sus$y <- 10+c(1:4,4)/2

df_plot_sus <- bind_rows(df_plot, df_sus)

P_617_country <- df_plot_sus %>% filter(!is.na(case_id)) %>% ggplot()+
	geom_histogram(aes(x = `Report date`, fill = `Country of importation`, color = id), binwidth = 1)+
	# scale_fill_discrete_qualitative(name = "Country",palette = "Dark 2")+
	scale_fill_manual(name = "Country", values = colors_countries_t, labels = labels_colors_countries_t)+
	scale_color_manual(values = rep("white", nrow(df_plot)), guide = "none")+
	scale_x_date(date_breaks = "1 week", date_labels = "%b-%d", expand = c(0, 0), limits = c(ymd("2021-03-25"), ymd("2021-07-16")), guide = guide_axis(n.dodge = 2))+
	# theme_minimal()+
	theme_bw()+
	scale_y_continuous(breaks = seq(0,15,5))+
	facet_wrap(vars(lineage_voc), ncol = 1)+
	# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	ylab("No. of imported cases")+	
	NULL

P_617_country_sus <- P_617_country + 
	new_scale_color()+
	geom_segment(aes(x = `Report date`, y = y, yend = y, color = color_country),
		xend = ymd("2021-07-16"),
		alpha = 0.9, 
		size = 2, data = df_sus %>% filter(`Flight policy`=="Ban"))+
	geom_segment(aes(x = `Report date`, y = y, yend = y, color = color_country),
		xend = ymd("2020-12-20"),
		alpha = 0.9, 
		size = 2, data = df_sus %>% filter(`Flight policy`!="Ban"))+
	# geom_segment(aes(x = `Report date`, xend = `Report date`, y = 0, yend = y, linetype=`Flight policy`, 
	# 	color = color_country),
	# 	# alpha = 0.8, 
	# 	size = 0.5, data = df_sus %>% select(-lineage_voc))+
	# scale_linetype_manual(name = , values = c(2,3))+
	geom_text_repel(aes(x = `Report date`, label = events, y = y, color = color_country),
		# alpha = 0.8,
		# nudge_x = -2,
		nudge_y = 15-df_sus$y,
		# hjust = 1,
		min.segment.length = 0.003, bg.color = "white", bg.r = 0.1, size = 3,
		# segment.curvature = 0.5,
		arrow = arrow(length = unit(0.015, "npc")),
		data = df_sus)+
	scale_color_identity()+
	theme(legend.position = c(0.8, 0.24), legend.background = element_rect(fill = "transparent", colour = "black"))+
	NULL

ggsave("../results/Imported_617_hist.pdf", width = 10, height = 8, plot = P_617_country_sus)


# chisq test
names(df_plot)
df_plot %>% group_by(lineage_raw) %>% summarise(num = n(), n_asymp = sum(`s/s`=="asymptomatic"), n_det_quarant = sum(!is.na(`Home/\r\nhotel Confinee`) | !is.na(`from Quarantine\r\ncentre/ camp`)))
df_plot %>% group_by(lineage_raw) %>% summarise(num = n(), prop_asymp = round(sum(`s/s`=="asymptomatic")/num,2), prop_det_quarant = round(sum(!is.na(`Home/\r\nhotel Confinee`) | !is.na(`from Quarantine\r\ncentre/ camp`))/num,2))

m1 <- matrix(c(c(70*0.53, 70*(1-0.53)), c(886*0.8, 886*(1-0.8))), ncol=2) # Delta
chisq.test(m1)
m2 <- matrix(c(c(42*0.81, 42*(1-0.81)), c(886*0.8, 886*(1-0.8))), ncol=2) # Kappa
chisq.test(m2)

