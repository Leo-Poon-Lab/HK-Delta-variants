library(tidyverse)
library(Biostrings)
library(ggplot2)
library(scico)
library(colorspace)
library(ggtree)
library(phangorn)
library(patchwork)
library(lubridate)

data_meta <- read_csv("../../2021-06-24_merge_metadata/results/cleaned_metadata.csv", guess_max = 20000)
data_delta <- readxl::read_excel("../results/metadata_hk_delta.xlsx")
seqs_hk_617 <- readDNAStringSet("../results/seqs_hk_delta.fasta")
df_lin <- read_csv("../results/lineage.csv")

tail(data_meta)
data_meta %>% filter(`Report date`>= ymd("2021-03-27")) %>% .$Classification %>% table()
t1 <- data_meta %>% filter(`Report date`>= ymd("2021-03-27")) %>% .$Classification %>% grepl(pattern = "ported") %>% sum()
t2 <- data_meta %>% filter(`Report date`>= ymd("2021-03-27")) %>% .$Classification %>% grepl(pattern = "ported") %>% length()
round((t1/t2)*100, 2) # HK adapted an elimination strategy to control COVID-19, and 85.3% (N=441) of all RT-PCR positive COVID-19 cases were imported cases 

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

data_meta_imported <- data_meta %>% filter(grepl("mport", Classification))
table(data_meta$Classification)
table(data_meta_imported$Classification)

data_meta_imported$sequenced <- data_meta_imported$case_id %in% df_lin$case_id
data_meta_imported$lineage_raw <- data_meta_imported$lineage
data_meta_imported$lineage[!data_meta_imported$sequenced] <- "Not sequenced"

# we sequenced 48.53% (N=214) of all imported/imported-related cases with Illumina NovaSeq/iSeq platform (Appendix Methods), and 53.27% (N=114) of them are Delta or Kappa variants.
num_seqed <- data_meta_imported %>% filter(`Report date`>= ymd("2021-03-27") & `Report date`<= ymd("2021-07-16")) %>% .$sequenced %>% table()
num_seqed
round(num_seqed[names(num_seqed)=="TRUE"]/sum(num_seqed)*100,2) # proportion of sequenced cases after 2021-03-27
num_lin <- data_meta_imported %>% filter(`Report date`>= ymd("2021-03-27") & `Report date`<= ymd("2021-07-16")) %>% filter(sequenced) %>% .$lineage %>% table()
num_lin/sum(num_lin)
sum(num_lin[grep("B.1.617", names(num_lin))])
round(sum(num_lin[grep("B.1.617", names(num_lin))])/sum(num_lin)*100, 2)

lineages_minor <- names(table(data_meta_imported$lineage)[table(data_meta_imported$lineage)<15])
data_meta_imported$lineage[data_meta_imported$lineage %in% lineages_minor] <- "Other lineages (N < 15)"

tmp <- table(data_meta_imported$lineage)

colors <- scico(length(tmp)-2, palette = "lapaz", direction = -1)
colors <- c(colors, "#636363", "#c7c7c7")
names(colors) <- names(tmp)
tmp <- tmp[order(tmp, decreasing = T)]
tmp <- tmp[c(3:length(tmp), 2, 1)]
labels <- paste0(names(tmp), " (N = ", tmp, ")")
labels[length(tmp)-1] <- "Other lineages (N < 15)"
data_meta_imported$lineage <- factor(data_meta_imported$lineage, levels = names(tmp), labels = labels)
colors <- colors[match(names(tmp), names(colors))]
names(colors) <- labels

color_lin <- c("#b2df8a", "#a6cee3", "#fb9a99")
names(color_lin) <- c("B.1.617.1", "B.1.617.2", "B.1.617.3")
colors[grepl("617.1", names(colors))] <- color_lin[["B.1.617.1"]]
colors[grepl("617.2", names(colors))] <- color_lin[["B.1.617.2"]]

# incidence of imported cases in 2021 colored by lineages
p_1 <- data_meta_imported %>% filter(`Report date` >= "2021-01-01") %>% ggplot()+
	geom_histogram(aes(x = `Report date`, fill = lineage), binwidth = 7)+
	scale_x_date(date_breaks = "1 week", date_labels = "%b-%d", expand = c(0, 0), limits = c(lubridate::ymd("2021-01-01"), NA))+
	theme_minimal()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	scale_fill_manual(name = "Lineage", values = colors)+
	# scale_color_manual(values = colors)+
	ylab("No. of imported cases")+
	NULL
ggsave("../results/Imported_2021_lineage.pdf", width = 10, height = 6, plot = p_1)


# Inbound cases by countries
df_plot <- data_meta_imported %>% filter(`Report date` >= "2021-01-01" & sequenced) %>% filter(grepl("617", lineage)
) %>% filter(Classification=="Imported")
df_plot$id <- as.character(seq_len(nrow(df_plot)))
countries_t <- sort(unique(df_plot$`Country of importation`))
colors_countries_t <- sequential_hcl(length(countries_t), "Batlow", rev = T)
names(colors_countries_t) <- countries_t
labels_colors_countries_t <- sapply(names(colors_countries_t), function(x) {
	cnt <- sum(df_plot$`Country of importation` == x)
	paste0(x, " (N=", cnt, ")")
})

P_617_country <- df_plot %>% ggplot()+
	geom_histogram(aes(x = `Report date`, fill = `Country of importation`, color = id), binwidth = 1)+
	# scale_fill_discrete_qualitative(name = "Country",palette = "Dark 2")+
	scale_fill_manual(name = "Country", values = colors_countries_t, guide = "none", labels = labels_colors_countries_t)+
	scale_color_manual(values = rep("white", nrow(df_plot)), guide = "none")+
	scale_x_date(date_breaks = "1 week", date_labels = "%b-%d", expand = c(0, 0), limits = c(ymd("2021-03-25"), ymd("2021-07-16")), guide = guide_axis(n.dodge = 2))+
	theme_minimal()+
	# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	ylab("No. of imported cases")+	
	NULL
	
# annotate control measures
# https://www.news.gov.hk/eng/2021/04/20210415/20210415_172942_093.html
# https://www.news.gov.hk/chi/2021/04/20210401/20210401_182214_863.html
# https://www.info.gov.hk/gia/general/202104/18/P2021041800768.htm?fontSize=1
place_suspension <- c("India", "Pakistan", "Philippines", "Nepal", "Indonesia", "United Kingdom")
date_start <- ymd(c("2021-04-20", "2021-04-20", "2021-04-20", "2021-05-01", "2021-06-25", "2021-07-01"))
df_sus <- tibble(place_suspension = place_suspension, date_start = date_start)
df_sus$y <- seq_len(nrow(df_sus))
df_sus$y_end <- seq_len(nrow(df_sus))
df_sus$date_stop <- ymd("2021-07-16")

color_t <- colors_countries_t
color_t_lack <- rep(NA, sum(!df_sus$place_suspension %in% names(color_t)))
color_t_lack <- sequential_hcl(length(color_t_lack), "Hawaii", rev = T)
names(color_t_lack) <- df_sus$place_suspension[!df_sus$place_suspension %in% names(color_t)]
color_t <- c(color_t, color_t_lack)

save(color_t, colors_countries_t, color_lin, file = "./colors.rdata")

p_suspension <- ggplot(df_sus)+
	geom_segment(aes(x = date_start, xend = date_stop, y = y, yend = y_end, color = place_suspension), size = 2)+
	scale_x_date(date_breaks = "1 week", date_labels = "%b-%d", expand = c(0, 0), limits = c(ymd("2021-03-25"), ymd("2021-07-16")))+
	scale_color_manual(name = "Country", values = color_t, labels = labels_colors_countries_t)+
	theme_minimal()+
	theme(axis.title.x=element_blank(),
		# axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())+
	ylab("Suspension")+
	NULL

color_lin <- c("#b2df8a", "#a6cee3")
names(color_lin) <- c("B.1.617.1", "B.1.617.2")
labels_color_lin <- sapply(names(color_lin), function(x) {
	cnt <- sum(df_plot$lineage_raw == x)
	paste0(x, " (N=", cnt, ")")
})
labels_color_lin <- gsub("B.1.617.1", "Kappa", labels_color_lin)
labels_color_lin <- gsub("B.1.617.2", "Delta", labels_color_lin)

P_617_lineage <- df_plot %>% ggplot()+
	geom_bar(aes(x = `Report date`, color = id, fill = lineage_raw), width = 1)+
	scale_color_manual(values = rep("white", nrow(df_plot)), guide = "none")+
	scale_fill_manual(name = "Variant", values = color_lin, labels = labels_color_lin)+
	scale_x_date(date_breaks = "1 week", date_labels = "%b-%d", expand = c(0, 0), limits = c(ymd("2021-03-25"), ymd("2021-07-16")), guide = guide_axis(n.dodge = 2))+
	theme_minimal()+
	# theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))+
	ylab("No. of imported cases")+
	theme(axis.title.x=element_blank(),
		# axis.title.y=element_blank(),
        # axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        # axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())+
	NULL

P_617_cl <- p_suspension + P_617_lineage + P_617_country + plot_layout(height = c(0.5, 2, 2), ncol = 1, guides = "collect") + plot_annotation(tag_levels = 'A')
ggsave("../results/Imported_617_hist.pdf", width = 10, height = 8, plot = P_617_cl)

names(df_plot)
df_plot %>% group_by(lineage_raw) %>% summarise(num = n(), n_asymp = sum(`s/s`=="asymptomatic"), n_det_quarant = sum(!is.na(`Home/\r\nhotel Confinee`) | !is.na(`from Quarantine\r\ncentre/ camp`)))
df_plot %>% group_by(lineage_raw) %>% summarise(num = n(), prop_asymp = round(sum(`s/s`=="asymptomatic")/num,2), prop_det_quarant = round(sum(!is.na(`Home/\r\nhotel Confinee`) | !is.na(`from Quarantine\r\ncentre/ camp`))/num,2))

m1 <- matrix(c(c(70*0.53, 70*(1-0.53)), c(886*0.8, 886*(1-0.8))), ncol=2) # Delta
chisq.test(m1)
m2 <- matrix(c(c(42*0.81, 42*(1-0.81)), c(886*0.8, 886*(1-0.8))), ncol=2) # Kappa
chisq.test(m2)

