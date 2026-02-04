#Library packages ####
#install.packages("reshape2")
library(reshape2)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("cowplot")
library(cowplot)
#install.packages("patchwork")
library(patchwork)
#install.packages("ggrepel")
library(ggrepel)
#install.packages("readxl")
library(readxl)
#install.packages("sf")
library(sf)
#install.packages("stringr")
library(stringr)
#install.packages("gt")
library(gt)
#install.packages("tidyr")
library(tidyr)
#install.packages("vegan")
library(vegan)
#install.packages("dplyr")
library(dplyr)
#install.packages("rnaturalearth")
library(rnaturalearth)
#install.packages("openxlsx")
library(openxlsx)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("ade4")
library(ade4)

# FIGURE 1 MPB DIV ####

  # Panel A, MAP

lat <- 47.875
lon <- -3.918

crs <- paste0("+proj=laea +lat_0=", lat, " +lon_0=", lon, " +datum=WGS84 +units=m +no_defs")

ctrys50m <- ne_countries(scale = 50, type = "countries", returnclass = "sf") %>%
  select(iso_a3, iso_n3, admin) %>%
  st_transform(crs = crs)

sphere <- st_graticule(ndiscr = 10000, margin = 10e-6) %>%
  st_transform(crs = crs) %>%
  st_convex_hull() %>%
  summarise(geometry = st_union(geometry))

concarneau <- data.frame(lon = lon, lat = lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = crs)

coords <- st_coordinates(concarneau)

marge <- 1e6 # Map limits 

xlim <- c(coords[1] - 7e5, coords[1] + 8.5e5)
ylim <- c(coords[2] - marge, coords[2] + marge)

ggplot() +
  geom_sf(data = sphere, fill = "#e2ebf7", alpha = 0.7) +
  geom_sf(data = ctrys50m, fill = "grey99", color = "black") + 
  geom_sf(data = concarneau, color = "white", size = 4) +
  coord_sf(crs = crs, xlim = xlim, ylim = ylim) +
  theme_test(base_family = "Arial") +
  theme(
    panel.grid.major = element_line(color = "black", size = 0.5)
  )+
  theme(
    axis.text.x = element_text(color = "black", face = "bold", size = 14),
    axis.text.y = element_text(color = "black", face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5))+
  plot_annotation(
    title = "A",
    caption = "") & theme(plot.caption = element_text(hjust = 0.5, size = 16, face = "bold"),
                                    plot.title = element_text(face="bold", size="20"))

  # Panel D, barplot

df_18 <- psmelt(final_biofilm18S)

# Phylum group
df_agg_18 <- df_18 %>%
  group_by(Sample, Conditions, Biofilm_type, Class) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop')

# Abundance means
df_agg_18 <- df_agg_18 %>%
  group_by(Conditions, Class) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = 'drop')

# Relative abundances
df_agg_18 <- df_agg_18 %>%
  group_by(Class, Conditions) %>%
  summarise(Mean_Abundance = sum(Mean_Abundance), .groups = 'drop') %>%
  group_by(Conditions) %>%
  mutate(Relative_Abundance = Mean_Abundance / sum(Mean_Abundance)*100)

# Threshold
df_agg_18 <- df_agg_18 %>%
  mutate(Class = ifelse(Relative_Abundance < 2, "Others", Class)) 

# Pool 'Others'$Conditions
df_agg_18 <- df_agg_18 %>%
  group_by(Conditions, Class) %>%
  summarise(
    Mean_Abundance = sum(Mean_Abundance, na.rm = TRUE),
    Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
    .groups = "drop")

# Cleaning names
df_agg_18$Class <- gsub("\\s+", "", df_agg_18$Class)
df_agg_18 <- df_agg_18 %>%rename(Phylum = Class)

# df subset
df_WSB_18 <- subset(df_agg_18, grepl("^WSB_", Conditions))
df_SB_18 <- subset(df_agg_18, grepl("^SB_", Conditions))

# Colors settings
all_Phylum_18 <- union(unique(df_SB_18$Phylum), unique(df_WSB_18$Phylum))
color_Phylum_18 <- c(
  "Annelida" = "#1a759f",
  "Arthropoda" = "#cae9ff",
  "Bacillariophyceae" = "#4ba898",
  "Dinophyceae" = "#99CC87",
  "Gregarinomorphea" = "grey55",
  "Mediophyceae" = "#257939",
  "Nematoda" = "#e5d4e7",
  "Platyhelminthes" = "#986eac",
  "IncertaeSedis" = "white",
  "Others" = "black",
  "Endomyxa" = "grey30"
)

Phylum_colors_18 <- data.frame(Phylum = all_Phylum_18, Color = color_Phylum_18)

# Phylum names
labels_italics <- setNames(
  lapply(all_Phylum_18, function(x) {
    if(x == "Others") {
      x  
    } else {
      bquote(italic(.(x)))  
    }
  }),
  all_Phylum_18
)

# Plot SB
plot_SB_phylum_18 <- ggplot(df_SB_18, aes(x = Conditions, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.8) +
  theme_test(base_family = "Arial") +
  labs(title = "Sediment biofilm", x = "", y = "Relative abundance (%)") +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = setNames(Phylum_colors_18$Color, Phylum_colors_18$Phylum),
                    labels = labels_italics)+
  theme(plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.ticks = element_line(color = "black"),
        legend.position = "bottom")+
  annotate("text", 
           x = 1, y = -4, 
           label = "DA", 
           size = 6, 
           colour = "black",
           family = "Arial")+
  annotate("text", 
           x = 2, y = -4, 
           label = "HL", 
           size = 6, 
           colour = "black",
           family = "Arial")+
  annotate("text", 
           x = 3, y = -4, 
           label = "HP", 
           size = 6, 
           colour = "black",
           family = "Arial")+
  annotate("text", 
           x = 4, y = -4, 
           label = "PL", 
           size = 6, 
           colour = "black",
           family = "Arial")

# PLot WSB
plot_WSB_phylum_18 <- ggplot(df_WSB_18, aes(x = Conditions, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.8) +
  theme_test(base_family = "Arial") +
  labs(title = "Without sediment biofilm", x = "", y = "") +
  theme(plot.title = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = setNames(Phylum_colors_18$Color, Phylum_colors_18$Phylum))+
  annotate("text", 
           x = 1, y = -4, 
           label = "DA", 
           size = 6, 
           colour = "black",
           family = "Arial")+
  annotate("text", 
           x = 2, y = -4, 
           label = "HL", 
           size = 6, 
           colour = "black",
           family = "Arial")+
  annotate("text", 
           x = 3, y = -4, 
           label = "HP", 
           size = 6, 
           colour = "black",
           family = "Arial")+
  annotate("text", 
           x = 4, y = -4, 
           label = "PL", 
           size = 6, 
           colour = "black",
           family = "Arial")

(plot_SB_phylum_18 + plot_WSB_phylum_18) +
  plot_annotation(
    title = "D",
    caption = "Treatments") & theme(plot.caption = element_text(hjust = 0.5, size = 16))

  # Pie chart (E)

df_WSB_18_pc <- df_WSB_18 %>%
  group_by(Phylum) %>%
  summarise(sum = sum(Mean_Abundance), .groups = 'drop') %>% 
  mutate(mean = 100 * sum / sum(sum))

# Label position 
df_WSB_18_pc <- df_WSB_18_pc %>%
  arrange(desc(Phylum)) %>%
  mutate(pos = cumsum(mean) - mean / 2) 

# Plot 1
ggplot(df_WSB_18_pc, aes(x = "", y = mean, fill = Phylum)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.9, position = "stack") +
  coord_polar(theta = "y", direction = -1) +  
  geom_label_repel(aes(label = paste0(round(mean, 2), " %"), y = pos),
                   color = "white", size = 10, fontface = "bold", 
                   family = "Arial",  
                   nudge_x = 0.73,              
                   force = 0.6,                 
                   box.padding = 0.2,
                   show.legend = FALSE,
                   segment.color = "black") +
  theme_void(base_family = "Arial") +
  labs(title = "E") +
  theme(plot.title = element_text(face = "bold", size = 16),
                                  legend.position = "none") +
  scale_fill_manual(values = setNames(Phylum_colors_18$Color, Phylum_colors_18$Phylum))

  # Pleurosigma domination, Pie chart n°2

dfeuk <-df_18

# Cleanning df
dfeuk[, 10:18] <- lapply(dfeuk[, 10:18], function(x) gsub("\\s+", "", x))

# Bacillariophyceae
dfeuk <- dfeuk %>% filter(Class %in% "Bacillariophyceae")

# Unidentified Bacillariophyceae correction (blastn ncbi)
correction_file <- "1_Data_formatting/data/Metabarcoding_data/unidentified_bacill_correction_blastn.xlsx"
correction_file <- read_excel(correction_file)
correction_file <- correction_file[, 1:2] # retrieve col 1 and 2
correction_file <- na.omit(correction_file) # unidentified taxon suppression

# Df formatting
correction_file <- correction_file %>%
  separate(
    blastn_id_order, 
    into = c("Order", "Family", "Genus", "Species"),
    sep = ";",
    remove = TRUE)

# Merge & correction
dfeuk <- dfeuk %>%
  left_join(correction_file, by = c("OTU" = "ID_ASV"))

dfeuk <- dfeuk %>%
  mutate(
    Order = coalesce(Order.y, Order.x),  
    Family = coalesce(Family.y, Family.x),
    Genus = coalesce(Genus.y, Genus.x),
    Species = coalesce(Species.y, Species.x)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# Subset
dfeuk <- subset(dfeuk, grepl("^WSB_", Conditions))

# Species
dfeuk_sp <- dfeuk[, c("Species", "Abundance")]
dfeuk_sp$Species <- sub("_n.*", "", dfeuk_sp$Species) # Tag suppr

# Calcul by species and % transformation
dfeuk_sp <- dfeuk_sp %>%
  group_by(Species) %>%
  summarise(sum = sum(Abundance), .groups = 'drop') %>%
  mutate(mean = 100 * sum / sum(sum))

# Name format
dfeuk_sp <- dfeuk_sp %>% mutate(Species = gsub("_", " ", Species))

# Threshold 
dfeuk_sp <- dfeuk_sp %>%mutate(Species = ifelse(mean < 2, "Others", Species))

# merge other
dfeuk_sp <- dfeuk_sp %>%
  group_by(Species) %>%
  summarise(sum = sum(mean), .groups = 'drop') %>%
  mutate(mean = 100 * sum / sum(sum))

# Color settings
ord_sp_bacill <- union(unique(dfeuk_sp$Species), unique(dfeuk_sp$Species))
col_bacill <- c(
  "IncertaeSedis" = "white", 
  "Other" = "#0b090a", 
  "Pleurosigma strigosum" = "#31572c")
ord_bacill <- data.frame(Species = ord_sp_bacill, Color = col_bacill) 

# Label position
dfeuk_sp <- dfeuk_sp %>%
  arrange(desc(Species)) %>%
  mutate(pos = cumsum(mean) - mean / 2)

# Species names
labels_italics <- setNames(
  lapply(ord_bacill$Species, function(x) bquote(italic(.(x)))),
  ord_bacill$Species)

ggplot(dfeuk_sp, aes(x = "", y = mean, fill = Species)) +
  geom_bar(stat = "identity", width = 1, position = "stack", color = "black", size = 1) +
  coord_polar(theta = "y", start = pi) +
  geom_label_repel(aes(label = paste0(round(mean, 2), " %"), y = pos),
                   color = "white", size = 10, fontface = "bold",
                   nudge_x = 0.7,              
                   force = 0.5,                 
                   box.padding = 0.2,
                   show.legend = FALSE,
                   segment.color = "black",
                   family = "Arial") +
  theme_void() +
  labs(title = "") +
  theme(legend.position = "right",
        legend.text = element_text(family = "Arial"),
        legend.title = element_text(family = "Arial")) +
  scale_fill_manual(values = setNames(ord_bacill$Color, ord_bacill$Species),
                    labels = labels_italics)

# FIGURE 2 MIGRATION ####

y_limits <- c(2000, 18900)

# Plot HL
plot_hl <- ggplot(File_HL_data, aes(x = minutes, y = `F0`, group = Conditions, color = Conditions)) + 
  annotate("rect", xmin = 13, xmax = 43, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(group = Conditions), method = "gam", se = TRUE, formula = y ~ s(x), span = 0.4) +
  geom_ribbon(aes(ymin = ..ymin.., ymax = ..ymax.., group = Conditions, fill = Conditions),
              stat = "smooth", method = "gam", se = TRUE, formula = y ~ s(x), alpha = 0.2, color = NA) + 
  scale_color_manual(values = c("#324851", "#E7B75F")) + 
  scale_fill_manual(values = c("#324851", "#E7B75F")) +  
  labs(x = "", y = "Surface microalgal biomass (F0 a.u.)", 
       title = "", 
       color = "Conditions", fill = "Conditions") + 
  theme_test(base_family = "Arial")+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),  
        axis.text.x = element_text(size = 18, color = "black"),                
        axis.text.y = element_text(size = 18, color = "black"),
        axis.ticks = element_line(color = "black")) +
  ylim(y_limits)

# Plot HP
plot_hp <- ggplot(File_HP_data, aes(x = minutes, y = `F0`, group = Conditions, color = Conditions)) + 
  annotate("rect", xmin = 25, xmax = 32.3, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(group = Conditions), method = "gam", se = TRUE, formula = y ~ s(x), span = 0.4) + 
  geom_ribbon(aes(ymin = ..ymin.., ymax = ..ymax.., group = Conditions, fill = Conditions),
              stat = "smooth", method = "gam", se = TRUE, formula = y ~ s(x), alpha = 0.2, color = NA) + 
  scale_color_manual(values = c("#324851", "#5C88C4")) + 
  scale_fill_manual(values = c("#324851", "#5C88C4")) +  
  labs(x = "Time (Minutes)", y = "", 
       title = "", 
       color = "Treatments", fill = "Treatments") + 
  theme_test(base_family = "Arial")+
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),  
        axis.text.x = element_text(size = 18, color = "black"),                
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black")) +
  ylim(y_limits)

# Plot PL
plot_pl <- ggplot(File_PL_data, aes(x = minutes, y = `F0`, group = Conditions, color = Conditions)) + 
  annotate("rect", xmin = 20, xmax = 31, ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.3) +
  geom_point(alpha = 0.5) +  
  stat_smooth(aes(group = Conditions), method = "gam", se = TRUE, formula = y ~ s(x), span = 0.4) + 
  geom_ribbon(aes(ymin = ..ymin.., ymax = ..ymax.., group = Conditions, fill = Conditions),
              stat = "smooth", method = "gam", se = TRUE, formula = y ~ s(x), alpha = 0.2, color = NA) + 
  scale_color_manual(values = c("#324851", "#937DC2")) + 
  scale_fill_manual(values = c("#324851", "#937DC2")) +  
  labs(x = "", y = "", 
       title = "", 
       color = "Conditions", fill = "Conditions") + 
  theme_test(base_family = "Arial")+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),  
        axis.text.x = element_text(size = 18, color = "black"),                
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black")) +
  ylim(y_limits)

# Final plot
(plot_hl | plot_hp | plot_pl) 

# FIGURE 3 HEATMAP ####

df_heatmap <- df_Pparameters %>%select(Time, Name, `α-slope`, `rETRm obs`, "Fv/Fm", NPQ1000, YNPQ1000, `YNOm obs`)

# Mobility subset
df_heatmap_wsb <- subset(df_heatmap, grepl("^WSB_", Name))
df_heatmap_sb <- subset(df_heatmap, grepl("^SB_", Name))

# Differences After-Before
df_heatmap_sb_diff <- df_heatmap_sb %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `α-slope` = `α-slope`[Time == "After"] - `α-slope`[Time == "Before"],
    `Fv/Fm` = `Fv/Fm`[Time == "After"] - `Fv/Fm`[Time == "Before"],
    `rETRm obs` = `rETRm obs`[Time == "After"] - `rETRm obs`[Time == "Before"],
    `NPQ 1000` = NPQ1000[Time == "After"] - NPQ1000[Time == "Before"],
    `Y(NPQ) 1000` = YNPQ1000[Time == "After"] - YNPQ1000[Time == "Before"],
    `Y(NO)m obs` = `YNOm obs`[Time == "After"] - `YNOm obs`[Time == "Before"])

df_heatmap_wsb_diff <- df_heatmap_wsb %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `α-slope` = `α-slope`[Time == "After"] - `α-slope`[Time == "Before"],
    `Fv/Fm` = `Fv/Fm`[Time == "After"] - `Fv/Fm`[Time == "Before"],
    `rETRm obs` = `rETRm obs`[Time == "After"] - `rETRm obs`[Time == "Before"],
    `NPQ 1000` = NPQ1000[Time == "After"] - NPQ1000[Time == "Before"],
    `Y(NPQ) 1000` = YNPQ1000[Time == "After"] - YNPQ1000[Time == "Before"],
    `Y(NO)m obs` = `YNOm obs`[Time == "After"] - `YNOm obs`[Time == "Before"])

# SB comparaison groups by day
SB_DA_1_3 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_DA_1", "SB_DA_2", "SB_DA_3")) # J1 experiment
SB_HL_1_3 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_HL_1", "SB_HL_2", "SB_HL_3")) # J1 experiment
SB_DA_7_9 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_DA_7", "SB_DA_8", "SB_DA_9")) # J3 experiment
SB_HP_1_3 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_HP_1", "SB_HP_2", "SB_HP_3")) # J3 experiment
SB_DA_4_6 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_DA_4", "SB_DA_5", "SB_DA_6")) # J2 experiment
SB_PL_1_3 <- df_heatmap_sb_diff %>% filter(Name %in% c("SB_PL_1", "SB_PL_2", "SB_PL_3")) # J2 experiment

# WSB comparaison groups by day
WSB_DA_1_3 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_DA_1", "WSB_DA_2", "WSB_DA_3"))
WSB_HL_1_3 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_HL_1", "WSB_HL_2", "WSB_HL_3"))
WSB_DA_7_9 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_DA_7", "WSB_DA_8", "WSB_DA_9"))
WSB_HP_1_3 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_HP_1", "WSB_HP_2", "WSB_HP_3"))
WSB_DA_4_6 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_DA_4", "WSB_DA_5", "WSB_DA_6"))
WSB_PL_1_3 <- df_heatmap_wsb_diff %>% filter(Name %in% c("WSB_PL_1", "WSB_PL_2", "WSB_PL_3"))

# P-values loop function
perform_t_tests_pv_SB <- function(group1, group2) { # For SB group
  p_values <- c()
  for (param in colnames(df_heatmap_sb_diff)[2:7]) {
    t_test <- t.test(group1[[param]], group2[[param]])
    p_values <- c(p_values, t_test$p.value)
  }
  return(p_values)
}

perform_t_tests_pv_WSB <- function(group1, group2) { # For WSB group
  p_values <- c()
  for (param in colnames(df_heatmap_wsb_diff)[2:7]) {
    t_test <- t.test(group1[[param]], group2[[param]])
    p_values <- c(p_values, t_test$p.value)
  }
  return(p_values)
}

# P-values calculations stress vs DA
p_values_HL_SB <- perform_t_tests_pv_SB(SB_DA_1_3, SB_HL_1_3)
p_values_HP_SB <- perform_t_tests_pv_SB(SB_DA_7_9, SB_HP_1_3)
p_values_PL_SB <- perform_t_tests_pv_SB(SB_DA_4_6, SB_PL_1_3)

p_values_HL_WSB <- perform_t_tests_pv_WSB(WSB_DA_1_3, WSB_HL_1_3)
p_values_HP_WSB <- perform_t_tests_pv_WSB(WSB_DA_7_9, WSB_HP_1_3)
p_values_PL_WSB <- perform_t_tests_pv_WSB(WSB_DA_4_6, WSB_PL_1_3)

# P-values df
df_heatmap_diff_pv <- data.frame(
  "HL" = p_values_HL_SB,
  "HP" = p_values_HP_SB,
  "PL" = p_values_PL_SB,
  "HL." = p_values_HL_WSB,
  "HP." = p_values_HP_WSB,
  "PL." = p_values_PL_WSB)

# Add photosynthetic parameters names metadata
df_heatmap_diff_pv <- t(df_heatmap_diff_pv)
colnames(df_heatmap_diff_pv) <- colnames(df_heatmap_sb_diff)[2:7]

# Df formatting
df_heatmap_diff_pv_long <- as.data.frame(as.table(as.matrix(df_heatmap_diff_pv)))

colnames(df_heatmap_diff_pv_long) <- c("Condition", "Variable", "Value")

# Add p-values categories
df_heatmap_diff_pv_long$Categorie <- cut(df_heatmap_diff_pv_long$Value,breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                         labels = c("< 0.001",
                                                    "0.001 < > 0.01",
                                                    "0.01 < > 0.05",
                                                    "> 0.05"), right = FALSE)

# % differences After-Before
df_heatmap_sb_diff_perc <- df_heatmap_sb %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `α-slope` = (((`α-slope`[Time == "After"]*100) / `α-slope`[Time == "Before"])-100),
    `Fv/Fm` = (((`Fv/Fm`[Time == "After"]*100) / `Fv/Fm`[Time == "Before"])-100),
    `rETRm obs` = (((`rETRm obs`[Time == "After"]*100) / `rETRm obs`[Time == "Before"])-100),
    `NPQ 1000` = (((NPQ1000[Time == "After"]*100) / NPQ1000[Time == "Before"])-100),
    `Y(NPQ) 1000` = (((YNPQ1000[Time == "After"]*100) / YNPQ1000[Time == "Before"])-100),
    `Y(NO)m obs` = (((`YNOm obs`[Time == "After"]*100) / `YNOm obs`[Time == "Before"])-100))

df_heatmap_wsb_diff_perc <- df_heatmap_wsb %>%
  group_by(Name) %>%
  filter(n() == 2) %>%
  summarise(
    `α-slope` = (((`α-slope`[Time == "After"]*100) / `α-slope`[Time == "Before"])-100),
    `Fv/Fm` = (((`Fv/Fm`[Time == "After"]*100) / `Fv/Fm`[Time == "Before"])-100),
    `rETRm obs` = (((`rETRm obs`[Time == "After"]*100) / `rETRm obs`[Time == "Before"])-100),
    `NPQ 1000` = (((NPQ1000[Time == "After"]*100) / NPQ1000[Time == "Before"])-100),
    `Y(NPQ) 1000` = (((YNPQ1000[Time == "After"]*100) / YNPQ1000[Time == "Before"])-100),
    `Y(NO)m obs` = (((`YNOm obs`[Time == "After"]*100) / `YNOm obs`[Time == "Before"])-100))

# Mean % loop function
perform_mean_diff_sb <- function(group1, group2) {
  mean_diff <- c()
  for (param in colnames(df_heatmap_sb_diff_perc)[2:7]) {
    diff <- mean(group2[[param]]) - mean(group1[[param]])
    mean_diff <- c(mean_diff, diff)
  }
  return(mean_diff)
}

perform_mean_diff_wsb <- function(group1, group2) {
  mean_diff <- c()
  for (param in colnames(df_heatmap_wsb_diff_perc)[2:7]) {
    diff <- mean(group2[[param]]) - mean(group1[[param]])
    mean_diff <- c(mean_diff, diff)
  }
  return(mean_diff)
}

# Groups
SB_DA_1_3_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_DA_1", "SB_DA_2", "SB_DA_3")) # J1 experiment
SB_HL_1_3_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_HL_1", "SB_HL_2", "SB_HL_3")) # J1 experiment
SB_DA_7_9_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_DA_7", "SB_DA_8", "SB_DA_9")) # J3 experiment
SB_HP_1_3_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_HP_1", "SB_HP_2", "SB_HP_3")) # J3 experiment
SB_DA_4_6_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_DA_4", "SB_DA_5", "SB_DA_6")) # J2 experiment
SB_PL_1_3_perc <- df_heatmap_sb_diff_perc %>% filter(Name %in% c("SB_PL_1", "SB_PL_2", "SB_PL_3")) # J2 experiment

WSB_DA_1_3_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_DA_1", "WSB_DA_2", "WSB_DA_3"))
WSB_HL_1_3_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_HL_1", "WSB_HL_2", "WSB_HL_3"))
WSB_DA_7_9_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_DA_7", "WSB_DA_8", "WSB_DA_9"))
WSB_HP_1_3_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_HP_1", "WSB_HP_2", "WSB_HP_3"))
WSB_DA_4_6_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_DA_4", "WSB_DA_5", "WSB_DA_6"))
WSB_PL_1_3_perc <- df_heatmap_wsb_diff_perc %>% filter(Name %in% c("WSB_PL_1", "WSB_PL_2", "WSB_PL_3"))

# Means calculations stress vs DA
mean_diff_HL_SB_perc <- perform_mean_diff_sb(SB_DA_1_3_perc, SB_HL_1_3_perc)
mean_diff_HP_SB_perc <- perform_mean_diff_sb(SB_DA_7_9_perc, SB_HP_1_3_perc)
mean_diff_PL_SB_perc <- perform_mean_diff_sb(SB_DA_4_6_perc, SB_PL_1_3_perc)

mean_diff_HL_WSB_perc <- perform_mean_diff_wsb(WSB_DA_1_3_perc, WSB_HL_1_3_perc)
mean_diff_HP_WSB_perc <- perform_mean_diff_wsb(WSB_DA_7_9_perc, WSB_HP_1_3_perc)
mean_diff_PL_WSB_perc <- perform_mean_diff_wsb(WSB_DA_4_6_perc, WSB_PL_1_3_perc)

# Means df
df_heatmap_mean_diff_perc <- data.frame(
  "HL" = mean_diff_HL_SB_perc,
  "HP" = mean_diff_HP_SB_perc,
  "PL" = mean_diff_PL_SB_perc,
  "HL." = mean_diff_HL_WSB_perc,
  "HP." = mean_diff_HP_WSB_perc,
  "PL." = mean_diff_PL_WSB_perc)

# Add photosynthetic parameters names metadata
df_heatmap_mean_diff_perc <- t(df_heatmap_mean_diff_perc)
colnames(df_heatmap_mean_diff_perc) <- colnames(df_heatmap_sb_diff)[2:7]

# Df formatting
df_heatmap_mean_diff_perc_long <- as.data.frame(as.table(as.matrix(df_heatmap_mean_diff_perc)))

colnames(df_heatmap_mean_diff_perc_long) <- c("Condition", "Variable", "Perc")

# Add arrows pictograms and colors
df_heatmap_mean_diff_perc_long$Arrow <- ifelse(df_heatmap_mean_diff_perc_long$Perc > 0, "↑", "↓")
df_heatmap_mean_diff_perc_long$Arrow_Color <- ifelse(df_heatmap_mean_diff_perc_long$Arrow == "↓", "firebrick", "black")

# Add 'Perc' % value column
df_heatmap_diff_pv_long <- cbind(df_heatmap_diff_pv_long, df_heatmap_mean_diff_perc_long[c("Arrow", "Arrow_Color", "Perc")])

# Converting 'Perc' (Absolute values)
df_heatmap_diff_pv_long$Perc <- abs(df_heatmap_diff_pv_long$Perc)

# Plot graph
ggplot(df_heatmap_diff_pv_long, aes(x = Condition, y = factor(Variable, levels = rev(unique(Variable))), fill = Categorie)) +
  geom_tile(color = "black", size = 0.7) + 
  scale_fill_manual(name = "P-values Welch Two Sample t-test", 
                    values = c("< 0.001" = "#216BA3",
                               "0.001 < > 0.01" = "#6096ba",
                               "0.01 < > 0.05" = "#a3cef1",
                               "> 0.05" = "white"))+
  labs(title = "", x = "Treatments", y = "") +
  geom_text(data = subset(df_heatmap_diff_pv_long, !Categorie %in% c("> 0.05")), 
            aes(label = Arrow, color = Arrow_Color), size = 12,  hjust = 0.5, vjust = 0.4) +
  scale_color_identity()+
  theme_test(base_family = "Arial")+
  theme(axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")+
  annotate("segment", 
           x = 0.51, xend = 3.48,  
           y = 6.8, yend = 6.8,  
           colour = "black", 
           size = 0.7) +
  annotate("segment", 
           x = 3.51, xend = 6.5,  
           y = 6.8, yend = 6.8,  
           colour = "black", 
           size = 0.7) +
  annotate("text", 
           x = 2, y = 7.3, 
           label = "Sediment biofilm", 
           size = 10, 
           colour = "black")+
  annotate("text", 
           x = 5, y = 7.3, 
           label = "Without sediment biofilm", 
           size = 10, 
           colour = "black")+
  scale_y_discrete(expand = expansion(mult = c(0.2, 0.37)))

# FIGURE 4 PPARAMETERS WSB ####

# Retrieve parameters
Pparameters <- merge(
  df_Pparameters,
  rETR[, c("No", names(rETR)[2:9])],
  by = "No",
  all.x = TRUE)

df_long_Pp <- Pparameters %>%
  mutate(RowID = row_number()) %>%
  pivot_longer(cols = c(23:50, 70:76),
               names_to = "Pparameter",
               values_to = "Value") %>%
  select(2:3, Pparameter, Value)

# Df formatting
df_long_Pp <- df_long_Pp %>%
  mutate(Irradiance = case_when(
    grepl("_l1", Pparameter) ~ "10",
    grepl("_l2", Pparameter) ~ "20",
    grepl("_l3", Pparameter) ~ "50",
    grepl("_l4", Pparameter) ~ "100",
    grepl("_l5", Pparameter) ~ "300",
    grepl("_l6", Pparameter) ~ "500",
    grepl("_l7", Pparameter) ~ "1000"))

df_long_Pp <- df_long_Pp %>%
  mutate(Pparameter = case_when(
    grepl("^qy", Pparameter)   ~ "Fv/Fm",
    grepl("^NPQ", Pparameter)  ~ "NPQ",
    grepl("^YNO", Pparameter)  ~ "YNO",
    grepl("^YNPQ", Pparameter) ~ "YNPQ",
    grepl("^rETR", Pparameter) ~ "rETR"))

df_long_Pp <- df_long_Pp %>% filter(str_starts(Name, "WSB"))

df_long_Pp <- df_long_Pp %>% mutate(Conditions = str_extract(Name, "(?<=_)[^_]+(?=_)"))

df_long_Pp <- df_long_Pp %>% mutate(Conditions = paste(Conditions, Time, sep = " "))

# Set Irradiance levels
df_long_Pp <- df_long_Pp %>%mutate(Irradiance = factor(Irradiance, levels = c(10, 20, 50, 100, 300, 500, 1000), ordered = TRUE))

# Set colors
Pp_colors <- c(
  "DA Before" = "#CBC9C9",
  "DA After" = "#040404",
  
  "HL Before" = "#EADD95",
  "HL After"  = "#E7B75F",
  
  "HP Before" = "#94B5FD",
  "HP After"  = "#5C88C4",
  
  "PL After"  = "#937DC2",
  "PL Before" = "#DAB7FE")

# Plot Fv/Fm
df_fvfm <- df_long_Pp %>%
  filter(Pparameter == "Fv/Fm") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

aPp <- ggplot(df_fvfm, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "",
       y = "PSII quantum yield (unitless)",
       color = "Treatments") +
  theme_test(base_family = "Arial") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),  
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks = element_line(color = "black"))

# Plot rETR
df_retr <- df_long_Pp %>%
  filter(Pparameter == "rETR") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

bPp <- ggplot(df_retr, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "",
       y = "rETR (arbitrary units)",
       color = "Treatments") +
  theme_test(base_family = "Arial") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),  
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks = element_line(color = "black"))

# Plot NPQ
df_npq <- df_long_Pp %>%
  filter(Pparameter == "NPQ") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

cPp <- ggplot(df_npq, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "Irradiance (µmol photons m-2 s-1)",
       y = "NPQ (unitless)",
       color = "Treatments") +
  theme_test(base_family = "Arial") +
  theme(axis.title.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),  
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks = element_line(color = "black")) +
  theme(legend.position = "bottom")


# Plot YNPQ
df_ynpq <- df_long_Pp %>%
  filter(Pparameter == "YNPQ") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

dPp <- ggplot(df_ynpq, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "Irradiance (µmol photons m-2 s-1)",
       y = "Y(NPQ) (unitless)",
       color = "Treatments") +
  theme_test(base_family = "Arial") +
  theme(axis.title.x = element_text(size = 18, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),  
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks = element_line(color = "black")) +
  theme(legend.position = "bottom")

# Plot YNO
df_yno <- df_long_Pp %>%
  filter(Pparameter == "YNO") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

ePp <- ggplot(df_yno, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "",
       y = "Y(NO) (unitless)",
       color = "Treatments") +
  theme_test(base_family = "Arial") +
  theme( legend.position = "none",
         axis.title.x = element_text(size = 18, color = "black"), 
         axis.title.y = element_text(size = 18, color = "black"),  
         axis.text.x = element_text(size = 14, color = "black"),
         axis.text.y = element_text(size = 14, color = "black"),
         axis.ticks = element_line(color = "black"))

(aPp | bPp | dPp | ePp) +
  plot_annotation(
    title = "",
    theme = theme(plot.title = element_text(size = 18)))

# FIGURE 5 BARPLOT RATIOS ####
gdf_pigment_ratio <- gdf_pigment

# Ratio calculations
gdf_pigment_ratio <- gdf_pigment_ratio %>%
  mutate(rChla_deriv = (`Chla-derivatives` / (`Chla-derivatives` + Chla + Pheopigments)) * 100)

gdf_pigment_ratio <- gdf_pigment_ratio %>%
  mutate(rPheo = (Pheopigments / (`Chla-derivatives` + Chla + Pheopigments)) * 100)

gdf_pigment_ratio <- gdf_pigment_ratio %>%
  mutate(DES = (Dt / (Dd + Dt)) * 100)


# Df barplot
barplot_ratios <- as.data.frame(gdf_pigment_ratio[, -c(5:34)]) 
barplot_ratios$Conditions <- factor(barplot_ratios$Conditions, levels = c("PL", "HP", "HL", "DA"))

# Standard deviation and mean calculations
barplot_ratios_rChla_deriv <- barplot_ratios %>%
  group_by(Conditions, Groups) %>%
  summarise(
    mean = mean(rChla_deriv),
    sd = sd(rChla_deriv),
    se = sd(rChla_deriv) / sqrt(n()))
barplot_ratios_rChla_deriv$Conditions <- factor(barplot_ratios_rChla_deriv$Conditions, levels = c("PL", "HP", "HL", "DA"))

barplot_ratios_rPheo <- barplot_ratios %>%
  group_by(Conditions, Groups) %>%
  summarise(
    mean = mean(rPheo),
    sd = sd(rPheo),
    se = sd(rPheo) / sqrt(n()))
barplot_ratios_rPheo$Conditions <- factor(barplot_ratios_rPheo$Conditions, levels = c("PL", "HP", "HL", "DA"))

barplot_ratios_DES <- barplot_ratios %>%
  group_by(Conditions, Groups) %>%
  summarise(
    mean = mean(DES),
    sd = sd(DES),
    se = sd(DES) / sqrt(n()))
barplot_ratios_DES$Conditions <- factor(barplot_ratios_DES$Conditions, levels = c("PL", "HP", "HL", "DA"))

# Reverse values
barplot_ratios <- barplot_ratios %>%
  mutate(across(5:7, ~ ifelse(Groups == "SB", -.x, .x)))

barplot_ratios_rChla_deriv <- barplot_ratios_rChla_deriv %>%
  mutate(mean_plot = ifelse(Groups == "SB", -mean, mean),
         se_plot = se)

barplot_ratios_rPheo <- barplot_ratios_rPheo %>%
  mutate(mean_plot = ifelse(Groups == "SB", -mean, mean),
         se_plot = se)

barplot_ratios_DES <- barplot_ratios_DES %>%
  mutate(mean_plot = ifelse(Groups == "SB", -mean, mean),
         se_plot = se)

# Colors
barplot_ratios_colors <- c("DA" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2")

# DES plot
p1 <- ggplot() +
  geom_bar(data = barplot_ratios_DES, 
           aes(x = Conditions, y = mean_plot, fill = Conditions),
           stat = "identity", 
           width = 0.75,
           color = "black",
           linewidth = 0.75) +
  geom_jitter(data = barplot_ratios, 
              aes(x = Conditions, y = DES, fill = Conditions),
              position = position_jitter(width = 0.15, height = 0),
              size = 3, shape = 21, fill = "white", color = "black") +
  
  geom_errorbar(data = barplot_ratios_DES,
                aes(x = Conditions, ymin = mean_plot - se_plot, ymax = mean_plot + se_plot),
                width = 0.2, 
                position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(values = barplot_ratios_colors, breaks = c("DA", "HL", "HP", "PL")) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  
  coord_flip() +
  
  theme_linedraw(base_family = "Arial") +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none")+
  
  scale_y_continuous(breaks = seq(-50, 50, 10),
                     labels = abs(seq(-50, 50, 10))) +
  
  labs(title = "")+
  ylab("Mean DES (%)")+
  

  annotate("text", 
           x = 5.5, y = -50,
           label = "", 
           size = 5)+
labs(title = "A") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0))


# Chla_deriv plot
p2 <- ggplot() +
  geom_bar(data = barplot_ratios_rChla_deriv, 
           aes(x = Conditions, y = mean_plot, fill = Conditions),
           stat = "identity", 
           width = 0.75,
           color = "black",
           linewidth = 0.75) +
  geom_jitter(data = barplot_ratios, 
              aes(x = Conditions, y = rChla_deriv, fill = Conditions),
              position = position_jitter(width = 0.15, height = 0),
              size = 3, shape = 21, fill = "white", color = "black") +
  
  geom_errorbar(data = barplot_ratios_rChla_deriv,
                aes(x = Conditions, ymin = mean_plot - se_plot, ymax = mean_plot + se_plot),
                width = 0.2, 
                position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(values = barplot_ratios_colors, breaks = c("DA", "HL", "HP", "PL")) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  
  coord_flip() +
  
  theme_linedraw(base_family = "Arial") +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none")+
  
  scale_y_continuous(breaks = seq(-30, 30, 10),
                     labels = abs(seq(-30, 30, 10))) +
  
  labs(title = "")+
  ylab("Mean Chlorophyll a-derivatives normalized to Chlorophyll a (%)")+
  

  annotate("text", 
           x =5.5, y = -26,
           label = "", 
           size = 5)+
  annotate("text", 
           x = 4, y = 25.2,
           label = "", 
           size = 5)+
labs(title = "B") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0))

# Pheo plot
p3 <- ggplot() +
  geom_bar(data = barplot_ratios_rPheo, 
           aes(x = Conditions, y = mean_plot, fill = Conditions),
           stat = "identity", 
           width = 0.75,
           color = "black",
           linewidth = 0.75) +
  geom_jitter(data = barplot_ratios, 
              aes(x = Conditions, y = rPheo, fill = Conditions),
              position = position_jitter(width = 0.15, height = 0),
              size = 3, shape = 21, fill = "white", color = "black") +
  
  geom_errorbar(data = barplot_ratios_rPheo,
                aes(x = Conditions, ymin = mean_plot - se_plot, ymax = mean_plot + se_plot),
                width = 0.2, 
                position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(name = "Treatments", values = barplot_ratios_colors, breaks = c("DA", "HL", "HP", "PL")) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  
  coord_flip() +
  
  theme_linedraw(base_family = "Arial") +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom")+
  
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     labels = abs(seq(-10, 10, 5))) +
  
  labs(title = "")+
  ylab("Mean Pheopigments normalized to Chlorophyll a (%)")+
  theme(plot.title = element_text(face = "bold"))+
  
  annotate("text", 
           x =5.5, y = -9,
           label = "", 
           size = 5)+
  annotate("text", 
           x = 4, y = 8.7,
           label = "", 
           size = 5)+
labs(title = "C") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0))
  
(p1) / (p2 | p3)
  

# P-values calculations
wilcox.test(DES ~ Conditions, data = barplot_ratios, 
            subset = Groups == "SB" & Conditions %in% c("DA", "PL"),
            exact = TRUE)

# FIGURE S1 QPHAR ####

# Panel A
ggplot(spectre, aes(x = Wavelength, y = Q_umol)) +
  geom_line(size = 1) +
  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "Irradiance (µmol photons m-2 s-1)",
    title = "Incident HL irradiance spectrum"
  ) +
  annotate(
    "text", 
    x = min(spectre$Wavelength),  
    y = max(spectre$Q_umol), 
    label = paste0("Total PPFD = ", round(total_PPFD, 1), " µmol photons m-2 s-1"),
    hjust = 0, vjust = 1,  
    size = 5)

# Panel B
ggplot(a_rec_un, aes(x = lambda)) +
  geom_line(aes(y = abs_rec_un, color = "Original")) +
  geom_line(aes(y = abs_rec_corrected_pack, color = "Package-corrected")) +  
  labs(x = "Wavelength (nm)", y = "Reconstitued spectrum (m-1)", color = "Reconstituted spectrum", title = "") +
  theme_test()+
  theme(legend.position = "bottom")

# Panel C
ggplot(a_rec_un, aes(x = lambda)) +
  geom_line(aes(y = decrease_package_effect)) +
  labs(x = "Wavelength (nm)", y = "Decrease pEffect (%)", title = "DA WSB") +
  theme_test()

# Panel D
ggplot(a_mpb_interp, aes(x = lambda, y = AbsorptionTotal, color = Condition, fill = Condition)) +
  
  stat_summary(
    fun = mean,
    geom = "line",
    size = 1) +
  
  stat_summary(
    fun.data = mean_sdl,       
    fun.args = list(mult = 1),  
    geom = "ribbon",
    alpha = 0.4,
    color = NA) +
  
  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "Surfacic absorption spectrum rec,un (dimensionless)",
    title = "Reconstructed spectrums without cell packing correction")

# Panel E
ggplot(a_mpb_interp, aes(x = lambda, y = AbsorptionTotal_corrected, color = Condition, fill = Condition)) +
  
  stat_summary(
    fun = mean,
    geom = "line",
    size = 1) +
  
  stat_summary(
    fun.data = mean_sdl,       
    fun.args = list(mult = 1),  
    geom = "ribbon",
    alpha = 0.4,
    color = NA) +
  
  theme_test() +
  labs(
    x = "Wavelength (nm)",
    y = "Surfacic absorption spectrum rec,pac (dimensionless)",
    title = "Reconstructed spectrums with cell packing correction")


# FIGURE S2 EUK DIV ####

# EUKARYOTES

richesse_spe18 <- estimate_richness(biofilm18S)

# Pielou calculation
richesse_spe18$Pielou <- richesse_spe18$Shannon / log(richesse_spe18$Observed) 

# Metadata
metadata_div<- df_community[, 1:3]  
rownames(metadata_div) <- metadata_div$File
richesse_spe18 <- merge(richesse_spe18, metadata_div, by = "row.names", all.x = TRUE)
richesse_spe18 <- richesse_spe18 %>%mutate(Treatment = paste(Conditions, Groups, sep = "_"))

# Df format
richesse_spe18_long <- richesse_spe18 %>%
  select(Observed, Shannon, InvSimpson, Pielou, Groups, Treatment) %>%
  gather(key = "Indice", value = "Valeur", Observed, Shannon, InvSimpson, Pielou)

# Order
richesse_spe18_long$Indice <- factor(richesse_spe18_long$Indice, levels = c("Observed", "Shannon", "InvSimpson", "Pielou"))

# Colors settings
colors_groups <- c("SB" = "#CC6677", "WSB" = "#81B29A")

# Stat comparisons 
my_comparisons <- list(c("WSB", "SB"))

# Plot
ggplot(richesse_spe18_long, aes(x = Groups, y = Valeur, fill = Groups, color = Groups)) +
  geom_boxplot(outlier.shape = NA, size = 0.5, alpha = 1) + 
  geom_jitter(width = 0.2, shape = 16, size = 2, alpha = 1) + 
  scale_fill_manual(values = colors_groups) + 
  scale_color_manual(values = c("SB" = "black", "WSB" = "black")) + 
  facet_wrap(~ Indice, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.10))) +
  labs(x = "", y = "Index values", fill = "Treatments", title = "Eukaryotic community") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif")+
  theme(axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(size = 18, hjust = 0.5))



# FIGURE S3 PPARAMETERS SB ####

# Retrieve parameters
Pparameters <- merge(
  df_Pparameters,
  rETR[, c("No", names(rETR)[2:9])],
  by = "No",
  all.x = TRUE)

df_long_Pp <- Pparameters %>%
  mutate(RowID = row_number()) %>%
  pivot_longer(cols = c(23:50, 70:76),
               names_to = "Pparameter",
               values_to = "Value") %>%
  select(2:3, Pparameter, Value)

# Df formatting
df_long_Pp <- df_long_Pp %>%
  mutate(Irradiance = case_when(
    grepl("_l1", Pparameter) ~ "10",
    grepl("_l2", Pparameter) ~ "20",
    grepl("_l3", Pparameter) ~ "50",
    grepl("_l4", Pparameter) ~ "100",
    grepl("_l5", Pparameter) ~ "300",
    grepl("_l6", Pparameter) ~ "500",
    grepl("_l7", Pparameter) ~ "1000"))

df_long_Pp <- df_long_Pp %>%
  mutate(Pparameter = case_when(
    grepl("^qy", Pparameter)   ~ "Fv/Fm",
    grepl("^NPQ", Pparameter)  ~ "NPQ",
    grepl("^YNO", Pparameter)  ~ "YNO",
    grepl("^YNPQ", Pparameter) ~ "YNPQ",
    grepl("^rETR", Pparameter) ~ "rETR"))

df_long_Pp <- df_long_Pp %>% filter(str_starts(Name, "SB"))

df_long_Pp <- df_long_Pp %>% mutate(Conditions = str_extract(Name, "(?<=_)[^_]+(?=_)"))

df_long_Pp <- df_long_Pp %>% mutate(Conditions = paste(Conditions, Time, sep = " "))

# Set Irradiance levels
df_long_Pp <- df_long_Pp %>%mutate(Irradiance = factor(Irradiance, levels = c(10, 20, 50, 100, 300, 500, 1000), ordered = TRUE))

# Set colors
Pp_colors <- c(
  "DA Before" = "#CBC9C9",
  "DA After" = "#040404",
  
  "HL Before" = "#EADD95",
  "HL After"  = "#E7B75F",
  
  "HP Before" = "#94B5FD",
  "HP After"  = "#5C88C4",
  
  "PL After"  = "#937DC2",
  "PL Before" = "#DAB7FE")

# Plot Fv/Fm
df_fvfm <- df_long_Pp %>%
  filter(Pparameter == "Fv/Fm") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

aPp <- ggplot(df_fvfm, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "",
       y = "PSII quantum yield (unitless)",
       color = "Treatments") +
  theme_test() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# Plot rETR
df_retr <- df_long_Pp %>%
  filter(Pparameter == "rETR") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

bPp <- ggplot(df_retr, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "",
       y = "rETR (arbitrary units)",
       color = "Treatments") +
  theme_test() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# Plot NPQ
df_npq <- df_long_Pp %>%
  filter(Pparameter == "NPQ") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

cPp <- ggplot(df_npq, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "Irradiance (µmol photons.m-2.s-1)",
       y = "NPQ (unitless)",
       color = "Treatments") +
  theme_test() +
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "bottom")


# Plot YNPQ
df_ynpq <- df_long_Pp %>%
  filter(Pparameter == "YNPQ") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

dPp <- ggplot(df_ynpq, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "Irradiance (µmol photons m-2 s-1)",
       y = "Y(NPQ) (unitless)",
       color = "Treatments") +
  theme_test() +
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "bottom")

# Plot YNO
df_yno <- df_long_Pp %>%
  filter(Pparameter == "YNO") %>%
  mutate(Irradiance = as.numeric(trimws(Irradiance)))

ePp <- ggplot(df_yno, aes(x = Irradiance, y = Value, color = Conditions, group = Conditions)) +
  geom_point(alpha = 0.5)+
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) + 
  scale_color_manual(values = Pp_colors) +
  labs(title = "",
       x = "",
       y = "Y(NO) (unitless)",
       color = "Treatments") +
  theme_test() +
  theme( legend.position = "none",
         axis.title.x = element_text(size = 18), 
         axis.title.y = element_text(size = 18),  
         axis.text.x = element_text(size = 14),
         axis.text.y = element_text(size = 14))

(aPp | bPp | dPp | ePp) +
  plot_annotation(
    title = "",
    theme = theme(plot.title = element_text(size = 18)))

# FIGURE S4 PIG DIV ####

df_barplot_pigment <- gdf_pigment_perc

#Carotenoid-like group
df_barplot_pigment <- df_barplot_pigment %>%
  mutate(`Unknown-carotenoids` = rowSums(select(., starts_with("Uc")))) %>%
  select(-starts_with("Uc"))

#fucoxanthin-like group
df_barplot_pigment<-df_barplot_pigment %>%
  mutate(`Fucoxanthin-likes` = rowSums(select(., starts_with("F-l")))) %>%
  select(-starts_with("F-l"))

colors_barplot_pigments <- c(
  "V" = "#774993", 
  "A" = "#1f77b4", 
  "Z" = "#93375A", 
  "Dt" = "#0A0701", 
  "Dd" = "#c5b0d5",
  "Chla" = "#4C8543", 
  "Chla-derivatives" = "#B3EE98", 
  "Pheopigments" = "darkgreen",
  "F" = "#81724F", 
  "Fucoxanthin-likes" = "#c49c94", 
  "Unknown-carotenoids" = "#ff9896", 
  "Chlc2" = "#ffbb78", 
  "ββ" = "#ff7f0e", 
  "βε" = "#aec7e8", 
  "L" = "#dbdb8d", 
  "N" = "#382469")

# Subset by groups
df_barplot_pigment_sb <- subset(df_barplot_pigment, Groups == "SB")
df_barplot_pigment_wsb <- subset(df_barplot_pigment, Groups == "WSB")

#SB
data_perc_sb <- df_barplot_pigment_sb %>%
  select(Conditions, 5:20) %>%
  pivot_longer(cols = -Conditions, names_to = "Pigments", values_to = "Percentage")

# Plot
data_perc_sb$Conditions <- factor(data_perc_sb$Conditions, levels = c("DA", "HL", "HP", "PL"))

barplot_pigment1 <- ggplot(data_perc_sb, aes(x = Conditions, y = Percentage, fill = Pigments)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) + 
  theme_test() +
  labs(title = "Sediment biofilm", x = "", y = "Relative abundance (%)") +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = colors_barplot_pigments)+
  theme(axis.title.y = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20),
        legend.position = "bottom")+
  annotate("text", 
           x = 1, y = -0.05, 
           label = "DA", 
           size = 6, 
           colour = "black")+
  annotate("text", 
           x = 2, y = -0.05, 
           label = "HL", 
           size = 6, 
           colour = "black")+
  annotate("text", 
           x = 3, y = -0.05, 
           label = "HP", 
           size = 6, 
           colour = "black")+
  annotate("text", 
           x = 4, y = -0.05, 
           label = "PL", 
           size = 6, 
           colour = "black")

#WSB
data_perc_wsb <- df_barplot_pigment_wsb %>%
  select(Conditions, 5:20) %>%
  pivot_longer(cols = -Conditions, names_to = "Pigments", values_to = "Percentage")

# Plot
data_perc_wsb$Conditions <- factor(data_perc_wsb$Conditions, levels = c("DA", "HL", "HP", "PL"))

barplot_pigment2 <- ggplot(data_perc_wsb, aes(x = Conditions, y = Percentage, fill = Pigments)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) + 
  theme_test() +
  labs(title = "Without sediment biofilm", x = "", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = colors_barplot_pigments)+
  theme(axis.title.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20),
        legend.position = "none")+
  annotate("text", 
           x = 1, y = -0.05, 
           label = "DA", 
           size = 6, 

           colour = "black")+
  annotate("text", 
           x = 2, y = -0.05, 
           label = "HL", 
           size = 6, 
           colour = "black")+
  annotate("text", 
           x = 3, y = -0.05, 
           label = "HP", 
           size = 6, 
           colour = "black")+
  annotate("text", 
           x = 4, y = -0.05, 
           label = "PL", 
           size = 6, 
           colour = "black")


(barplot_pigment1 + barplot_pigment2) + 
  plot_annotation(
    title = "",
    caption = "Treatments") & theme(plot.caption = element_text(hjust = 0.5, size = 20),
                                    plot.title = element_text(size="20"))

# Average pigment composition across all samples
barplot_pig_all <- rbind(data_perc_sb, data_perc_wsb)
aggregate(Percentage ~ Pigments, data = barplot_pig_all, mean)

# Average pigment composition only WSB DA
data_perc_wsb_DA <- data_perc_wsb[data_perc_wsb$Conditions == "DA", ]
aggregate(Percentage ~ Pigments, data = data_perc_wsb_DA, mean)

# FIGURE S5 PIG BIPLOTS ####

# Pannel A
pigments_all <- df_pigment_perc %>% select(5:ncol(df_pigment_perc))
pigments_hell_all <- decostand(pigments_all, method = "hellinger")

# PCA
pca_pigments <- prcomp(pigments_hell_all, scale = T) 

# Add PCA results
data <- cbind(df_pigment_perc, pca_pigments$x)

# Graphical rescale
n<-45 

# Convex function
source("https://raw.githubusercontent.com/cmartin/ggConvexHull/refs/heads/master/R/geom_convexhull.R") 


biplot_all <- ggplot(data, aes(x = -PC1, y = -PC2, color = Conditions, shape = Groups)) +
  geom_point(size = 3) +
  
  scale_color_manual(values = c("DA" = "#324851", "HP" = "#5C88C4", "HL"="#E7B800", "PL"="#937DC2")) +
  scale_fill_manual(values = c("DA" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2"))+
  geom_convexhull(aes(fill = Conditions), alpha=0.6)+
  stat_ellipse(aes(fill = Groups), type = "norm", level = 0.95, alpha = 0.2, geom = "polygon", linetype = "dashed", color = "black") + 
  annotate("segment",
           x = 0, y = 0,
           xend = -pca_pigments$rotation[,1] * n,
           yend = -pca_pigments$rotation[,2] * n,
           arrow = arrow(length = unit(0.2, "cm")),
           color = "darkgrey",
           linewidth = 0.2) +
  geom_text_repel(
    data = as.data.frame(pca_pigments$rotation),
    aes(x = -PC1 * n, y = -PC2 * n, label = colnames(pigments_all)),
    inherit.aes = FALSE,
    color = "black",
    size = 5,
    segment.color = "grey30",
    segment.size = 0.3,
    max.overlaps = 100) +
  labs(
    x = paste0("PC 1 (", round(100 * pca_pigments$sdev[1]^2 / sum(pca_pigments$sdev^2), 1), "%)"),
    y = paste0("PC 2 (", round(100 * pca_pigments$sdev[2]^2 / sum(pca_pigments$sdev^2), 1), "%)")) +
  guides(
    color = guide_legend(override.aes = list(shape = NA)),
    shape = guide_legend(),
    fill = guide_legend(override.aes = list(shape = c(19, 17), color = "black", linetype = "dashed"))) +
  theme_test() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "bottom")+
  labs(title = "A") +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0))


# Anosim, diff between SB and WSB in their pigment compositions
coord_axes <- pca_pigments$x[,1]

# metadata factor
group_fac <- df_pigment_perc$Groups 

coord_axes <- data.frame(coord_axes, group = group_fac)

# Euclidean distance matrix
distance_matrix_pigment <- dist(coord_axes$coord_axes, method = "euclidean")
anosim(distance_matrix_pigment, coord_axes$group, permutations = 9999)

#Panel B SB

#Table
df_pigment_perc_sb <- subset(df_pigment_perc, Groups == "SB")

simper_sb <- df_pigment_perc_sb[, -c(1:4)]
simper_sb_var <- df_pigment_perc_sb[, -c(5:48)]
sim_sb <- simper(simper_sb, simper_sb_var$Conditions, permutations = 9999)
summary(sim_sb)

sim_sb_df <- bind_rows(
  lapply(names(sim_sb), function(name) {
    sim <- sim_sb[[name]]
    df <- as.data.frame(sim)
    df$species <- rownames(df)
    df <- df[, c("species", "average", "sd", "p")]
    df$contrast <- name
    df[, c("contrast", "species", "average", "sd", "p")]
  })
)

# Show only significant values
sim_sb_df <- sim_sb_df %>%filter(p > 1e-7 & p < 0.05) 

sim_sb_df %>%
  gt() %>%
  tab_header(
    title = "SIMPER SB results - Stress Contrast"
  ) %>%
  cols_label(
    contrast = "Contrast",
    species = "Species",
    average = "Mean difference",
    sd = "sd",
    p = "p-value"
  ) %>%
  fmt_number(
    columns = c(average, sd, p),
    decimals = 4
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything()))%>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title"))

# Pigment number impacted
sim_sb_df_contrast <- sim_sb_df %>%
  count(contrast)
sim_sb_df_contrast

#Plot
pigments_sb <- df_pigment_perc_sb %>% select(5:ncol(df_pigment_perc_sb))
pigments_hell_sb <- decostand(pigments_sb, method = "hellinger")
pca_pigments_sb <- prcomp(pigments_hell_sb, scale = TRUE)

data_sb <- cbind(df_pigment_perc_sb, pca_pigments_sb$x)

n<-45

source("https://raw.githubusercontent.com/cmartin/ggConvexHull/refs/heads/master/R/geom_convexhull.R") 


biplot_sb <- ggplot(data_sb, aes(x = -PC1, y = -PC2, color = Conditions)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("DA" = "#324851", "HP" = "#5C88C4", "HL"="#E7B800", "PL"="#937DC2")) +
  scale_fill_manual(values = c("DA" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2"))+
  geom_convexhull(aes(fill = Conditions), alpha=0.6)+
  annotate("segment",
           x = 0, y = 0,
           xend = -pca_pigments_sb$rotation[,1] * n,
           yend = -pca_pigments_sb$rotation[,2] * n,
           arrow = arrow(length = unit(0.2, "cm")),
           color = "darkgrey",
           linewidth = 0.2) +
  geom_text_repel(
    data = as.data.frame(pca_pigments_sb$rotation),
    aes(x = -PC1 * n, y = -PC2 * n, label = colnames(pigments_sb)),
    inherit.aes = FALSE,
    color = "black",
    size = 5,
    segment.color = "grey30",
    segment.size = 0.3,
    max.overlaps = 100) +
  labs(
    x = paste0("PC 1 (", round(100 * pca_pigments_sb$sdev[1]^2 / sum(pca_pigments_sb$sdev^2), 1), "%)"),
    y = paste0("PC 2 (", round(100 * pca_pigments_sb$sdev[2]^2 / sum(pca_pigments_sb$sdev^2), 1), "%)")) +
  guides(
    color = guide_legend(override.aes = list(shape = NA)),
    shape = guide_legend(),
    fill = guide_legend(override.aes = list(shape = c(19, 17), color = "black", linetype = "dashed"))) +
  theme_test() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none")+
  labs(title = "B") +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0))

coord_axes <- pca_pigments_sb$x[,1]
group_fac <- df_pigment_perc_sb$Conditions

coord_axes <- data.frame(coord_axes, group = group_fac)

distance_matrix_pigment <- dist(coord_axes$coord_axes, method = "euclidean")
anosim(distance_matrix_pigment, coord_axes$group, permutations = 9999)


#Pannel C WSB

#Table
df_pigment_perc_wsb <- subset(df_pigment_perc, Groups == "WSB")

simper_wsb <- df_pigment_perc_wsb[, -c(1:4)]
simper_wsb_var <- df_pigment_perc_wsb[, -c(5:48)]
sim_wsb <- simper(simper_wsb, simper_wsb_var$Conditions, permutations = 9999)
summary(sim_wsb)

sim_wsb_df <- bind_rows(
  lapply(names(sim_wsb), function(name) {
    sim <- sim_wsb[[name]]
    df <- as.data.frame(sim)
    df$species <- rownames(df)
    df <- df[, c("species", "average", "sd", "p")]
    df$contrast <- name
    df[, c("contrast", "species", "average", "sd", "p")]
  })
)

# Show only significant values
sim_wsb_df <- sim_wsb_df %>%filter(p > 1e-7 & p < 0.05) 

sim_wsb_df %>%
  gt() %>%
  tab_header(
    title = "SIMPER WSB results - Stress Contrast"
  ) %>%
  cols_label(
    contrast = "Contrast",
    species = "Species",
    average = "Mean difference",
    sd = "sd",
    p = "p-value"
  ) %>%
  fmt_number(
    columns = c(average, sd, p),
    decimals = 4
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything()))%>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title"))

# Pigment number impacted
sim_wsb_df_contrast <-sim_wsb_df %>%
  count(contrast)
sim_wsb_df_contrast

# Plot
pigments_wsb <- df_pigment_perc_wsb %>% select(5:ncol(df_pigment_perc_wsb))
pigments_hell_wsb <- decostand(pigments_wsb, method = "hellinger")
pca_pigments_wsb <- prcomp(pigments_hell_wsb, scale = TRUE) 

data_wsb <- cbind(df_pigment_perc_wsb, pca_pigments_wsb$x)

n<-25

biplot_wsb <- ggplot(data_wsb, aes(x = -PC1, y = -PC2, color = Conditions)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("DA" = "#324851", "HP" = "#5C88C4", "HL"="#E7B800", "PL"="#937DC2")) +
  scale_fill_manual(values = c("DA" = "#324851", "HL" = "#E7B75F", "HP" = "#5C88C4", "PL" = "#937DC2"))+
  geom_convexhull(aes(fill = Conditions), alpha=0.6)+
  annotate("segment",
           x = 0, y = 0,
           xend = -pca_pigments_wsb$rotation[,1] * n,
           yend = -pca_pigments_wsb$rotation[,2] * n,
           arrow = arrow(length = unit(0.2, "cm")),
           color = "darkgrey",
           linewidth = 0.2) +
  geom_text_repel(
    data = as.data.frame(pca_pigments_wsb$rotation),
    aes(x = -PC1 * n, y = -PC2 * n, label = colnames(pigments_wsb)),
    inherit.aes = FALSE,
    color = "black",
    size = 5,
    segment.color = "grey30",
    segment.size = 0.3,
    max.overlaps = 100) +
  labs(
    x = paste0("PC 1 (", round(100 * pca_pigments_wsb$sdev[1]^2 / sum(pca_pigments_wsb$sdev^2), 1), "%)"),
    y = paste0("PC 2 (", round(100 * pca_pigments_wsb$sdev[2]^2 / sum(pca_pigments_wsb$sdev^2), 1), "%)")) +
  guides(
    color = guide_legend(override.aes = list(shape = NA)),
    shape = guide_legend(),
    fill = guide_legend(override.aes = list(shape = c(19, 17), color = "black", linetype = "dashed"))) +
  theme_test() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "right")+
  labs(title = "C",
       color = "Treatments",
       fill  = "Treatments") +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0))

# Final graph
sim_wsb_df_contrast$Group <- "WSB"
sim_sb_df_contrast$Group <- "SB"
values_df <- rbind(sim_sb_df_contrast, sim_wsb_df_contrast)
values_df <- values_df[values_df$contrast %in% c("DA_HL", "DA_HP", "DA_PL"), ]
values_df$Percentage <- (values_df$n / 44) * 100
values_df <- values_df %>%rename(Contrast = contrast)

pig_impact <- ggplot(values_df, aes(x = Contrast, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("SB" = "#CC6677", "WSB" = "#81B29A")) +
  theme_test(base_size = 16) +
  labs(
    title = "C",
    x = "Contrast",
    y = "Percentage",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "bottom")

# Anosim, diff between stresses in WSB on their pigment compositions
coord_axes <- pca_pigments_wsb$x[,1]
group_fac <- df_pigment_perc_wsb$Conditions # metadata factor

coord_axes <- data.frame(coord_axes, group = group_fac)

# Euclidean distance matrix
distance_matrix_pigment <- dist(coord_axes$coord_axes, method = "euclidean")
anosim(distance_matrix_pigment, coord_axes$group, permutations = 9999)


(biplot_all) / (biplot_sb | biplot_wsb)


# FIGURE S6 BCA ####

bca_sb <- subset(df_pigment, Groups == "SB")
bca_sb_hell <- decostand(bca_sb[, -c(1:4)], method = "hellinger")

bca_wsb <- subset(df_pigment, Groups == "WSB")
bca_wsb_hell <- decostand(bca_wsb[, -c(1:4)], method = "hellinger")

source("https://raw.githubusercontent.com/cmartin/ggConvexHull/refs/heads/master/R/geom_convexhull.R")

# SB BCA
res.pca_sb <- dudi.pca(bca_sb_hell, scannf = FALSE, nf = 2) # PCA

res.bca.pigment_sb <- bca(res.pca_sb, as.factor(bca_sb$Conditions), scannf = FALSE, nf = 2)
res.bca.pigment_sb$ratio 

# Eig values
varexp2<-res.bca.pigment_sb$eig*100/sum(res.bca.pigment_sb$eig)

fac<-as.factor(bca_sb$Conditions)

# Plot
bca1 <- ggplot(res.bca.pigment_sb$ls,aes(x=res.bca.pigment_sb$ls[,1],y=res.bca.pigment_sb$ls[,2],col=fac))+
  geom_point()+
  scale_fill_manual(values=c("DA" = "#324851", "HP" = "#5C88C4", "HL" = "#E7B75F", "PL" = "#937DC2")) +
  scale_color_manual(values=c("DA" = "#324851", "HP" = "#5C88C4", "HL" = "#E7B75F", "PL" = "#937DC2")) +
  geom_convexhull(alpha = 0.3,aes(fill = fac))+
  xlab(paste("PC",1," : ",round(varexp2[1],1),"%"))+
  ylab(paste("PC",2," : ",round(varexp2[2],1),"%"))+
  labs(fill="Explanatory factors",col="Explanatory factors")+
  theme_test() +
  theme(axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom")+
  annotate("text", 
           x = -4.4, y = 5, 
           label = "Sediment biofilm (84%)", 
           size = 4, 
           colour = "black")


# WSB BCA
res.pca.ade4 <- dudi.pca(bca_wsb_hell, scannf = FALSE, nf = 2) # PCA
res.bca.pigment_wsb <- bca(res.pca.ade4, as.factor(bca_wsb$Conditions), scannf = FALSE, nf = 2)
res.bca.pigment_wsb$ratio 

varexp2<-res.bca.pigment_wsb$eig*100/sum(res.bca.pigment_wsb$eig)

fac<-as.factor(bca_wsb$Conditions)

# Plot
bca2 <- ggplot(res.bca.pigment_wsb$ls,aes(x=res.bca.pigment_wsb$ls[,1],y=res.bca.pigment_wsb$ls[,2],col=fac))+
  geom_point()+
  scale_fill_manual(values=c("DA" = "#324851", "HP" = "#5C88C4", "HL" = "#E7B75F", "PL" = "#937DC2")) +
  scale_color_manual(values=c("DA" = "#324851", "HP" = "#5C88C4", "HL" = "#E7B75F", "PL" = "#937DC2")) +
  geom_convexhull(alpha = 0.3,aes(fill = fac))+
  xlab(paste("PC",1," : ",round(varexp2[1],1),"%"))+
  ylab(paste("PC",2," : ",round(varexp2[2],1),"%"))+
  labs(fill="Explanatory factors",col="Explanatory factors")+
  theme_test() +
  
  theme(axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none")+
  annotate("text", 
           x = -3, y = 5, 
           label = "Without sediment biofilm (67%)", 
           size = 4, 
           colour = "black")

(bca1 + bca2) + 
  plot_annotation(
    title = "",
    caption = "(Treatments)") & theme(plot.caption = element_text(hjust = 0.5, size = 20),
                                    plot.title = element_text(size="20"))




# TABLE 1 MAT MET ####

mm <- read_excel("1_Data_formatting/data/mat_met.xlsx")

mm %>%
  gt() %>%
  
  fmt_missing(
    columns = everything(),
    missing_text = ""
  ) %>%
  
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = 1)
  ) %>%
  
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      rows = (nrow(mm)-1):nrow(mm),
      columns = -1  
    )
  ) %>%
  
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_fill(color = "lightgrey")
    ),
    locations = cells_column_labels(columns = everything()))

# TABLE 2 MIGRATION ####

# HL values
hl_var <- File_HL_data %>%
  select(-c(1:4)) %>%
  filter((minutes >= 11 & minutes <= 13) | (minutes >= 44 & minutes <= 47)) %>%  
  mutate(Time = case_when(
    minutes >= 11 & minutes <= 13 ~ "Before",
    minutes >= 44 & minutes <= 47 ~ "After"))

# F0 values mean before and after
hl_sum <- hl_var %>%
  group_by(Conditions, Time) %>%
  summarise(
    mean = mean(`F0`, na.rm = TRUE),
    sd = sd(`F0`, na.rm = TRUE),
    .groups = "drop")

# Format
hl_sum_wide <- hl_sum %>%
  pivot_wider(names_from = Time, values_from = c(mean, sd))

hl_sum <- hl_sum_wide %>%
  mutate(
    diff_After_Before = mean_After - mean_Before,
    percent_change = (diff_After_Before / mean_Before) * 100,
    diff_After_Before = sprintf("%+.0f", diff_After_Before),
    percent_change = sprintf("%+.0f%%", percent_change),
    Before = sprintf("%.0f (%.0f)", mean_Before, sd_Before),
    After = sprintf("%.0f (%.0f)", mean_After, sd_After)
  ) %>%
  select(Conditions, Before, After, diff_After_Before, percent_change)

# Add DA (HL) tag
hl_sum <- hl_sum %>%
  mutate(Conditions = if_else(str_detect(Conditions, "DA"), 
                          paste0(Conditions, " (1)"), 
                          Conditions))

# HP values
hp_var <- File_HP_data %>%
  select(-c(1:4)) %>%
  filter(minutes >= 24 & minutes <= 33) %>%
  mutate(Time = case_when(
    minutes >= 24 & minutes <= 26 ~ "Before",
    minutes >= 30 & minutes <= 33 ~ "After"))

hp_sum <- hp_var %>%
  group_by(Conditions, Time) %>%
  summarise(
    mean = mean(`F0`, na.rm = TRUE),
    sd = sd(`F0`, na.rm = TRUE),
    .groups = "drop")

hp_sum_wide <- hp_sum %>%
  pivot_wider(names_from = Time, values_from = c(mean, sd))

hp_sum <- hp_sum_wide %>%
  mutate(
    diff_After_Before = mean_After - mean_Before,
    percent_change = (diff_After_Before / mean_Before) * 100,
    diff_After_Before = sprintf("%+.0f", diff_After_Before),
    percent_change = sprintf("%+.0f%%", percent_change),
    Before = sprintf("%.0f (%.0f)", mean_Before, sd_Before),
    After = sprintf("%.0f (%.0f)", mean_After, sd_After)
  ) %>%
  select(Conditions, Before, After, diff_After_Before, percent_change)

hp_sum <- hp_sum %>%
  mutate(Conditions = if_else(str_detect(Conditions, "DA"), 
                          paste0(Conditions, " (2)"), 
                          Conditions))

# PL values
pl_var <- File_PL_data %>%
  select(-c(1:4)) %>%
  filter(minutes >= 17 & minutes <= 34) %>%
  mutate(Time = case_when(
    minutes >= 17 & minutes <= 19 ~ "Before",
    minutes >= 30 & minutes <= 34 ~ "After"))

pl_sum <- pl_var %>%
  group_by(Conditions, Time) %>%
  summarise(
    mean = mean(`F0`, na.rm = TRUE),
    sd = sd(`F0`, na.rm = TRUE),
    .groups = "drop")

pl_sum_wide <- pl_sum %>%
  pivot_wider(names_from = Time, values_from = c(mean, sd))

pl_sum <- pl_sum_wide %>%
  mutate(
    diff_After_Before = mean_After - mean_Before,
    percent_change = (diff_After_Before / mean_Before) * 100,
    diff_After_Before = sprintf("%+.0f", diff_After_Before),
    percent_change = sprintf("%+.0f%%", percent_change),
    Before = sprintf("%.0f (%.0f)", mean_Before, sd_Before),
    After = sprintf("%.0f (%.0f)", mean_After, sd_After)
  ) %>%
  select(Conditions, Before, After, diff_After_Before, percent_change)

pl_sum <- pl_sum %>%
  mutate(Conditions = if_else(str_detect(Conditions, "DA"), 
                          paste0(Conditions, " (3)"), 
                          Conditions))

# Post-stress resilience (60 mins), slope

# HL 1 hour recovery
HL_slope <- File_HL_data[File_HL_data$minutes >= 43 & File_HL_data$minutes <= 105, ]

# Slope linear model
df_HL_slope <- HL_slope %>%
  group_by(Conditions) %>%
  do(modele = lm(`F0` ~ minutes, data = .)) %>%
  mutate(Slope = sprintf("%+.0f", coef(modele)[["minutes"]]))

ggplot(HL_slope, aes(x = minutes, y = F0, color = Conditions)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

# HP
HP_slope <- File_HP_data[File_HP_data$minutes >= 32 & File_HP_data$minutes <= 90, ]

df_HP_slope <- HP_slope %>%
  group_by(Conditions) %>%
  do(modele = lm(`F0` ~ minutes, data = .)) %>%
  mutate(Slope = sprintf("%+.0f", coef(modele)[["minutes"]]))

ggplot(HP_slope, aes(x = minutes, y = F0, color = Conditions)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

# PL
PL_slope <- File_PL_data[File_PL_data$minutes >= 30 & File_PL_data$minutes <= 93, ]

df_PL_slope <- PL_slope %>%
  group_by(Conditions) %>%
  do(modele = lm(`F0` ~ minutes, data = .)) %>%
  mutate(Slope = sprintf("%+.0f", coef(modele)[["minutes"]]))

ggplot(PL_slope, aes(x = minutes, y = F0, color = Conditions)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

# merge df_slope
df_slope <- bind_rows(df_HL_slope, df_HP_slope, df_PL_slope)

# Edit table
compiled_table <- bind_rows(hl_sum, hp_sum, pl_sum)
compiled_table$Slope <- df_slope$Slope # Add slope
compiled_table <- compiled_table %>% rename(Treatments = Conditions)

compiled_table <- compiled_table %>%
  rename(
    `F0 After (a.u.)` = `After`,
    `F0 Before (a.u.)` = `Before`,
    `Δ After-Before (a.u.)` = `diff_After_Before`,
    `Changes (%)` = `percent_change`,
    `Recovery slope (F0.min-1)` = `Slope`)

compiled_table %>% 
  gt() %>%
  
  cols_width(
    everything() ~ px(200)
  ) %>%
  
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = c(last_col(), last_col(offset = 1)))
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = 1)
  ) %>%
  
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_fill(color = "lightgrey")
    ),
    locations = cells_column_labels(columns = everything()))

# Export table
write.xlsx(compiled_table, file = "/Users/alexandre/Downloads/compiled_table.xlsx")

# TABLE S2 PIGMENT DIV ####

pig_tab <- read_excel("1_Data_formatting/data/Pigments_data/Pigment_diversity.xlsx")

# Table
pig_tab %>%
  gt() %>%
  
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%

  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_fill(color = "lightgrey")
    ),
    locations = cells_column_labels(columns = everything())
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = 1))


# TABLE S3 RAW PPARAMETERS ####
df_Pp_raw <- df_Pparameters %>%select(Time, Treatment, Day, `α-slope`, `Fv/Fm`, `rETRm obs`, NPQ1000, YNPQ1000, `YNOm obs`, F0, Eopt, Ek)

df_Pp_raw$Treatment <- sapply(strsplit(df_Pp_raw$Treatment, "_"), function(x) paste(rev(x), collapse = "_"))

df_Pp_raw$Treatment <- ifelse(
  grepl("DA", df_Pp_raw$Treatment),
  paste0(df_Pp_raw$Treatment, 
         ifelse(df_Pp_raw$Day == 1, "_⁽¹⁾",
                ifelse(df_Pp_raw$Day == 2, "_⁽³⁾",
                       ifelse(df_Pp_raw$Day == 3, "_⁽²⁾", "")))),
df_Pp_raw$Treatment)

# Treatment name
df_Pp_raw$Treatment <- gsub("_", " ", df_Pp_raw$Treatment)

df_Pp_raw$Treatment <- paste(df_Pp_raw$Treatment, df_Pp_raw$Time)

# Df structure
treatm_order <- c("DA SB ⁽¹⁾ Before", "DA SB ⁽¹⁾ After", "HL SB Before", "HL SB After", 
                  "DA WSB ⁽¹⁾ Before", "DA WSB ⁽¹⁾ After", "HL WSB Before", "HL WSB After",
                  "DA SB ⁽²⁾ Before", "DA SB ⁽²⁾ After", "HP SB Before", "HP SB After", 
                  "DA WSB ⁽²⁾ Before", "DA WSB ⁽²⁾ After", "HP WSB Before", "HP WSB After", 
                  "DA SB ⁽³⁾ Before", "DA SB ⁽³⁾ After", "PL SB Before", "PL SB After", 
                  "DA WSB ⁽³⁾ Before", "DA WSB ⁽³⁾ After", "PL WSB Before", "PL WSB After")

# Mean sd Chlalculation
Pp_table <- df_Pp_raw %>%
  mutate(Treatment = factor(Treatment, levels = treatm_order)) %>%
  group_by(Treatment) %>%
  summarise(
    `α-slope` = sprintf("%.2f (%.2f)", mean(`α-slope`), sd(`α-slope`)),
    `Fv/Fm` = sprintf("%.2f (%.2f)", mean(`Fv/Fm`), sd(`Fv/Fm`)),
    `rETRm obs` = sprintf("%.2f (%.2f)", mean(`rETRm obs`), sd(`rETRm obs`)),
    NPQ1000 = sprintf("%.2f (%.2f)", mean(NPQ1000), sd(NPQ1000)),
    YNPQ1000 = sprintf("%.2f (%.2f)", mean(YNPQ1000), sd(YNPQ1000)),
    `YNOm obs` = sprintf("%.2f (%.2f)", mean(`YNOm obs`), sd(`YNOm obs`)),
    F0 = sprintf("%.2f (%.2f)", mean(F0), sd(F0)),
    Eopt = sprintf("%.2f (%.2f)", mean(Eopt), sd(Eopt)),
    Ek = sprintf("%.2f (%.2f)", mean(Ek), sd(Ek)))

# Table
Pp_table %>%
  gt() %>%
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_fill(color = "lightgrey")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = 1)
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_fill(color = "grey92"),
    locations = cells_body(
      rows = grepl("WSB", Treatment)
    )
  )

# TABLE S4 PPARAMETERS CHANGES AF-BF ####

Pp_table_av <- Pp_table  

for (col in names(Pp_table_av)) {
  Pp_table_av[[col]] <- gsub(" \\([^)]*\\)", "", Pp_table_av[[col]])
}

idx_before <- seq(1, nrow(Pp_table_av), by = 2)
idx_after  <- seq(2, nrow(Pp_table_av), by = 2)

Pp_table_av_changes <- Pp_table_av[idx_after, "Treatment", drop = FALSE]

for (j in 2:10) {
  colname <- names(Pp_table_av)[j]
  
  before_vals <- as.numeric(sapply(Pp_table_av[idx_before, j], as.character))
  after_vals  <- as.numeric(sapply(Pp_table_av[idx_after,  j], as.character))
  
  Pp_table_av_changes[[colname]] <- (after_vals * 100 / before_vals) - 100
}

# Formatting
Pp_table_av_changes$Treatment <- gsub(" After$", "", Pp_table_av_changes$Treatment)

Pp_table_av_changes <- Pp_table_av_changes %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Table
Pp_table_av_changes %>%
  gt() %>%
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_fill(color = "lightgrey")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = 1)
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_fill(color = "grey92"),
    locations = cells_body(
      rows = grepl("WSB", Treatment))) 



# TABLE S5 XANTHOPHYLLS DATA ####
gdf_pigment_concentration <-gdf_pigment

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(cDd = (Dd / (`Chla-derivatives` + Chla + Pheopigments)))

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(cDt = (Dt / (`Chla-derivatives` + Chla + Pheopigments)))

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(`cDd+Dt` = ((Dd / (`Chla-derivatives` + Chla + Pheopigments)) + (Dt / (`Chla-derivatives` + Chla + Pheopigments))))

gdf_pigment_concentration <- gdf_pigment_concentration %>%
  mutate(DES = (Dt / (Dd + Dt)) * 100)

gdf_pigment_concentration <- gdf_pigment_concentration[, c("Treatment", "cDd", "cDt", "cDd+Dt", "DES")]

# Treatment metadata
gdf_pigment_concentration$Treatment <- gsub("_", " ", gdf_pigment_concentration$Treatment)

# Df structure
treatm_order <- c("DA SB", "HL SB", "HP SB", "PL SB", "DA WSB", "HL WSB", "HP WSB", "PL WSB")

# Mean sd Chlalculation
xc_table <- gdf_pigment_concentration %>%
  mutate(Treatment = factor(Treatment, levels = treatm_order)) %>%
  group_by(Treatment) %>%
  summarise(
    `[Diadinoxanthin]` = sprintf("%.1e (%.1e)", mean(cDd), sd(cDd)),
    `[Diatoxanthin]` = sprintf("%.1e (%.1e)", mean(cDt), sd(cDt)),
    `[Diadinoxanthin + Diatoxanthin]` = sprintf("%.1e (%.1e)", mean(`cDd+Dt`), sd(`cDd+Dt`)),
    `DES%` = sprintf("%.1f (%.1f)", mean(DES), sd(DES)))

# Table
xc_table %>%
  gt() %>%
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_fill(color = "lightgrey")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = 1)
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_column_labels(columns = everything()))