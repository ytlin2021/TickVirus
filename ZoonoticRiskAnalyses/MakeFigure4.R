## =================================================================================================
## Plot known tick viruses with complete genomes
## =================================================================================================
library(readr)
library(ggplot2)
library(cowplot)
library(tidyverse)

colour_palatte <- c("#D63163", "#FB9A29", "#FEC44F", "#019c64")
#colour_palatte <- c("#CC4C02", "#FB9A29", "#FEC44F", "#FFF7BC")
setwd("~/Project-TickVirus/SharedDataWithDavid/tickvirus-github/ZoonoticRiskAnalyses")

taxonomy <- read_csv('TBVCompleteGenome.csv') %>% 
	select("Virus", "Acronym", "Family", "Know to infect human?") %>% 
	rename("Human" = "Know to infect human?")

## Load predictions (known TBVs with complete gemones)
pred_col_types <- cols(Name = 'c', 
                       Acronym = 'c',
                       bagged_prediction = 'c', 
                       priority_category = 'c', 
                       .default = 'd')

predictions_TBVCompleteGenome <- read_csv(file.path('TBVCompleteGenome.predictions.csv'),
                              col_types = pred_col_types)

## Tick vectors of viruses
hosts <- read_csv("TickVector.csv") %>% 
  distinct() 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Ranks --------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Rank
TBVCompleteGenome_ranks <- predictions_TBVCompleteGenome %>% 
  mutate(Rank = rank(-.data$calibrated_score_mean, ties.method = 'min')) %>% 
  arrange(-.data$Rank) %>% 
  mutate(Acronym = factor(.data$Acronym, levels = .data$Acronym),
         Priority = factor(.data$priority_category, 
                           levels = rev(c('Low', 'Medium', 'High', 'Very high'))))

host_data <- TBVCompleteGenome_ranks %>% 
  left_join(hosts, by = c("Name" = "Virus")) %>% 
  mutate(Acronym = factor(.data$Acronym, levels = levels(TBVCompleteGenome_ranks$Acronym)))

# Determine zoom region:
# - Zoomed plot shows only non-human-associated viruses, so need to find the area covered by 
#   the top 25 viruses sampled from other species
Human <- host_data %>% 
  mutate(Acronym = as.character(.data$Acronym)) %>%   # Needed to avoid losing factor ordering during join to taxonomy
  left_join(taxonomy, by = "Acronym") %>% 
  filter(.data$Human == "Y")  
  
zoom_cutoff <- 25
zoom_start <- TBVCompleteGenome_ranks %>% 
  filter(!.data$Acronym %in% Human$Acronym) %>% 
  top_n(n = 1, wt = -.data$Rank) %>%  # The highest-ranked virus among the top 25
  pull(.data$Acronym)

zoom_stop <- TBVCompleteGenome_ranks %>% 
  filter(!.data$Acronym %in% Human$Acronym) %>% 
  top_n(n = zoom_cutoff, wt = -.data$Rank) %>% 
  top_n(n = 1, wt = .data$Rank) %>%  # The lowest-ranked virus among the top 25
  pull(.data$Acronym)

## Plots
rank_plot <- ggplot(TBVCompleteGenome_ranks, aes(x = Acronym, y = calibrated_score_mean, colour = Priority)) +
	geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf), colour = NA, fill = 'grey90') +
	geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.6) +
	geom_step(group = 1, colour = 'grey10', size = 0.5) +
  geom_rug(sides = 't', colour = 'grey10', data = Human) +
	geom_hline(yintercept = 0.293, colour = 'grey10', linetype = 2) +
	scale_y_continuous(limits = c(0, 1)) +
	scale_colour_manual(values = colour_palatte, guide = "none") +
	labs(x = NULL, y = "Predicted probability") +
	theme_linedraw() +
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.text.x = element_text(angle = 90, hjust = 0, size = 6),
				axis.ticks.x = element_blank(),
				axis.title.y = element_text(size = 10),
				panel.grid.major.x = element_blank(),
				panel.grid.minor.x = element_blank(),
				plot.margin = margin(t = 5.5, r = 5.5, b = 1, l = 5.5))+
  scale_x_discrete(position = "top") 

indicator_plot <- ggplot(TBVCompleteGenome_ranks, aes(x = Acronym, y = Priority, fill = Priority)) +
	geom_blank() +
	geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf), fill = 'grey90', colour = NA) +
	geom_tile() +
	scale_x_discrete(expand = expansion(add = 0)) +
	scale_fill_manual(values = colour_palatte, guide = "none") +
	labs(x = NULL, y = "Zoonotic\npotential") +
	theme_linedraw() +
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				axis.title.y = element_text(size = 10),
				panel.grid = element_blank(),
				plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5))

# Show hosts as another indicator plot:
source_data <- host_data %>% 
  group_by(.data$TickGenera) %>% 
  mutate(N = n()) %>% 
  ungroup() %>% 
  mutate(facet = case_when(is.na(.data$TickFamily) ~ "D",
                           .data$TickFamily == "Ixodidae" ~ "A",
                           .data$TickFamily == "Argasidae" ~ "B",
                           TRUE ~ "C"),
         source = case_when(TRUE ~ .data$TickGenera)) %>% 
  select(.data$Acronym, .data$source, .data$facet, .data$Rank, .data$Priority)

source_name_genera <- sort(unique(source_data$source), decreasing = TRUE)
source_name_genera <- c("Unknown", 
                       source_name_genera[source_name_genera != "Unknown"])

source_data <- source_data %>% 
  mutate(source = factor(.data$source, levels = source_name_genera))

source_plot <- ggplot(source_data, aes(x = Acronym, y = source, fill = Priority)) +
  geom_blank() +
  geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf), fill = 'grey90', colour = NA) +
  geom_tile() +
  facet_grid(rows = vars(facet), scales = "free", space = "free") + 
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0)) +
  scale_fill_manual(values = colour_palatte, guide = "none") +
  labs(x = NULL, y = "Tick genera") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = 'italic'),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(-0.5, 'pt'),
        plot.margin = margin(t = 0, r = 5.5, b = -2, l = 5.5))

# Show families as another indicator plot:
family_data <- host_data %>% 
  mutate(Acronym = as.character(.data$Acronym)) %>%   # Needed to avoid losing factor ordering during join to taxonomy
  left_join(taxonomy, by = "Acronym") %>% 
  group_by(.data$Family) %>% 
  mutate(N = n()) %>% 
  ungroup() %>% 
  filter(!duplicated(Acronym)) %>%
  select(.data$Acronym, .data$Family, .data$Rank, .data$Priority)

family_name <- sort(unique(family_data$Family), decreasing = TRUE)

family_data <- family_data %>% 
  mutate(Family = factor(.data$Family, levels = family_name),
         Acronym = factor(.data$Acronym, levels = levels(TBVCompleteGenome_ranks$Acronym)))


family_plot <- ggplot(family_data, aes(x = Acronym, y = Family, fill = Priority)) +
  geom_blank() +
  geom_rect(aes(xmin = 94, xmax = 130, ymin = -Inf, ymax = Inf), fill = 'grey90', colour = NA) +
  geom_tile() +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0)) +
  scale_fill_manual(values = colour_palatte, guide = "none") +
  labs(x = NULL, y = "Virus family") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = 'italic'),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(-0.5, 'pt'),
        plot.margin = margin(t = 0, r = 5.5, b = -2, l = 5.5))

# Combine
rank_plot_combined <- plot_grid(rank_plot, indicator_plot, family_plot, source_plot, 
																ncol = 1, rel_heights = c(3, 1, 2.5, 2),
																align = 'v', axis = 'lr')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Zoom in and label top non-human-associated viruses ---------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
top_viruses <- TBVCompleteGenome_ranks %>% 
  filter(!.data$Name %in% Human$Name) %>% 
  top_n(n = zoom_cutoff, wt = -.data$Rank) %>%  
	mutate(Virus = as.character(.data$Name)) %>% # Needed to avoid losing factor ordering during join to taxonomy
	left_join(taxonomy, by = 'Virus') 

top_host_genera <- host_data %>% 
  filter(.data$Name %in% top_viruses$Name) %>% 
  mutate(facet = case_when(.data$TickFamily == "Ixodidae" ~ "Hard ticks",
                           .data$TickFamily == "Argasidae" ~ "Soft ticks"),
         source = if_else(is.na(.data$TickGenera), "Unknown Aegasidae", .data$TickGenera),
         facet = factor(.data$facet, levels = c("Hard ticks", "Soft ticks")),
         source = factor(.data$source, levels = sort(unique(.data$source))))

top_virus_plot <-	ggplot(top_viruses, aes(x = reorder(Name, desc(Rank)), y = calibrated_score_mean, fill = Priority)) +
  geom_hline(yintercept = 0.293, colour = 'grey60', linetype = 2) +
  geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.5, colour = 'grey10') +
  geom_point(shape = 22, size = 2) +
  scale_fill_manual(values = colour_palatte, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = c(0.02, 0.02))) +
  labs(x = NULL, y = "Predicted probability") +
  coord_flip() +
	theme_linedraw() +
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.title.x = element_text(size = 10),
				panel.grid.major.y = element_line(colour = 'grey92'),
				plot.margin = margin(t = 5.5, r = 1, b = 5.5, l = 0))

family_indicator_plot <- ggplot(top_viruses, aes(x = Family, y = reorder(Name, desc(Rank)), fill = Priority)) +
	geom_tile(colour = NA) + 
	scale_fill_manual(values = colour_palatte, guide = "none") +
	xlab("Virus family") +
	theme_linedraw() +
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.x = element_text(size = 10),
				axis.text.x = element_text(face = 'italic', angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.65),
				panel.grid.major.y = element_line(colour = 'grey92'),
				plot.margin = margin(t = 5.5, r = 1, b = 5.5, l = 0))

host_indicator_plot <- ggplot(top_host_genera, aes(x = source, y = reorder(Name, desc(Rank)), fill = Priority)) +
  geom_tile(colour = NA) + 
  scale_fill_manual(values = colour_palatte, guide = "none") +
  facet_grid(cols = vars(facet), scales = "free", space = "free") + 
  xlab("Tick genera") +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = 'italic', angle = 90, hjust = 1, vjust = 0.5, lineheight = 0.65),
        panel.grid.major.y = element_line(colour = 'grey92'),
        strip.background = element_blank(),
        panel.spacing = unit(-0.5, 'pt'),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0))

# Combine
bottom_row <- plot_grid(top_virus_plot, family_indicator_plot, host_indicator_plot, 
												ncol = 3, rel_widths = c(2, 1.3, 1.3),
												align = 'h', axis = 'tb',
												labels = c('B', '', ''))


final_plot <- plot_grid(rank_plot_combined, bottom_row,
												nrow = 2, rel_heights = c(1.5, 1.25),
												labels = c('A', ''))



ggsave2(file.path('Figure4TBVZoonoticPotential_ExcludeHuman.png'), final_plot, width = 11, height = 11)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Values mentioned in text -------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cat('\nNovel viruses found in humans:\n')
deduped.human <- unique(Human[ , 1:8, 11:12 ] )
human_tab <- table(deduped.human$Priority)

print(human_tab)
print(human_tab/sum(human_tab))

cat('\nNovel viruses with no link to humans in the very high priority category:\n')
TBVCompleteGenome_ranks %>%
  filter(!.data$Name %in% Human$Name) %>%
  pull(.data$Priority) %>%
  table() %>%
  print()

high_priority <- TBVCompleteGenome_ranks %>%
  filter(!.data$Name %in% Human$Name) %>%
  filter(.data$Priority %in% c('Very high', 'High')) %>% 
  left_join(taxonomy, by = 'Name')

high_priority_tab <- table(high_priority$Family) 

print(high_priority_tab)
print(high_priority_tab/sum(high_priority_tab))
