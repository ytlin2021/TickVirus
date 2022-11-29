## =================================================================================================
## Plot known tick viruses with complete genomes
## =================================================================================================
library(readr)
library(ggplot2)
library(cowplot)
library(tidyverse)

setwd("~/TickVirus/ZoonoticRiskAnalyses")

taxonomy <- read_csv('TBVCompleteGenome.csv') %>% 
	select("Virus", "Family", "Know to infect human?") %>% 
	rename("Human" = "Know to infect human?")

## Load predictions (known TBVs with complete gemones)
pred_col_types <- cols(Name = 'c', 
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
  mutate(Name = factor(.data$Name, levels = .data$Name),
         Priority = factor(.data$priority_category, 
                           levels = rev(c('Low', 'Medium', 'High', 'Very high'))))

host_data <- TBVCompleteGenome_ranks %>% 
  left_join(hosts, by = c("Name" = "Virus")) %>% 
  mutate(Name = factor(.data$Name, levels = levels(TBVCompleteGenome_ranks$Name)))

# Determine zoom region:
# - Zoomed plot shows only non-human-associated viruses, so need to find the area covered by 
#   the top 25 viruses sampled from other species
Human <- host_data %>% 
  mutate(Name = as.character(.data$Name)) %>% # Needed to avoid losing factor ordering during join to taxonomy
  left_join(taxonomy, by = c("Name" = "Virus")) %>% 
  filter(.data$Human == "Y")

zoom_cutoff <- 25
zoom_start <- TBVCompleteGenome_ranks %>% 
  filter(!.data$Name %in% Human$Name) %>% 
  top_n(n = 1, wt = -.data$Rank) %>%  # The highest-ranked virus among the top 25
  pull(.data$Name)

zoom_stop <- TBVCompleteGenome_ranks %>% 
  filter(!.data$Name %in% Human$Name) %>% 
  top_n(n = zoom_cutoff, wt = -.data$Rank) %>% 
  top_n(n = 1, wt = .data$Rank) %>%  # The lowest-ranked virus among the top 25
  pull(.data$Name)

## Plots
rank_plot <- ggplot(TBVCompleteGenome_ranks, aes(x = Name, y = calibrated_score_mean, colour = Priority)) +
	geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf), colour = NA, fill = 'grey90') +
	geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.6) +
	geom_step(group = 1, colour = 'grey10', size = 0.5) +
  geom_rug(sides = 't', colour = 'grey10', data = Human) +
	geom_hline(yintercept = 0.293, colour = 'grey10', linetype = 2) +
	scale_y_continuous(limits = c(0, 1)) +
	scale_colour_manual(values = c("#CC4C02", "#FB9A29", "#FEC44F", "#FFF7BC"), guide = "none") +
	labs(x = NULL, y = "Predicted probability") +
	theme_linedraw() +
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				axis.title.y = element_text(size = 10),
				panel.grid.major.x = element_blank(),
				panel.grid.minor.x = element_blank(),
				plot.margin = margin(t = 5.5, r = 5.5, b = 1, l = 5.5))

indicator_plot <- ggplot(TBVCompleteGenome_ranks, aes(x = Name, y = Priority, fill = Priority)) +
	geom_blank() +
	geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf), fill = 'grey90', colour = NA) +
	geom_tile() +
	scale_x_discrete(expand = expansion(add = 0)) +
	scale_fill_manual(values = c("#CC4C02", "#FB9A29", "#FEC44F", "#FFF7BC"), guide = "none") +
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
  select(.data$Name, .data$source, .data$facet, .data$Rank, .data$Priority)

source_name_genera <- sort(unique(source_data$source), decreasing = TRUE)
source_name_genera <- c("Unknown", 
                       source_name_genera[source_name_genera != "Unknown"])

source_data <- source_data %>% 
  mutate(source = factor(.data$source, levels = source_name_genera))


source_plot <- ggplot(source_data, aes(x = Name, y = source, fill = Priority)) +
  geom_blank() +
  geom_rect(aes(xmin = zoom_start, xmax = zoom_stop, ymin = -Inf, ymax = Inf), fill = 'grey90', colour = NA) +
  geom_tile() +
  
  facet_grid(rows = vars(facet), scales = "free", space = "free") + 
  
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0)) +
  scale_fill_manual(values = c("#CC4C02", "#FB9A29", "#FEC44F", "#FFF7BC"), guide = "none") +
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

# Combine
rank_plot_combined <- plot_grid(rank_plot, indicator_plot, source_plot, 
																ncol = 1, rel_heights = c(3, 1.15, 3),
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

top_virus_plot <-	ggplot(top_viruses, aes(x = Name, y = calibrated_score_mean, fill = Priority)) +
  geom_hline(yintercept = 0.293, colour = 'grey60', linetype = 2) +
  geom_errorbar(aes(ymin = calibrated_score_lower, ymax = calibrated_score_upper), width = 0.5, colour = 'grey10') +
  geom_point(shape = 22, size = 2) +
  scale_fill_manual(values = c("#CC4C02", "#FB9A29", "#FEC44F", "#FFF7BC"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = c(0.02, 0.02))) +
  labs(x = NULL, y = "Predicted probability") +
  coord_flip() +
	theme_linedraw() +
	theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.title.x = element_text(size = 10),
				panel.grid.major.y = element_line(colour = 'grey92'),
				plot.margin = margin(t = 5.5, r = 1, b = 5.5, l = 0))

family_indicator_plot <- ggplot(top_viruses, aes(x = Family, y = Name, fill = Priority)) +
	geom_tile(colour = NA) + 
	scale_fill_manual(values = c("#CC4C02", "#FB9A29", "#FEC44F", "#FFF7BC"), guide = "none") +
	xlab("Family") +
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

host_indicator_plot <- ggplot(top_host_genera, aes(x = source, y = Name, fill = Priority)) +
  geom_tile(colour = NA) + 
  scale_fill_manual(values = c("#CC4C02", "#FB9A29", "#FEC44F", "#FFF7BC"), guide = "none") +
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
												nrow = 2, rel_heights = c(1, 1.25),
												labels = c('A', ''))


ggsave2(file.path('Figure4TBVZoonoticPotential_ExcludeHuman.pdf'), final_plot, width = 8.5, height = 9)

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
  left_join(taxonomy, by = c('Name' = 'Virus'))

high_priority_tab <- table(high_priority$Family) 

print(high_priority_tab)
print(high_priority_tab/sum(high_priority_tab))
