adjImm_threshold <- 0.93 ## threshold for contact-adjusted immunity

scenarios <- scenarios_adjImm %>%
    mutate(country=factor(country),
           scenario=factor(scenario)) %>%
    arrange(country) %>%
    gather(model, value, ends_with("_immunity")) %>%
    mutate(model=recode_factor(model,
                               adjusted_immunity="age-specific",
                               mean_immunity="homogeneous"))

countries <-
  grep("\\(", invert=TRUE, unique(as.character(scenarios$country)), value=TRUE)

change_order <- c("decrease", "increase", "none")

age_levels <- c("0-4", "5-9", "10-14", "15-19", "20+")
target_levels <- tibble(age = factor(age_levels, age_levels),
                        target = c(0.85, 0.9, 0.95, 0.95, 0.95),
                        scenario = "current",
                        change = factor("none", change_order))

change_levels <-
  change_levels <- list(increase_infant = c(0.05, 0, 0, 0, 0),
                        catchup_kids = c(0, 0.05, 0, 0, 0),
                        catchup_kids_less_teenagers = c(0, 0.05, -0.05, 0, 0),
                        catchup_kids_less_adolescents = c(0, 0.05, 0, -0.05, 0))

levels <- list()
for (changed in names(change_levels))
{
  subtract_levels <- change_levels[[changed]]
  subtract_levels[subtract_levels > 0] <- 0
  add_levels <- change_levels[[changed]]
  add_levels[add_levels < 0] <- 0
  changed_target_levels <- target_levels$target + subtract_levels
  levels[[changed]] <-
    rbind(tibble(age = target_levels$age,
                 target = -subtract_levels,
                 scenario = changed,
                 change = factor("decrease", change_order)),
          tibble(age = target_levels$age,
                 target = add_levels,
                 scenario = changed,
                 change = factor("increase", change_order)))
  levels[[changed]] <-
    rbind(target_levels %>% mutate(scenario = changed,
                                   target = changed_target_levels),
          levels[[changed]])
}

all_levels <- rbind(target_levels, bind_rows(levels)) %>%
  mutate(scenario = factor(scenario,
                           levels = c("current", names(change_levels))))

age_specific_scenarios <- scenarios %>%
  filter(model == "age-specific",
         scenario %in% unique(all_levels$scenario),
         settings=="all") %>%
  mutate(scenario = factor(scenario, levels(all_levels$scenario)),
         country = factor(country, levels=countries))

p_imm <- ggplot(all_levels %>%
                mutate(scenario =
                         factor(scenario,
                                labels = letters[1:length(unique(scenario))])),
                aes(x = age, y = target, fill = change)) +
  geom_bar(stat = "identity", color = "black", position="stack") +
  scale_y_continuous("immunity", label=scales::percent_format(accuracy=1),
                     limits = c(0, 1)) +
  facet_grid( ~ scenario) +
  coord_cartesian(ylim = c(0.8, 1)) +
  theme(strip.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 1.94), "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_x_discrete("") +
  scale_fill_manual(values = c("white", "lightblue", "black"))

eI_breaks <- c(0.9, 0.95)
eI_limits <- eI_breaks[c(1, length(eI_breaks))] * c(0.99, 1.01)
p_eI <- ggplot(age_specific_scenarios) +
  facet_grid( ~ scenario) +
  geom_boxplot(aes(x = factor(country, rev(levels(country))), y = value)) +
  geom_hline(yintercept=adjImm_threshold, linetype="dashed") +
  xlab("") +
  scale_y_continuous("Contact-adjusted immunity",
                     breaks = eI_breaks,
                     labels = scales::percent_format(accuracy=1),
                     limits = eI_limits) +
  ## geom_vline(xintercept=length(regions)+0.5, linetype="dashed") +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  coord_flip()

p <- plot_grid(p_imm, p_eI, nrow = 2, rel_heights=c(0.7, length(countries) / 10))
ggsave("figure_3.pdf", p, width = 15, height = length(countries) / 2)

## in the text: countries with >10% of outbreaks under current levels
current_elim_prob <- age_specific_scenarios %>%
    filter(scenario=="current") %>%
    group_by(country) %>%
    summarise(prob=mean(value < 0.93)) %>%
    filter(prob > 0.1) %>%
    mutate(prob=formattable::percent(prob, digits=0))

increase_infant_elim_prob <- age_specific_scenarios %>%
    filter(scenario=="increase_infant") %>%
    group_by(country) %>%
    summarise(prob=mean(value < 0.93)) %>%
    filter(prob > 0.1) %>%
    mutate(prob=formattable::percent(prob, digits=0))

catchup_kids_elim_prob <- age_specific_scenarios %>%
    filter(scenario=="catchup_kids") %>%
    group_by(country) %>%
    summarise(prob=mean(value < 0.93)) %>%
    filter(prob > 0.1) %>%
    mutate(prob=formattable::percent(prob, digits=1))
