data(ms_world) ## measles case data
data(contact_countries) ## mapping of countries in ESEN2 to countries in POLYMOD

## country-specific R0 scaling
data(polymod)

contact_r0 <- list()
for (country in survey_countries(polymod))
{
  m <- contact_matrix(survey=polymod, countries=country, n=100, symmetric=TRUE,
                      estimated.contact.age = "sample", sample.all.age.groups = TRUE,
                      age.limits=seq(0, 65, by=5))
  contact_r0[[country]] <- vapply(m$matrices, function(x) {
    as.numeric(eigen(x$matrix, only.values = TRUE)$values[1])
  }, 0)
}

contact_r0_means <- vapply(contact_r0, mean, 0)
ccr0 <- data.frame(country=names(contact_r0_means), contact_r0=contact_r0_means) %>%
  remove_rownames %>%
  mutate(rel_r0 = contact_r0 / mean(contact_r0))

cr0 <- tibble(country=factor(names(contact_countries)), polymod=contact_countries) %>%
    remove_rownames %>%
    left_join(ccr0 %>%
              dplyr::rename(polymod=country) %>%
              mutate(polymod=as.character(polymod)), by="polymod") %>%
    mutate(polymod=factor(polymod),
           country=countrycode(country, "country.name", "eurostat.name"))

model_adjImm <- sero_adjImm %>%
    gather(variable, immunity, ends_with("_immunity")) %>%
    separate(variable, c("model", "dump")) %>%
    mutate(model=recode_factor(model, adjusted="contact-adjusted", mean="plain"),
           country=countrycode(country, "country.name", "eurostat.name")) %>%
    select(-dump) %>%
    left_join(cr0, by=c("country")) %>%
    rename(scaled=rel_r0) %>%
    mutate(fixed=1) %>%
    gather(r0_model, rel_r0, fixed, scaled) %>%
    mutate(susceptibility=(1-immunity)*rel_r0) %>%
    select(-contact_r0, -rel_r0)

survey_years <- model_adjImm %>%
    group_by(country) %>%
    summarise(survey_yr=as.integer(max(year, na.rm=TRUE)))

ms_cases <- ms_world %>%
    tbl_df() %>%
    inner_join(survey_years, by=("country"))

horizon <- 10

horizon_col <- paste("years", horizon, sep=".")
case_col <- paste("cases", horizon, sep=".")
ms_cases %<>%
  replace_na(list(cases=0)) %>%
  mutate(!!horizon_col := if_else(year > survey_yr & year <= survey_yr + horizon, 1, 0),
         !!case_col := !!sym(horizon_col) * cases)

msw <- ms_cases %>%
    select(-year) %>%
    rename(year=survey_yr) %>%
    select(-starts_with("years.")) %>%
    gather(variable, cases, starts_with("cases.")) %>%
    separate(variable, c("dummy", "window"), sep="\\.") %>%
    mutate(window=as.integer(window)) %>%
    select(-dummy) %>%
    group_by(country, window, year) %>%
    summarise(sum.cases=sum(cases),
              max.cases=max(cases)) %>%
    rowwise %>%
    mutate(pop_df=list(wpp_age(country, year)),
           population=sum(pop_df$population)) %>%
    select(-pop_df) %>%
    ungroup %>%
    mutate(cases.per.million.per.year = signif(sum.cases / (population*window)*1e6),
           max.cases.per.million = signif(max.cases / (population)*1e6))

threshold_test <- seq(0.8, 1, by=0.01)

measure <- "cases.per.million.per.year"
measure_threshold <- 5
method <- "spearman"

corr_coeff <- function(df)
{
    corr <- cor.test(1 - df$susceptibility, log(df[[measure]]), method=method)
    correlation <- data.frame(corr[c("estimate", "p.value")])

    outbreaks <- df[[measure]] > measure_threshold
    immunity <- 1 - df$susceptibility
    correct_threshold <-
        lapply(threshold_test, function(t) {sum(outbreaks == (immunity < t))})
    correct <-
        data.frame(threshold=threshold_test, correct=unlist(correct_threshold))

    n <- nrow(df)
    return(tibble(correlation=list(correlation), correct=list(correct), n))
}

cc <- model_adjImm %>%
    left_join(msw, by=c("country", "year")) %>%
    group_by(sample, eqi, r0_model, model, vaccination, window) %>%
    nest %>%
    mutate(corr=map(data, corr_coeff)) %>%
    select(-data) %>%
    unnest(corr)

scorr <- cc %>%
    unnest(correlation) %>%
    group_by(eqi, r0_model, model, vaccination, window) %>%
    summarise(p.value=mean(p.value),
              median=median(estimate),
              min.estimate=quantile(estimate, 0.05),
              max.estimate=quantile(estimate, 0.95),
              n=n()) %>%
    ungroup() %>%
    arrange(median) %>%
    mutate_at(vars(5:ncol(.)), signif, digits=2)

## correlation summary table
opts_current$set(label = "correlations")
scorr %>%
    mutate(mean_ci=paste0(median, "~(", min.estimate, ", ", max.estimate, ")")) %>%
    arrange(model, r0_model, eqi) %>%
    select(`Immunity model`=model, `$R_0$ model`=r0_model, `Vaccination model`=vaccination, `Equivocal samples`=eqi,
           `Correlation~(90\\% CI)`=mean_ci) %>%
    kable("latex",
          escape=FALSE, linesep="", booktabs=TRUE,
          caption="Spearman's rank correlation between immunity estimated from nationwide serology and~(if contact-adjusted) contact studies on the one hand and the mean number of cases in the 10 years following the studies on the other. The model with the greatest absolute correlation is highlighted in bold.") %>%
    kable_styling(latex_options=c("hold_position")) %>%
    write("table_2.tex")

## misclassification error
smce <- cc %>%
    unnest(correct) %>%
    group_by(eqi, r0_model, model, vaccination, window, threshold) %>%
    dplyr::summarise(median=median(1-correct/n),
                     mean=mean(1-correct/n),
                     min=mean-sd(1-correct/n),
                     max=mean+sd(1-correct/n)) %>%
    ungroup

## figures for paper
p <- ggplot(msw %>%
            arrange(cases.per.million.per.year) %>%
            mutate(country=factor(country, levels=unique(country))),
            aes(x=country, y=100*cases.per.million.per.year)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=scales::comma(max.cases)), position=position_dodge(width=0.9), vjust=-0.3) +
    expand_limits(y=50000) +
    theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_x_discrete("") +
    geom_vline(xintercept=9.5, linetype="dotted", lwd=0.7) +
    geom_hline(yintercept=500, linetype="dotted", lwd=0.7) +
    scale_y_log10("Cases per year per million",
                  breaks=10^(0:4), labels=10^((-2):2))
save_plot("figure_1.pdf", p, base_aspect_ratio = 1.9)

smce_select <- smce %>%
    filter(eqi == "positive", r0_model=="fixed", vaccination=="projected", threshold < 1)

p <- ggplot(smce_select, aes(x=threshold, y=mean, ymin=min, ymax=max)) +
    facet_grid(~model) +
    geom_point()  +
    geom_line() +
    geom_ribbon(alpha=0.25) +
    scale_x_continuous("Threshold level", label=scales::percent) +
    scale_y_continuous("Misclassification error")

save_plot("figure_2.pdf", p, base_aspect_ratio = 1.9, base_height = 3)

number_formatter <- function(x)
{
    return(as.character(if_else(x >= 10, round(x), if_else(x >= 1, round(x, 1), signif(x, 2)))))
}

esen2_adjImm <- model_adjImm %>%
    filter(eqi == "positive", vaccination=="projected") %>%
    group_by(country, model) %>%
    summarise(immunity=mean(immunity)) %>%
    mutate(immunity=signif(immunity, 2)) %>%
    spread(model, immunity)

## pretty format summary table
opts_current$set(label = "outbreaks")
msw %>%
    left_join(esen2_adjImm) %>%
    arrange(max.cases.per.million) %>%
    dplyr::select(country, sum.cases, max.cases, cases.per.million.per.year, max.cases.per.million,
           `contact-adjusted`, plain) %>%
    mutate(country=as.character(country),
           sum.cases=as.integer(sum.cases),
           max.cases=as.integer(max.cases),
           cases.per.million.per.year=formattable(cases.per.million.per.year, formatter=number_formatter),
           max.cases.per.million=formattable(max.cases.per.million, formatter=number_formatter)) %>%
    mutate_all(linebreak) %>%
    kable("latex",
          col.names=linebreak(c("Country",
                                "Total\n(10 years)",
                                "Maximum\nannual",
                                "Mean annual\n(per million)",
                                "Maximum\nannual\n(per million)",
                                "contact-\nadjusted",
                                "plain"),
                              align="c"),
          escape=FALSE, linesep="", booktabs=TRUE,
          caption="Measles cases in the 10 years following the ESEN2 serological study, and mean estimated population immunity~(contact-adjusted or not, with fixed $R_0$ and equivocal samples interpreted as positive) based on the study and adjusted for vaccination uptake.") %>%
    add_header_above(c(" ", Cases=4, Immunity=2)) %>%
    kable_styling(latex_options=c("hold_position")) %>%
    ## kable_as_image("esen2") %>%
    write("table_1.tex")
