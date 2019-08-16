if (!exists("nsamples") || !(is.numeric(nsamples))) stop("nsamples must be set")

age_limits <- seq(0, 20, by = 5)
targets <- list(current = c(0.85, 0.9, 0.95, 0.95, 0.95),
                increase_infant = c(0.9, 0.9, 0.95, 0.95, 0.95),
                fewer_infants_catchup_kids = c(0.8, 0.9, 0.95, 0.95, 0.95),
                catchup_kids = c(0.85, 0.95, 0.95, 0.95, 0.95),
                increase_infant_and_catchup_kids = c(0.9, 0.95, 0.95, 0.95, 0.95),
                all_95 = c(0.95, 0.95, 0.95, 0.95, 0.95),
                current_less_teenagers = c(0.85, 0.9, 0.9, 0.95, 0.95),
                catchup_kids_less_teenagers = c(0.85, 0.95, 0.9, 0.95, 0.95),
                current_even_less_teenagers = c(0.85, 0.9, 0.85, 0.95, 0.95),
                catchup_kids_even_less_teenagers = c(0.85, 0.95, 0.85, 0.95, 0.95),
                current_less_adolescents = c(0.85, 0.9, 0.95, 0.90, 0.95),
                catchup_kids_less_adolescents = c(0.85, 0.95, 0.95, 0.90, 0.95),
                current_even_less_adolescents = c(0.85, 0.9, 0.95, 0.85, 0.95),
                catchup_kids_even_less_adolescents = c(0.85, 0.95, 0.95, 0.85, 0.95),
                current_less_older = c(0.85, 0.9, 0.95, 0.95, 0.90),
                catchup_kids_less_older = c(0.85, 0.95, 0.95, 0.95, 0.90),
                current_even_less_older = c(0.85, 0.9, 0.95, 0.95, 0.85),
                catchup_kids_even_less_older = c(0.85, 0.95, 0.95, 0.95, 0.85))
targets <- lapply(targets, function(x) { names(x) = age_limits; x })

data(polymod)

## Uganda and SMILI contact data available upon request from the original authors
zimbabwe <- get_survey("https://doi.org/10.5281/zenodo.1127693")
peru <- get_survey("https://doi.org/10.5281/zenodo.1095664")
france <- get_survey("https://doi.org/10.5281/zenodo.1157918")
hongkong <- get_survey("https://doi.org/10.5281/zenodo.1165561")

surveys <- list(polymod=polymod,
                zimbabwe=zimbabwe,
                peru=peru,
                france=france,
                hongkong=hongkong)

datasets <-
  tibble(survey=names(surveys)) %>%
  rowwise %>%
  mutate(settings=
             list(c("all", unique(surveys[[survey]][["participants"]][["setting"]]))),
         country=list(c(survey_countries(surveys[[survey]], quiet=TRUE)))) %>%
  unnest(settings, .preserve = c(survey, country)) %>%
  unnest(country, .preserve=survey)

base_parameters <-
    list(age.limits=age_limits, symmetric=TRUE, weigh.dayofweek=TRUE,
         sample.all.age.groups=TRUE, estimated.contact.age="sample",
         missing.participant.age="remove", missing.contact.age="sample",
         quiet=TRUE, n=nsamples)

R <- tibble(scenario=names(targets)) %>%
    rowwise() %>%
    mutate(immunity=list(targets[[scenario]])) %>%
    ungroup() %>%
    crossing(datasets) %>%
    arrange(survey, settings, country, scenario)

matrix_options <- R %>%
    group_by(settings, survey, country) %>%
    summarise() %>%
    rowwise() %>%
    mutate(country_parameters = list(list(countries=country)),
           filter_parameters =
               if_else(settings=="all",
                       list(list()),
                       list(list(filter=list(setting=settings)))),
           options=list(c(base_parameters, country_parameters, filter_parameters,
                          list(survey=surveys[[survey]]))))

country_matrices <- matrix_options %>%
    rowwise %>%
    mutate(matrix_ret=list(do.call(contact_matrix, options)),
           matrix=list(matrix_ret$matrices),
           demography=list(matrix_ret$demography$proportion)) %>%
    ungroup

scenarios_adjImm <- R %>%
    left_join(country_matrices, by = c("settings", "survey", "country")) %>%
    unnest(matrix, .preserve=c(demography, immunity)) %>%
    rowwise %>%
    mutate(adjusted_immunity=
               adjust_immunity(mixing=matrix, immunity=unlist(immunity)),
           mean_immunity= sum(demography * immunity) / sum(demography)) %>%
    ungroup %>%
    select(-survey, -matrix, -demography)
