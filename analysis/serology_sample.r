nsamples <- 1000

data(mcv_world_estimate) ## measles coverage estimates
data(esen2_schedules) ## vaccination schedules
data(contact_countries) ## ESEN2-to-POLYMOD country mapping

wms <- ms.sero %>%
  tbl_df %>%
  mutate(survey.yr = as.integer(sub("[^0-9].*$", "", survey.yr)),
         country = if_else(country == "Czech", "Czech Republic", as.character(country)),
         country = countrycode(country, "country.name", "eurostat.name"),
         country = factor(country),
         lower.age.limit = as.integer(sub("[-+ ].*$", "", age1)))

## work out immunity levels by country and age group
wms_sero <- list()

base_parameters <-
    list(symmetric=TRUE, weigh.dayofweek=TRUE, sample.all.age.groups=TRUE,
         estimated.contact.age="sample", missing.participant.age="remove", quiet=TRUE)

survey_parameters <- list(survey=polymod)

R_countries <- list()

eqi_scenarios <- c("positive", "negative", "removed")

sample_sero <- function(df, nsamples = 1)
{
  lapply(1:nsamples, function (x)
  {
    temp <- df[sample(1:nrow(df), replace=TRUE), ]
    temp$sample <- x
    temp
  }) %>%
    bind_rows %>%
    group_by(stdres, sample, lower.age.limit) %>%
    summarise(n = n()) %>%
    ungroup %>%
    complete(stdres, sample, lower.age.limit, fill=list(n=0)) %>%
    spread(stdres, n) %>%
    mutate(positive = (EQI + POS) / (EQI + POS + NEG),
           negative = POS / (EQI + POS + NEG),
           removed = POS / (POS + NEG)) %>%
    select(-EQI, -NEG, -POS) %>%
    gather(eqi, immunity, -sample, -lower.age.limit) %>%
    complete(sample, lower.age.limit, eqi, fill=list(immunity=1)) %>%
    spread(lower.age.limit, immunity) %>%
    nest(-sample, -eqi) %>%
    rowwise() %>%
    mutate(immunity=list(unlist(data[1,]))) %>%
    select(-data)
}

maternal.immunity <- wms %>%
  filter(lower.age.limit==0) %>%
  nest() %>%
  mutate(sero_sample=map(data, sample_sero, nsamples = nsamples)) %>%
  unnest(sero_sample) %>%
  mutate(immunity=unlist(immunity))

project_country <- function(country, year, eqi, immunity, schedule, vaccination=c("projected", "ignored"), years=10)
{
  vaccination <- match.arg(vaccination)
  vacc_country <- country
  sero_year <- year
  sero_eqi <- eqi
  schedule <- schedule[!is.na(schedule)]

  uptake <- mcv_world_estimate %>%
    filter(country==vacc_country) %>%
    mutate(uptake=uptake/100) %>%
    spread(vaccine, uptake) %>%
    replace_na(list(MCV1=0, MCV2=0))

  if (country=="Lithuania") {
    uptake %<>% mutate(MCV3=MCV2)
  }

  uptake %<>%
    column_to_rownames("year") %>%
    select(starts_with("MCV")) %>%
    as.matrix %>%
    t

  ## fill missing dose information
  while (nrow(uptake) < length(schedule))
  {
    uptake <- rbind(uptake, rep(0, ncol(uptake)))
  }

  if (vaccination=="ignored") years <- 1

  immunities <- list()
  mat.imm <- maternal.immunity %>%
    filter(sample==sample(seq_len(max(sample)), 1),
           eqi==sero_eqi) %>%
    .$immunity

  for (i in seq_len(years)) {
    projected <-
      project_immunity(immunity, year, year + i, uptake, schedule, mat.imm, 0.95)
    immunities[[i]] <-
      tibble(lower.age.limit=as.integer(names(projected)), immunity=projected,
                              year=year+i)
  }

  immunities <- immunities %>%
    bind_rows %>%
    group_by(lower.age.limit) %>%
    summarise(immunity=mean(immunity))

  ret <- immunities$immunity
  names(ret) <- immunities$lower.age.limit

  return(ret)
}

## create samples form serology
swms <- wms %>%
  filter(!is.na(survey.yr)) %>%
  unite(country_yr, country, survey.yr, by="country") %>%
  nest(-country_yr) %>%
  mutate(sero_sample=map(data, sample_sero, nsamples = nsamples)) %>%
  unnest(sero_sample) %>%
  separate(country_yr, c("country", "year"), sep="_") %>%
  mutate(year=as.integer(year)) %>%
  left_join(esen2_schedules, by="country") %>%
  crossing(vaccination=c("projected", "ignored")) %>%
  rowwise() %>%
  mutate(immunity=list(project_country(country, year, eqi, immunity, schedule, vaccination))) %>%
  ungroup()

country_matrices <- swms %>%
  group_by(country, vaccination) %>%
  summarise(n=n(), age.groups=list(as.integer(unique(names(unlist(immunity)))))) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(options=
           list(c(base_parameters, survey_parameters,
                  list(n=n,
                       age.limits=as.integer(grep("^[0-9]", age.groups,
                                                  value=TRUE)),
                       countries=contact_countries[country]))),
         matrix=list(do.call(contact_matrix, options)),
         matrices=list(matrix$matrices),
         demography=list(matrix$demography$proportion)) %>%
  unnest(matrices, .preserve=demography) %>%
  rename(matrix=matrices) %>%
  group_by(country, vaccination) %>%
  mutate(id=1:n()) %>%
  select(-n)

sero_adjImm <- swms %>%
  group_by(country, vaccination) %>%
  mutate(id=1:n()) %>%
  ungroup %>%
  left_join(country_matrices, by=c("country", "vaccination", "id")) %>%
  select(-id) %>%
  rowwise() %>%
  mutate(adjusted_immunity=adjust_immunity(mixing=matrix, immunity=unlist(immunity)),
         mean_immunity=sum(demography * immunity)) %>%
  ungroup %>%
  select(country, year, sample, eqi, vaccination,
         adjusted_immunity, mean_immunity)
