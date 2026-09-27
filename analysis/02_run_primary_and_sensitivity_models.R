suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(survey)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "analysis/02_run_primary_and_sensitivity_models.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
data_dir <- Sys.getenv("CHARLS_DATA_DIR", unset = file.path(repo_root, "data"))
out_dir <- file.path(repo_root, "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
norm_id <- function(x) {
  x <- trimws(as.character(x))
  ifelse(nchar(x) == 11, paste0(substr(x, 1, 9), "0", substr(x, 10, 11)), x)
}

waves <- c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018)
intervals <- list(c(1, 2, "2011-2013"), c(2, 3, "2013-2015"), c(3, 4, "2015-2018"))
diseases <- c("hibpe", "diabe", "cancre", "lunge", "hearte", "stroke", "psyche",
              "arthre", "dyslipe", "livere", "kidneye", "asthmae", "memrye")

profile_specs <- list(
  `1` = list("w1_health.dta", c("da010_10_s1", "da010_10_s2", "da010_10_s3", "da010_10_s4")),
  `2` = list("w2_health.dta", c("da010_10_s1", "da010_10_s2", "da010_10_s3", "da010_10_s4")),
  `3` = list("w3_health.dta", c("da010_10_s1", "da010_10_s2", "da010_10_s3", "da010_10_s4")),
  `4` = list("w4_health.dta", c("da010_w4_10__s1", "da010_w4_10__s2", "da010_w4_10__s3", "da010_w4_10__s4"))
)

build_profile <- function(w) {
  spec <- profile_specs[[as.character(w)]]
  d <- read_dta(file.path(data_dir, spec[[1]]), col_select = c("ID", all_of(spec[[2]])))
  d$ID <- norm_id(d$ID)
  s <- lapply(seq_along(spec[[2]]), function(i) {
    out <- num(d[[spec[[2]][i]]]) == i
    out[is.na(out)] <- FALSE
    out
  })
  observed <- rowSums(!is.na(as.data.frame(d[spec[[2]]]))) > 0
  contradiction <- s[[4]] & (s[[1]] | s[[2]] | s[[3]])
  profile <- rep(NA_character_, nrow(d))
  profile[observed & !contradiction & s[[4]]] <- "None listed"
  profile[observed & !contradiction & s[[1]] & !s[[2]]] <- "TCM without Western"
  profile[observed & !contradiction & s[[2]] & !s[[1]]] <- "Western without TCM"
  profile[observed & !contradiction & s[[1]] & s[[2]]] <- "Combined TCM and Western"
  profile[observed & !contradiction & s[[3]] & !s[[1]] & !s[[2]]] <- "Other only"
  tibble(ID = d$ID, profile = profile) |> distinct(ID, .keep_all = TRUE)
}

needed <- c("ID", "communityID", "ragender", "raeduc_c")
for (w in 1:4) {
  needed <- c(needed, paste0("r", w, c("iwstat", "agey", "digeste", "digestf", "rural2",
                                       "mstat", "higov", "adla_c", "wtrespb", "rxdigest", "rxdigest_c")),
              paste0("h", w, "rural"), paste0("r", w, diseases))
}
h <- read_dta(file.path(data_dir, "harmonized.dta"), col_select = all_of(needed))
h$ID <- norm_id(h$ID)

pw <- bind_rows(lapply(1:4, function(w) {
  p <- build_profile(w)
  z <- h |> select(ID, communityID, ragender, raeduc_c,
                   all_of(paste0("r", w, c("iwstat", "agey", "digeste", "digestf", "rural2", "mstat", "higov", "adla_c", "wtrespb"))),
                   all_of(paste0("h", w, "rural")), all_of(paste0("r", w, diseases))) |>
    rename(iwstat = all_of(paste0("r", w, "iwstat")), age = all_of(paste0("r", w, "agey")),
           digest = all_of(paste0("r", w, "digeste")), dispute = all_of(paste0("r", w, "digestf")),
           rural_hukou = all_of(paste0("r", w, "rural2")), rural_residence = all_of(paste0("h", w, "rural")),
           mstat = all_of(paste0("r", w, "mstat")), public_insurance = all_of(paste0("r", w, "higov")),
           adl_count = all_of(paste0("r", w, "adla_c")), weight = all_of(paste0("r", w, "wtrespb"))) |>
    left_join(p, by = "ID")
  dv <- as.data.frame(lapply(paste0("r", w, diseases), function(nm) num(z[[nm]])))
  z$comorbidity_count <- ifelse(rowSums(!is.na(dv)) >= 10, rowSums(dv == 1, na.rm = TRUE), NA_real_)
  edu <- num(z$raeduc_c)
  z |> transmute(
    ID, communityID = as.character(communityID), wave = w, year = waves[as.character(w)],
    iwstat = num(iwstat), age = num(age), digest = num(digest), dispute = num(dispute),
    rural_hukou = num(rural_hukou), rural_residence = num(rural_residence), weight = num(weight), profile,
    female = ifelse(is.na(num(ragender)), NA_real_, as.numeric(num(ragender) == 2)),
    education = factor(case_when(edu %in% 1:5 ~ "Lower", edu %in% 6:7 ~ "Upper secondary/vocational",
                                 edu %in% 8:10 ~ "Tertiary", TRUE ~ NA_character_),
                       levels = c("Lower", "Upper secondary/vocational", "Tertiary")),
    partnered = ifelse(is.na(num(mstat)), NA_real_, as.numeric(num(mstat) %in% c(1, 3))),
    public_insurance = ifelse(is.na(num(public_insurance)), NA_real_, as.numeric(num(public_insurance) == 1)),
    adl_difficulty = ifelse(is.na(num(adl_count)), NA_real_, as.numeric(num(adl_count) > 0)),
    comorbidity_count,
    joint_group = factor(case_when(
      rural_hukou == 0 & rural_residence == 0 ~ "Urban hukou/urban residence",
      rural_hukou == 1 & rural_residence == 0 ~ "Rural hukou/urban residence",
      rural_hukou == 0 & rural_residence == 1 ~ "Urban hukou/rural residence",
      rural_hukou == 1 & rural_residence == 1 ~ "Rural hukou/rural residence",
      TRUE ~ NA_character_),
      levels = c("Urban hukou/urban residence", "Rural hukou/urban residence",
                 "Urban hukou/rural residence", "Rural hukou/rural residence")),
    none_listed = profile == "None listed", any_listed = !is.na(profile) & profile != "None listed"
  )
}))

iv <- bind_rows(lapply(intervals, function(it) {
  w0 <- as.numeric(it[1]); w1 <- as.numeric(it[2]); lab <- it[3]
  a <- pw |> filter(wave == w0) |> rename_with(~paste0(.x, "_base"), -ID)
  b <- pw |> filter(wave == w1) |> select(ID, digest, dispute, profile, none_listed, any_listed, iwstat) |>
    rename_with(~paste0(.x, "_follow"), -ID)
  a |> left_join(b, by = "ID") |> mutate(
    interval = lab,
    eligible_base = digest_base == 1 & !is.na(profile_base) & age_base >= 45 & weight_base > 0,
    observed_follow = digest_follow == 1 & !is.na(profile_follow)
  )
}))

analytic <- iv |> filter(eligible_base, observed_follow) |>
  group_by(interval) |> mutate(weight_norm = weight_base / mean(weight_base, na.rm = TRUE)) |> ungroup() |>
  mutate(to_listed = as.numeric(none_listed_base & any_listed_follow),
         to_none = as.numeric(any_listed_base & none_listed_follow))


# Number of distinct respondents contributing one, two, or three observed intervals.
contrib <- analytic |> count(ID, name = "intervals_contributed") |> count(intervals_contributed, name = "respondents")
write.csv(contrib, file.path(out_dir, "respondent_interval_contributions.csv"), row.names = FALSE)

# Weight distribution and covariate completeness.
weights <- analytic |> group_by(interval) |> summarise(
  n = n(), min = min(weight_norm), p1 = quantile(weight_norm, .01), p50 = median(weight_norm),
  mean = mean(weight_norm), p99 = quantile(weight_norm, .99), max = max(weight_norm), .groups = "drop")
write.csv(weights, file.path(out_dir, "normalized_weight_distribution.csv"), row.names = FALSE)

covars <- c("age_base", "female_base", "education_base", "partnered_base", "public_insurance_base",
            "adl_difficulty_base", "comorbidity_count_base", "communityID_base")
miss <- analytic |> group_by(joint_group_base) |> summarise(n = n(),
  across(all_of(covars), ~sum(is.na(.x)), .names = "missing_{.col}"), .groups = "drop")
write.csv(miss, file.path(out_dir, "baseline_covariate_missingness.csv"), row.names = FALSE)

wmean <- function(x, w) weighted.mean(x, w, na.rm = TRUE)
baseline <- analytic |> group_by(joint_group_base) |> summarise(
  records = n(), respondents = n_distinct(ID),
  age_mean = wmean(age_base, weight_norm), female_pct = 100*wmean(female_base, weight_norm),
  lower_education_pct = 100*wmean(as.numeric(education_base == "Lower"), weight_norm),
  partnered_pct = 100*wmean(partnered_base, weight_norm),
  public_insurance_pct = 100*wmean(public_insurance_base, weight_norm),
  adl_difficulty_pct = 100*wmean(adl_difficulty_base, weight_norm),
  comorbidity_mean = wmean(comorbidity_count_base, weight_norm),
  baseline_none_pct = 100*wmean(as.numeric(none_listed_base), weight_norm),
  int_2011_2013_pct = 100*wmean(as.numeric(interval == "2011-2013"), weight_norm),
  int_2013_2015_pct = 100*wmean(as.numeric(interval == "2013-2015"), weight_norm),
  int_2015_2018_pct = 100*wmean(as.numeric(interval == "2015-2018"), weight_norm), .groups = "drop")
write.csv(baseline, file.path(out_dir, "baseline_by_hukou_residence.csv"), row.names = FALSE)

# Weighted standardized differences versus the urban-hukou/urban-residence reference.
smd_data <- analytic |> filter(!is.na(joint_group_base)) |> mutate(
  lower_education = as.numeric(education_base == "Lower"),
  baseline_none = as.numeric(none_listed_base),
  int_2011_2013 = as.numeric(interval == "2011-2013"),
  int_2013_2015 = as.numeric(interval == "2013-2015"),
  int_2015_2018 = as.numeric(interval == "2015-2018"))
smd_vars <- c("age_base", "female_base", "lower_education", "partnered_base", "public_insurance_base",
              "adl_difficulty_base", "comorbidity_count_base", "baseline_none",
              "int_2011_2013", "int_2013_2015", "int_2015_2018")
smd_rows <- bind_rows(lapply(smd_vars, function(v) {
  stats <- smd_data |> group_by(joint_group_base) |> summarise(
    mean = weighted.mean(.data[[v]], weight_norm, na.rm = TRUE),
    variance = weighted.mean((.data[[v]] - weighted.mean(.data[[v]], weight_norm, na.rm = TRUE))^2,
                             weight_norm, na.rm = TRUE), .groups = "drop")
  ref <- stats |> filter(joint_group_base == "Urban hukou/urban residence")
  stats |> mutate(variable = v,
                  smd_vs_urban_urban = (mean - ref$mean[1]) / sqrt((variance + ref$variance[1])/2))
}))
write.csv(smd_rows, file.path(out_dir, "baseline_standardized_differences.csv"), row.names = FALSE)

# Unadjusted group-specific transition probabilities with respondent-cluster robust CIs.
unadjusted_group <- function(dat, outcome, label) {
  dat <- dat |> filter(!is.na(joint_group_base), !is.na(.data[[outcome]]), !is.na(weight_norm))
  des <- svydesign(ids = ~ID, weights = ~weight_norm, data = dat, nest = TRUE)
  est <- svyby(as.formula(paste0("~", outcome)), ~joint_group_base, des, svymean,
               vartype = c("se", "ci"), na.rm = TRUE, level = 0.95, keep.names = FALSE)
  counts <- dat |> group_by(joint_group_base) |> summarise(
    unweighted_n = n(), events = sum(.data[[outcome]] == 1), .groups = "drop")
  names(est)[2:5] <- c("weighted_probability", "se", "ci_l", "ci_u")
  as_tibble(est) |> left_join(counts, by = "joint_group_base") |>
    transmute(outcome = label, joint_group = joint_group_base, unweighted_n, events,
              weighted_probability, ci_l, ci_u)
}
unadjusted <- bind_rows(
  unadjusted_group(filter(analytic, none_listed_base), "to_listed", "None-listed to listed"),
  unadjusted_group(filter(analytic, any_listed_base), "to_none", "Listed to none-listed")
)
write.csv(unadjusted, file.path(out_dir, "unadjusted_transitions_by_hukou_residence.csv"), row.names = FALSE)

formula_a <- to_listed ~ joint_group_base + factor(interval) + age_base + female_base + education_base +
  partnered_base + public_insurance_base + adl_difficulty_base + comorbidity_count_base
formula_b <- to_none ~ joint_group_base + factor(interval) + age_base + female_base + education_base +
  partnered_base + public_insurance_base + adl_difficulty_base + comorbidity_count_base

fit_extract <- function(dat, formula, cluster, label) {
  keep <- complete.cases(model.frame(formula, data = dat, na.action = na.pass)) & !is.na(dat[[cluster]])
  dat <- dat[keep, , drop = FALSE]
  des <- svydesign(ids = as.formula(paste0("~", cluster)), weights = ~weight_norm, data = dat, nest = TRUE)
  fit <- svyglm(formula, design = des, family = quasipoisson(link = "log"))
  co <- summary(fit)$coefficients
  tibble(analysis = label, cluster = cluster, n = nrow(dat), events = sum(model.response(model.frame(fit))),
         term = rownames(co), estimate = co[,1], se = co[,2], p = co[,4],
         rr = exp(co[,1]), lcl = exp(co[,1] - 1.96*co[,2]), ucl = exp(co[,1] + 1.96*co[,2]),
         converged = isTRUE(fit$converged))
}

standardized_primary <- function(dat, formula, outcome_label) {
  keep <- complete.cases(model.frame(formula, data = dat, na.action = na.pass)) & !is.na(dat$ID)
  dat <- droplevels(dat[keep, , drop = FALSE])
  des <- svydesign(ids = ~ID, weights = ~weight_norm, data = dat, nest = TRUE)
  fit <- svyglm(formula, design = des, family = quasipoisson(link = "log"))
  beta <- coef(fit); vv <- vcov(fit); tt <- delete.response(terms(fit))
  groups <- levels(dat$joint_group_base)
  calc <- function(g) {
    nd <- dat; nd$joint_group_base <- factor(g, levels = groups)
    mm <- model.matrix(tt, nd, contrasts.arg = fit$contrasts)
    mu <- as.numeric(exp(mm %*% beta))
    grad <- colMeans(mu * mm)
    p <- mean(mu); se <- sqrt(as.numeric(t(grad) %*% vv %*% grad))
    list(p = p, se = se, grad = grad)
  }
  vals <- lapply(groups, calc); ref <- vals[[1]]
  bind_rows(lapply(seq_along(groups), function(i) {
    x <- vals[[i]]
    if (i == 1) return(tibble(outcome = outcome_label, joint_group = groups[i], n = nrow(dat),
                              events = sum(model.response(model.frame(formula, dat))),
                              probability = x$p, probability_lcl = max(0, x$p-1.96*x$se), probability_ucl = min(1, x$p+1.96*x$se),
                              risk_difference = 0, rd_lcl = NA_real_, rd_ucl = NA_real_, risk_ratio = 1, rr_lcl = NA_real_, rr_ucl = NA_real_))
    gd <- x$grad - ref$grad; sed <- sqrt(as.numeric(t(gd) %*% vv %*% gd)); rd <- x$p-ref$p
    gl <- x$grad/x$p - ref$grad/ref$p; sel <- sqrt(as.numeric(t(gl) %*% vv %*% gl)); rr <- x$p/ref$p
    tibble(outcome = outcome_label, joint_group = groups[i], n = nrow(dat),
           events = sum(model.response(model.frame(formula, dat))),
           probability = x$p, probability_lcl = max(0, x$p-1.96*x$se), probability_ucl = min(1, x$p+1.96*x$se),
           risk_difference = rd, rd_lcl = rd-1.96*sed, rd_ucl = rd+1.96*sed,
           risk_ratio = rr, rr_lcl = exp(log(rr)-1.96*sel), rr_ucl = exp(log(rr)+1.96*sel))
  }))
}

model_rows <- bind_rows(
  fit_extract(filter(analytic, none_listed_base), formula_a, "ID", "None-listed to listed, all intervals"),
  fit_extract(filter(analytic, any_listed_base), formula_b, "ID", "Listed to none-listed, all intervals"),
  fit_extract(filter(analytic, none_listed_base), formula_a, "communityID_base", "None-listed to listed, community-cluster"),
  fit_extract(filter(analytic, any_listed_base), formula_b, "communityID_base", "Listed to none-listed, community-cluster"),
  fit_extract(filter(analytic, none_listed_base, interval != "2015-2018"), formula_a, "ID", "None-listed to listed, 2011-2015 only"),
  fit_extract(filter(analytic, any_listed_base, interval != "2015-2018"), formula_b, "ID", "Listed to none-listed, 2011-2015 only")
)
write.csv(model_rows, file.path(out_dir, "model_coefficients_and_sensitivity.csv"), row.names = FALSE)

primary_standardized <- bind_rows(
  standardized_primary(filter(analytic, none_listed_base), formula_a, "None-listed to listed"),
  standardized_primary(filter(analytic, any_listed_base), formula_b, "Listed to none-listed")
)
write.csv(primary_standardized, file.path(out_dir, "final_primary_standardized.csv"), row.names = FALSE)

standardized_intervals <- function(dat, outcome, outcome_label) {
  dat$interval <- factor(dat$interval, levels = c("2011-2013", "2013-2015", "2015-2018"))
  form <- as.formula(paste0(outcome, " ~ joint_group_base * interval + age_base + female_base + education_base + partnered_base + public_insurance_base + adl_difficulty_base + comorbidity_count_base"))
  keep <- complete.cases(model.frame(form, data = dat, na.action = na.pass)) & !is.na(dat$ID)
  dat <- droplevels(dat[keep, , drop = FALSE]); des <- svydesign(ids=~ID, weights=~weight_norm, data=dat, nest=TRUE)
  fit <- svyglm(form, design=des, family=quasipoisson(link="log")); beta<-coef(fit); vv<-vcov(fit); tt<-delete.response(terms(fit))
  groups <- levels(dat$joint_group_base)
  bind_rows(lapply(levels(dat$interval), function(intv) {
    eval <- dat[dat$interval == intv, , drop=FALSE]
    calc <- function(g) { nd<-eval; nd$joint_group_base<-factor(g,levels=groups); mm<-model.matrix(tt,nd,contrasts.arg=fit$contrasts); mu<-as.numeric(exp(mm%*%beta)); gr<-colMeans(mu*mm); list(p=mean(mu),gr=gr,se=sqrt(as.numeric(t(gr)%*%vv%*%gr))) }
    vals<-lapply(groups,calc); ref<-vals[[1]]
    bind_rows(lapply(seq_along(groups), function(i){x<-vals[[i]]; if(i==1) return(tibble(outcome=outcome_label,interval=intv,joint_group=groups[i],probability=x$p,probability_lcl=max(0,x$p-1.96*x$se),probability_ucl=min(1,x$p+1.96*x$se),risk_difference=0,rd_lcl=NA_real_,rd_ucl=NA_real_,risk_ratio=1,rr_lcl=NA_real_,rr_ucl=NA_real_)); gd<-x$gr-ref$gr; sed<-sqrt(as.numeric(t(gd)%*%vv%*%gd)); rd<-x$p-ref$p; gl<-x$gr/x$p-ref$gr/ref$p; sel<-sqrt(as.numeric(t(gl)%*%vv%*%gl)); rr<-x$p/ref$p; tibble(outcome=outcome_label,interval=intv,joint_group=groups[i],probability=x$p,probability_lcl=max(0,x$p-1.96*x$se),probability_ucl=min(1,x$p+1.96*x$se),risk_difference=rd,rd_lcl=rd-1.96*sed,rd_ucl=rd+1.96*sed,risk_ratio=rr,rr_lcl=exp(log(rr)-1.96*sel),rr_ucl=exp(log(rr)+1.96*sel)) }))
  }))
}
interval_standardized <- bind_rows(
  standardized_intervals(filter(analytic, none_listed_base), "to_listed", "None-listed to listed"),
  standardized_intervals(filter(analytic, any_listed_base), "to_none", "Listed to none-listed"))
write.csv(interval_standardized, file.path(out_dir, "final_interval_standardized.csv"), row.names=FALSE)

cat("Completed primary and sensitivity analyses.\n")
print(contrib)
print(model_rows |> filter(grepl("Rural hukou/rural residence", term)) |>
        select(analysis, cluster, n, events, rr, lcl, ucl, p, converged))
