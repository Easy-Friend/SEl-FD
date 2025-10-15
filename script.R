install.packages("renv")
renv::init()
renv::install(c("tidyverse","readxl","janitor",
                "cobalt","survey", "pscl"
))
renv::snapshot()

library(tidyverse)
library(readxl)
library(janitor)
library(cobalt)
library(survey)
library(pscl)

raw <- readxl::read_excel("C:/Users/kkw10/OneDrive/바탕 화면/강우/공부/논문/pf한현진 SURPASS DWI/data/Final_1.1.1age.xlsx") |> janitor::clean_names()

names(raw)
raw <- raw %>%
  rename(
    group          = group,                
    procedure_time = procedure_time,
    sex            = sex_m_1_f_2,
    ais_hx = ais_hx,
    sah_hx         = sah_hx,
    fhx            = f_hx,
    smoking_current       = smoking_current_0,   
    other_mhx      = other_m_hx,           
    anti_thrombotics_prev = anti_thrombotics_ijeonbuteo_bog_yong, 
    aru            = pre_op_anti_plt_verify_now_aru,
    anti_plt_resistance = anti_plt_resistance,
    tx_location = tx_location_parent_artery_location,
    stent_length   = stent_gil_i,
    stent_no       = stent,
    immediate_dwi  = immediate_dwi,        
    dwi_count      = dwi_lesion_gaesu_reader_1
  )

vars <- c("procedure_time", "sex", "age" , "bmi", "htn", "dm", "dyslipidemia", "cva_hx",
          "fhx", "smoking_current", "alcohol", "anti_thrombotics_prev", "anti_plt_resistance",
          "tx_location", "an_morphology", "multiplicity", "balloon_angioplasty", "adjuvant_coil", 
          "tirofiban", "immediate_dwi", "width", "dwi_count", "stent_no", 
          "stent_diameter", "stent_length")

cont_vars <- c("procedure_time", "dwi_count", "stent_diameter", "stent_length" )
cat_vars  <- setdiff(vars, cont_vars)

cleaned <- raw |>
  filter(if_all(vars, ~ !is.na(.)))

cleaned_inversed <- cleaned |>
  mutate(across(all_of(setdiff(cat_vars, "sex")), ~ 1-.))

#####propensity score calculation & overlap loveplot 확인 & SMD adjustment 확인
#sah_hx는 perfect separation 된 관계로 confounder에서 제외
confounders <- c("sex", "age" , "bmi", "htn", "dm", "dyslipidemia", "cva_hx",
          "smoking_current", "anti_thrombotics_prev", "anti_plt_resistance",
          "tx_location", "an_morphology")

ps_model <- glm(
  formula = as.formula (
    paste("group ~", paste(confounders, collapse = "+")),
  ),
  data = cleaned_inversed,
  family = binomial()
  )

cleaned_inversed <- cleaned_inversed |>
  mutate(
    ps     = predict(ps_model, type = "response"),
    ps_cap = pmin(pmax(ps, 0.01), 0.99)
  )

ps_overlap <- cleaned_inversed |>
  mutate(group=factor(group, levels = c(0, 1), labels = c("EVOLVE", "ELITE")))

ggplot(ps_overlap, aes(x = ps_cap, fill = group)) + 
  geom_density(alpha = 0.4, adjust = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Propensity Score Overlap",
       x = "Propensity Score (ps_cap)", y = "Density", fill = "Group") +
  theme_minimal()

p_elite <- mean(cleaned_inversed$group ==1)
cleaned_inversed_weighted <- cleaned_inversed |>
  mutate(
    sw = ifelse(group ==1, p_elite / ps_cap, (1-p_elite) /(1-ps_cap))
  )

bal <- bal.tab(
  as.formula(paste("group ~", paste(confounders, collapse = "+"))),
  data = cleaned_inversed_weighted,
  weights = cleaned_inversed_weighted$sw,
  un = TRUE
)

love.plot(
  bal,
  un = TRUE,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = .1),
  colors = c("grey60", "steelblue"),
  sample.names = c("Unadjusted","Adjusted"),
  title = "SMD"
)

#####IPTW-weighted Logistic Regression
design_iptw <- svydesign(ids = ~1, weights = ~sw, data = cleaned_inversed_weighted)
fit_lr <- svyglm(immediate_dwi ~ group, design = design_iptw, family = quasibinomial())
summary(fit_lr)

#####IPTW-weighted Quasi-Poisson analysis
fit_qp <- svyglm(dwi_count ~ group, design = design_iptw, family = quasipoisson(link = "log"))
summary(fit_qp)

#####IPTW-weighted ZINB(혹시나... 통계적 의미 잘 모르겠음)
fit_zinb2 <- zeroinfl(
  dwi_count ~ group | 1,
  dist = "negbin",
  data = cleaned_inversed_weighted,
  weights = sw,
  link = "logit"
)
summary(fit_zinb2)

#####단순 두 그룹 통계 비교
#1. Welch t-test for cont_var
mean_ci <- function(x, conf.level= 0.95){
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  se <- s/sqrt(n)
  tcrit <- qt((1 + conf.level)/2, df = n - 1)
  c(mean = m, low = m - tcrit * se, high = m + tcrit * se, n = n)
}

fmt_mean_ci <- function(m, lo, hi, digits = 2) {
  sprintf("%.*f (95%% CI %.*f–%.*f)", digits, m, digits, lo, digits, hi)
}

base_cont <- function(data, var, group = "group", conf.level = 0.95) {
  gsum <- data |>
    group_by(.data[[group]]) |> 
    summarise(stats = list(mean_ci(.data[[var]], conf.level = conf.level)),
              .groups = "drop") |>
    mutate(mean = map_dbl(stats, ~ .x["mean"]),
           lo   = map_dbl(stats, ~ .x["low"]),
           hi   = map_dbl(stats, ~ .x["high"]))
  
  g0 <- gsum |> filter(.data[[group]] == 0)
  g1 <- gsum |> filter(.data[[group]] == 1)
  
  tt <- t.test(as.formula(paste(var, "~", group)), data = data, var.equal = FALSE, conf.level = conf.level)
  
  tibble(
    Variable = var,
    Group0   = fmt_mean_ci(g0$mean, g0$lo, g0$hi),
    Group1   = fmt_mean_ci(g1$mean, g1$lo, g1$hi),
    p_value  = signif(tt$p.value, 3)
  )
}

base_welch <- map_dfr(cont_vars, ~ base_cont(cleaned_inversed, .x))

#2. chi-square for cat_vars
prop_ci <- function(event, total, conf.level = 0.95) {
  pt <- prop.test(event, total, conf.level = conf.level)
  c(p = event/total, lo = pt$conf.int[1], hi = pt$conf.int[2])
}

fmt_prop_ci <- function(p, lo, hi, digits = 2){
  sprintf("%.1f%% (95%% CI %.1f–%.1f%%)", 100*p, 100*lo, 100*hi)
}

base_cat <- function(data, var, group = "group", event_level = 1, conf.level = 0.95){
  g <- data[[group]]
  x <- data[[var]]
  levs <- sort(unique(g))
  g0 <- levs[1]; g1 <- levs[2]
  
  tot0 <- sum(!is.na(x[g == g0]))
  tot1 <- sum(!is.na(x[g == g1]))
  ev0 <- sum(x[g == g0] == event_level)
  ev1 <- sum(x[g == g1] == event_level)
  
  ci0 <- prop_ci(ev0, tot0, conf.level)
  ci1 <- prop_ci(ev1, tot1, conf.level)
  g0_txt <- fmt_prop_ci(ci0["p"], ci0["lo"], ci0["hi"])
  g1_txt <- fmt_prop_ci(ci1["p"], ci1["lo"], ci1["hi"])
  
  mat <- matrix(c(ev0, tot0 - ev0, ev1, tot1 - ev1), nrow = 2, byrow = TRUE)
  pval <- suppressWarnings(chisq.test(mat))$p.value
  
  tibble(
    Variable = var,
    Group0   = g0_txt,
    Group1   = g1_txt,
    p_value  = signif(pval, 3)
  )
}

event_for <- function(v) if (v == "sex") 2 else 1
base_chi <- map_dfr(cat_vars, ~ base_cat(cleaned_inversed, .x, event_level = event_for(.x)))

#####univariate logistic regression
cat_vars2 <- setdiff(cat_vars, c("immediate_dwi"))
cat_vars2 <- union(cat_vars2, "group")
cont_vars2 <- setdiff(cont_vars, "dwi_count")

or_ci <- function(fit, term){
  b  <- coef(fit)[term]
  ci <- confint(fit, parm = term)
  tibble(OR = unname(exp(b)),
         CI_low  = unname(exp(ci[1])),
         CI_high = unname(exp(ci[2])))
}

uni_cont <- function(df, var){
  f   <- as.formula(paste("immediate_dwi ~", var))
  fit <- glm(f, data = df, family = binomial())
  s   <- summary(fit)$coefficients
  out <- or_ci(fit, term = var) %>%
    mutate(Variable = var,
           Type = "continuous",
           p_value = s[2,4],
           Note = NA_character_)
  out[, c("Variable", "OR","CI_low","CI_high","p_value")]
}

uni_cat <- function(df, var){
  x <- df[[var]]
  k <- length(unique(x))
  
  f   <- as.formula(paste("immediate_dwi ~", var))
  fit <- glm(f, data = df, family = binomial())
  term_name <- colnames(model.matrix(f, df))[2]
  s <- summary(fit)$coefficients
  out <- or_ci(fit, term = term_name) %>%
    mutate(Variable = var,
           Type = "categorical(2)",
           p_value = s[2,4],
           Note = paste0("ref = ", levels(factor(x))[1]))
  out[, c("Variable","OR","CI_low","CI_high","p_value")]
}

uni_cont_res <- map_dfr(cont_vars2, ~ uni_cont(cleaned_inversed, .x))
uni_cat_res  <- map_dfr(cat_vars2,  ~ uni_cat (cleaned_inversed, .x))

#####multivariate logistic regression
multi_vars <- c("smoking_current", "group", "sex", "tx_location")
multi_model <- glm(
  formula = as.formula(paste("immediate_dwi ~ ", paste(multi_vars, collapse = "+"))),
  data = cleaned_inversed,
  family = binomial
)

summary(multi_model)

###### table 1
library(dplyr)
library(tibble)
library(purrr)

confs <- c("sex", "age" , "bmi", "htn", "dm", "dyslipidemia", "cva_hx",
           "smoking_current", "anti_thrombotics_prev", "anti_plt_resistance",
           "tx_location", "an_morphology")

cont_set <- c("age")  # age만 연속형

# event 레벨 규칙: sex=2, 그 외=1
event_level_for <- function(v) if (v == "sex") 2 else 1

cleaned_inversed_weighted <- cleaned_inversed_weighted |>
  mutate(sex = factor(sex, levels = c(1, 2)))

get_smds <- function(bal_obj, var) {
  B  <- bal_obj$Balance
  rn <- rownames(B)
  
  if (var %in% cont_set) {
    idx <- which(rn == var)
    if (length(idx) == 1) return(c(B[idx, "Diff.Un"], B[idx, "Diff.Adj"]))
    else return(c(NA_real_, NA_real_))
  }
  ev  <- as.character(event_level_for(var))
  key <- paste0(var, "_", ev)
  idx <- which(rn == key)
  
  # 없으면 var로 시작하는 행 중 첫 번째 사용(백업)
  if (length(idx) == 0) {
    cand <- which(grepl(paste0("^", var, "([_:].*)?$"), rn))
    if (length(cand) == 0) return(c(NA_real_, NA_real_))
    # 여러 개면 |Diff.Un| 가장 큰 거
    idx <- cand[which.max(abs(B[cand, "Diff.Un"]))]
  }
  
  c(B[idx, "Diff.Un"], B[idx, "Diff.Adj"])
}

one_row <- function(df, var, bal_obj) {
  g  <- df$group
  w  <- df$sw
  x  <- df[[var]]
  
  smds <- get_smds(bal_obj, var)
  smd_b <- smds[1]
  smd_a <- smds[2]
  
  if (var %in% cont_set) {
    # age: 평균/가중평균
    g0_b <- mean(x[g == 0], na.rm = TRUE)
    g1_b <- mean(x[g == 1], na.rm = TRUE)
    
    g0_a <- weighted.mean(x[g == 0], w = w[g == 0], na.rm = TRUE)
    g1_a <- weighted.mean(x[g == 1], w = w[g == 1], na.rm = TRUE)
    
    tibble(
      variables        = var,
      `group 0 (before)` = round(g0_b, 2),
      `group 1 (before)` = round(g1_b, 2),
      `SMD (before)`     = round(smd_b, 3),
      `group 0 (after)`  = round(g0_a, 2),
      `group 1 (after)`  = round(g1_a, 2),
      `SMD (after)`      = round(smd_a, 3)
    )
  } else {
    # 범주형: event 수 / 가중 event 수
    ev  <- event_level_for(var)
    n0_b <- sum(x[g == 0] == ev, na.rm = TRUE)
    n1_b <- sum(x[g == 1] == ev, na.rm = TRUE)
    
    n0_a <- sum(w[g == 0 & x == ev], na.rm = TRUE)
    n1_a <- sum(w[g == 1 & x == ev], na.rm = TRUE)
    
    tibble(
      variables        = var,
      `group 0 (before)` = n0_b,
      `group 1 (before)` = n1_b,
      `SMD (before)`     = round(smd_b, 3),
      `group 0 (after)`  = round(n0_a, 1),  # 가중 event "수" (소수점 발생 가능)
      `group 1 (after)`  = round(n1_a, 1),
      `SMD (after)`      = round(smd_a, 3)
    )
  }
}

table1_ipwt <- map_dfr(confs, ~ one_row(cleaned_inversed_weighted, .x, bal)) |>
  arrange(variables)

table1_ipwt

#####table 2 
inc_by  <- survey::svyby(~immediate_dwi, ~group, design_iptw, survey::svymean) %>% as.data.frame()
mean_by <- survey::svyby(~dwi_count,    ~group, design_iptw, survey::svymean) %>% as.data.frame()

inc_g0  <- 100 * inc_by$immediate_dwi[inc_by$group == 0]
inc_g1  <- 100 * inc_by$immediate_dwi[inc_by$group == 1]
cnt_g0  <-      mean_by$dwi_count     [mean_by$group == 0]
cnt_g1  <-      mean_by$dwi_count     [mean_by$group == 1]

## 2) 효과크기(OR/RR), 95% CI, p-value 추출
# (a) IPTW-Logistic (quasibinomial): OR
b_lr   <- coef(fit_lr)["group"]
ci_lr  <- confint(fit_lr)["group", ]             # 로그 오즈 스케일
OR     <-  exp(b_lr)
OR_lo  <-  exp(ci_lr[1])
OR_hi  <-  exp(ci_lr[2])
p_lr   <- summary(fit_lr)$coefficients["group", "Pr(>|t|)"]

# (b) IPTW-Quasi-Poisson: RR
b_qp   <- coef(fit_qp)["group"]
ci_qp  <- confint(fit_qp)["group", ]             # 로그 rate 스케일
RR     <-  exp(b_qp)
RR_lo  <-  exp(ci_qp[1])
RR_hi  <-  exp(ci_qp[2])
p_qp   <- summary(fit_qp)$coefficients["group", "Pr(>|t|)"]

## 3) 표 생성 (요청 포맷)
table2 <- tribble(
  ~result,                    ~`result group0`,                  ~`result group1`,                  ~`OR or RR (95% CI)`,                                  ~`p-value`,
  "immediate_dwi incidence",  sprintf("%.1f%%", inc_g0),         sprintf("%.1f%%", inc_g1),         sprintf("OR %.2f (95%% CI %.2f–%.2f)", OR, OR_lo, OR_hi), sprintf("%.3f", p_lr),
  "dwi_count",                sprintf("%.2f",  cnt_g0),          sprintf("%.2f",  cnt_g1),          sprintf("RR %.2f (95%% CI %.2f–%.2f)", RR, RR_lo, RR_hi), sprintf("%.3f", p_qp)
)

table2

#####table 3 
uni_all <- rbind(uni_cont_res, uni_cat_res) |> 
  arrange(p_value)

multi_keep <- c("smoking_current","group","sex","tx_location")

multi_sum <- summary(multi_model)$coefficients
multi_ci  <- confint(multi_model)

multi_res <- tibble(
  Variable  = rownames(multi_sum),
  OR_m      = exp(multi_sum[, "Estimate"]),
  CI_low_m  = exp(multi_ci[, 1]),
  CI_high_m = exp(multi_ci[, 2]),
  p_value_m = multi_sum[, "Pr(>|z|)"]
) |> 
  filter(Variable %in% multi_keep)

table3_raw <- uni_all |> 
  left_join(multi_res, by = c("Variable" = "Variable"))

table3 <- table3_raw |> 
  transmute(
    Variable,
    `Odds ratio (uni)`   = sprintf("%.2f", OR),
    `95%% CI (uni)`      = sprintf("%.2f–%.2f", CI_low, CI_high),
    `p value (uni)`      = sprintf("%.3f", p_value),
    `Odds ratio (multi)` = ifelse(is.na(OR_m), "", sprintf("%.2f", OR_m)),
    `95%% CI (multi)`    = ifelse(is.na(CI_low_m), "", sprintf("%.2f–%.2f", CI_low_m, CI_high_m)),
    `p value (multi)`    = ifelse(is.na(p_value_m), "", sprintf("%.3f", p_value_m))
  ) |> 
  arrange(as.numeric(`p value (uni)`))

# 출력
table3