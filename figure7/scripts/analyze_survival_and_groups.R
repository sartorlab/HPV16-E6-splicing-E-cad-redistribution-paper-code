# Figure 7B-E: survival, molecular subtype, and HPV integration analyses.
# Statistical analysis: Bailey F. Garb; portable adaptation: Shaomiao Xia.
#
# Usage: Rscript analyze_survival_and_groups.R MODULE INPUTS RESULTS PANEL
# PANEL is all, 7B, 7C, 7D, or 7E. The runner creates the output directory.
# Inputs are patient-level scores and the original endpoint/annotation snapshots.
# Numerical results include Kaplan-Meier steps, risk tables, adjusted Cox models,
# and the two-group comparisons needed to draw the selected figure panels.
options(stringsAsFactors = FALSE, digits = 17)
suppressPackageStartupMessages(library(survival))
args <- commandArgs(TRUE)
stopifnot(length(args) == 4)
input <- normalizePath(args[2])
out <- normalizePath(args[3])
selection <- args[4]

write_result <- function(x, name) {
    write.csv(x, file.path(out, name), row.names = FALSE, na = "")
}

if (selection %in% c("all", "7B", "7E")) {
    d <- read.csv(file.path(input, "survival_n230.csv"), check.names = TRUE)
    stopifnot(nrow(d) == 230, !anyDuplicated(d$samples))
    d$e6fl <- as.numeric(d$e6.FL_influence_scores)
    d$Age <- as.numeric(d$Age)

    # Preserve the original stage pooling and five-cohort adjustment. The four
    # UM67 sequencing batches form one FFPE clinical cohort for survival models.
    d$cat_stage <- ifelse(
        d$ajcc_updated %in% c("Stage III", "Stage IV"),
        "Stage III/IV", d$ajcc_updated
    )
    d$cohort <- as.character(d$batch)
    d$cohort[d$cohort %in% c("2536-KZ", "3687-SS", "3860-SS", "4165-SS")] <- "FFPE"
    cohort_names <- c(FF = "UM18", "9067-MS" = "UM20", FFPE = "UM67", HVC = "HVC", TCGA = "TCGA")
    stopifnot(!anyNA(unname(cohort_names[d$cohort])))

    genotypes <- if (selection == "7B") {
        "All HPV"
    } else if (selection == "7E") {
        "HPV16"
    } else {
        c("All HPV", "HPV16")
    }
    models <- risk_tables <- steps <- list()
    for (genotype in genotypes) {
        population <- if (genotype == "HPV16") d[d$HPV_status == "HPV16", ] else d

        # The high/low threshold is the genotype-population mean, calculated
        # BEFORE excluding patients without a usable endpoint. It is not the
        # endpoint-specific mean. Cox models below use the continuous score.
        threshold <- mean(population$e6fl, na.rm = TRUE)
        population$g_km <- factor(
            ifelse(population$e6fl >= threshold, "E6FL High", "E6FL Low"),
            levels = c("E6FL Low", "E6FL High")
        )
        for (endpoint in c("OS", "RFS")) {
            time_col <- if (endpoint == "OS") "OS.time" else "Rtime"
            event_col <- if (endpoint == "OS") "OS" else "recurstatus"
            eligible <- !is.na(population[[time_col]]) &
                !is.na(population[[event_col]]) & population[[time_col]] > 0 &
                !is.na(population$g_km) & !is.na(population$e6fl)
            z <- population[eligible, ]

            # Administrative censoring at 60 months: later events are censored,
            # whereas an event exactly at month 60 remains an event.
            z$t_km <- pmin(z[[time_col]], 60)
            z$e_km <- ifelse(z[[time_col]] > 60, 0L, as.integer(z[[event_col]]))
            km <- survfit(Surv(t_km, e_km) ~ g_km, data = z)
            lr <- survdiff(Surv(t_km, e_km) ~ g_km, data = z)
            lr_p <- 1 - pchisq(lr$chisq, df = 1)

            # Preserve coxph's original defaults (including Efron ties and
            # complete-case handling of age). Report HR per 50 score points
            # from the unrounded coefficient and its confidence interval.
            cox_df <- z[!is.na(z$cat_stage) & !is.na(z$cohort), ]
            cox_df$cohort <- droplevels(factor(cox_df$cohort))
            cox <- coxph(Surv(t_km, e_km) ~ e6fl + Age + cat_stage + cohort, data = cox_df)
            cf <- summary(cox)$coefficients["e6fl", ]
            ci <- confint(cox)["e6fl", ]
            panel <- if (genotype == "All HPV") "Figure7B" else "Figure7E"
            key <- paste(panel, endpoint, sep = "_")
            models[[key]] <- data.frame(
                panel = panel, genotype = genotype, endpoint = endpoint,
                population_n = nrow(population), split_mean = threshold,
                km_n = nrow(z), km_low = sum(z$g_km == "E6FL Low"),
                km_high = sum(z$g_km == "E6FL High"), km_events = sum(z$e_km),
                cox_n = cox$n, cox_events = cox$nevent,
                logrank_p = lr_p, cox_p = unname(cf["Pr(>|z|)"]),
                cox_beta_per_point = unname(cf["coef"]),
                cox_se_per_point = unname(cf["se(coef)"]),
                HR_per_50_points = exp(50 * unname(cf["coef"])),
                HR50_ci_low = exp(50 * ci[1]), HR50_ci_high = exp(50 * ci[2]),
                cox_formula = paste0(
                    "Surv(pmin(time,60), event_at_or_before_60) ~ ",
                    "continuous_E6_score + age + pooled_stage + clinical_cohort"
                ), row.names = NULL
            )
            sk <- summary(km, times = seq(0, 60, 12), extend = TRUE)
            risk_tables[[key]] <- data.frame(
                panel = panel, genotype = genotype, endpoint = endpoint,
                group = sub("g_km=", "", as.character(sk$strata)),
                month = sk$time, n_risk = sk$n.risk
            )
            ss <- summary(km, censored = TRUE)
            steps[[key]] <- data.frame(
                panel = panel, genotype = genotype, endpoint = endpoint,
                group = sub("g_km=", "", as.character(ss$strata)),
                month = ss$time, n_risk = ss$n.risk, n_event = ss$n.event,
                n_censor = ss$n.censor, survival = ss$surv,
                lower = ss$lower, upper = ss$upper
            )
        }
    }
    models <- do.call(rbind, models)
    rownames(models) <- NULL
    risks <- do.call(rbind, risk_tables)
    rownames(risks) <- NULL
    write_result(models, "Figure7BE_survival_statistics.csv")
    write_result(risks, "Figure7BE_risk_tables.csv")
    write_result(do.call(rbind, steps), "Figure7BE_KM_steps.csv")
}

if (selection %in% c("all", "7C", "7D")) {
    # Preserve the score-master subtype/integration annotations. Each comparison
    # has 204 observed labels; later clinical-table enrichment is not substituted.
    cd <- read.csv(file.path(input, "subtype_integration_n230.csv"), check.names = TRUE)
    stopifnot(nrow(cd) == 230, !anyDuplicated(cd$samples))
    cd$score <- as.numeric(cd$e6.FL_influence_scores)
    cd$subtype <- ifelse(is.na(cd$subtype_IMU), NA, ifelse(cd$subtype_IMU, "IMU", "KRT"))
    cd$integration <- ifelse(
        is.na(cd$HPV_integration), NA,
        ifelse(cd$HPV_integration == "neg", "Negative", "Positive")
    )
    panels <- if (selection == "all") c("Figure7C", "Figure7D") else paste0("Figure", selection)
    tests <- plot_data <- descriptives <- list()
    for (panel in panels) {
        variable <- if (panel == "Figure7C") "subtype" else "integration"
        levels <- if (panel == "Figure7C") c("IMU", "KRT") else c("Negative", "Positive")
        z <- cd[!is.na(cd[[variable]]) & !is.na(cd$score), ]
        z$group <- factor(z[[variable]], levels = levels)
        stopifnot(nrow(z) == 204)
        plot_data[[panel]] <- data.frame(
            panel = panel, sample_id = z$samples,
            group = as.character(z$group), E6FL_E6ALL_score = z$score
        )
        for (group in levels) {
            values <- z$score[z$group == group]
            descriptives[[length(descriptives) + 1]] <- data.frame(
                panel = panel, group = group, n = length(values),
                mean = mean(values), median = median(values), sd = sd(values)
            )
        }

        # The displayed comparison is an unpaired, two-sided Student t test
        # with equal variances. Group order determines the mean-difference sign.
        test <- t.test(score ~ group, data = z, var.equal = TRUE, alternative = "two.sided")
        tests[[length(tests) + 1]] <- data.frame(
            panel = panel, test = "Student unpaired t-test",
            variance_assumption = "equal", two_sided = TRUE,
            group1 = levels[1], group2 = levels[2],
            n1 = sum(z$group == levels[1]), n2 = sum(z$group == levels[2]),
            mean_difference = unname(test$estimate[1] - test$estimate[2]),
            ci95_low = test$conf.int[1], ci95_high = test$conf.int[2],
            t = unname(test$statistic), df = unname(test$parameter), p_value = test$p.value
        )
    }
    write_result(do.call(rbind, tests), "Figure7CD_ttests.csv")
    write_result(do.call(rbind, plot_data), "Figure7CD_plot_data.csv")
    write_result(do.call(rbind, descriptives), "Figure7CD_group_statistics.csv")
}
cat("Calculated the selected survival or molecular-group analysis.\n")
