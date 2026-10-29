# Figure 7F: batch-adjusted associations for the 18 specified display candidates.
# Statistical definitions: Bailey F. Garb; portable adaptation: Shaomiao Xia.
#
# Usage: Rscript analyze_network.R MODULE INPUTS RESULTS PANEL
# All 18 candidates are fitted before applying the nominal p < 0.05 display rule.
# The E6 influence score is continuous. Eight sequencing batches are included
# as a categorical adjustment; this differs from the five-cohort survival model.
# HPV integration uses a binomial logit model (negative = 0, positive = 1).
# Other candidates use ordinary least squares in their original numeric units.
options(stringsAsFactors = FALSE, digits = 17)
args <- commandArgs(TRUE)
stopifnot(length(args) == 4)
root <- normalizePath(args[1])
input <- normalizePath(args[2])
out <- normalizePath(args[3])
d <- read.csv(file.path(input, "network_n230.csv"), check.names = TRUE)
candidates <- read.csv(file.path(root, "data/network_candidates.csv"), check.names = FALSE)
stopifnot(nrow(d) == 230, !anyDuplicated(d$samples), nrow(candidates) == 18)
d$batch <- factor(d$batch)
d$HPV_integration <- factor(d$HPV_integration)
stopifnot(nlevels(d$batch) == 8)
score <- "e6.FL_influence_scores"

models <- lapply(seq_len(nrow(candidates)), function(index) {
    variable <- candidates$R_variable[index]
    logistic <- variable == "HPV_integration"
    z <- d
    if (logistic) {
        stopifnot(identical(levels(z[[variable]]), c("neg", "pos")))
        z[[variable]] <- ifelse(
            is.na(z[[variable]]), NA_real_,
            as.numeric(as.character(z[[variable]]) == "pos")
        )
    }

    # Use candidate-specific complete cases and remove unused batch levels.
    # This retains 204 patients for integration and 230 for numeric candidates.
    z <- droplevels(z[complete.cases(z[, c(variable, score, "batch")]), , drop = FALSE])
    formula <- reformulate(c(score, "batch"), response = variable)
    fit <- if (logistic) glm(formula, data = z, family = binomial()) else lm(formula, data = z)
    if (logistic) stopifnot(fit$converged)
    coefficient <- coef(summary(fit))[score, ]
    beta <- unname(coefficient[1])
    se <- unname(coefficient[2])
    p <- unname(coefficient[4])

    # Preserve the original Wald normal interval for the logistic model and
    # residual-df Student interval for linear models. Edge selection uses the
    # unadjusted model p value; no new multiple-testing correction is applied.
    critical <- if (logistic) qnorm(.975) else qt(.975, df.residual(fit))
    data.frame(
        variable = variable, Variable = candidates$original_metadata_column[index],
        Label = candidates$figure_label[index], N = nrow(z),
        model = if (logistic) "binomial_logit" else "ordinary_linear_model",
        formula = paste(deparse(formula), collapse = " "),
        Beta = beta, SE = se, CI95_lower = beta - critical * se,
        CI95_upper = beta + critical * se, P_value = p,
        direction = if (beta > 0) "positive" else "negative", displayed = p < .05
    )
})
models <- do.call(rbind, models)
write.csv(models, file.path(out, "Figure7F_candidate_models.csv"), row.names = FALSE, na = "")
cat("Calculated 18 candidate associations;", sum(models$displayed), "meet p < 0.05.\n")
