# Illustration : two-stage analysis of generation and detection structures


library(dplyr)
library(tidyr)
library(knitr)

n <- 1000
B <- 100
n_random <- 20
seed <- 123
output_dir <- "numerical_illustrations_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

generate_data <- function(
  n = 1000, mechanism = c("CRCR", "CRNAR"), crnar_rho = -0.5, j = 10, response_location = 3,
  item_thresholds = c(1.5, 2.5, 3.5, 4.5),
  lambda_min = 0.60, lambda_max = 0.85,
  careless_prevalence = 0.30, cr_type = c("straightlining", "random"),
  sens = 0.60, spec = 0.80,
  time_careless_effect = -0.70, time_c_effect = -0.40,
  time_error_sd = 1, time_cutoff = 200, pattern_cutoff = 4,
  beta_c_y = 0.50, y_intercept = 0, y_error_sd = 1
) {
  cr_type <- match.arg(cr_type)
  mechanism <- match.arg(mechanism)
  stopifnot(abs(crnar_rho) < 1)
  stopifnot(
    n > 0, j > 1, n == as.integer(n), j == as.integer(j),
    careless_prevalence > 0, careless_prevalence < 1,
    lambda_min > 0, lambda_max < 1, lambda_min <= lambda_max,
    sens >= 0, sens <= 1, spec >= 0, spec <= 1,
    length(item_thresholds) == 4, all(is.finite(item_thresholds)),
    all(diff(item_thresholds) > 0), time_error_sd >= 0,
    y_error_sd > 0
  )
  C <- rnorm(n)
  carelessness_error <- rnorm(n)
  Carelessness <- if (mechanism == "CRCR") carelessness_error else
    crnar_rho * C + sqrt(1 - crnar_rho^2) * carelessness_error
  # In CRNAR, unobserved C causes both clean responses and response state.
  # Carelessness remains marginally N(0,1), preserving population prevalence.
  # C is simulation truth, NOT an observed adjustment variable.
  R_CR <- as.integer(Carelessness < qnorm(1 - careless_prevalence))
  lambda <- runif(j, lambda_min, lambda_max)
  item_error_sd <- sqrt(1 - lambda^2)
  V_star_continuous <- matrix(NA_real_, n, j)
  for (item in seq_len(j)) {
    V_star_continuous[, item] <- response_location + lambda[item] * C +
      rnorm(n, sd = item_error_sd[item])
  }
  V_star <- matrix(NA_integer_, n, j)
  for (item in seq_len(j)) {
    V_star[, item] <- findInterval(V_star_continuous[, item], item_thresholds) + 1L
  }
  if (cr_type == "straightlining") {
    U <- matrix(sample(1:5, n, replace = TRUE), n, j)
  } else {
    U <- matrix(sample(1:5, n * j, replace = TRUE), n, j)
  }
  storage.mode(U) <- "integer"
  V_obs <- V_star
  V_obs[R_CR == 0, ] <- U[R_CR == 0, ]
  CR_binary <- 1L - R_CR
  bogus_flag <- rbinom(n, 1, ifelse(CR_binary == 1, sens, 1 - spec))
  D_bogus <- 1L - bogus_flag
  log_total_time <- log(300) + time_careless_effect * Carelessness +
    time_c_effect * C + rnorm(n, sd = time_error_sd)
  total_time <- exp(log_total_time)
  D_time <- as.integer(total_time >= time_cutoff)
  Q_longstring <- apply(V_obs, 1, function(x) max(rle(x)$lengths))
  D_pattern <- as.integer(Q_longstring < pattern_cutoff)
  Y <- y_intercept + beta_c_y * C + rnorm(n, sd = y_error_sd)
  list(
    C = C, Carelessness = Carelessness, R_CR = R_CR,
    CR_binary = CR_binary, bogus_flag = bogus_flag,
    V_star = V_star, V_star_continuous = V_star_continuous,
    V_obs = V_obs, U = U, lambda = lambda, item_error_sd = item_error_sd,
    item_thresholds = item_thresholds, Y = Y, beta_c_y = beta_c_y,
    total_time = total_time, log_total_time = log_total_time,
    Q_longstring = Q_longstring, D_oracle = R_CR,
    D_pattern = D_pattern, D_bogus = D_bogus, D_time = D_time,
    response_location = response_location,
    careless_prevalence = careless_prevalence, cr_type = cr_type
  )
}

# Uniform random deletion at the EXACT observed retained counts in each stratum.
# sample.int handles singleton strata correctly.
matched_random <- function(r, keep, clean_score, observed_score, draws) {
  ids_a <- which(r == 1)
  ids_c <- which(r == 0)
  ka <- sum(keep[ids_a])
  kc <- sum(keep[ids_c])
  pick <- function(ids, k) if (k == 0) integer(0) else ids[sample.int(length(ids), k)]
  values <- replicate(draws, {
    ia <- pick(ids_a, ka)
    ic <- pick(ids_c, kc)
    stopifnot(length(ia) == ka, length(ic) == kc,
              all(r[ia] == 1), all(r[ic] == 0), !anyDuplicated(c(ia, ic)))
    selected <- c(ia, ic)
    c(clean = mean(clean_score[selected]), observed = mean(observed_score[selected]))
  })
  c(clean = mean(values["clean", ]), observed = mean(values["observed", ]),
    sensitivity = 1 - kc / length(ids_c), specificity = ka / length(ids_a))
}

analyze_dataset <- function(dat, mechanism, n_random = 20) {
  r <- dat$R_CR
  a <- r == 1
  stopifnot(any(a), any(!a), n_random >= 1)
  clean_score <- rowMeans(dat$V_star)
  observed_score <- rowMeans(dat$V_obs)
  full <- mean(clean_score)
  true_state <- mean(observed_score[a])
  generation_bias <- true_state - full
  naive_bias <- mean(observed_score) - full
  replacement <- mean(!a) * mean(observed_score[!a] - clean_score[!a])

  # STAGE 1: true response state, before using any detector.
  generation <- tibble(
    full = full, naive = mean(observed_score), true_state = true_state,
    naive_bias = naive_bias, generation_bias = generation_bias,
    deletion_change = true_state - mean(observed_score),
    contamination_replacement = replacement,
    prop_attentive = mean(a),
    check_replacement = naive_bias - replacement,
    check_generation = generation_bias - naive_bias -
      (true_state - mean(observed_score)),
    check_true_state = mean(observed_score[a] - clean_score[a])
  )

  # Stop after generation analysis for CRNAR.
  # Detector variables remain in generate_data to preserve its RNG sequence.
  if (mechanism != "CRCR") {
    return(list(generation = generation, detection = NULL))
  }

  # STAGE 2: detector structure is evaluated ONLY under CRCR.
  # Total bias includes remaining contamination; excess_total isolates the
  # contrast with random deletion at the same realized accuracy, not a
  # universal causal effect independent of the response distribution.
  strategies <- list(Bogus = dat$D_bogus == 1, Time = dat$D_time == 1,
                     Pattern = dat$D_pattern == 1)
  detection <- bind_rows(lapply(names(strategies), function(method) {
    keep <- strategies[[method]]
    if (!any(keep)) stop("No retained respondents: ", method)
    clean <- mean(clean_score[keep])
    observed <- mean(observed_score[keep])
    random <- matched_random(r, keep, clean_score, observed_score, n_random)
    selection <- clean - full
    contamination <- observed - clean
    total <- observed - full
    detection_selection <- clean - true_state
    detection_change <- observed - true_state
    random_selection <- unname(random["clean"]) - full
    random_contamination <- unname(random["observed"] - random["clean"])
    random_total <- unname(random["observed"]) - full
    excess_selection <- selection - random_selection
    excess_contamination <- contamination - random_contamination
    excess_total <- total - random_total
    sensitivity <- mean(!keep[!a])
    specificity <- mean(keep[a])
    tibble(
      method = method, full = full, true_state = true_state,
      clean = clean, observed = observed,
      generation_bias = generation_bias,
      detection_selection = detection_selection,
      contamination = contamination, detection_change = detection_change,
      selection = selection, total = total,
      random_clean = unname(random["clean"]),
      random_observed = unname(random["observed"]),
      random_selection = random_selection, random_contamination = random_contamination,
      random_total = random_total,
      excess_selection = excess_selection, excess_contamination = excess_contamination,
      excess_total = excess_total,
      sensitivity = sensitivity, specificity = specificity,
      n_retained = sum(keep), n_attentive_retained = sum(keep & a),
      n_careless_retained = sum(keep & !a), retained_careless = mean(!a[keep]),
      check_total = total - selection - contamination,
      check_detection = detection_change - detection_selection - contamination,
      check_bridge = total - generation_bias - detection_change,
      check_selection = selection - generation_bias - detection_selection,
      check_random = random_total - random_selection - random_contamination,
      check_excess = excess_total - excess_selection - excess_contamination,
      check_structure = total - random_total - excess_total,
      check_sensitivity = sensitivity - unname(random["sensitivity"]),
      check_specificity = specificity - unname(random["specificity"])
    )
  }))
  list(generation = generation, detection = detection)
}

conditions <- expand.grid(mechanism = c("CRCR", "CRNAR"),
                          response_location = c(3, 4),
                          cr_type = c("straightlining", "random"),
                          stringsAsFactors = FALSE)
generation_list <- detection_list <- vector("list", nrow(conditions) * B)
index <- 0L
# Independent reproducible streams prevent comparator draws from changing data.
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
data_rng <- .Random.seed
comparator_rng <- parallel::nextRNGStream(data_rng)
for (condition in seq_len(nrow(conditions))) {
  cfg <- conditions[condition, ]
  message("Condition ", condition, "/", nrow(conditions), ": ",
          cfg$mechanism, ", location=", cfg$response_location, ", ", cfg$cr_type)
  for (b in seq_len(B)) {
    .Random.seed <- data_rng
    dat <- generate_data(n = n, mechanism = cfg$mechanism,
                         response_location = cfg$response_location, cr_type = cfg$cr_type)
    data_rng <- .Random.seed
    .Random.seed <- comparator_rng
    result <- analyze_dataset(dat, mechanism = cfg$mechanism, n_random = n_random)
    comparator_rng <- .Random.seed
    index <- index + 1L
    add_keys <- function(x) mutate(x, mechanism = cfg$mechanism,
      response_location = cfg$response_location, cr_type = cfg$cr_type,
      replication = b, sample_size = n)
    generation_list[[index]] <- add_keys(result$generation)
    if (!is.null(result$detection)) {
      detection_list[[index]] <- add_keys(result$detection)
    }
  }
}
generation_results <- bind_rows(generation_list)
detection_results <- bind_rows(detection_list)
# Scope and row-count checks prevent CRNAR from entering Stage 2 outputs.
stopifnot(setequal(unique(generation_results$mechanism), c("CRCR", "CRNAR")),
          identical(unique(detection_results$mechanism), "CRCR"),
          nrow(generation_results) == 8L * B,
          nrow(detection_results) == 4L * B * 3L)
keys <- c("mechanism", "response_location", "cr_type")

component_summary <- function(data, components, groups) {
  data |>
    select(all_of(groups), all_of(components)) |>
    pivot_longer(all_of(components), names_to = "component", values_to = "value") |>
    group_by(across(all_of(groups)), component) |>
    summarise(replications = n(), mean = mean(value), sd = sd(value),
              mc_se = sd(value) / sqrt(n()), ci_lo = mean - 1.96 * mc_se,
              ci_hi = mean + 1.96 * mc_se, .groups = "drop")
}
stage1_summary <- component_summary(generation_results,
  c("naive_bias", "generation_bias", "deletion_change"), keys)
stage2_summary <- component_summary(detection_results,
  c("generation_bias", "detection_selection", "contamination",
    "detection_change", "selection", "total"), c(keys, "method"))
structure_summary <- component_summary(detection_results,
  c("random_selection", "random_contamination", "random_total",
    "excess_selection", "excess_contamination", "excess_total"), c(keys, "method"))
accuracy_summary <- detection_results |>
  group_by(across(all_of(c(keys, "method")))) |>
  summarise(across(c(sensitivity, specificity, n_retained, n_attentive_retained,
                     n_careless_retained, retained_careless), mean), .groups = "drop")

# Separate estimator summaries. Random is the average of 20 matched draws,
# so its SD/RMSE are NOT the performance of a single random-deletion procedure.
estimator_results <- bind_rows(
  generation_results |> select(all_of(keys), replication, full, naive, true_state) |>
    pivot_longer(c(naive, true_state), names_to = "method", values_to = "estimate") |>
    mutate(comparator = "Reference"),
  detection_results |> select(all_of(keys), replication, full, method, observed,
                              random_observed) |>
    pivot_longer(c(observed, random_observed), names_to = "comparator",
                 values_to = "estimate")
) |> mutate(error = estimate - full)
table_estimators <- estimator_results |>
  group_by(across(all_of(c(keys, "method", "comparator")))) |>
  summarise(replications = n(), mean_target = mean(full), mean_estimate = mean(estimate),
            bias = mean(error), empirical_sd = sd(estimate),
            rmse_paired_reference = sqrt(mean(error^2)),
            mc_se_bias = sd(error) / sqrt(n()), .groups = "drop")

stage1_estimators <- filter(table_estimators, comparator == "Reference")
stage2_estimators <- filter(table_estimators, mechanism == "CRCR")
stopifnot(all(stage2_estimators$mechanism == "CRCR"))

check_summary <- function(data) {
  data |> group_by(across(all_of(keys))) |>
    summarise(across(starts_with("check_"), ~ max(abs(.x))), .groups = "drop")
}
generation_checks <- check_summary(generation_results)
detection_checks <- check_summary(detection_results)
for (checks in list(generation_checks, detection_checks)) {
  values <- as.matrix(select(checks, starts_with("check_")))
  stopifnot(all(is.finite(values)), all(values < 1e-10))
}
# Pattern removes all straightliners; no residual contamination can remain.
stopifnot(all(abs(filter(detection_results, cr_type == "straightlining",
                         method == "Pattern")$contamination) < 1e-10))

summaries <- list(stage1_generation = stage1_summary,
                  stage2_detection = stage2_summary,
                  stage2_matched_structure = structure_summary,
                  detector_accuracy = accuracy_summary,
                  stage1_estimators = stage1_estimators,
                  stage2_estimators = stage2_estimators,
                  generation_checks = generation_checks,
                  detection_checks = detection_checks)
for (name in names(summaries)) {
  write.csv(summaries[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  writeLines(kable(summaries[[name]], format = "latex", booktabs = TRUE,
                   digits = 4, escape = TRUE, caption = gsub("_", " ", name)),
             file.path(output_dir, paste0(name, ".tex")))
}
write.csv(generation_results, file.path(output_dir, "generation_replications.csv"), row.names = FALSE)
write.csv(detection_results, file.path(output_dir, "detection_replications.csv"), row.names = FALSE)
saveRDS(list(generation = generation_results, detection = detection_results,
             summaries = summaries, generation_conditions = conditions,
             detection_conditions = subset(conditions, mechanism == "CRCR"),
             settings = list(n = n, B = B, n_random = n_random, seed = seed,
                             crnar_rho = -0.5)), file.path(output_dir, "all_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
print(stage1_summary, n = Inf, width = Inf)
print(stage2_summary, n = Inf, width = Inf)
print(structure_summary, n = Inf, width = Inf)
print(accuracy_summary, n = Inf, width = Inf)
message("All identity and exact matching checks passed. Results: ", output_dir)



# Figures: horizontal component bars, following the reference layout.
# Positive and negative components stack separately around zero.
# The black diamond is their ALGEBRAIC SUM; opposite-sign bar lengths must
# not be interpreted as the total. These figures show Monte Carlo means;
# confidence intervals remain available in the exported tables.
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("ragg", quietly = TRUE)) {
  stop("Tables saved. Install plotting packages: install.packages(c('ggplot2','ragg')).")
}
library(ggplot2)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Actual decomposition: include the two reference strategies under CRCR.
actual_plot <- bind_rows(
  generation_results |>
    filter(mechanism == "CRCR") |>
    transmute(mechanism, response_location, cr_type, replication,
              method = "Retain all", selection = 0,
              contamination = naive_bias, total = naive_bias),
  detection_results |>
    select(mechanism, response_location, cr_type, replication, method,
           selection, contamination, total),
  generation_results |>
    filter(mechanism == "CRCR") |>
    transmute(mechanism, response_location, cr_type, replication,
              method = "True-state deletion", selection = generation_bias,
              contamination = 0, total = generation_bias)
) |>
  group_by(mechanism, response_location, cr_type, method) |>
  summarise(across(c(selection, contamination, total), mean), .groups = "drop")

# Generation plot retains both mechanisms, without any detector.
generation_plot <- bind_rows(
  generation_results |>
    transmute(mechanism, response_location, cr_type, replication,
              method = "Retain all", selection = 0,
              contamination = naive_bias, total = naive_bias),
  generation_results |>
    transmute(mechanism, response_location, cr_type, replication,
              method = "True-state deletion", selection = generation_bias,
              contamination = 0, total = generation_bias)
) |>
  group_by(mechanism, response_location, cr_type, method) |>
  summarise(across(c(selection, contamination, total), mean), .groups = "drop")

# Difference plot: compare each detector with its OWN matched random baseline.
structure_plot <- detection_results |>
  group_by(mechanism, response_location, cr_type, method) |>
  summarise(selection = mean(excess_selection),
            contamination = mean(excess_contamination),
            total = mean(excess_total), .groups = "drop")



make_bar_plot <- function(data, difference = FALSE,
                          include_mechanism = FALSE) {
  data <- data |>
    mutate(
      method = factor(method, levels = rev(c(
        "Retain all", "Bogus", "Time", "Pattern", "True-state deletion"
      ))),
      Location = factor(
        response_location, levels = c(3, 4),
        labels = c("Response location = 3", "Response location = 4")
      ),
      Pattern = factor(
        cr_type, levels = c("straightlining", "random"),
        labels = c("Straightlining", "Random responses")
      ),
      mechanism = factor(mechanism, levels = c("CRCR", "CRNAR"))
    )
  
  bars <- data |>
    pivot_longer(
      c(selection, contamination),
      names_to = "Component",
      values_to = "Change"
    ) |>
    mutate(
      Component = factor(
        Component, levels = c("selection", "contamination")
      )
    )
  
  legend_labels <- if (difference) {
    c("Selection difference", "Contamination difference")
  } else {
    c("Selection", "Remaining contamination")
  }
  
  p <- ggplot(bars, aes(x = Change, y = method, fill = Component)) +
    geom_vline(
      xintercept = 0, color = "grey55", linewidth = .45
    ) +
    geom_col(
      width = .55,
      position = position_stack(reverse = TRUE),
      orientation = "y"
    ) +
    geom_point(
      data = data,
      aes(x = total, y = method, shape = "Total"),
      inherit.aes = FALSE,
      color = "black", size = 2.7
    ) +
    scale_fill_manual(
      values = c(selection = "#0072B2", contamination = "#E69F00"),
      labels = legend_labels,
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(Total = 18),
      labels = if (difference) "Total difference" else "Total"
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(.05, .08))
    ) +
    labs(
      x = if (difference) {
        "Change in bias: actual minus matched random deletion"
      } else {
        "Change from the full clean sample"
      },
      y = NULL, fill = NULL, shape = NULL,
      caption = NULL
    ) +
    theme_minimal(base_size = 18) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      strip.text = element_text(size = 17, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      plot.title = element_text(size = 22),
      plot.subtitle = element_text(size = 17),
      plot.tag = element_text(size = 22),
      plot.caption = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.spacing = grid::unit(1.2, "lines"),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = margin(12, 16, 12, 12)
    ) +
    guides(
      fill = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    )
  
  if (include_mechanism) {
    p + facet_grid(mechanism + Pattern ~ Location)
  } else {
    p + facet_grid(Pattern ~ Location)
  }
}

plots <- list(
  figure1_generation = make_bar_plot(
    generation_plot,
    include_mechanism = TRUE
  ),
  figure2_detection_decomposition = make_bar_plot(actual_plot)
)

# Uses your existing figure_dir.
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

save_tikz <- function(plot, filename, width, height,
                      standalone = FALSE) {
  tikzDevice::tikz(
    file = filename,
    width = width,
    height = height,
    standAlone = standalone,
    engine = "pdftex",
    sanitize = TRUE,
    bg = "white"
  )
  
  # Close the device even if printing fails.
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

for (name in names(plots)) {
  figure_height <- if (name == "figure1_generation") 9 else 7
  
  # For inclusion in your manuscript with \input{}.
  save_tikz(
    plot = plots[[name]],
    filename = file.path(figure_dir, paste0(name, ".tex")),
    width = 11,
    height = figure_height
  )
  
  # Complete LaTeX document for separate compilation.
  save_tikz(
    plot = plots[[name]],
    filename = file.path(
      figure_dir, paste0(name, "_standalone.tex")
    ),
    width = 11,
    height = figure_height,
    standalone = TRUE
  )
}

write.csv(
  generation_plot,
  file.path(figure_dir, "figure1_data.csv"),
  row.names = FALSE
)
write.csv(
  actual_plot,
  file.path(figure_dir, "figure2_data.csv"),
  row.names = FALSE
)
saveRDS(plots, file.path(figure_dir, "plot_objects.rds"))

message(
  "Two figures saved as includable and standalone TikZ files to: ",
  normalizePath(figure_dir)
)

      # Record session information for reproducibility
sessionInfo()

# Save session information as a text file
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)
