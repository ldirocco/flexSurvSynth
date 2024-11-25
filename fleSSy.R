source("R/extended_syn.r")


fleSSY <- function(input_data, time_column, event_column, covariates, k, seq, pred, n, s,syn_methods=NULL) {

  library("survival")
  library("flexsurv")
  library("stringr")
  library("rpart")

  data <- input_data

  colnames(data)[colnames(data) == time_column] <- "tsurv"
  colnames(data)[colnames(data) == event_column] <- "surv"

  # Fitting the Survival Model on Observed Covariates
  formula <- as.formula(paste("Surv(tsurv, surv) ~ ", paste(covariates, collapse = " + ")))
  surv_model <- flexsurvspline(formula,
                               data = data,
                               scale = "hazard",
                               k = k)

  # Generate the synthesis specification for covariates
  if (is.null(syn_methods)) {
    syn.spec <- list()

    for (var in colnames(data)) {
      if (is.numeric(data[[var]]) && length(unique(data[[var]])) > 10) {
        syn.spec[[var]] <- list("method" = "normrank")  # Continuous variable
      } else if (is.factor(data[[var]])) {
        levels_count <- length(levels(data[[var]]))
        if (levels_count == 2) {
          syn.spec[[var]] <- list("method" = "logistic")  # Binary factor
        } else if (is.ordered(data[[var]])) {
          syn.spec[[var]] <- list("method" = "polr")  # Ordered factor
        } else {
          syn.spec[[var]] <- list("method" = "polyreg")  # Categorical factor
        }
      } else if (is.logical(data[[var]])) {
        syn.spec[[var]] <- list("method" = "logistic")  # Binary logical
      } else {
        syn.spec[[var]] <- list("method" = "normrank")  # Default to normrank
      }
    }
  } else {
    # Use user-provided preferred_methods as syn.spec
    syn.spec <- syn_methods
  }

  n_vars <- ncol(data)
  pred_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
  colnames(pred_matrix) <- rownames(pred_matrix) <- colnames(data)
  for (var in colnames(data)) {
    pred_matrix[var, pred] <- 1
  }

  # Perform the synthesis
  synth.obj <- extended_syn(
    data = data,
    method = syn.spec,
    visit.sequence = seq,
    pred = pred_matrix,
    m = n,
    seed = s
  )

  db_syn = synth.obj$syn
  db_syn = db_syn[,seq]

  db_syn_surv_complete = complete(mice(db_syn_surv))
  #seed:161
  synthetic_surv = simulate(surv_model, newdata = db_syn_surv_complete, nsim = 1,
                            tidy = TRUE, seed = 161, censtime = max(original_data$tsurv))

  synthetic_surv$time = round(synthetic_surv$time,0)


  synthetic_surv$tsurv = synthetic_surv$tsurv

  synthetic_surv$surv = synthetic_surv$surv

  synth_db = cbind(synthetic_surv[,(ncol(synthetic_surv) - 1):ncol(synthetic_surv)], db_syn)

  # merging the two datasets
  original_data$type = "Original"
  synth_db$type = "Synthetic"

  db = full_join(original_data,synth_db)

  return(db)
}
