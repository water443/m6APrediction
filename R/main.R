
#' Predict m6A status for multiple candidate sites
#'
#' Use a trained random forest model to predict m6A status
#' for multiple RNA sites.
#'
#' @param ml_fit A trained random forest classification model.
#' @param feature_df A data frame containing the seven required
#'   RNA-related features, with one row per site.
#' @param positive_threshold Numeric probability cutoff between
#'   0 and 1. Scores strictly greater than this value are Positive.
#'
#' @return The input data frame with two additional columns:
#'   predicted_m6A_prob and predicted_m6A_status.
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' features <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' predictions <- prediction_multiple(
#'   ml_fit, head(features, 3), positive_threshold = 0.6
#' )
#' predictions[, c("predicted_m6A_prob", "predicted_m6A_status")]

prediction_multiple <- function(
    ml_fit,
    feature_df,
    positive_threshold = 0.5
) {
  
  # Check input data
  if (!is.data.frame(feature_df)) {
    stop("feature_df must be a data frame.")
  }
  
  if (nrow(feature_df) == 0L) {
    stop("feature_df must contain at least one row.")
  }
  
  if (anyDuplicated(names(feature_df))) {
    stop("feature_df must not contain duplicated column names.")
  }
  
  required_features <- c(
    "gc_content",
    "RNA_type",
    "RNA_region",
    "exon_length",
    "distance_to_junction",
    "evolutionary_conservation",
    "DNA_5mer"
  )
  
  missing_features <- setdiff(
    required_features,
    names(feature_df)
  )
  
  if (length(missing_features) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  # Check prediction threshold
  if (
    !is.numeric(positive_threshold) ||
    length(positive_threshold) != 1L ||
    is.na(positive_threshold) ||
    !is.finite(positive_threshold) ||
    positive_threshold < 0 ||
    positive_threshold > 1
  ) {
    stop("positive_threshold must be a single finite number between 0 and 1.")
  }
  
  # Check numerical features
  numeric_features <- c(
    "gc_content",
    "exon_length",
    "distance_to_junction",
    "evolutionary_conservation"
  )
  
  for (nm in numeric_features) {
    
    if (!is.numeric(feature_df[[nm]])) {
      stop(nm, " must be numeric.")
    }
    
    if (
      anyNA(feature_df[[nm]]) ||
      any(!is.finite(feature_df[[nm]]))
    ) {
      stop(nm, " must not contain NA or non-finite values.")
    }
  }
  
  # GC content must be between 0 and 1
  if (any(
    feature_df$gc_content < 0 |
    feature_df$gc_content > 1
  )) {
    stop("gc_content must be between 0 and 1.")
  }
  
  # Check random forest model
  if (
    !inherits(ml_fit, "randomForest") ||
    !identical(ml_fit$type, "classification")
  ) {
    stop("ml_fit must be a random forest classification model.")
  }
  
  if (
    is.null(ml_fit$classes) ||
    !("Positive" %in% ml_fit$classes)
  ) {
    stop("The model must contain a Positive class.")
  }
  
  # Obtain factor levels from the fitted model
  
  # Obtain factor levels and predictors from the fitted model
  model_xlevels <- ml_fit$forest$xlevels
  
  model_predictors <- names(ml_fit$forest$ncat)
  
  if (
    is.null(model_xlevels) ||
    is.null(model_predictors)
  ) {
    stop("The model does not contain the expected predictor information.")
  }
  
  # Check required biological features
  required_base_predictors <- c(
    "gc_content",
    "RNA_type",
    "RNA_region",
    "exon_length",
    "distance_to_junction",
    "evolutionary_conservation"
  )
  
  missing_base_predictors <- setdiff(
    required_base_predictors,
    model_predictors
  )
  
  if (length(missing_base_predictors) > 0L) {
    stop(
      "The model is missing required predictors: ",
      paste(missing_base_predictors, collapse = ", ")
    )
  }
  
  # Identify nucleotide positions used by the model
  nucleotide_columns <- intersect(
    paste0("nt_pos", 1:5),
    model_predictors
  )
  
  if (length(nucleotide_columns) == 0L) {
    stop("The model does not contain DNA nucleotide predictors.")
  }
  
  
  # Check categorical features
  allowed_categories <- list(
    RNA_type = c(
      "mRNA", "lincRNA", "lncRNA", "pseudogene"
    ),
    RNA_region = c(
      "CDS", "intron", "3'UTR", "5'UTR"
    )
  )
  
  output_df <- feature_df
  
  for (nm in names(allowed_categories)) {
    
    if (!(
      is.character(feature_df[[nm]]) ||
      is.factor(feature_df[[nm]])
    )) {
      stop(nm, " must be a character or factor variable.")
    }
    
    model_levels <- as.character(model_xlevels[[nm]])
    
    if (!setequal(
      model_levels,
      allowed_categories[[nm]]
    )) {
      stop(
        "The model factor levels for ", nm,
        " do not match the expected categories."
      )
    }
    
    values <- as.character(feature_df[[nm]])
    
    if (anyNA(values)) {
      stop(nm, " must not contain missing values.")
    }
    
    if (any(!values %in% model_levels)) {
      stop(
        "Invalid values in ", nm,
        ". Please check the allowed categories."
      )
    }
    
    output_df[[nm]] <- factor(
      values,
      levels = model_levels
    )
  }
  
  # Check DNA 5-mer sequences
  if (!(
    is.character(feature_df$DNA_5mer) ||
    is.factor(feature_df$DNA_5mer)
  )) {
    stop("DNA_5mer must be a character or factor variable.")
  }
  
  dna_strings <- as.character(feature_df$DNA_5mer)
  
  if (anyNA(dna_strings)) {
    stop("DNA_5mer must not contain missing values.")
  }
  
  if (any(!grepl(
    "^[ATCG]{2}A[ATCG]{2}$",
    dna_strings
  ))) {
    stop(
      "DNA_5mer must contain five uppercase A/T/C/G bases ",
      "with A at the central position."
    )
  }
  
  # Encode DNA sequences
  dna_df <- dna_encoding(dna_strings)
  
  # Match nucleotide factor levels to the model
  
  # Match nucleotide factor levels to the model
  for (nm in nucleotide_columns) {
    
    model_levels <- model_xlevels[[nm]]
    
    if (is.null(model_levels) || length(model_levels) == 0L) {
      stop("Missing nucleotide factor levels for ", nm)
    }
    
    model_levels <- as.character(model_levels)
    
    if (!setequal(
      model_levels,
      c("A", "T", "C", "G")
    )) {
      stop("Unexpected nucleotide levels for ", nm)
    }
    
    dna_df[[nm]] <- factor(
      as.character(dna_df[[nm]]),
      levels = model_levels
    )
  }
  
  
  # Avoid duplicated encoded feature columns
  if (any(nucleotide_columns %in% names(output_df))) {
    stop(
      "Please provide the original seven features only. ",
      "DNA position columns are generated automatically."
    )
  }
  
  # Prepare prediction data
  prediction_data <- cbind(
    output_df,
    dna_df
  )
  
  # Keep only predictors used by the fitted model
  missing_model_columns <- setdiff(
    model_predictors,
    names(prediction_data)
  )
  
  if (length(missing_model_columns) > 0L) {
    stop(
      "Missing model predictors: ",
      paste(missing_model_columns, collapse = ", ")
    )
  }
  
  prediction_data <- prediction_data[
    ,
    model_predictors,
    drop = FALSE
  ]
  
  
  # Predict m6A probabilities
  prediction_result <- predict(
    ml_fit,
    newdata = prediction_data,
    type = "prob"
  )
  
  if (
    is.null(dim(prediction_result)) ||
    is.null(colnames(prediction_result)) ||
    !("Positive" %in% colnames(prediction_result))
  ) {
    stop("The model did not return a Positive-class probability.")
  }
  
  positive_probability <- as.numeric(
    prediction_result[, "Positive"]
  )
  
  if (
    length(positive_probability) != nrow(feature_df) ||
    anyNA(positive_probability) ||
    any(!is.finite(positive_probability))
  ) {
    stop("The model returned invalid prediction probabilities.")
  }
  
  # Add prediction results
  output_df$predicted_m6A_prob <- positive_probability
  
  output_df$predicted_m6A_status <- ifelse(
    positive_probability > positive_threshold,
    "Positive",
    "Negative"
  )
  
  return(output_df)
}




#' Convert DNA sequences to nucleotide-position factors
#'
#' Split equal-length DNA sequences into individual nucleotide
#' positions and convert them into categorical factors.
#'
#' @param dna_strings A non-empty character or factor vector
#'   containing equal-length uppercase DNA sequences (A, T, C, G).
#' @return A data frame containing one factor column for each
#'   nucleotide position, with levels A, T, C and G.
#'
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "AGACT"))

dna_encoding <- function(dna_strings) {
  
  # Check input type
  if (!(
    is.character(dna_strings) ||
    is.factor(dna_strings)
  )) {
    stop("dna_strings must be a character or factor vector.")
  }
  
  dna_strings <- as.character(dna_strings)
  
  if (length(dna_strings) == 0L) {
    stop("dna_strings must not be empty.")
  }
  
  if (anyNA(dna_strings)) {
    stop("dna_strings must not contain missing values.")
  }
  
  # Check nucleotide characters
  if (any(!grepl("^[ATCG]+$", dna_strings))) {
    stop("DNA sequences must contain uppercase A, T, C or G only.")
  }
  
  # Check sequence lengths
  sequence_lengths <- nchar(dna_strings)
  
  if (length(unique(sequence_lengths)) != 1L) {
    stop("All DNA sequences must have the same length.")
  }
  
  # Split DNA sequences into nucleotide positions
  nn <- sequence_lengths[1]
  
  seq_m <- matrix(
    unlist(strsplit(dna_strings, split = "", fixed = TRUE)),
    ncol = nn,
    byrow = TRUE
  )
  
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  
  seq_df <- as.data.frame(
    seq_m,
    stringsAsFactors = FALSE
  )
  
  # Convert nucleotide positions to factors
  seq_df[] <- lapply(
    seq_df,
    factor,
    levels = c("A", "T", "C", "G")
  )
  
  return(seq_df)
}





#' Predict m6A status for a single candidate site
#'
#' Predict the m6A probability and classification status
#' for one RNA site using a trained random forest model.
#'
#' @param ml_fit A trained random forest classification model.
#' @param gc_content GC proportion between 0 and 1.
#' @param RNA_type RNA type: mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region RNA region: CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length Log2-transformed exon length.
#' @param distance_to_junction Log2-transformed distance to splice junction.
#' @param evolutionary_conservation Numerical conservation value.
#' @param DNA_5mer Five-base DNA sequence with A at the central position.
#' @param positive_threshold Probability cutoff between 0 and 1.
#'   Scores strictly greater than this value are Positive.
#'
#' @return A named character vector containing predicted_m6A_prob
#'   and predicted_m6A_status.
#'
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' prediction_single(
#'   ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )

prediction_single <- function(
    ml_fit,
    gc_content,
    RNA_type,
    RNA_region,
    exon_length,
    distance_to_junction,
    evolutionary_conservation,
    DNA_5mer,
    positive_threshold = 0.5
) {
  
  # Check that each feature contains one value
  single_inputs <- list(
    gc_content,
    RNA_type,
    RNA_region,
    exon_length,
    distance_to_junction,
    evolutionary_conservation,
    DNA_5mer
  )
  
  input_lengths <- vapply(
    single_inputs,
    length,
    integer(1)
  )
  
  if (any(input_lengths != 1L)) {
    stop("Each input feature must contain exactly one value.")
  }
  
  # Create a single-row feature data frame
  feature_df <- data.frame(
    
    gc_content = gc_content,
    
    RNA_type = RNA_type,
    
    RNA_region = RNA_region,
    
    exon_length = exon_length,
    
    distance_to_junction = distance_to_junction,
    
    evolutionary_conservation = evolutionary_conservation,
    
    DNA_5mer = DNA_5mer,
    
    stringsAsFactors = FALSE
  )
  
  # Reuse the multiple prediction function
  prediction_result <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = feature_df,
    positive_threshold = positive_threshold
  )
  
  # Return a named character vector
  returned_vector <- c(
    
    predicted_m6A_prob = as.character(
      prediction_result$predicted_m6A_prob[1]
    ),
    
    predicted_m6A_status = as.character(
      prediction_result$predicted_m6A_status[1]
    )
    
  )
  
  return(returned_vector)
}

