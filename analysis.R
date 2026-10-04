# Lung cancer survival-status classification
# Educational portfolio project - not for clinical use.

set.seed(123)

required_packages <- c(
  "caret", "corrplot", "dplyr", "ggplot2", "naniar", "pROC",
  "randomForest", "rpart", "rpart.plot", "smotefamily"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Install the missing packages first: ",
    paste(missing_packages, collapse = ", "),
    ". Run source('packages.R') if you want the helper to install them."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

normalise_names <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

normalise_target <- function(x) {
  value <- tolower(trimws(as.character(x)))
  mapped <- ifelse(
    value %in% c("1", "yes", "true", "survived"), "Yes",
    ifelse(value %in% c("0", "no", "false", "did_not_survive", "not_survived"), "No", NA)
  )
  factor(mapped, levels = c("No", "Yes"))
}

mode_value <- function(x) {
  observed <- x[!is.na(x)]
  if (length(observed) == 0) return(NA)
  names(sort(table(observed), decreasing = TRUE))[1]
}

fit_imputer <- function(frame) {
  list(
    numeric = vapply(frame[vapply(frame, is.numeric, logical(1))], median,
                     numeric(1), na.rm = TRUE),
    categorical = vapply(frame[!vapply(frame, is.numeric, logical(1))], mode_value,
                         character(1))
  )
}

apply_imputer <- function(frame, imputer) {
  for (name in names(imputer$numeric)) {
    frame[[name]][is.na(frame[[name]])] <- imputer$numeric[[name]]
  }
  for (name in names(imputer$categorical)) {
    replacement <- imputer$categorical[[name]]
    frame[[name]][is.na(frame[[name]])] <- replacement
  }
  frame
}

align_factor_levels <- function(train_frame, test_frame) {
  factor_names <- names(train_frame)[vapply(train_frame, is.factor, logical(1))]
  for (name in factor_names) {
    train_levels <- levels(train_frame[[name]])
    test_values <- as.character(test_frame[[name]])
    test_values[!test_values %in% train_levels] <- NA_character_
    test_frame[[name]] <- factor(test_values, levels = train_levels)
  }
  list(train = train_frame, test = test_frame)
}

to_class <- function(probability, threshold = 0.5) {
  factor(ifelse(probability >= threshold, "Yes", "No"), levels = c("No", "Yes"))
}

evaluate_model <- function(name, probability, truth) {
  prediction <- to_class(probability)
  matrix <- caret::confusionMatrix(prediction, truth, positive = "Yes")
  roc_object <- pROC::roc(
    response = truth,
    predictor = as.numeric(probability),
    levels = c("No", "Yes"),
    quiet = TRUE
  )

  list(
    name = name,
    probability = probability,
    confusion = matrix,
    roc = roc_object,
    metrics = data.frame(
      model = name,
      accuracy = unname(matrix$overall[["Accuracy"]]),
      kappa = unname(matrix$overall[["Kappa"]]),
      sensitivity = unname(matrix$byClass[["Sensitivity"]]),
      specificity = unname(matrix$byClass[["Specificity"]]),
      auc = as.numeric(pROC::auc(roc_object)),
      stringsAsFactors = FALSE
    )
  )
}

args <- commandArgs(trailingOnly = TRUE)
data_path <- if (length(args) >= 1) args[[1]] else file.path("data", "raw", "Lung Cancer.csv")

if (!file.exists(data_path)) {
  stop(
    "Dataset not found at '", data_path,
    "'. See data/README.md for the required location and schema."
  )
}

data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
names(data) <- normalise_names(names(data))

data <- dplyr::distinct(data)
data <- data[rowSums(is.na(data)) != ncol(data), , drop = FALSE]
data <- data[, colSums(is.na(data)) < nrow(data), drop = FALSE]

if (!"survived" %in% names(data)) stop("The dataset must contain a 'survived' target column.")
data$survived <- normalise_target(data$survived)
if (anyNA(data$survived)) stop("The 'survived' column contains unsupported or missing values.")
if (length(unique(data$survived)) != 2) stop("The target must contain both No and Yes classes.")

date_columns <- intersect(c("diagnosis_date", "end_treatment_date"), names(data))
if (length(date_columns) == 2) {
  diagnosis_date <- as.Date(data$diagnosis_date)
  end_treatment_date <- as.Date(data$end_treatment_date)
  data$treatment_duration <- as.numeric(difftime(end_treatment_date, diagnosis_date, units = "days"))
}

if ("age" %in% names(data)) {
  data$age_group <- cut(
    as.numeric(data$age), c(0, 30, 50, 70, Inf),
    labels = c("Young", "Middle-aged", "Senior", "Elderly"),
    include.lowest = TRUE
  )
}
if ("bmi" %in% names(data)) {
  data$bmi_category <- cut(
    as.numeric(data$bmi), c(0, 18.5, 24.9, 29.9, Inf),
    labels = c("Underweight", "Normal", "Overweight", "Obese"),
    right = FALSE
  )
}
if ("cholesterol_level" %in% names(data)) {
  data$cholesterol_category <- cut(
    as.numeric(data$cholesterol_level), c(0, 200, 240, Inf),
    labels = c("Desirable", "Borderline High", "High"),
    include.lowest = TRUE
  )
}
if ("treatment_duration" %in% names(data)) {
  data$treatment_group <- cut(
    data$treatment_duration, c(-Inf, 0, 200, 400, 600, Inf),
    labels = c("Invalid", "Short", "Medium", "Long", "Very Long")
  )
}

comorbidities <- intersect(c("asthma", "cirrhosis", "other_cancer"), names(data))
if (length(comorbidities) > 0) {
  comorbidity_frame <- data[comorbidities]
  for (name in names(comorbidity_frame)) {
    if (!is.numeric(comorbidity_frame[[name]])) {
      values <- tolower(as.character(comorbidity_frame[[name]]))
      comorbidity_frame[[name]] <- ifelse(values %in% c("1", "yes", "true"), 1, 0)
    }
  }
  comorbidity_frame[is.na(comorbidity_frame)] <- 0
  data$risk_score <- rowSums(comorbidity_frame)
}

character_columns <- names(data)[vapply(data, is.character, logical(1))]
data[character_columns] <- lapply(data[character_columns], factor)

dir.create("outputs", showWarnings = FALSE)

numeric_for_eda <- data[vapply(data, is.numeric, logical(1))]
numeric_for_eda <- numeric_for_eda[vapply(numeric_for_eda, function(x) sd(x, na.rm = TRUE) > 0, logical(1))]
if (ncol(numeric_for_eda) >= 2) {
  png(file.path("outputs", "correlation_matrix.png"), width = 1100, height = 900)
  corrplot::corrplot(
    cor(numeric_for_eda, use = "pairwise.complete.obs"),
    method = "circle", title = "Correlation Matrix", mar = c(0, 0, 2, 0)
  )
  dev.off()
}

if ("bmi" %in% names(data)) {
  bmi_plot <- ggplot2::ggplot(data, ggplot2::aes(x = survived, y = as.numeric(bmi))) +
    ggplot2::geom_boxplot(fill = "#4F81BD") +
    ggplot2::labs(title = "BMI by Survival Status", x = "Survived", y = "BMI") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path("outputs", "bmi_by_survival.png"), bmi_plot, width = 8, height = 5)
}

# Split before learning imputation, encoding, or scaling parameters.
split_index <- caret::createDataPartition(data$survived, p = 0.8, list = FALSE)
train <- data[split_index, , drop = FALSE]
test <- data[-split_index, , drop = FALSE]

drop_columns <- intersect(c("id", "diagnosis_date", "end_treatment_date"), names(train))
train <- train[, setdiff(names(train), drop_columns), drop = FALSE]
test <- test[, setdiff(names(test), drop_columns), drop = FALSE]

y_train <- factor(train$survived, levels = c("No", "Yes"))
y_test <- factor(test$survived, levels = c("No", "Yes"))
x_train <- train[, setdiff(names(train), "survived"), drop = FALSE]
x_test <- test[, setdiff(names(test), "survived"), drop = FALSE]

aligned <- align_factor_levels(x_train, x_test)
x_train <- aligned$train
x_test <- aligned$test

imputer <- fit_imputer(x_train)
x_train <- apply_imputer(x_train, imputer)
x_test <- apply_imputer(x_test, imputer)
if (anyNA(x_test)) stop("Test data contains values or categories that could not be imputed.")

encoder <- caret::dummyVars(~ ., data = x_train, fullRank = TRUE)
x_train_numeric <- as.data.frame(predict(encoder, newdata = x_train))
x_test_numeric <- as.data.frame(predict(encoder, newdata = x_test))

preprocessor <- caret::preProcess(x_train_numeric, method = c("center", "scale"))
x_train_scaled <- as.data.frame(predict(preprocessor, x_train_numeric))
x_test_scaled <- as.data.frame(predict(preprocessor, x_test_numeric))

# Apply SMOTE only to the training set.
minority_count <- min(table(y_train))
if (minority_count >= 2 && length(unique(y_train)) == 2) {
  neighbours <- min(5, minority_count - 1)
  smote_result <- smotefamily::SMOTE(
    X = x_train_scaled,
    target = y_train,
    K = neighbours,
    dup_size = 0
  )
  smote_frame <- smote_result$data
  smote_target <- smote_frame[[ncol(smote_frame)]]
  smote_predictors <- smote_frame[, -ncol(smote_frame), drop = FALSE]
  names(smote_predictors) <- names(x_train_scaled)
  smote_target <- normalise_target(smote_target)
  if (anyNA(smote_target)) stop("SMOTE returned unexpected class labels.")
  train_model <- cbind(smote_predictors, survived = smote_target)
} else {
  warning("SMOTE skipped because the minority class has fewer than two records.")
  train_model <- cbind(x_train_scaled, survived = y_train)
}

test_model <- cbind(x_test_scaled, survived = y_test)

control <- caret::trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 2,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  savePredictions = "final"
)

set.seed(123)
logistic_model <- caret::train(
  survived ~ ., data = train_model,
  method = "glm", family = binomial(),
  trControl = control, metric = "ROC"
)

set.seed(123)
tree_model <- rpart::rpart(
  survived ~ ., data = train_model, method = "class",
  control = rpart::rpart.control(cp = 0.001, minsplit = 20, xval = 10)
)
best_cp <- tree_model$cptable[which.min(tree_model$cptable[, "xerror"]), "CP"]
tree_model <- rpart::prune(tree_model, cp = best_cp)

png(file.path("outputs", "decision_tree.png"), width = 1200, height = 800)
rpart.plot::rpart.plot(tree_model, main = "Pruned CART Decision Tree")
dev.off()

set.seed(123)
random_forest_model <- randomForest::randomForest(
  survived ~ ., data = train_model,
  ntree = 300,
  mtry = max(1, floor(sqrt(ncol(train_model) - 1))),
  importance = TRUE
)

png(file.path("outputs", "random_forest_importance.png"), width = 1000, height = 800)
randomForest::varImpPlot(random_forest_model, main = "Random Forest - Variable Importance")
dev.off()

probability_logistic <- predict(logistic_model, newdata = test_model, type = "prob")[, "Yes"]
probability_tree <- predict(tree_model, newdata = test_model, type = "prob")[, "Yes"]
probability_forest <- predict(random_forest_model, newdata = test_model, type = "prob")[, "Yes"]

results <- list(
  evaluate_model("Logistic Regression", probability_logistic, y_test),
  evaluate_model("Decision Tree (CART)", probability_tree, y_test),
  evaluate_model("Random Forest", probability_forest, y_test)
)

performance <- dplyr::bind_rows(lapply(results, `[[`, "metrics"))
numeric_metrics <- names(performance)[vapply(performance, is.numeric, logical(1))]
performance[numeric_metrics] <- lapply(performance[numeric_metrics], round, digits = 4)
write.csv(performance, file.path("outputs", "model_performance.csv"), row.names = FALSE)

png(file.path("outputs", "roc_curves.png"), width = 900, height = 700)
plot(results[[1]]$roc, legacy.axes = TRUE, main = "ROC Curves - Held-out Test Set")
plot(results[[2]]$roc, add = TRUE)
plot(results[[3]]$roc, add = TRUE)
legend(
  "bottomright",
  legend = sprintf("%s (AUC = %.3f)", performance$model, performance$auc),
  lty = 1
)
dev.off()

capture.output(sessionInfo(), file = file.path("outputs", "session_info.txt"))
print(performance)

