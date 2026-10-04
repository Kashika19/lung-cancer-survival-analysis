required_packages <- c(
  "caret", "corrplot", "dplyr", "ggplot2", "naniar", "pROC",
  "randomForest", "rpart", "rpart.plot", "smotefamily"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) == 0) {
  message("All required packages are already installed.")
} else {
  install.packages(missing_packages, dependencies = TRUE)
}

