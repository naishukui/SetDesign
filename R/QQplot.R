#this code generates QQ plot for  STAAR pipeline results under different configurations
library(ggplot2)
library(dplyr)

# Define the qqvals function to calculate expected and observed values
qqvals <- function(observedPValues) {
  x <- -log10(1:length(observedPValues) / length(observedPValues))
  y <- -log10(sort(observedPValues))
  return(data.frame(Expected = x, Observed = y))
}
# Load datasets into a list
datasets <- list(
  setting1 = get(load("output_file1.R")),
  setting2 = get(load("output_file2.R")),
  setting3 = get(load("output_file3.R")),
  setting4 = get(load("output_file4.R")),
  setting5 = get(load("output_file5.R"))
)

# Initialize an empty data frame to store QQ plot values for all datasets
all_qq_data <- data.frame()

# Loop through each dataset to calculate QQ values and add to all_qq_data
for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  STAAR_P <- dataset[, 6]
  qq_data <- qqvals(STAAR_P)
  qq_data$Dataset <- dataset_name  # Add a column for the dataset name
  all_qq_data <- rbind(all_qq_data, qq_data)
}

# Plot all datasets on the same QQ plot
ggplot(all_qq_data, aes(x = Expected, y = Observed, color = Dataset)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_classic() +
  xlab("Expected -log10(p)") +
  ylab("Observed -log10(p)") +
  ggtitle("QQ Plot Across Five STAAR Pipeline Configurations") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
