#' ---
#' title: "CHANGEME Seidr threshold plot analysis"
#' author: "CHANGEME and Elena van Zalen"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' 
# Load the necessary library
library(ggplot2)
library(tidyr)

# Replace 'your_data.tsv' with the actual file path to your TSV file
data <- read.table("/mnt/ada/projects/spruce/nstreet/conifer-networks/Sara/noRootsMetabolome/thresholded2.tsv", 
                   header = FALSE, sep = "\t",)

# Assign column names manually
colnames(data) <- c("Threshold", "Number of nodes", "Number of edges", "Scale free fit (R^2)", "Average clustering coefficient")

# Reshape data into long format for ggplot
data_long <- pivot_longer(data, cols = -Threshold, names_to = "Variable", values_to = "Value")

# Create a ggplot object to plot each column against the 'Threshold' column
# and facet it into a 2x2 grid
pdf("~/Git/ConiferNetworks/results/FullMetabolome_threshold.pdf", width = 8, height = 6)
plot <- ggplot(data_long, aes(x = Threshold, y = Value, color = Variable)) +
  geom_line(size = 1) +
  labs(
    x = "Threshold",
    y = "Values",
    title = "Seidr threshold values for Metabolome no root data"
  ) +
  theme_minimal() +
  facet_wrap(~Variable, scales = "free") +
  theme(legend.position = "none")
dev.off()

# Set individual axis limits for specific facets
zoomed_plot <- plot +
  facet_wrap(~Variable, scales = "free") +
  scale_y_continuous(
    limits = ifelse(data_long$Variable == "Number of edges", c(0, 1000000), NULL)
  ) +
  scale_x_continuous(
    limits = ifelse(data_long$Variable == "Number of edges", c(0.1, 0.9), NULL)
  )

# Display the zoomed plot
print(zoomed_plot)


