##09182024
##wllab
##SHAN Jiang
##NHF1 Targeted Damage-seq
##Dipyrimidine counting from -2 to -1


# Load necessary libraries
library(ggplot2)
library(scales)

# File paths for input data
file_path_CPD1 <- "CPD1.dinucleotide_counts.txt"
file_path_CPD2 <- "CPD2.dinucleotide_counts.txt"
file_path_CIS1 <- "CIS1.dinucleotide_counts.txt"
file_path_CIS2 <- "CIS2.dinucleotide_counts.txt"
file_path_CTR1 <- "CTR1.dinucleotide_counts.txt"
file_path_CTR2 <- "CTR2.dinucleotide_counts.txt"

# Read datasets
data_CPD1 <- read.table(file_path_CPD1, header = TRUE, sep = "\t")
data_CPD2 <- read.table(file_path_CPD2, header = TRUE, sep = "\t")
data_CIS1 <- read.table(file_path_CIS1, header = TRUE, sep = "\t")
data_CIS2 <- read.table(file_path_CIS2, header = TRUE, sep = "\t")
data_CTR1 <- read.table(file_path_CTR1, header = TRUE, sep = "\t")
data_CTR2 <- read.table(file_path_CTR2, header = TRUE, sep = "\t")

# Summing dinucleotide counts
sum_CPD <- data_CPD1$dinucleotide_counts + data_CPD2$dinucleotide_counts
sum_CIS <- data_CIS1$dinucleotide_counts + data_CIS2$dinucleotide_counts
sum_CTR <- data_CTR1$dinucleotide_counts + data_CTR2$dinucleotide_counts

# Create a combined dataframe
row_names <- data_CPD1$X  # Assuming the dinucleotides in column 'X' are the same across datasets
combined_data <- data.frame(
  row_names = rep(row_names, 3),
  Position_Sum = c(sum_CPD, sum_CIS, sum_CTR),
  Condition = rep(c("CPD", "CIS", "CTR"), each = length(row_names))
)

# Define the custom order of dinucleotides
custom_order <- c("TT", "TC", "CT", "CC", 
                  "TA", "CA", 
                  "AA", "AC", "AT", "CG", "GC",
                  "GT", "TG", "GA", "AG", "GG")

# Set row_names as a factor with the custom order
combined_data$row_names <- factor(combined_data$row_names, levels = custom_order)

# Ensure Condition factor levels are in the correct order
combined_data$Condition <- factor(combined_data$Condition, levels = c("CPD", "CIS", "CTR"))

# Create the bar plot
ggplot(combined_data, aes(x = row_names, y = Position_Sum, fill = Condition)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge()) + 
  labs(
    title = "Comparison of Dinucleotides Between CPD, CIS, and CTR", 
    x = "Dinucleotides", 
    y = "Summed Frequency"
  ) + 
  scale_fill_manual(
    values = c("CPD" = "#F28C9D", "CIS" = "#F3E065", "CTR" = "#9EAAE0")
  ) + 
  scale_y_continuous(labels = scales::comma) + 
  theme_classic(base_size = 14) +  # White background
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered bold title
    axis.text.x = element_text(angle = 45, hjust = 1),                 # Rotate x-axis labels
    axis.text.y = element_text(size = 10),                             # Adjust y-axis text size
    legend.title = element_blank(),                                    # Remove legend title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add black outer box
  )

