



##########################################################

library(dplyr)
library(ggplot2)
library(reshape2)

# Read the two CSV files
CPD1_df <- read.csv("CPD1_127.csv", stringsAsFactors = FALSE, header = FALSE)
CPD2_df <- read.csv("CPD2_127.csv", stringsAsFactors = FALSE, header = FALSE)
CIS1_df <- read.csv("CIS1_127.csv", stringsAsFactors = FALSE, header = FALSE)
CIS2_df <- read.csv("CIS2_127.csv", stringsAsFactors = FALSE, header = FALSE)
CTR1_df <- read.csv("CTR1_127.csv", stringsAsFactors = FALSE, header = FALSE)
CTR2_df <- read.csv("CTR2_127.csv", stringsAsFactors = FALSE, header = FALSE)

# Rename columns
colnames(CPD1_df) <- c("Gene", "Count")
colnames(CPD2_df) <- c("Gene", "Count")
colnames(CIS1_df) <- c("Gene", "Count")
colnames(CIS2_df) <- c("Gene", "Count")
colnames(CTR1_df) <- c("Gene", "Count")
colnames(CTR2_df) <- c("Gene", "Count")

# Sum duplicated genes
CPD1_df <- CPD1_df %>% group_by(Gene) %>% summarise(Count = sum(Count), .groups = "drop")
CPD2_df <- CPD2_df %>% group_by(Gene) %>% summarise(Count = sum(Count), .groups = "drop")
CIS1_df <- CIS1_df %>% group_by(Gene) %>% summarise(Count = sum(Count), .groups = "drop")
CIS2_df <- CIS2_df %>% group_by(Gene) %>% summarise(Count = sum(Count), .groups = "drop")
CTR1_df <- CTR1_df %>% group_by(Gene) %>% summarise(Count = sum(Count), .groups = "drop")
CTR2_df <- CTR2_df %>% group_by(Gene) %>% summarise(Count = sum(Count), .groups = "drop")

# Merge
merged_CPD <- merge(CPD1_df, CPD2_df, by = "Gene")
merged_CIS <- merge(CIS1_df, CIS2_df, by = "Gene")
merged_CTR <- merge(CTR1_df, CTR2_df, by = "Gene")

colnames(merged_CPD)[2:3] <- c("CPD1", "CPD2")
colnames(merged_CIS)[2:3] <- c("CIS1", "CIS2")
colnames(merged_CTR)[2:3] <- c("CTR1", "CTR2")

# -----------------------------------
# Pearson correlation + p-value
# -----------------------------------
cpd_test <- cor.test(merged_CPD$CPD1, merged_CPD$CPD2, method = "pearson")
cis_test <- cor.test(merged_CIS$CIS1, merged_CIS$CIS2, method = "pearson")
ctr_test <- cor.test(merged_CTR$CTR1, merged_CTR$CTR2, method = "pearson")

# Print results
cat("CPD: r =", cpd_test$estimate, " p =", cpd_test$p.value, "\n")
cat("CIS: r =", cis_test$estimate, " p =", cis_test$p.value, "\n")
cat("CTR: r =", ctr_test$estimate, " p =", ctr_test$p.value, "\n")
#CPD: r = 0.9927621  p = 6.274549e-117 
#CIS: r = 0.99547  p = 1.298404e-129 
#CTR: r = 0.9997705  p = 1.643885e-210 


# Make a small summary table
cor_summary <- data.frame(
  Sample = c("CPD", "CIS", "CTR"),
  Pearson_r = c(unname(cpd_test$estimate), unname(cis_test$estimate), unname(ctr_test$estimate)),
  P_value = c(cpd_test$p.value, cis_test$p.value, ctr_test$p.value)
)

write.csv(cor_summary, "Pearson_correlation_summary.csv", row.names = FALSE)

# -----------------------------------
# Correlation matrices for heatmap
# -----------------------------------
cormat_CPD <- matrix(c(1, unname(cpd_test$estimate),
                       unname(cpd_test$estimate), 1),
                     nrow = 2, byrow = TRUE)
rownames(cormat_CPD) <- colnames(cormat_CPD) <- c("CPD1", "CPD2")

cormat_CIS <- matrix(c(1, unname(cis_test$estimate),
                       unname(cis_test$estimate), 1),
                     nrow = 2, byrow = TRUE)
rownames(cormat_CIS) <- colnames(cormat_CIS) <- c("CIS1", "CIS2")

cormat_CTR <- matrix(c(1, unname(ctr_test$estimate),
                       unname(ctr_test$estimate), 1),
                     nrow = 2, byrow = TRUE)
rownames(cormat_CTR) <- colnames(cormat_CTR) <- c("CTR1", "CTR2")

# p-value matrices
pmat_CPD <- matrix(c(NA, cpd_test$p.value,
                     cpd_test$p.value, NA),
                   nrow = 2, byrow = TRUE)
rownames(pmat_CPD) <- colnames(pmat_CPD) <- c("CPD1", "CPD2")

pmat_CIS <- matrix(c(NA, cis_test$p.value,
                     cis_test$p.value, NA),
                   nrow = 2, byrow = TRUE)
rownames(pmat_CIS) <- colnames(pmat_CIS) <- c("CIS1", "CIS2")

pmat_CTR <- matrix(c(NA, ctr_test$p.value,
                     ctr_test$p.value, NA),
                   nrow = 2, byrow = TRUE)
rownames(pmat_CTR) <- colnames(pmat_CTR) <- c("CTR1", "CTR2")

# Keep upper triangle
get_upper_tri <- function(mat) {
  mat[lower.tri(mat)] <- NA
  mat
}

upper_tri_CPD <- get_upper_tri(cormat_CPD)
upper_tri_CIS <- get_upper_tri(cormat_CIS)
upper_tri_CTR <- get_upper_tri(cormat_CTR)

upper_p_CPD <- get_upper_tri(pmat_CPD)
upper_p_CIS <- get_upper_tri(pmat_CIS)
upper_p_CTR <- get_upper_tri(pmat_CTR)

# Melt
melted_cormat_CPD <- melt(upper_tri_CPD, na.rm = TRUE)
melted_cormat_CIS <- melt(upper_tri_CIS, na.rm = TRUE)
melted_cormat_CTR <- melt(upper_tri_CTR, na.rm = TRUE)

melted_p_CPD <- melt(upper_p_CPD, na.rm = TRUE)
melted_p_CIS <- melt(upper_p_CIS, na.rm = TRUE)
melted_p_CTR <- melt(upper_p_CTR, na.rm = TRUE)

# Add p-values into heatmap table
melted_cormat_CPD$p_value <- melted_p_CPD$value
melted_cormat_CIS$p_value <- melted_p_CIS$value
melted_cormat_CTR$p_value <- melted_p_CTR$value

# Label with r and p
melted_cormat_CPD$label <- paste0("r = ", round(melted_cormat_CPD$value, 5),
                                  "\nP = ", signif(melted_cormat_CPD$p_value, 3))
melted_cormat_CIS$label <- paste0("r = ", round(melted_cormat_CIS$value, 5),
                                  "\nP = ", signif(melted_cormat_CIS$p_value, 3))
melted_cormat_CTR$label <- paste0("r = ", round(melted_cormat_CTR$value, 5),
                                  "\nP = ", signif(melted_cormat_CTR$p_value, 3))

# -----------------------------------
# Heatmaps
# -----------------------------------
plot_heatmap <- function(df, title_text) {
  ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), color = "black", size = 4) +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0.5, limit = c(0, 1),
      name = "Pearson\nCorrelation"
    ) +
    labs(title = title_text) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    coord_fixed() +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
}

plot_heatmap(melted_cormat_CPD, "CPD Correlation")
plot_heatmap(melted_cormat_CIS, "CIS Correlation")
plot_heatmap(melted_cormat_CTR, "CTR Correlation")

# -----------------------------------
# Scatter plots with r and p
# -----------------------------------
ggplot(merged_CPD, aes(x = CPD1, y = CPD2)) +
  geom_point(color = "#9EAAE0", size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#F3E065", linewidth = 1) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("r = ", round(cpd_test$estimate, 5),
                          "\nP = ", signif(cpd_test$p.value, 3)),
           hjust = 1.1, vjust = 1.5, size = 5) +
  labs(title = "CPD Correlation", x = "CPD1/Gene", y = "CPD2/Gene") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggplot(merged_CIS, aes(x = CIS1, y = CIS2)) +
  geom_point(color = "#9EAAE0", size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#F3E065", linewidth = 1) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("r = ", round(cis_test$estimate, 5),
                          "\nP = ", signif(cis_test$p.value, 3)),
           hjust = 1.1, vjust = 1.5, size = 5) +
  labs(title = "CIS Correlation", x = "CIS1/Gene", y = "CIS2/Gene") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggplot(merged_CTR, aes(x = CTR1, y = CTR2)) +
  geom_point(color = "#9EAAE0", size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#F3E065", linewidth = 1) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("r = ", round(ctr_test$estimate, 5),
                          "\nP = ", signif(ctr_test$p.value, 3)),
           hjust = 1.1, vjust = 1.5, size = 5) +
  labs(title = "CTR Correlation", x = "CTR1/Gene", y = "CTR2/Gene") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


























