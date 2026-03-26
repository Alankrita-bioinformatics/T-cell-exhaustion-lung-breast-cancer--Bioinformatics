# Load and install the specific library needed for the Sina plot
library(ggforce)
library(ggplot2)

# Ensure Subtype is a factor for consistent coloring
tex_df2$Subtype <- factor(tex_df2$Subtype, levels = c("Normal", "Tumor"))

# Generate the plot
sina_plot <- ggplot(tex_df2, aes(x = Subtype, y = TEX_Score, color = Subtype)) +
  # Add the boxplot layer
  geom_boxplot(width = 0.3, outlier.shape = NA, fatten = 2) +
  # Add the Sina layer (points jittered based on density)
  geom_sina(alpha = 0.6, size = 2.5) +
  # Match the colors from your image (Blue for Normal, Red for Tumor)
  scale_color_manual(values = c("Normal" = "#3182bd", "Tumor" = "#de2d26")) +
  # Apply a clean theme
  theme_bw() +
  labs(
    title = "T Cell Exhaustion Scores: Tumor vs Normal",
    x = "Sample Group",
    y = "TEX Score",
    color = "Group"
  ) +
  # Refine the visual theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_line(color = "#f0f0f0"),
    panel.grid.major = element_line(color = "#f0f0f0")
  )

# Save the plot
ggsave("TEX_Sina_Plot_Final.png", plot = sina_plot, width = 10, height = 8, dpi = 600)

# Display the plot
print(sina_plot)
