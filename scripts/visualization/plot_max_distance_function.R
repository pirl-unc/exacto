library(ggplot2)
library(scales)


OUTPUT.PDF <- "/Users/leework/Documents/Research/projects/project_exacto/exacto/misc/exacto_breakpoint_max_distance_plot.pdf"
max.distance <- 1000
tau <- 2000
variant.size <- seq(0, 20000, by = 1)
y <- max.distance * (1 - exp(-variant.size / tau))
data <- data.frame(variant.size = variant.size, distance = y)

x.values <- c(10, 100, 1000, 2000, 4000, 6000, 8000, 10000)
y.values <- max.distance * (1 - exp(-x.values / tau))

annotations <- data.frame(x = x.values, y = y.values, label = sprintf("(%d, %.0f)", x.values, y.values))

p <- ggplot(data, aes(x = variant.size, y = distance)) +
  geom_line(color = "blue", size = 1) +
  labs(
    x = "Variant Size",
    y = "Breakpoint Max Distance",
    title = "Breakpoint Max Distance vs Variant Size"
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  geom_hline(yintercept = max.distance, linetype = "dashed", color = "red", size = 1) +
  geom_point(data = annotations, aes(x = x, y = y), color = "blue", size = 3) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), color = "black", hjust = -0.1, size = 6) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 24)
  )

print(p)

ggsave(plot = p, filename = OUTPUT.PDF, width = 16, height = 9, dpi = 300)
