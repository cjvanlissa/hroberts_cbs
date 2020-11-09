library(motley)
library(ggplot2)
library(ggforce)

layout <- get_layout(read.csv("layout.csv", header = FALSE, stringsAsFactors = FALSE))

df_nodes <-  data.frame(node_id = 1:length(layout$param), param = layout$param, stringsAsFactors = FALSE)
df_nodes$shape <- "oval"
df_nodes$shape[grepl("(phys)", df_nodes$param)] <- "rect"
labels <- list("ne" = "Natural environment",
               "plea" = "Pleasantness",
               "dist" = "Disturbance",
               "saf"  = "Safety",
               "coh" = "Cohesion",
               "phys" = "Physical activity",
               "stress" = "Stress",
               "dep" ="Depression"
)
df_nodes$label <- unlist(labels[match(df_nodes$param, names(labels))])
#df_nodes$param <- df_nodes$node_id
df_edges <- data.frame(matrix(c(
  1, 6, "last", "+", 
  2, 6, "last", "+", 
  3, 6, "last", "-", 
  4, 6, "last", "+", 
  5, 6, "last", "+", 
  
  1, 7, "last", "-", 
  2, 7, "last", "-", 
  3, 7, "last", "+", 
  4, 7, "last", "-", 
  5, 7, "last", "-", 
  6, 8, "last", "-", 
  7, 8, "last", "+"), ncol = 4, byrow = TRUE), stringsAsFactors = FALSE)
names(df_edges) <- c("from", "to", "arrow", "label")

df_edges$from <- df_nodes$param[as.numeric(df_edges$from)]
df_edges$to <- df_nodes$param[as.numeric(df_edges$to)]
#df_edges$label[grepl("^ri", df_edges$from)] <- "1"
prep <- prepare_plot_sem(nodes = df_nodes, layout = layout, edges = df_edges)
prep$edges$connect_from <- "right"
prep$edges$connect_to[c(5,6)] <- c("bottom", "top")
prep$edges$connect_to <- c(rep("left", 10), "top", "bottom")
p <- plot(prep)

ggsave("Figure1.png", p, width = 10, height = 3)
