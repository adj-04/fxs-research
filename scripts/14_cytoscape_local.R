# ============================================================
# FXS PIPELINE v5 — STEP 14 (local): Cytoscape visualisation
#
# Colab has no display, so Cytoscape runs on your own machine.
# Colab writes the network tables; this script loads them into a
# running Cytoscape via RCy3, applies the FXS style, lays out, and
# exports the figure, XGMML, the style file and a .cys session.
#
# Before running:
#   1. Download the results bundle from the Colab notebook (cell 12)
#      and unzip it, or sync the Drive folder (fxs-research-run).
#   2. Start Cytoscape (3.10+). Cap its memory in Cytoscape.vmoptions,
#      e.g. -Xmx3g, so it does not compete with R for RAM.
#   3. In RStudio: setwd("~/fxs-research") (v5 results are committed there),
#      or to the unzipped results folder, then source this file.
#
# Nodes in ALPACA modules significant under the label-permutation
# test (part 5, family-wise p < 0.05) get a thick gold border.
# ============================================================

library(RCy3)
library(dplyr)

ROOT <- Sys.getenv("FXS_ROOT", unset = getwd())
setwd(ROOT)
dir.create("networks", showWarnings = FALSE)
dir.create(file.path("figures", "ALPACA"), recursive = TRUE, showWarnings = FALSE)

# v5.1: prefer the corrected ALPACA modules from part 5 (control as reference)
nodes_f <- file.path("results", "robust_nulls", "ALPACA_corrected_cyto_nodes.csv")
edges_f <- file.path("results", "robust_nulls", "ALPACA_corrected_cyto_edges.csv")
if (!file.exists(nodes_f) || !file.exists(edges_f)) {
  cat("Part 5 tables not found — using the v5 ALPACA tables.\n")
  nodes_f <- file.path("results", "ALPACA", "FXS_ALPACA_nodes.csv")
  edges_f <- file.path("results", "ALPACA", "FXS_ALPACA_edges.csv")
}
if (!file.exists(nodes_f) || !file.exists(edges_f))
  stop("Run part2 on Colab first and download results/ALPACA/FXS_ALPACA_{nodes,edges}.csv")

node_table <- read.csv(nodes_f, stringsAsFactors = FALSE) %>%
  mutate(id = as.character(id), label = id, Module = as.character(Module))
edge_all <- read.csv(edges_f, stringsAsFactors = FALSE) %>%
  rename(source = TF, target = Gene) %>%
  filter(source %in% node_table$id, target %in% node_table$id) %>%
  mutate(interaction = "regulates")

# The top modules hold ~300k edges — too many to lay out or read. Keep:
#   * every edge touching a module significant in the label-permutation
#     test (family-wise p < 0.05; part 5), and
#   * the strongest remaining edges, up to MAX_EDGES in total.
MAX_EDGES <- as.numeric(Sys.getenv("FXS_CYTO_MAX_EDGES", 3000))
sig_f <- file.path("results", "robust_nulls", "ALPACA_labelperm_modules.csv")
sig_modules <- if (file.exists(sig_f)) {
  m <- read.csv(sig_f); as.character(m$Module[m$p_global_maxQ < 0.05])
} else character(0)
sig_nodes <- node_table$id[node_table$Module %in% sig_modules]
keep_sig  <- edge_all %>% filter(source %in% sig_nodes | target %in% sig_nodes)
keep_top  <- edge_all %>% anti_join(keep_sig, by = c("source", "target")) %>%
  arrange(desc(Score)) %>% head(max(0, MAX_EDGES - nrow(keep_sig)))
edge_table <- bind_rows(keep_sig, keep_top)
node_table <- node_table %>% filter(id %in% c(edge_table$source, edge_table$target))
node_table$Significant_module <- ifelse(node_table$Module %in% sig_modules, "yes", "no")
cat("Significant module(s):", if (length(sig_modules)) paste(sig_modules, collapse = ", ") else "none",
    "\nEdges kept:", nrow(edge_table), "of", nrow(edge_all),
    "| nodes:", nrow(node_table), "\n")

cytoscapePing()
cat("Cytoscape version:", cytoscapeVersionInfo()[["cytoscapeVersion"]], "\n")

net_suid <- createNetworkFromDataFrames(
  nodes = node_table, edges = edge_table,
  title = "FXS_ALPACA_TopModules", collection = "FXS_PANDA_ALPACA")

style_name <- "FXS_ALPACA_Style"
if (style_name %in% getVisualStyleNames()) deleteVisualStyle(style_name)
createVisualStyle(style_name)

setNodeColorMapping("Type", c("TF", "Target"), c("#E74C3C", "#3498DB"),
                    mapping.type = "d", style.name = style_name)
setNodeShapeMapping("Type", c("TF", "Target"), c("ELLIPSE", "ROUND_RECTANGLE"),
                    style.name = style_name)
setNodeSizeMapping("Score", c(0, 0.1, 0.3, 1.0), c(35, 50, 70, 100),
                   mapping.type = "c", style.name = style_name)
setNodeLabelMapping("label", style.name = style_name)
setEdgeLineWidthMapping("Score",
                        table.column.values = range(edge_table$Score, na.rm = TRUE),
                        widths = c(0.3, 3), mapping.type = "c", style.name = style_name)

setNodeBorderColorDefault("#FFFFFF", style.name = style_name)
setNodeBorderWidthDefault(2, style.name = style_name)
setNodeBorderColorMapping("Significant_module", c("yes", "no"), c("#F1C40F", "#FFFFFF"),
                          mapping.type = "d", style.name = style_name)
setNodeBorderWidthMapping("Significant_module", c("yes", "no"), c(8, 2),
                          mapping.type = "d", style.name = style_name)
setNodeFontSizeDefault(20, style.name = style_name)
setNodeLabelColorDefault("#000000", style.name = style_name)
setEdgeOpacityDefault(60, style.name = style_name)
setEdgeColorDefault("#AAAAAA", style.name = style_name)
setEdgeTargetArrowShapeDefault("DELTA", style.name = style_name)
setBackgroundColorDefault("#FFFFFF", style.name = style_name)
setVisualStyle(style_name)

layoutNetwork("kamada-kawai")
Sys.sleep(3)
fitContent()

# One child network per module, handy for module-by-module figures
for (m in unique(node_table$Module)) {
  clearSelection()
  selectNodes(node_table$id[node_table$Module == m], by.col = "id")
  if (length(getSelectedNodes()) >= 5) {
    createSubnetwork(subnetwork.name = paste0("Module_", m))
    setVisualStyle(style_name); layoutNetwork("kamada-kawai"); fitContent()
    exportImage(file.path("figures", "ALPACA", paste0("FXS_ALPACA_module_", m, ".png")),
                type = "PNG", resolution = 300, zoom = 200, overwriteFile = TRUE)
  }
  setCurrentNetwork(net_suid)
}
clearSelection()

exportImage(file.path("figures", "ALPACA", "FXS_ALPACA_network_final.png"),
            type = "PNG", resolution = 300, zoom = 200, overwriteFile = TRUE)
exportImage(file.path("figures", "ALPACA", "FXS_ALPACA_network_final.pdf"),
            type = "PDF", overwriteFile = TRUE)
exportNetwork(file.path(getwd(), "networks", "FXS_ALPACA_network"), type = "xgmml",
              overwriteFile = TRUE)
exportVisualStyles(file.path(getwd(), "networks", "FXS_ALPACA_Style"),
                   styles = style_name, overwriteFile = TRUE)
saveSession(file.path(getwd(), "networks", "FXS_ALPACA_session"))

cat("\nWritten:\n",
    " figures/ALPACA/FXS_ALPACA_network_final.png / .pdf\n",
    " figures/ALPACA/FXS_ALPACA_module_*.png\n",
    " networks/FXS_ALPACA_network.xgmml\n",
    " networks/FXS_ALPACA_Style.xml\n",
    " networks/FXS_ALPACA_session.cys\n")
