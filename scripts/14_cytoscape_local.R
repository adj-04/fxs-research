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
#   3. In RStudio: setwd() to the unzipped folder, then source this file.
# ============================================================

library(RCy3)
library(dplyr)

ROOT <- Sys.getenv("FXS_ROOT", unset = getwd())
setwd(ROOT)
dir.create("networks", showWarnings = FALSE)
dir.create(file.path("figures", "ALPACA"), recursive = TRUE, showWarnings = FALSE)

nodes_f <- file.path("results", "ALPACA", "FXS_ALPACA_nodes.csv")
edges_f <- file.path("results", "ALPACA", "FXS_ALPACA_edges.csv")
if (!file.exists(nodes_f) || !file.exists(edges_f))
  stop("Run part2 on Colab first and download results/ALPACA/FXS_ALPACA_{nodes,edges}.csv")

node_table <- read.csv(nodes_f, stringsAsFactors = FALSE) %>%
  mutate(id = as.character(id), label = id, Module = as.character(Module))
edge_table <- read.csv(edges_f, stringsAsFactors = FALSE) %>%
  rename(source = TF, target = Gene) %>%
  filter(source %in% node_table$id, target %in% node_table$id) %>%
  mutate(interaction = "regulates")
cat("Nodes:", nrow(node_table), "| edges:", nrow(edge_table), "\n")

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
