# Importazione delle librerie utilizzate
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggkegg)
library(tidygraph)
library(ggraph)


#' Plot PCA con simboli per MouseID e colori per trattamento
#'
#' Esegue un'analisi PCA su una matrice di espressione e produce un grafico
#' con punti colorati in base al trattamento e con forme diverse per ciascun MouseID.
#'
#' @param expr_matrix Una matrice numerica di espressione genica (geni x campioni).
#' @param sample_data Data frame contenente metadati dei campioni, incluso il trattamento (`treatment:ch1`).
#' @param map_names Vettore con nomi dei campioni come chiavi e nomi descrittivi come valori.
#' @param title Titolo del grafico.
#'
#' @return Un oggetto ggplot che rappresenta il PCA plot.
#'
#' @import ggplot2
#' @importFrom stringr str_extract
#' @examples
#' # plotPCA_mouse(expr_matrix, sample_data, map_names, "Titolo del PCA")
plotPCA_mouse <- function(expr_matrix, sample_data, map_names, title) {
  # PCA
  pca <- prcomp(t(expr_matrix), scale. = TRUE)
  
  # Percentuale di varianza spiegata
  percent_var <- round(100 * summary(pca)$importance[2, 1:2], 1)
  
  # Mappa inversa
  map_names_inv <- setNames(names(map_names), map_names)
  
  # Costruzione dataframe PCA
  pca_df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    Treatment = factor(sample_data$`treatment:ch1`),
    Sample = rownames(pca$x)
  )
  pca_df$Sample_name <- map_names_inv[pca_df$Sample]
  pca_df$MouseID <- factor(str_extract(pca_df$Sample_name, "^M_\\d+"))
  
  # Plot
  ggplot(pca_df, aes(PC1, PC2, color = Treatment, shape = MouseID)) +
    geom_point(size = 4) +
    theme_minimal() +
    labs(
      title = title,
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)")
    ) +
    scale_color_manual(values = c("CART19-vehicle" = "blue", "CART19-tazemetostat" = "red")) +
    theme(legend.position = "right")
}


#' Crea un vettore dei logFC associati a ID Entrez
#'
#' Estrae i valori di logFC da un data frame di risultati di DE e converte
#' gli ID genici da ENSEMBL a ENTREZ.
#'
#' @param DEGenes Vettore di ENSEMBL ID dei geni differenzialmente espressi.
#' @param top_noGSM_df Data frame contenente le statistiche di espressione differenziale (colonna `logFC`).
#'
#' @return Un vettore numerico di logFC, con nomi corrispondenti a ENTREZ ID.
#'
#' @importFrom clusterProfiler bitr
#' @import org.Hs.eg.db
#' @examples
#' # prepare_logfc(DEGenes, top_noGSM_df)
prepare_logfc <- function(DEGenes, top_noGSM_df) {
  # Estrae i valori di logFC per i geni differenzialmente espressi (DEGenes)
  logfc_map <- top_noGSM_df[DEGenes, "logFC"]
  names(logfc_map) <- DEGenes  # assegna i nomi (ENSEMBL ID)
  
  # Converte gli ID da ENSEMBL a ENTREZ usando il database org.Hs.eg.db
  id_map <- bitr(names(logfc_map),
                 fromType = "ENSEMBL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  
  # Crea un nuovo vettore logFC con nomi ENTREZID
  logfc_entrez <- logfc_map[id_map$ENSEMBL]
  names(logfc_entrez) <- id_map$ENTREZID
  
  # Rimuove eventuali valori NA (geni senza mappatura)
  na.omit(logfc_entrez)
}

#' Cnetplot per termini arricchiti e logFC
#'
#' Disegna una rete (cnetplot) che collega geni a termini GO o pathway arricchiti.
#' I nodi sono colorati in base ai valori di logFC.
#'
#' @param enrich_obj Oggetto di arricchimento (es. da enrichGO o enrichKEGG).
#' @param logfc_vector Vettore numerico con valori di logFC, con nomi ENTREZ ID.
#' @param title Titolo del grafico.
#'
#' @return Un oggetto ggplot con la rete dei geni e termini arricchiti.
#'
#' @import enrichplot
#' @importFrom ggplot2 ggtitle scale_color_gradient2
#' @examples
#' # plot_cnet(enrich_obj, logfc_vector, "Titolo del plot")
plot_cnet <- function(enrich_obj, logfc_vector, title) {
  # Calcola il massimo valore assoluto del logFC per impostare la scala dei colori
  max_abs <- max(abs(logfc_vector))
  
  # plot
  p <- cnetplot(enrich_obj,
                foldChange = logfc_vector, # vettore logFC usato per colorare i nodi
                circular = TRUE,           # rete circolare
                colorEdge = TRUE,
                showCategory = 5,
                node_label = "all") +      # etichette su tutti i nodi
    # Imposta la scala di colori del plot
    scale_color_gradient2(low = "#ADD8E6", mid = "#FFFFFF", high = "#00008B",
                          midpoint = 0, limits = c(-max_abs, max_abs)) +
    # Titolo del grafico
    ggtitle(title)
}


#' Analizza e visualizza un pathway KEGG
#'
#' Scarica e visualizza la rete di un pathway KEGG, evidenziando i geni
#' più connessi secondo centralità.
#'
#' @param pathway_id Stringa con l'ID del pathway KEGG (es. "hsa04620").
#'
#' @return Un oggetto ggplot che rappresenta il grafo KEGG.
#'
#' @import ggkegg
#' @import tidygraph
#' @import ggraph
#' @importFrom stringr strsplit
#' @examples
#' # analyze_pathway("hsa04620")
analyze_pathway <- function(pathway_id) {
  # Scarica la rappresentazione KEGG
  g <- ggkegg(pathway_id, return_tbl_graph = TRUE)
  
  # Calcola centralità e showname
  g <- g |>
    mutate(
      degree = centrality_degree(mode = "all"),
      betweenness = centrality_betweenness(),
      showname = strsplit(graphics_name, ",") |> 
        vapply(`[`, 1, FUN.VALUE = "a")
    )
  
  # Plot del grafo
  g |> 
    ggraph(layout = "manual", x = x, y = y) +
    geom_node_rect(aes(fill = degree, filter = type == "gene")) +
    overlay_raw_map() +
    scale_fill_viridis_c() +
    theme_void()
 
}

#' Volcano plot per analisi di espressione differenziale
#'
#' Crea un volcano plot da una tabella di risultati contenente logFC, PValue e FDR.
#' I geni up-regolati sono in rosso, quelli down-regolati in blu, e i non significativi in grigio.
#' Applica una doppia soglia: FD5 <0.05 e abs(logFC) >=1
#'
#' @param top_table Data frame con i risultati (deve contenere logFC, PValue, FDR).
#' @param titolo Titolo del grafico.
#'
#' @return Un grafico volcano plot.
#' #' @examples
#' # volcano_plot(tof_full,"Volcano plot")
volcano_plot <- function(top_table, titolo) {
  sig <- top_table$FDR <= 0.05 & abs(top_table$logFC) >= 1
  cols <- ifelse(sig, ifelse(top_table$logFC >= 1, "red", "blue"), "grey70")
  plot(top_table$logFC, -log10(top_table$PValue),
       pch = 20, col = cols,
       main = titolo,
       xlab = "log₂FC", ylab = "-log₁₀(p-value)")
  abline(v = c(-1, 1), lty = 2, col = "darkgrey")
  abline(h = -log10(0.05), lty = 2, col = "darkgrey")
  legend("topright", 
         legend = c("Up", "Down", "NotSig"),
         col = c("red", "blue", "grey70"), pch = 20, cex = 0.8)
}





