if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)



library(shiny)
library(plotly)
library(readxl)
library(dplyr)
library(DT)
library(BiocManager) 
library(limma)


# =====================================================
# LOAD DATA
# =====================================================

getwd()

df <- read_excel(
  "mer_proteins.xlsx",
  sheet = 4
)

sample_info <- df[,1:2]

protein_data0 <- df[,3:ncol(df)]



col2remove <- which(apply(protein_data0, 2, sd) == 0)
View(col2remove)

# remove cols with no variation

protein_data <- protein_data0[, -col2remove]
View(protein_data)


# log transform
protein_data <- log1p(protein_data)
View(protein_data)


sample_info$treatment <- factor(
  sample_info$treatment,
  levels = c(0,1),
  labels = c("Control","Treated")
)

# =====================================================
# PCA
# =====================================================

pca <- prcomp(
  protein_data,
  center = TRUE,
  scale. = TRUE
)

scores <- as.data.frame(pca$x)

scores$sample_id <- sample_info$sample_nr

scores$treatment <- sample_info$treatment

# =====================================================
# PCA LOADINGS
# =====================================================

loadings <- as.data.frame(
  pca$rotation
)

loadings$protein <- rownames(loadings)

rownames(loadings) <- NULL

# =====================================================
# LIMMA DIFFERENTIAL ANALYSIS
# =====================================================

design <- model.matrix(
  ~ treatment,
  data = sample_info
)

protein_matrix <- t(
  as.matrix(protein_data)
)

fit <- lmFit(
  protein_matrix,
  design
)

fit <- eBayes(fit)

limma_results <- topTable(
  fit,
  coef = 2,
  number = Inf,
  adjust.method = "BH"
)

limma_results$protein <- rownames(
  limma_results
)

limma_results$Significance <- case_when(
  limma_results$adj.P.Val < 0.001 ~ "***",
  limma_results$adj.P.Val < 0.01  ~ "**",
  limma_results$adj.P.Val < 0.05  ~ "*",
  TRUE ~ "ns"
)


limma_results <- limma_results %>%
  rename(
    `P-value` = P.Value,
    `FDR` = adj.P.Val
  )

rownames(limma_results) <- NULL

# =====================================================
# MERGE PCA LOADINGS + LIMMA RESULTS
# =====================================================

combined_results <- loadings %>%
  left_join(
    limma_results,
    by = "protein"
  )

# =====================================================
# UI
# =====================================================

ui <- fluidPage(
  
  titlePanel(
    "Proteomics PCA Explorer"
  ),
  
  wellPanel(
    
    h4("About this application"),
    
    p(
      "Interactive exploration of principal component analysis (PCA) and treatment-associated protein abundance differences."
    ),
    
    p(
      "The PCA plot visualises sample clustering while the accompanying tables identify proteins contributing to principal components and proteins associated with treatment effects."
    ),
    
    tags$ul(
      tags$li("Loading: contribution of a protein to the selected principal component."),
      tags$li("logFC: log2 fold change between treated and control samples."),
      tags$li("FDR: Benjamini-Hochberg adjusted p-value."),
      tags$li("Significance: ns = not significant, * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001.")
    ) 
    
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "xpc",
        "X-axis Principal Component",
        choices = paste0("PC",1:5),
        selected = "PC1"
      ),
      
      selectInput(
        "ypc",
        "Y-axis Principal Component",
        choices = paste0("PC",1:5),
        selected = "PC2"
      ),
      
      sliderInput(
        "nproteins",
        "Proteins displayed",
        min = 5,
        max = 50,
        value = 15
      )
      
    ),
    
    mainPanel(
      
      plotlyOutput(
        "pca_plot",
        height = "600px"
      ),
      
      br(),
      
      tabsetPanel(
        
        tabPanel(
          "Proteins Driving PCA",
          DTOutput("loading_table")
        ),
        
        tabPanel(
          "Treatment Effect (limma)",
          DTOutput("limma_table")
        )
        
      )
      
    )
    
  )
  
)

# =====================================================
# SERVER
# =====================================================

server <- function(input, output){
  
  output$pca_plot <- renderPlotly({
    
    plot_ly(
      data = scores,
      x = ~get(input$xpc),
      y = ~get(input$ypc),
      color = ~treatment,
      type = "scatter",
      mode = "markers+text",
      text = ~paste(
        "Sample ID:", sample_id,
        "<br>Treatment:", treatment
      ),
      hoverinfo = "text"
    )
    
  })
  
  # ---------------------------------------------
  # PCA DRIVER PROTEINS
  # ---------------------------------------------
  
  output$loading_table <- renderDT({
    
    pc <- input$xpc
    
    combined_results %>%
      arrange(
        desc(
          abs(.data[[pc]])
        )
      ) %>%
      select(
        Protein = protein,
        Loading = all_of(pc),
        logFC,
        FDR,
        Significance
      ) %>%
      head(input$nproteins)
    
  },
  rownames = FALSE)
  
  # ---------------------------------------------
  # LIMMA RESULTS
  # ---------------------------------------------
  
  output$limma_table <- renderDT({
    
    limma_results %>%
      select(
        Protein = protein,
        logFC,
        `P-value`,
        FDR,
        Significance
      ) %>%
      head(input$nproteins)
    
  },
  rownames = FALSE)

  
}



# =====================================================
# RUN APP
# =====================================================


shinyApp(ui, server)

