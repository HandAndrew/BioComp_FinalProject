# SNP PCA Analyzer App version 8.2 
# Andrew Hand 
# 5/01/2027

# version 8.2 of the shiny app
# conversion of FASTA files and STRUCTURE files
# SNP matrix generation
# PCA computation and plotting 
# Download buttons for both PCA results and SNP matrix 
# SNP matrix and PCA preview in the output


# Increase file upload limit (up to 400MB)
options(shiny.maxRequestSize = 400 * 1024^2)

# Load required libraries
library(shiny)
library(Biostrings)
library(ggplot2)

# ---------- Converting Functions ----------
# Converts FASTA sequences into a numeric SNP matrix
create_snp_matrix <- function(fasta) {
  seqs <- as.character(fasta)  # Convert DNAStringSet to character
  seq_matrix <- do.call(rbind, strsplit(seqs, split = ""))  # Split each sequence into bases
  
  # Keep only columns with variation (polymorphic sites)
  snp_cols <- apply(seq_matrix, 2, function(col) length(unique(col)) > 1)
  snp_matrix <- seq_matrix[, snp_cols, drop = FALSE]
  
  # Convert bases to numeric values (A=1, C=2, G=3, T=4, N=5; NA gets set to 0)
  base_to_num <- function(base) match(base, c("A", "C", "G", "T", "N"))
  snp_numeric <- apply(snp_matrix, c(1, 2), base_to_num)
  snp_numeric[is.na(snp_numeric)] <- 0
  
  return(snp_numeric)
}

# ---------- User Interface ----------
ui <- fluidPage(
  titlePanel("SNP PCA Analyzer"),
  
  tabsetPanel(
    # Tab for analyzing FASTA files
    tabPanel("FASTA File",
             sidebarLayout(
               sidebarPanel(
                 fileInput("fasta_file", "Upload .fasta File", accept = c(".fasta", ".fa")),
                 helpText("Upload aligned DNA sequences in FASTA format."),
                 downloadButton("download_pca", "Download PCA Results"),
                 downloadButton("download_snp", "Download SNP Matrix")
               ),
               mainPanel(
                 plotOutput("pca_plot"),  # PCA plot output
                 verbatimTextOutput("snp_summary")  # SNP matrix preview and summary
               )
             )
    ),
    
    # Tab for analyzing STRU files
    tabPanel("STRU File",
             sidebarLayout(
               sidebarPanel(
                 fileInput("stru_file", "Upload .stru File", accept = c(".stru", ".txt")),
                 helpText("Upload a STRUCTURE format SNP file (numeric)."),
                 downloadButton("download_pca_stru", "Download PCA Results"),
                 downloadButton("download_snp_stru", "Download SNP Matrix")
               ),
               mainPanel(
                 plotOutput("pca_plot_stru"),  # PCA plot for STRU
                 verbatimTextOutput("snp_summary_stru")  # STRU matrix preview
               )
             )
    )
  )
)

# ---------- Server ----------
server <- function(input, output) {
  
  # Read uploaded FASTA file reactively
  fasta_data <- reactive({
    req(input$fasta_file)
    readDNAStringSet(input$fasta_file$datapath)
  })
  
  # Create numeric SNP matrix from FASTA
  snp_matrix <- reactive({
    fasta <- fasta_data()
    withProgress(message = "Processing FASTA file...", value = 0.5, {
      create_snp_matrix(fasta)
    })
  })
  
  # Display dimensions and preview of SNP matrix
  output$snp_summary <- renderPrint({
    matrix <- snp_matrix()
    cat("Original SNP matrix dimensions:\n")
    print(dim(matrix))
    
    # Filter out invariant or NA-containing columns
    matrix <- matrix[, apply(matrix, 2, function(x) !any(is.na(x)) && var(x) > 0), drop = FALSE]
    
    cat("\nAfter filtering:\n")
    print(dim(matrix))
    
    # Preview first few rows/columns
    if (nrow(matrix) == 0 || ncol(matrix) == 0) {
      cat("\nNo usable SNP data remains after filtering.\n")
    } else {
      cat("\nPreview:\n")
      print(matrix[1:min(5, nrow(matrix)), 1:min(5, ncol(matrix)), drop = FALSE])
    }
  })
  
  # Render PCA plot from filtered SNP matrix
  output$pca_plot <- renderPlot({
    matrix <- snp_matrix()
    matrix <- matrix[, apply(matrix, 2, function(x) !any(is.na(x)) && var(x) > 0), drop = FALSE]
    
    # Ensure enough data for PCA
    validate(
      need(ncol(matrix) >= 2, "Too few valid SNPs for PCA."),
      need(nrow(matrix) >= 2, "Too few sequences for PCA.")
    )
    
    # Perform PCA plot
    pca <- prcomp(matrix, scale. = TRUE)
    df <- as.data.frame(pca$x[, 1:2])
    df$sample <- rownames(df)
    
    ggplot(df, aes(x = PC1, y = PC2)) +
      geom_point(color = "blue", size = 3) +
      labs(title = "PCA of SNP Matrix (FASTA)", x = "PC1", y = "PC2") +
      theme_minimal()
  })
  
  # Download PCA results as CSV file
  output$download_pca <- downloadHandler(
    filename = function() paste0("PCA_results_FASTA_", Sys.Date(), ".csv"),
    content = function(file) {
      matrix <- snp_matrix()
      matrix <- matrix[, apply(matrix, 2, function(x) !any(is.na(x)) && var(x) > 0), drop = FALSE]
      pca <- prcomp(matrix, scale. = TRUE)
      write.csv(pca$x, file)
    }
  )
  
  # Download SNP matrix as CSV
  output$download_snp <- downloadHandler(
    filename = function() paste0("SNP_matrix_FASTA_", Sys.Date(), ".csv"),
    content = function(file) {
      matrix <- snp_matrix()
      rownames(matrix) <- names(fasta_data())  # Set row names to sequence names
      write.csv(matrix, file)
    }
  )
  
  # ---------- STRUCTURE File ----------
  # Process uploaded STRU file and generate numeric SNP matrix
  stru_matrix <- reactive({
    req(input$stru_file)
    
    withProgress(message = "Processing STRU file...", value = 0.5, {
      raw <- read.table(
        input$stru_file$datapath,
        header = FALSE,
        stringsAsFactors = FALSE,
        fill = TRUE,
        strip.white = TRUE,
        sep = ""
      )
      
      # Check for even number of rows (2 per individual)
      if (nrow(raw) %% 2 != 0) {
        stop("STRU file must have an even number of rows (2 per individual).")
      }
      
      meta_cols <- 2  # First 2 columns are metadata
      data_only <- raw[, -(1:meta_cols)]  # Exclude metadata
      
      # Combine pairs of rows (diploid genotypes)
      combined <- matrix(NA, nrow = nrow(data_only) / 2, ncol = ncol(data_only))
      for (i in seq(1, nrow(data_only), by = 2)) {
        row1 <- data_only[i, ]
        row2 <- data_only[i + 1, ]
        combined[(i + 1) / 2, ] <- paste0(row1, row2)
      }
      
      snp_numeric <- apply(combined, 2, as.numeric)  # Convert to numeric
      return(snp_numeric)
    })
  })
  
  # Display summary of STRU matrix
  output$snp_summary_stru <- renderPrint({
    matrix <- stru_matrix()
    cat("Original STRU matrix dimensions:\n")
    print(dim(matrix))
    
    matrix <- matrix[, apply(matrix, 2, function(x) !any(is.na(x)) && var(x) > 0), drop = FALSE]
    
    cat("\nAfter filtering:\n")
    print(dim(matrix))
    
    if (nrow(matrix) == 0 || ncol(matrix) == 0) {
      cat("\nNo usable SNP data remains after filtering.\n")
    } else {
      cat("\nPreview:\n")
      print(matrix[1:min(5, nrow(matrix)), 1:min(5, ncol(matrix)), drop = FALSE])
    }
  })
  
  # Render PCA plot from STRU data
  output$pca_plot_stru <- renderPlot({
    matrix <- stru_matrix()
    matrix <- matrix[, apply(matrix, 2, function(x) !any(is.na(x)) && var(x) > 0), drop = FALSE]
    
    validate(
      need(ncol(matrix) >= 2, "Too few valid SNPs for PCA."),
      need(nrow(matrix) >= 2, "Too few individuals for PCA.")
    )
    
    pca <- prcomp(matrix, scale. = TRUE)
    df <- as.data.frame(pca$x[, 1:2])
    df$sample <- rownames(df)
    
    ggplot(df, aes(x = PC1, y = PC2)) +
      geom_point(color = "darkgreen", size = 3) +
      labs(title = "PCA of SNP Matrix (.stru)", x = "PC1", y = "PC2") +
      theme_minimal()
  })
  
  # Download PCA results from STRU file
  output$download_pca_stru <- downloadHandler(
    filename = function() paste0("PCA_results_STRU_", Sys.Date(), ".csv"),
    content = function(file) {
      matrix <- stru_matrix()
      matrix <- matrix[, apply(matrix, 2, function(x) !any(is.na(x)) && var(x) > 0), drop = FALSE]
      pca <- prcomp(matrix, scale. = TRUE)
      write.csv(pca$x, file)
    }
  )
  
  # Download SNP matrix from STRU file
  output$download_snp_stru <- downloadHandler(
    filename = function() paste0("SNP_matrix_STRU_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(stru_matrix(), file)
    }
  )
}

# ---------- Run App ----------
# Launch the Shiny app
shinyApp(ui = ui, server = server)
