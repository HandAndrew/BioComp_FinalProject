# R Shiny Application for SNP Matrix Generation and PCA Visualization from FASTA and STRUCTURE Files

## Author: Andrew Hand 


## Abstract

The rapid analysis of Single Nucleotide Polymorphisms (SNPs) is critical in modern genetic research. To improve and make accessible analysis tools for SNP data, I developed an interactive R Shiny application capable of converting aligned DNA FASTA sequence files and STRUCTURE files into numeric SNP matrices and Principal Component Analysis (PCA) visualization. The application enables users to upload sequence data, perform SNP extraction, conduct PCA, and download results with minimal programming knowledge. This tool was developed to address the need for user-friendly, flexible software that can handle large sequence datasets for population genetics, genomics, and evolutionary studies.

## Introduction 

Single-nucleotide polymorphisms (SNPs) are crucial genetic variation markers and are used in various biological and ecological studies. Efficiently processing raw sequence data into formats amenable to statistical analysis, such as PCA, is a common requirement in these fields.

The transformation of sequence data (e.g., FASTA and STRUCTURE formats) into numeric matrices suitable for statistical computation required a combination of command-line or script-heavy coding workflows, which present barriers for researchers unfamiliar with programming.

To overcome these challenges, I developed a user-centric R Shiny application capable of automatically parsing FASTA and STRUCTURE-formatted files. The application extracts SNPs, constructs clean numeric matrices, performs PCA analysis, and visualizes the results, all within a web browser user interface.

By combining the accessibility of R Shiny with bioinformatics pipelines, this tool makes SNP-based analyses available to a broader range of usages.

## Methods

The full R code can be found under `Full_R_Code.r` file in the main repository or [this link](https://github.com/HandAndrew/BioComp_FinalProject/blob/main/Full_R_Code.R) will take you there. 

### Software and Libraries 

This R Shiny Application was built using R studio version 4.3.0 . This Project used three libraries. Shiny was used for the user interface, Biostrings for handling the sequences, and ggplot2 for data visualization. This app has an upload limit of 400MB ` (options(shiny.maxRequestSize = 400 * 1024^2)` to allow the processing of large datasets. 

### Data Processing 

#### FASTA

Upon the upload of a FASTA file, DNA sequences were read into R using `readDNAStringSet`. Sequences were split into individual nucleotide columns. Only polymorphic sites (sites with more than one observed nucleotide variant) were retained. Nucleotides were numerically encoded (A=1, C=2, G=3, T=4, N=5), with missing or unrecognized characters replaced by 0. Invariant columns and columns containing NA values were filtered out before PCA analysis.

#### STRUCTURE

For STRUCTURE (.stru) formatted files, Data is parsed into a numeric matrix, excluding the first two columns reserved for metadata. Each individual was represented by two rows (diploid genotypes), which were merged into a single row. Combined genotypes were then converted into numeric form. As with FASTA, invariant or incomplete columns were filtered before PCA analysis.

### Principal Component Analysis 

PCA is performed using `prcomp()` function. Two seperate tabs were created for processing FASTA or STRU files. Results, including both the PCA scores and the numeric SNP matrices, were made available for download in CSV format via the `downloadHandler()` function.

## Results

After uploading your FASTA or STRU file, two outputs will be made. A numeric SNP matrix showing the transformation from raw nucleotide or genotype data to a numeric format. A PCA plot seperating individuals or populations based on SNP variation. SNP matrixes and PCA plots are available to be downlaoded as a CSV. An example plot of a FASTA file (Figure 1, Figure 2.) shows the user output. The data used is available on the `Figure1.png` and `Figure2.png` files in the main repository, or here [Figure 1.](https://github.com/HandAndrew/BioComp_FinalProject/blob/main/Figure1.png). [Figure 2.](https://github.com/HandAndrew/BioComp_FinalProject/blob/main/Figure2.png). The data used for this example can be here in the main repository under `Example_Data.fasta` or [here](https://github.com/HandAndrew/BioComp_FinalProject/blob/main/Example_Data.fasta). 

## Discussion and Conclusion

The developed R Shiny application effectively reduces the need for scientists to have programming knowledge to generate SNP matrices and PCA visualization. Users with minimal programming expertise can perform complex genetic data transformations and analyses by using this interface.

A key strength of the application is its flexibility. The app accommodates raw sequence data (FASTA) and already-processed STRUCTURE outputs, offering usages across various genetic studies.

Some limitations persist. The application assumes properly aligned FASTA files and correctly formatted STRUCTURE files. Poor input formatting can cause failures or inaccurate outputs. Moreover, while the PCA implementation suffices for initial explorations, more advanced applications have not been integrated yet. Future iterations could include enhanced error checking and expand the app’s support for additional file formats. 

Overall, this tool represents a significant step toward accessible, reproducible, and scalable SNP analyses. 

## Literature Cited 

> Alboukadel Kassambara. (n.d.). ggplot2: Guide to Create Beautiful Graphics in R. Alboukadel
>	KASSAMBARA. Used for learning how to create the plots with ggplot and other basic R functions.
>
>Alboukadel Kassambara. (2019). GGPlot2 essentials : great data visualization in R. Usa Datanovia. Used for
>	learning how to use ggplot and other basic R functions.
>
>Beeley, C. (2016). Web Application Development with R Using Shiny - Second Edition. Packt Publishing,
>	Limited. https://books.google.com/books?
>	hl=en&lr=&id=FW0dDAAAQBAJ&oi=fnd&pg=PP1&dq=shiny:+Web+Application+Framework+for
> 	+R&ots=kG8gfbxsEO&sig=geyKYQKwtoH6t_DBOsuJp2pv1to#v=onepage&q=shiny%3A%20Web
>	%20Application%20Framework%20for%20R&f=false. Used for learning how to create a shiny app.
>
>freeCodeCamp.org. (2021). R Shiny for Data Science Tutorial – Build Interactive Data-Driven Web Apps. In
>	YouTube. https://www.youtube.com/watch?v=9uFQECk30kA. Beginner tutorial for using R shiny.
>
>Li, Z., & Buck, M. (2021). Beyond history and “on a roll”: The list of the most well‐studied human protein
>	structures and overall trends in the protein data bank. Protein Science, 30(4), 745–760.
>	https://doi.org/10.1002/pro.4038. Table 1 gives the names of large and highly studied human protein
>	sequences used for testing the shiny app.	
