#!/usr/bin/env Rscript

## =============================================================================
## plot_signature.R
##
## Uses Christie Haskell's polarBarPlot to create a "Mutation Signature Fingerprint"
##
## USAGE: Rscript plot_signature.R <input path> <output path> <number of signatures> <number of columns>
##
## EXAMPLE: Rscript plot_signature.R ../example-files/signature_probabilities.txt signature_probabilities.svg 30 6
##
## The above line will produce an SVG file visualizing the 30 signatures
## contained in the example file. There will be 5 rows of 6 signatures each.
##
## This code developed by Eric Zhao, 2016
## as part of work supported by the BC Genome Sciences Centre
## See www.bcgsc.ca for more information
##
## This code is hosted on github at https://github.com/eyzhao/mutation-signature-tools
## =============================================================================

library(reshape2)
library(ggplot2)
library(plyr)
library(Cairo)
library(RSvgDevice)
library(stringr)
source("scripts/polarBarPlot.r")
source("scripts/multiplot.r")

args <- commandArgs(TRUE)
inputPath = args[1]
outputPath = args[2]
numSignatures = as.numeric(args[3])
numCols = as.numeric(args[4])

table <- read.table(inputPath, sep='\t', header=TRUE) # Input

outputFileType = substr(outputPath, sapply(gregexpr("\\.", outputPath), tail, 1), nchar(outputPath))
print(outputFileType)

if (outputFileType == '.pdf') {
    CairoPDF(file = outputPath, width = 10 * numCols, height = ceiling(numSignatures/numCols) * 10, onefile = TRUE, bg = 'transparent') # Output
} else if (outputFileType == '.png') {
    CairoPNG(file = outputPath, width = 800 * numCols, height = ceiling(numSignatures/numCols) * 800, onefile = TRUE, bg = 'white')
} else if (outputFileType == '.svg') {
    print('exporting SVG')
    devSVG(file = outputPath, width = 6 * numCols, height = ceiling(numSignatures/numCols) * 6,       
 onefile = TRUE, bg = 'transparent')
} else {
    print('Output file types supported: PDF, SVG, PNG')
    stop()
}

figureVector = list()

for (i in numSignatures : 1) {
factorCat <- function (a, b) {f = factor(c(as.character(a), as.character(b))); return(f);}
# For some reason, there's an error when there are not two 

  targetColumn = dim(table)[2] + 1 - i

  df<-data.frame(
    family=factorCat(table$Substitution.Type, table$Substitution.Type),
    item=factorCat(table$Somatic.Mutation.Type, table$Somatic.Mutation.Type),
    value=c(table[[targetColumn]], rep(0, 96)),
    score=c(rep(c('Mutations'), 96), rep(c('mutations'), 96)))

  p<-polarBarPlot(df,familyLabel=TRUE, legTitle="Legend");

  figureVector[[length(figureVector) + 1]] <- p
}

multiplot(plotlist = figureVector, cols = numCols)
 
dev.off()
