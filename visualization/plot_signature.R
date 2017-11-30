' plot_signature.R

Usage: plot_signature.R -i INPUT -o OUTPUT [ --ncol COLS --format FORMAT --ncores NCORES ]

Options:
    -i --input INPUT        Path to input mutation catalog table. For an example,
                                see example-files/signatures_probabilities.txt

    -o --output OUTPUT      Path to output figure

    --ncol COLS             Number of columns to arrange into (defaults to the number of signatures)

    --format FORMAT         Options: pdf, svg, png. If not provided, will try to guess based
                                on the file extension. Otherwise, defaults to PDF.

    --ncores NCORES         Specifies number of cores. If > 1, then runs in parallel mode.
' -> doc

library(docopt)
args <- docopt(doc)

library(plyr)
library(tidyverse)
library(cowplot)
source("scripts/polarBarPlot.r")
source("scripts/multiplot.r")
library(doParallel)

table <- read_tsv(args[['input']])

factorCat <- function (a, b) {f = factor(c(as.character(a), as.character(b))); return(f);}

figureVector = list()
n_signatures = table %>% dplyr::select(-category, -mutation_type) %>% dim %>% .[2]
n_mutation_types = dim(table)[1]

if (is.null(args[['ncol']])) {
    n_col = n_signatures
} else {
    n_col = as.numeric(args[['ncol']])
}

if (! is.null(args[['format']])) {
    stopIfNot(args[['format']] %in% c('pdf', 'png', 'svg'))
    output_format <- args[['format']]
} else if (endsWith(args[['output']], '.png')) {
    output_format <- 'png'
} else if (endsWith(args[['output']], '.svg')) {
    output_format <- 'svg'
} else {
    output_format <- 'pdf'
}

message(sprintf('Outputting file of type: %s', output_format))

if (output_format == 'png') {
    png(
        args[['output']], 
        width = 800 * n_col, 
        height = ceiling(n_signatures/n_col) * 800, 
        onefile = TRUE, 
        bg = 'white'
    )
} else if (output_format == 'svg') {
    svg(
        args[['output']], 
        width = 6 * n_col, 
        height = ceiling(n_signatures/n_col) * 6,       
        onefile = TRUE, 
        bg = 'transparent'
    )
} else if (output_format == 'pdf') {
    pdf(
        args[['output']], 
        width = 10 * n_col, 
        height = ceiling(n_signatures/n_col)* 10, 
        onefile = TRUE, 
        bg = 'transparent'
    )
}

catalogs <- table %>%
    tidyr::gather(signature, count, -category, -mutation_type)

if (!is.null(args[['ncores']])) {
    if (as.numeric(args[['ncores']]) != 1) {
        n_cores = as.numeric(args[['ncores']])
        message(sprintf("Running in parallel on %s cores.", n_cores))
        run_in_parallel = TRUE
        registerDoParallel(n_cores)
    } else {
        message("Running in serial on one core.")
        run_in_parallel = FALSE
    }
}

figure_list <- catalogs %>% dlply('signature', function(sig_table) {
  df <- data.frame(
    family=factorCat(sig_table$category, sig_table$category),
    item=factorCat(sig_table$mutation_type, sig_table$mutation_type),
    value=c(sig_table$count, rep(0, n_mutation_types)),
    score=c(rep(c('Mutations'), n_mutation_types), rep(c('mutations'), n_mutation_types))
  )

  p <- polarBarPlot(df, familyLabel=TRUE, legTitle="Legend") + theme(legend.position = 'none')
  message(sprintf('Created figure for %s', sig_table$signature[1]))
  return(p)
}, .parallel = run_in_parallel)

colour_df <- data.frame(
  five_prime = c(rep('A', 4), rep('C', 4), rep('G', 4), rep('T', 4)), 
  three_prime = rep(c('A', 'C', 'G', 'T'), 4),
  status = factor(1:16, levels=1:16),
  colour = hsv(h=sort(rep(c(0, 0.3, 0.6, 0.9), 4)), s=rep(c(0.1, 0.2, 0.5, 1), 4), v=rep(c(0.8, 0.8, 0.8, 0.8), 4))
)

legend_plot <- ggplot(colour_df) +
  labs(x = "5' context", y = "3' context") +
  geom_tile(aes(x = five_prime, y = three_prime, fill = status), height = 0.9, width = 0.9) +
  scale_fill_manual(values = as.character(colour_df$colour)) +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_text(size = 10)
  )

legend_plot <- plot_grid(legend_plot, NULL, nrow = 2, rel_heights = c(1, ceiling(n_signatures / n_col)))

plot_grid(plotlist = figure_list, ncol = n_col) %>% plot_grid(legend_plot, ncol = 2, rel_widths = c(n_col, 1))
 
dev.off()
