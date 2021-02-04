# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

gene.manhattan <- function(df, threshold, hlight, col, ylims, title, sugg){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(START)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, START) %>%
    mutate( BPcum=START+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(GENE %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(FDR < sugg, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=sign(COEF)*-log10(FDR))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
#    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(threshold), color = "red", size = 1) +
    geom_hline(yintercept = -log10(sugg), color = "red", size = 1, linetype="dashed") +
    geom_hline(yintercept = log10(threshold), color = "blue", size = 1) +
    geom_hline(yintercept = log10(sugg), color = "blue", size = 1, linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(GENE)), alpha=1, size=7, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 18) +
      theme( 
      plot.title = element_text(hjust = 0.0),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x = element_text(size = 12)
    )
}
