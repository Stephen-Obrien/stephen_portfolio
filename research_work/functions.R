library(tidyverse)
library(ggbreak)


#Reference Factors
rf_4methylguaiacol <- 1.187514325
rf_4ethylphenol <- 1.398982903
rf_guaiacol <- 1.303681539
rf_4ethylguaiacol <- 1.182673218


##FUNCTIONS

#importfunction
import_gc <- function(path){
  data <- read_csv(path,
                   col_select= c(1,2),
                   skip = 2,
                   col_names = FALSE)
  data <- data |> rename("time" = "X1", "abundance" = "X2") |> 
    mutate(time = as.numeric(time),
           abundance = as.numeric(abundance))
  return(data)
}

#  
GC_peak <- function(data, time1 = 11.8, time2 = 12.1) {
  library(pracma)
  library(knitr)
  ## EGDE area
  
  subset <- data |> 
    filter(time >= 8.55,
           time <= 8.75)
  
  time <- subset$time
  abundance <- subset$abundance
  
  EDGE <- round(trapz(time, abundance), digits = 2)
  
  ##phenol area
  subset2 <- data |> 
    filter(time >= time1,
           time <= time2)
  
  time <- subset2$time
  abundance <- subset2$abundance
  
  phenol<- round(trapz(time, abundance), digits = 2)
  
  #Summary table
  table<- tibble(
    Peak = c("EGDE", "phenol"),
    Area = c(EDGE, phenol))
  return(kable(table))
}

GC_plot <- function(data, time1 = 11.9, time2 = 12, title = NULL, ggbreak = NULL) {
  plot <- data |> ggplot(mapping = aes(x = time, y = abundance)) +
    geom_area(fill = "lightblue", color = "black") +
    xlim(8, time2 + .25) +
    ylim(0, 8e5) +
    labs(x = "Time (min)", y = "Abundance",
         title  = title) +
    theme(plot.title = element_text(hjust = 0.5, size = 18))
  
  return(plot)
}

GC_plot_break <- function(data, time1 = 11.9, time2 =12, title = NULL,
                          nudge = .15, break1 = 9, break2 = 11.5) {
  library(ggbreak)
  plot <- GC_plot(data, time1, time2, title)
  plot <- plot +
    scale_x_break(c(break1, break2))
  
  y <- data |> 
    mutate(type = if_else(time < 9, 1, 2)) |> 
    filter(time < 14)
  
  y <- y |> group_by(type) |> 
    slice_max(abundance)
  
  plot <- plot +
    geom_point(y, mapping = aes(x = time, y = abundance)) +
    geom_text(data = y, mapping = aes(label = time), nudge_x = nudge)
  
  return(plot)
}

GC_anal <- function(data, time1, time2, title = NULL) {
  plot <- GC_plot(data, time1, time2, title)
  print(plot)
  return(GC_peak(data, time1, time2))
}

fourMG_plot <- function(data, time1 = 12.85, time2 = 13, title = NULL) {
  GC_plot(data, time1, time2, title) +
    ggbreak::scale_x_break(c(9,12.5), scales = 1)
}

