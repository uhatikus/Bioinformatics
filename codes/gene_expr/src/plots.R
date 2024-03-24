library(ggplot2)

plot_exp_vec <- function(exp_vec){
  p <- ggplot(melt(as.matrix(-2*exp_vec+1)), aes(x = "", y = Var2, fill = as.factor(value))) +
    geom_tile(color = "white",
              lwd = 0.2,
              linetype = 1) +
    labs(title = "Pattern",
         x = "",
         y = "Genes") + 
    scale_fill_manual(values = c("-1" = "#00FF00", "1" = "#FF0000"), labels=c("-1" = "LogFC > 0", "1" = "LogFC < 0"), name="") +
    coord_fixed(1) +
    theme(axis.text.x=element_blank(), legend.position="left", plot.title= element_text(hjust = 0.5)) 
  return(p)
}

plot_d_exp <- function(d_exp){
  p <- ggplot(melt(as.matrix(- 2*d_exp + 1)), aes(x = Var1, y = Var2, fill = as.factor(value))) +
    geom_tile(color = "white",
              lwd = 0.05,
              linetype = 1) +
    labs(title = "Patients",
         x = "",
         y = "Genes")  +
    scale_fill_manual(values = c("-1" = "#00C957", "1" = "#DC143C"), labels=c("-1" = "Gene Expression > Median", "1" = "Gene Expression < Median"), name="") +
    coord_fixed(10) + 
    theme(axis.text.x = element_blank(), legend.position="right", plot.title= element_text(hjust = 0.5), axis.text.y=element_blank(), axis.title.y=element_blank())
  return(p)
}

plot_arrow_with_text_boxes <- function(xlim = c(-0.05, 1.05)){
  # Add the arrow
  arrow_df <- data.frame(
    x = c(0.0, 1.0),
    y = c(0, 0)
  )
  p <- ggplot(arrow_df, aes(x, y)) + theme_void()
  p <- p + geom_path(arrow = arrow(length = unit(0.3, "cm"))) +
    coord_cartesian(xlim = xlim, ylim = c(-0.1, 0.25))
  
  # add Matched rectangle
  p <- p + geom_rect(aes(xmin = 0.0, xmax = 0.2, ymin = 0.05, ymax = 0.15),
                     fill = "#00A7E1", color = "black", alpha = 0.7)
  p <- p + geom_text(aes(x = 0.1, y = 0.1, label = "Matched"))
  
  
  # add Unmatched rectangle
  p <- p + geom_rect(aes(xmin = 0.8, xmax = 1.0, ymin = 0.05, ymax = 0.15),
                     fill = "#F17720", color = "black", alpha = 0.7)
  p <- p + geom_text(aes(x = 0.9, y = 0.1, label = "Unmatched"))
  
  # add custom cutoff values
  p <- p + geom_text(aes(x = 0.05, y = 0.25, label = "0.2"))
  p <- p + geom_text(aes(x = 0.5, y = 0.25, label = "0.5"))
  
  
  return(p)
}

plot_survival <- function(tcga_data_for_coxph){
  gsp <- ggsurvplot(survfit(Surv(time, alive) ~ group, data=tcga_data_for_coxph),
                          legend.title = "",
                          legend.labs = c("Matched", "Unmatched"),
                          ggtheme = theme_classic(), palette = c("#00A7E1", "#F17720"))
  
  p <- gsp$plot
  
  p <- p + theme(plot.margin = unit(c(0,1, 1, 2), "cm"))+ theme(legend.position = c(0.9, 0.9))
  
  return(p)
}

plot_ggforest <- function(tcga_data_for_coxph){
  cox_model <-coxph(Surv(time, alive) ~ age_at_initial_pathologic_diagnosis + gender + group, data=tcga_data_for_coxph)
  p <- ggforest(cox_model, data = tcga_data_for_coxph)
  return(p)
}

plot_cutoff_vs_p_value <- function(cutoff_values, p_values, n_patient, coeff){
  p <- ggplot(data = data.frame(x = 1-cutoff_values, y1 = p_values, y2 = n_patient/coeff/10), aes(x = x)) +
    geom_point(aes(y = y1, color = "Y1")) +
    geom_point(aes(y = y2, color = "Y2")) +
    geom_line(linetype = "dashed", aes(y = y1), color = "#00A7E1") +
    geom_line(linetype = "dashed", aes(y = y2), color = "#F17720") +
    labs(title = "", x = "level of dissimilarity", y = "COXPH p-values") + 
    scale_y_continuous(name = "COXPH p-values", sec.axis = sec_axis(~.*coeff/10, name="Number of patients in analysis")) +
    scale_color_manual(name = "",
                       values = c("Y1" = "#00A7E1", "Y2" = "#F17720"),
                       labels = c("p-value", "Number of patients")) + 
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + 
    # geom_text(aes(x = 0.22, y = 0.02, label = "p-value = 0.05")) + 
    theme(legend.position = c(0.7, 0.9), plot.margin = unit(c(0,3, 1, 1), "cm")) +
    theme_bw() + 
    theme_classic()
  return(p)
}
