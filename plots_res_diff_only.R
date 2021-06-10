#Erin's preferences

comp_df_diff_plot <- comp_df_diff %>% filter(!str_detect(source, "cces")) %>% 
    mutate(method =  factor(case_when(str_detect(source, "unweighted") ~ "Raw Survey", 
                               str_detect(source, "rake") ~ "Mean Calibration",
                               str_detect(source, "post") ~ "Post-Strat",
                               str_detect(source, "kpop") ~ "Kpop", 
                               TRUE ~ NA_character_),
                            levels = c("Raw Survey", "Mean Calibration", "Post-Strat", "Kpop")),
           zero = 0,
           natl_diff =natl_diff*100,
           cces_target = comp_df_diff$est[comp_df_diff$source_name == "CCES Modeled"])

text_for_gg <- data.frame(source_name =  c(4.50, 4.50, 5),
                          method =factor("Kpop", levels(comp_df_diff_plot$method)),
                          est = c(comp_df_diff_plot$natl_diff[1], 
                                  comp_df_diff_plot$cces_target[1], 0),
                          label = c(" True\n National\n Diff", 
                                    " Weighted\n CCES\n Target", " "))


ggplot(comp_df_diff_plot) +
    geom_pointrange(aes(x = source_name, y = est, ymin = est - 1.96*SE, ymax = est + 1.96*SE)) +
    facet_grid(cols = vars(method),  scales = "free_x", space = "free_x") +
    geom_hline(aes(yintercept = zero)) +
    geom_hline(aes(yintercept = natl_diff), linetype = 2, color = "grey60") +
    geom_hline(aes(yintercept = cces_target), linetype = 2) +
   
    scale_y_continuous(breaks = c(natl_diff*100, seq(-5, 10, 5)),
                       minor_breaks = NULL,
                       labels = scales::percent_format(scale=1, accuracy = 0.1)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = NULL, y = "Estimated Modeled Vote Difference (95% CI)") +
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    geom_text(data = text_for_gg, aes(x = source_name, 
                                      y = est,
                                      label = label), 
              angle = c(-90,90, 90), color = c("grey60", "black", "purple"),
              hjust = 0, size = 2.6) +
    theme()
ggsave("./plots/bmaxvar_votediff_alt2.pdf", width = 8.45, height = 6)



ggplot(comp_df_diff_plot) +
    geom_pointrange(aes(x = source_name, y = est, ymin = est - 1.96*SE, ymax = est + 1.96*SE)) +
    facet_grid(cols = vars(method),  scales = "free_x", space = "free_x") +
    geom_hline(aes(yintercept = zero)) +
    geom_hline(aes(yintercept = natl_diff), linetype = 2, color = "grey60") +
    geom_hline(aes(yintercept = cces_target), linetype = 2) +
    
    scale_y_continuous(breaks = c(natl_diff*100, seq(-5, 10, 5)),
                       minor_breaks = NULL,
                       labels = scales::percent_format(scale=1, accuracy = 0.1)) +
    #theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = NULL, y = "Estimated Modeled Vote Difference (95% CI)") +
    #theme_classic() + #need this to get rid of the lines inbetween the facets
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    
    geom_text(data = text_for_gg, aes(x = source_name, 
                                      y = est,
                                      label = label), 
              angle = c(-90,90, 90), color = c("grey60", "black", "purple"),
              hjust = 0, size = 2.6) +
    theme( axis.title.y.right = element_blank(),                # 
           axis.text.y.right = element_blank(),                 # 
           axis.ticks.y = element_blank(),                      # 
           axis.text.y = element_text(margin = margin(r = 0)),  # 
           panel.spacing = unit(0, "mm") ,
           strip.background = element_rect(size = 0.5)) 
ggsave("./plots/bmaxvar_votediff_alt3.pdf", width = 8.45, height = 6)




comp_df_diff_plot <- rbind(comp_df_diff_plot[1,], blankrow,
                           comp_df_diff_plot[c(2:4),], blankrow,
                           comp_df_diff_plot[5,], blankrow,
                           comp_df_diff_plot[c(7:10),]) %>% 
    mutate(source_name = as.character(source_name) )

comp_df_diff_plot[c(2,6,8),"source_name"] <- c("null1", "null2", "null3")
comp_df_diff_plot <- comp_df_diff_plot %>%
    mutate(source_name = factor(source_name, 
                                levels = as.character(comp_df_diff_plot$source_name)))


#Chad's preferences

comp_df_diff_plot <- comp_df_diff_plot[c(1:4,6,5,7:10),]
comp_df_diff_plot <- comp_df_diff_plot %>% mutate(x_new = c(1,3,4,5,6,8,10,11,12,13))
#comp_df_plot_diff <-
    ggplot(comp_df_diff_plot) +
    geom_hline(yintercept = c(0, 
                              if(pop_weights) {natl_diff*100},
                              #if(pop_weights) {comp_df_diff$est[comp_df_diff$source_name == "CCES Orig"]},
                              comp_df_diff$est[comp_df_diff$source_name == "CCES Modeled"]),
               linetype = c("solid", 
                            if(pop_weights) {"dashed"}, 
                            #if(pop_weights) {"longdash"},
                            "dashed"),
               color = c("black", 
                         if(pop_weights) {"gray60"}, 
                         #if(pop_weights) {"black"},
                         "black")) +
    
    geom_pointrange(aes(x = x_new, y = est, ymin = est - 1.96*SE, ymax = est + 1.96*SE)) +
    scale_y_continuous(breaks = c(natl_diff*100, seq(-5, 10, 5)),
                       minor_breaks = NULL,
                       labels = scales::percent_format(scale=1, accuracy = 0.1)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    
    labs(x = NULL, y = "Estimated Modeled Vote Difference (95% CI)") +
    scale_x_continuous(breaks = comp_df_diff_plot$x_new, 
                       labels = as.character(comp_df_diff_plot$source_name),
                       limits = c(1,13.5)) +
   theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    
    annotate(geom = "text",
             x = 13.5,
             y = if(pop_weights) {natl_diff * 100} else {comp_df_diff$est[comp_df_diff$source_name == "CCES Modeled"]},
             label = if(pop_weights) {" True\n National\n Diff"} else {""},
              angle = -90,hjust = 0,
             color = "gray60", size = 3) +
    annotate(geom = "text",
             x = 13.5,
             y = comp_df_diff$est[comp_df_diff$source_name == "CCES Modeled"],
             label = " CCES\n Modeled Target",
             color = "black",hjust = 0,
             angle = 90, size = 3) +
    annotate(geom = "text",
             x = 1,
             y = 9.2,
             label = "Raw\nSurvey",
             color = "black", size = 4) +
        annotate(geom = "text",
                 x = 4.5,
                 y = 9.2,
                 label = "Mean\nCalibration",
                 color = "black", size = 4) +
        annotate(geom = "text",
                 x = 8,
                 y = 9.2,
                 label = "Post\nStratification",
                 color = "black", size = 4) +
        annotate(geom = "text",
                 x = 11.5,
                 y = 9.2,
                 label = "Kpop\n",
                 color = "black", size = 4) +
    theme(legend.title = element_blank())

ggsave("./plots/bmaxvar_alt_votediff.pdf", width = 8.4, height = 6)    
