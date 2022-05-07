rm(list = ls())
library ( car )
library(grDevices)
library(ggplot2)
library(scales)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(rms)
library(dplyr)

graphics.off()

setwd("~/Desktop/MDS_Project/only_test_8_cl")
MDS_data = read.csv("final_MDS.csv")
attach(MDS_data)
n = dim(MDS_data)[1]

clusters_df = read.table("clusters.txt")
test_idx_df = read.table("test_idx_real.txt")         #
test_idx_df = test_idx_df + 1
predicted_times_df = read.table("predicted_times.txt")
empirical_dead_df = read.table("empirical_dead.txt")
km_dead_dist_sca_df= read.table("km_dead_dist_sca.txt")

predicted_times = predicted_times_df[, 1]
empirical_dead = empirical_dead_df[, 1]
km_dead_dist_sca= km_dead_dist_sca_df[, 1]




fit = lm ( empirical_dead ~ km_dead_dist_sca)
summary(fit)

linearHypothesis ( fit, c(0, 1), 1 )



columns = c("gender", "cohort", "IPSSR", "genomic_group", "hgb",
            "plt", "bmblasts", "fup_status", "age", "time",
            # no cluster
            "leu", "neut")

df = MDS_data[test_idx_df[,1], columns]
df = cbind(df, clusters = clusters_df[,1])
head(df)
detach(MDS_data)
attach(df)

# #hemoglobin
# boxplot(hgb~clusters)
# hgb_mean<- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# for(i in 1:length(hgb_mean)) {
#   hgb_mean[i] <- mean(hgb[which(clusters==i-1)])
# }
# hgb_mean_wd=mean(hgb)
# hgb_mean_wd
# hgb_mean
# 
# 
# #bomboblast
# boxplot(bmblasts~clusters)
# bmblasts_mean<- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# for(i in 1:length(bmblasts_mean)) {
#   bmblasts_mean[i] <- mean(bmblasts[which(clusters==i-1)])
# }
# 
# bmblasts_mean_wd=mean(bmblasts)
# bmblasts_mean_wd
# bmblasts_mean
# 
# 
# boxplot(age~clusters)
# age_mean<- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# for(i in 1:length(age_mean)) {
#   age_mean[i] <- mean(age[which(clusters==i-1)])
# }
# age_mean_wd=mean(age)
# age_mean_wd
# age_mean
# 
# #Not used for clustering
# boxplot(leu~clusters)
# leu_mean<- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# for(i in 1:length(leu_mean)) {
#   leu_mean[i] <- mean(leu[which(clusters==i-1)])
# }
# leu_mean_wd=mean(leu)
# leu_mean_wd
# leu_mean


# boxplot(neut~clusters)
# boxplot(plt~clusters)
   


############## PLOT FROM HERE -> GG plot (barplots)
cohort.table = table(cohort, clusters)
cohort.table
# barplot(cohort.table, beside=TRUE, legend.text=rownames(cohort.table), ylab="absolute frequency")
plotdata <- df %>%
  group_by(clusters,cohort) %>%
  summarize(n = n()) %>% 
  mutate(pct = n/sum(n),
         lbl = scales::percent(pct))

graph = ggplot(plotdata, 
               aes(x = factor(clusters),
                   y = pct,
                   fill = factor(cohort))) +
  geom_bar(stat = "identity",
           position = "fill",
           size = 5,
           alpha = 1) +
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     label = percent) +
  geom_text(aes(label = lbl), 
            size = 4, 
            position = position_stack(vjust = 0.5), color="#f8f5ec") +
  labs(y = "",
       fill = "Cohort",
       x = "Cluster",
       title = "Patients by cohort and cluster",
  ) +
  theme(text = element_text(size = 12, color="#252e47"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.title = element_text(colour = "#252e47"),
        axis.title.x = element_text(colour = "#252e47"),
        axis.title.y = element_text(colour = "#252e47"),
        axis.text.x = element_text(colour = "#252e47"),
        axis.text.y = element_text(colour = "#252e47"),
        plot.title = element_text(colour = "#252e47"),
        title = element_text(colour = "#252e47")
        
  ) + scale_fill_manual(values=c("#90e0ef",      #COLORS
                                 "#0077b6"))
#x11()
graph
ggsave(graph, filename = "cohort_clusters.png",  bg = "transparent", dpi=600)     #SAVING




#Gender, not working
# graph = ggplot(plotdata, 
#                aes(x = factor(clusters),
#                    y = pct,
#                    fill = factor(gender))) +
#   geom_bar(stat = "identity",
#            position = "fill",
#            size = 5,
#            alpha = 0.92) +
#   scale_y_continuous(breaks = seq(0, 1, .2), 
#                      label = percent) +
#   geom_text(aes(label = lbl), 
#             size = 5, 
#             position = position_stack(vjust = 0.5)) +
#   labs(y = "",
#        fill = "Cluster",
#        x = "Gender",
#        title = "Patients by gender and cluster",
#   ) +
#   theme(text = element_text(size = 12)) +
#   theme_minimal()
# #x11()
# graph



# genomic_group.table = table(genomic_group, clusters)
# genomic_group.table
# plotdata <- df %>%
#   group_by(clusters,genomic_group) %>%
#   summarize(n = n()) %>% 
#   mutate(pct = n/sum(n),
#          lbl = scales::percent(pct))
# 
# graph = ggplot(plotdata, 
#                aes(x = factor(clusters),
#                    y = pct,
#                    fill = factor(genomic_group))) +
#   geom_bar(stat = "identity",
#            position = "fill",
#            size = 5,
#            alpha = 0.92) +
#   scale_y_continuous(breaks = seq(0, 1, .2), 
#                      label = percent) +
#   geom_text(aes(label = lbl), 
#             size = 5, 
#             position = position_stack(vjust = 0.5)) +
#   labs(y = "",
#        fill = "genomic_group",
#        x = "Cluster",
#        title = "Patients by genomic_group and cluster",
#   ) +
#   theme(text = element_text(size = 12)) +
#   theme_minimal()
# #x11()
# graph
# 
# 
# 
# 
# 
# genomic_group.table = table(genomic_group, clusters)
# genomic_group.table
# plotdata <- df %>%
#   group_by(genomic_group,clusters) %>%
#   summarize(n = n()) %>% 
#   mutate(pct = n/sum(n),
#          lbl = scales::percent(pct))
# 
# 
# graph = ggplot(plotdata, 
#                aes(x = factor(genomic_group),
#                    y = pct,
#                    fill = factor(clusters))) +
#   geom_bar(stat = "identity",
#            position = "fill",
#            size = 5,
#            alpha = 0.92) +
#   scale_y_continuous(breaks = seq(0, 1, .2), 
#                      label = percent) +
#   geom_text(aes(label = lbl), 
#             size = 5, 
#             position = position_stack(vjust = 0.5)) +
#   labs(y = "",
#        fill = "clusters",
#        x = "genomic group",
#        title = "Patients by genomic_group and cluster",
#   ) +
#   theme(text = element_text(size = 12)) +
#   theme_minimal()
# #x11()
# graph


#ANOTHER BARPLOT

IPSSR.table = table(IPSSR, clusters)
IPSSR.table
plotdata <- df %>%
  group_by(clusters,IPSSR) %>%
  summarize(n = n()) %>% 
  mutate(pct = n/sum(n),
         lbl = scales::percent(pct))

graph = ggplot(plotdata, 
               aes(x = factor(clusters),
                   y = pct,
                   fill = factor(IPSSR))) +
  geom_bar(stat = "identity",
           position = "fill",
           size = 5,
           alpha = 1) +
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     label = percent) +
  geom_text(aes(label = lbl), 
            size = 4, 
            position = position_stack(vjust = 0.5), color="#f8f5ec") +
  labs(y = "",
       fill = "IPSSR",
       x = "Cluster",
       title = "Patients by IPSSR and cluster",
  ) +
  theme(text = element_text(size = 12, color="#252e47"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.title = element_text(colour = "#252e47"),
    axis.title.x = element_text(colour = "#252e47"),
    axis.title.y = element_text(colour = "#252e47"),
    axis.text.x = element_text(colour = "#252e47"),
    axis.text.y = element_text(colour = "#252e47"),
    plot.title = element_text(colour = "#252e47"),
    title = element_text(colour = "#252e47")
    
  ) + scale_fill_manual(values=c("#caf0f8", "#90e0ef", "#00b4d8",
                                               "#0077b6", "#03045e"))
#x11()
graph
ggsave(graph, filename = "IPSSR.png",  bg = "transparent", dpi=600)


# fup_status.table = table(fup_status, clusters)
# fup_status.table
# # barplot(cohort.table, beside=TRUE, legend.text=rownames(cohort.table), ylab="absolute frequency")
# plotdata <- df %>%
#   group_by(clusters,fup_status) %>%
#   summarize(n = n()) %>% 
#   mutate(pct = n/sum(n),
#          lbl = scales::percent(pct))
# 
# graph = ggplot(plotdata, 
#                aes(x = factor(clusters),
#                    y = pct,
#                    fill = factor(fup_status))) +
#   geom_bar(stat = "identity",
#            position = "fill",
#            size = 5,
#            alpha = 0.92) +
#   scale_y_continuous(breaks = seq(0, 1, .2), 
#                      label = percent) +
#   geom_text(aes(label = lbl), 
#             size = 5, 
#             position = position_stack(vjust = 0.5)) +
#   labs(y = "",
#        fill = "Cluster",
#        x = "fup_status",
#        title = "Patients by fup_status and cluster",
#   ) +
#   theme(text = element_text(size = 12)) +
#   theme_minimal()
# x11()
# graph





#par(mfrow=c(3,3))
# for(i in 1:length(unique(clusters))) {
# pie(prop.table(table(genomic_group[which(clusters == i-1)])))
# }




# graph = ggplot(df,
#                aes(x = factor(clusters),
#                    y = hgb, fill =factor(clusters))) +
#   geom_boxplot(size=0.5,
#                outlier.size = 1,
#                alpha = 0.8) +
#   labs(title="hgb distribution by clusters",
#        x = "clusters",
#        y = "hgb") +
#   theme_minimal() +
#   theme(legend.position = "right")
# graph
# 
# 
# graph = ggplot(df,
#                aes(x = factor(clusters),
#                    y = bmblasts, fill =factor(clusters))) +
#   geom_boxplot(size=0.5,
#                outlier.size = 1,
#                alpha = 0.8) +
#   labs(title="bmblasts distribution by clusters",
#        x = "clusters",
#        y = "bmblasts") +
#   theme_minimal() +
#   theme(legend.position = "right")
# graph
# 
# 
# graph = ggplot(df,
#                aes(x = factor(clusters),
#                    y = age, fill =factor(clusters))) +
#   geom_boxplot(size=0.5,
#                outlier.size = 1,
#                alpha = 0.8) +
#   labs(title="age distribution by clusters",
#        x = "clusters",
#        y = " age") +
#   theme_minimal() +
#   theme(legend.position = "right")
# graph
# 
# 
# graph = ggplot(df,
#                aes(x = factor(clusters),
#                    y = hgb, fill =factor(clusters))) +
#   geom_boxplot(size=0.5,
#                outlier.size = 1,
#                alpha = 0.8) +
#   labs(title="hgb distribution by clusters",
#        x = "clusters",
#        y = "hgb") +
#   theme_minimal() +
#   theme(legend.position = "right")
# graph



####BOXPLOTS

graph = ggplot(df,
               aes(x = factor(clusters),
                   y = bmblasts, fill =factor(clusters))) +
  geom_boxplot(size=0.5,
               outlier.size = 1,
               alpha = 0.8) +
  labs(title="Blast cells [%] distribution by clusters",
       x = "Clusters",
       y = "Blast cells [%]", fill="Clusters") +
  theme_minimal() +
  theme(text = element_text(size = 12, color="#252e47"),legend.position = "right",
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.title = element_text(colour = "#252e47"),
        axis.title.x = element_text(colour = "#252e47"),
        axis.title.y = element_text(colour = "#252e47"),
        axis.text.x = element_text(colour = "#252e47"),
        axis.text.y = element_text(colour = "#252e47"),
        plot.title = element_text(colour = "#252e47"),
        title = element_text(colour = "#252e47"))+
  scale_fill_manual(values=c("#f77189", "#c69432", "#82a931",
                             "#34af8a", "#37aabb", 
                             "#8197f4", "#f45deb"))
graph
ggsave(graph, filename = "bmblasts_clusters.png",  bg = "transparent", dpi=600)


graph = ggplot(df,
               aes(x = factor(clusters),
                   y = age, fill =factor(clusters))) +
  geom_boxplot(size=0.5,
               outlier.size = 1,
               alpha = 0.8) +
  labs(title="Age distribution by clusters",
       x = "Clusters",
       y = " Age", fill="Clusters") +
  theme_minimal() +
  theme(text = element_text(size = 12, color="#252e47"),legend.position = "right",
                                          panel.background = element_rect(fill = "transparent"), # bg of the panel
                                          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                                          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                                          legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
                                          axis.title = element_text(colour = "#252e47"),
                                          axis.title.x = element_text(colour = "#252e47"),
                                          axis.title.y = element_text(colour = "#252e47"),
                                          axis.text.x = element_text(colour = "#252e47"),
                                          axis.text.y = element_text(colour = "#252e47"),
                                          plot.title = element_text(colour = "#252e47"),
                                          title = element_text(colour = "#252e47"))+
         scale_fill_manual(values=c("#f77189", "#c69432", "#82a931",
                                                               "#34af8a", "#37aabb", 
                                                            "#8197f4", "#f45deb"))
graph
ggsave(graph, filename = "age_clusters.png",  bg = "transparent", dpi=600)








df[3, ]
predicted_times[3]

diff=abs(df[, 10]-predicted_times)
err=diff/df[, 10]
summary(fit)
shapiro.test(fit$residuals)

