library("survival")
library("survminer")
library(ggplot2)

rm(list=ls())

MDS_data = read.csv("MDS_data.csv")
n = length(MDS_data)
attach(MDS_data)

dodx = pmin(dodx_ngs, dodx_hsct, na.rm = TRUE)
time = as.integer(difftime(dofup, dodx, units="weeks") / 30)

MDS_data$IPSSR <- as.character(MDS_data$IPSSR)
IPSSR <- as.character(IPSSR)

MDS_data$IPSSR[which(MDS_data$IPSSR == "very low")] = 1
MDS_data$IPSSR[which(MDS_data$IPSSR == "low")] = 2
MDS_data$IPSSR[which(MDS_data$IPSSR == "intermediate")] = 3
MDS_data$IPSSR[which(MDS_data$IPSSR == "high")] = 4
MDS_data$IPSSR[which(MDS_data$IPSSR == "very high")] = 5

IPSSR[which(IPSSR == "very low")] = 1
IPSSR[which(IPSSR == "low")] = 2
IPSSR[which(IPSSR == "intermediate")] = 3
IPSSR[which(IPSSR == "high")] = 4
IPSSR[which(IPSSR == "very high")] = 5

MDS_data$IPSSR = as.numeric(MDS_data$IPSSR)
IPSSR = as.numeric(IPSSR)

MDS_data <- cbind(MDS_data, dodx, time, IPSSR)
MDS_data <- data.frame(MDS_data, dodx, time, IPSSR)
detach(MDS_data)
attach(MDS_data)

plots <- list()

km.model <- survfit(Surv(time, fup_status)~1, type='kaplan-meier')
Surv(time, fup_status)
km.model
summary(km.model)

# plot(km.model, conf.int = F, xlab="Time (weeks)", ylab = "ALive = S(t)", main="KM Model")
# plot(km.model, conf.int = T, xlab="Time (weeks)", ylab = "ALive = S(t)", main="KM Model")

plots[[1]] <- ggsurvplot(km.model,
           data= MDS_data,
           pval = FALSE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           title="KM simple model")

km.model2 <- survfit(Surv(time, fup_status)~cohort)
summary(km.model2)
# plot(km.model2, conf.int = F, xlab="Time (weeks)", ylab = "ALive = S(t)", main="KM Model wrt Cohort", col=c('red', 'blue'))
# plot(km.model2, conf.int = T, xlab="Time (weeks)", ylab = "ALive = S(t)", main="KM Model wrt Cohort", col=c('red', 'blue'))
plots[[2]] <- ggsurvplot(km.model2,
           data= MDS_data,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           title = "KM wrt cohort")

arrange_ggsurvplots(plots)
graphics.off()

### KM vs genomic_groups
x11()
km.genomic.model <- survfit(Surv(time, fup_status)~genomic_group)
ggsurvplot(km.genomic.model,
                         data= MDS_data,
                         pval = TRUE, conf.int = FALSE,
                         risk.table = TRUE, # Add risk table
                         risk.table.col = "strata", # Change risk table color by groups
                         linetype = "strata", # Change line type by groups
                         surv.median.line = "hv", # Specify median survival
                         ggtheme = theme_bw(), # Change ggplot2 theme
                         title = "KM wrt genomic group")
graphics.off()

### KM vs IPSSR
x11()
km.IPSSR.model <- survfit(Surv(time, fup_status)~IPSSR)
ggsurvplot(km.IPSSR.model,
           data= MDS_data,
           pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           title = "KM wrt IPSSR")
graphics.off()

## KM (DNH vs HSCT) for each genomic group
plots <- list()
km.genomic.model <- list()
km.genomic.model.diff <- list()

for (i in 0:(length(levels(as.factor(genomic_group)))-1)) {
  print(paste("Group:", i))
  
  km.genomic.model[[i + 1]] <- survfit(Surv(time[which(genomic_group == i)], fup_status[which(genomic_group == i)])
                              ~cohort[which(genomic_group == i)], type = "kaplan-meier")
  
  km.genomic.model.diff[[i + 1]] <- survdiff(Surv(time[which(genomic_group == i)], fup_status[which(genomic_group == i)])
           ~cohort[which(genomic_group == i)])
  print(km.genomic.model.diff[[i + 1]])
  plots[[i + 1]] = ggsurvplot(km.genomic.model[[i + 1]],
                           data = MDS_data,
                           pval = TRUE, conf.int = FALSE,
                           risk.table = TRUE, # Add risk table
                           risk.table.col = "strata", # Change risk table color by groups
                           linetype = "strata", # Change line type by groups
                           surv.median.line = "hv", # Specify median survival
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           legend.title = "Cohort",
                           legend.labs = c("DNH", "HSCT"),
                           title = paste("KM of genomic group", toString(i), "vs cohort")
                            )
}

arrange_ggsurvplots(plots)
graphics.off()

## KM (DNH vs HSCT) for each IPSSR
plots <- list()
km.IPSSR.model <- list()
km.IPSSR.model.diff <- list()

for (i in 1:length(levels(as.factor(IPSSR)))) {
  print(paste("IPSSR:", i))
  
  km.IPSSR.model[[i]] <- survfit(Surv(time[which(IPSSR == i)], fup_status[which(IPSSR == i)])
                                    ~cohort[which(IPSSR == i)], type = "kaplan-meier")
  
  km.IPSSR.model.diff[[i]] <- survdiff(Surv(time[which(IPSSR == i)], fup_status[which(IPSSR == i)])
                                          ~cohort[which(IPSSR == i)])
  print(km.IPSSR.model.diff[[i]])
  plots[[i]] = ggsurvplot(km.IPSSR.model[[i]],
                              data = MDS_data,
                              pval = TRUE, conf.int = FALSE,
                              risk.table = TRUE, # Add risk table
                              risk.table.col = "strata", # Change risk table color by groups
                              linetype = "strata", # Change line type by groups
                              surv.median.line = "hv", # Specify median survival
                              ggtheme = theme_bw(), # Change ggplot2 theme
                              legend.title = "Cohort",
                              legend.labs = c("DNH", "HSCT"),
                              title = paste("KM of IPSSR", toString(i), "vs cohort")
  )
}

arrange_ggsurvplots(plots)
graphics.off()

x11()
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(as.factor(IPSSR) ~ cohort, main="IPSSR vs cohort")
plot(as.factor(IPSSR[which(cohort == "HSCT")]), main="IPSSR for HSCT")
plot(as.factor(IPSSR[which(cohort == "DNH")]), main="IPSSR for DNH")
print(paste("IPSSR for HSCT mean:", mean(IPSSR[which(cohort == "HSCT")], na.rm = TRUE)))
print(paste("IPSSR for DNH mean:", mean(IPSSR[which(cohort == "DNH")], na.rm = TRUE)))
graphics.off()

## KM (DNH vs HSCT) group 2 for each IPSSR
plots <- list()
km.IPSSR.model <- list()
km.IPSSR.model.diff <- list()

for (i in 3:length(levels(as.factor(IPSSR)))) {
  print(paste("IPSSR:", i))
  
  km.IPSSR.model[[i]] <- survfit(Surv(time[which(IPSSR == i & genomic_group == 2)], fup_status[which(IPSSR == i & genomic_group == 2)])
                                 ~cohort[which(IPSSR == i & genomic_group == 2)], type = "kaplan-meier")
  
  km.IPSSR.model.diff[[i]] <- survdiff(Surv(time[which(IPSSR == i & genomic_group == 2)], fup_status[which(IPSSR == i & genomic_group == 2)])
                                       ~cohort[which(IPSSR == i & genomic_group == 2)])
  print(km.IPSSR.model.diff[[i]])
  plots[[i-2]] = ggsurvplot(km.IPSSR.model[[i]],
                          data = MDS_data,
                          pval = TRUE, conf.int = FALSE,
                          risk.table = TRUE, # Add risk table
                          risk.table.col = "strata", # Change risk table color by groups
                          linetype = "strata", # Change line type by groups
                          surv.median.line = "hv", # Specify median survival
                          ggtheme = theme_bw(), # Change ggplot2 theme
                          legend.title = "Cohort",
                          legend.labs = c("DNH", "HSCT"),
                          title = paste("KM of IPSSR", toString(i), " for group 2 vs cohort")
  )
}

arrange_ggsurvplots(plots)
graphics.off()

x11()
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(as.factor(IPSSR) ~ cohort, main="IPSSR vs cohort")
plot(as.factor(IPSSR[which(cohort == "HSCT")]), main="IPSSR for HSCT")
plot(as.factor(IPSSR[which(cohort == "DNH")]), main="IPSSR for DNH")
print(paste("IPSSR for HSCT mean:", mean(IPSSR[which(cohort == "HSCT")], na.rm = TRUE)))
print(paste("IPSSR for DNH mean:", mean(IPSSR[which(cohort == "DNH")], na.rm = TRUE)))
graphics.off()
