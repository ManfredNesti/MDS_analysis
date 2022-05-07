library(survival)
library(survminer)

cox.model <- coxph(Surv(time, fup_status)~cohort, MDS_data)
cox.model
summary(cox.model)
# colnames(MDS_data)

covariates <- c("gender",  "cohort", "genomic_group", "IPSSR")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, fup_status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = MDS_data)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
# As seen from the results, p.value for cohort is very low 1.7e-09 thus indicating that it is not significant, 
# which is strange. Recall that with KM we have seen significant diff.

# Now let's do multivariate cox regression

res.cox <- coxph(Surv(time, fup_status) ~ gender + cohort + genomic_group, data = MDS_data)
summary(res.cox)
ggsurvplot(survfit(res.cox), data=MDS_data, palette = "#2E9FDF",
           ggtheme = theme_minimal())
# graphics.off()
plots <- list()
cohort_df_F <- with(MDS_data,
               data.frame(cohort = c("HSCT", "DNH"),
                          genomic_group = 0,
                          gender = "M",
                          # IPSSR = "very low" 
               )
)
cohort_df_M <- with(MDS_data,
                    data.frame(cohort = c("HSCT", "DNH"),
                               genomic_group = 2,
                               gender = "M",
                               # IPSSR = "very high" 
                    )
)
fit_F <- survfit(res.cox, newdata = cohort_df_F, data=MDS_data)
plots[[1]] <- ggsurvplot(fit_F, conf.int = TRUE, legend.labs=c("HSCT", "DNH"),
           ggtheme = theme_minimal())
fit_M <- survfit(res.cox, newdata = cohort_df_M, data=MDS_data)
plots[[2]] <- ggsurvplot(fit_M, conf.int = TRUE, legend.labs=c("HSCT", "DNH"),
           ggtheme = theme_minimal())
print(arrange_ggsurvplots(plots))

# test for proportional hazard

test.ph <- cox.zph(res.cox)
test.ph
ggcoxzph(test.ph)

