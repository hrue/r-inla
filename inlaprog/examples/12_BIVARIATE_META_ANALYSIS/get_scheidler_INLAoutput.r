library(xtable) 

invlogit <- function(x){
    return(exp(x)/(1+exp(x)))
}


# rho <- as.numeric(read.table("hyperparameter-user-scale-rho-intern-for-2diid/quantiles.dat")[c(5, 3,7)])
# sigma_sens <- as.numeric(read.table("hyperparameter-user-scale-log-precision-for-2diid-second-component/quantiles.dat")[c(5, 3,7)])
# sigma_spec <- as.numeric(read.table("hyperparameter-user-scale-log-precision-for-2diid-first-component/quantiles.dat")[c(5, 3,7)])
ct_tn <- as.numeric(read.table("fixed.effect_ct_tn/quantiles.dat")[c(5, 3,7)])
ct_tp <- as.numeric(read.table("fixed.effect_ct_tp/quantiles.dat")[c(5, 3,7)])
mr_tn <- as.numeric(read.table("fixed.effect_mr_tn/quantiles.dat")[c(5, 3,7)])
mr_tp <- as.numeric(read.table("fixed.effect_mr_tp/quantiles.dat")[c(5, 3,7)])
lag_tn <- as.numeric(read.table("fixed.effect_lag_tn/quantiles.dat")[c(5, 3,7)])
lag_tp <- as.numeric(read.table("fixed.effect_lag_tp/quantiles.dat")[c(5, 3,7)])


table <- rbind(c(invlogit(lag_tp), invlogit(lag_tn)), 
                c(invlogit(ct_tp), invlogit(ct_tn)),
                c(invlogit(mr_tp), invlogit(mr_tn)))
rownames(table) <- c("LAG", "CT", "MR")
                
xtable_scheidler <- xtable(table)

print(xtable_scheidler, digits=2, only.contents=T, include.colnames=F, include.rownames=T, hline.after=-1)
