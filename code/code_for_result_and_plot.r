install.packages(c("pROC", "PropCIs","dplyr"))
#--------------- Figure S1 ------------------------
#######ROC curve plotting#######
library(pROC)
# The classifier was first validated on an independent set of 75 FNA samples with known pathological results
# Load the data
load(file = "../data/data.Rda")
ran_roc <- roc(data$True_Label, data$Prediction_Score)
# Calculate AUC and 95% CI
roc_ci <- ci.auc(ran_roc, conf.level = 0.95)
auc_text <- paste0("AUC: ", round(roc_ci[2]*100, 2), "%",
                   " (", round(roc_ci[1]*100, 2),"%", "-", round(roc_ci[3]*100, 2),"%", ")")
pdf("../results/validation_set_ROC.pdf", width = 6, height = 6)
plot(
  ran_roc,
  print.auc=F, print.auc.cex=1.5,
  auc.polygon=F, 
  grid = T,
  cex.lab=2, cex.axis=1.5,
  xlim=c(1,0),ylim = c(0, 1), asp = 1,
  max.auc.polygon=F,
)
text(x = 0.3, y = 0.3, labels = auc_text, cex = 1.2, col = "black")
dev.off()


#######calculating accuracy, sensitivity, specificity, PPV, and NPV with 95% Wilson confidence intervals
#tp is true positive, fp is false positive, fn is false negative, and tn is true negative
# in the entire cohort, the performance of ThyroScan is:
tp <- 152
fp <- 7
fn <- 6
tn <- 52
accuracy <- (tp + tn) / (tp + fp + fn + tn)
sensitivity <- tp / (tp + fn)  
specificity <- tn / (tn + fp)  
ppv <- tp / (tp + fp)      
npv <- tn / (tn + fn)         
# using prop.tes to calculate 95% Wilson confidence intervals
ci_acc <- prop.test(tp + tn, tp + fp + fn + tn, conf.level = 0.95)$conf.int
ci_sens <- prop.test(tp, tp + fn, conf.level = 0.95)$conf.int
ci_spec <- prop.test(tn, tn + fp, conf.level = 0.95)$conf.int
ci_ppv <- prop.test(tp, tp + fp, conf.level = 0.95)$conf.int
ci_npv <- prop.test(tn, tn + fn, conf.level = 0.95)$conf.int

# print results
cat("accuracy (95% CI):", round(accuracy*100,1), "[", round(ci_acc[1]*100,1), ",", round(ci_acc[2]*100,1), "]\n",
    "sensitivity (95% CI):", round(sensitivity*100,1), "[", round(ci_sens[1]*100,1), ",", round(ci_sens[2]*100,1), "]\n",
    "specificity (95% CI):", round(specificity*100,1), "[", round(ci_spec[1]*100,1), ",", round(ci_spec[2]*100,1), "]\n",
    "PPV (95% CI):", round(ppv*100,1), "[", round(ci_ppv[1]*100,1), ",", round(ci_ppv[2]*100,1), "]\n",
    "NPV (95% CI):", round(npv*100,1), "[", round(ci_npv[1]*100,1), ",", round(ci_npv[2]*100,1), "]\n")


#--------------- Figure 2 ------------------------
#######Given sensitivity and specificity, predict performance of Multigene Genomic Classifier across varying cancer prevalence rates
library(PropCIs)
library(dplyr)
# Se is sensitivity and Sp is specificity
# in the entire cohort, the performance of ThyroScan is:
Se <- 0.962
Sp <- 0.8813
FNR=1-Se
FPR=1-Sp
# Create a sequence of prevalence values from 0 to 1
prev <- seq(0, 1, by = 0.01)
# Compute PPV/NPV 95% CI for each prevalence
ppv_ci_list <- list()
npv_ci_list <- list()
for (prevalence in prev) {
  total_pos <- round(1000 * prevalence)  # Given a total sample size of 1000, calculate the number of malignant (pathologically positive) samples.
  total_neg <- round(1000 * (1-prevalence))  # Given a total sample size of 1000, calculate the number of benign (pathologically negative) samples.
  TP <- round(total_pos * Se)
  FP <- round(total_neg * FPR)
  TN <- round(total_neg * Sp)
  FN <- round(total_pos * FNR)
  # Using the exact Clopper-Pearson method
  ci <- exactci(TN, TN + FN, 0.95)
  npv_ci_list[[as.character(prevalence)]] <- c(
    prevalence = prevalence,
    NPV = TN / (TN + FN),
    lower = ci$conf.int[1],
    upper = ci$conf.int[2]
  )
  ci <- exactci(TP, TP + FP, 0.95)
  ppv_ci_list[[as.character(prevalence)]] <- c(
    prevalence = prevalence,
    PPV = TP / (TP + FP),
    lower = ci$conf.int[1],
    upper = ci$conf.int[2]
  )
}
# turn the list into a data frame
npv_ci_df <- do.call(rbind, npv_ci_list) %>% data.frame()
ppv_ci_df <- do.call(rbind, ppv_ci_list) %>% data.frame()
pdf("../results/predicted_curves_across_prevalance.pdf", width = 6, height = 6)
# Plot the base graph and PPV curve
plot(
  prev*100, ppv_ci_df$PPV*100, 
  type = "l", col = "red", lwd = 2, 
  xlab = "Prevalence of Malignancy (%)", ylab = "NPV, PPV (%)",
  # main = "PPV and NPV vs. Prevalence"
)
# Add NPV curve
lines(prev*100, npv_ci_df$NPV*100, col = "blue", lwd = 2)
# Draw grid lines confined within plot box
x_vals <- seq(0, 100, by = 10)  # Vertical grid lines at intervals of 0.1
y_vals <- seq(0, 100, by = 10)  # Horizontal grid lines at intervals of 0.2
# Add vertical grid lines
for (x in x_vals) {
  segments(x0 = x, y0 = 0, x1 = x, y1 = 100, col = "lightgray", lty = "dotted")
}
# Add horizontal grid lines
for (y in y_vals) {
  segments(x0 = 0, y0 = y, x1 = 100, y1 = y, col = "lightgray", lty = "dotted")
}
# Add 95% CI
lines(prev*100, ppv_ci_df$lower*100, col = "red", lwd = 2, lty =2)
lines(prev*100, ppv_ci_df$upper*100, col = "red", lwd = 2, lty =2)
lines(prev*100, npv_ci_df$lower*100, col = "blue", lwd = 2, lty =2)
lines(prev*100, npv_ci_df$upper*100, col = "blue", lwd = 2, lty =2)
# Add the legend at the bottom center
legend("bottom", inset = c(0, 0), legend = c("PPV_ThyroScan", "NPV_ThyroScan"), col = c("red", "blue"), lty = 1, lwd = 2, horiz = F, bty = "n", xpd = TRUE)
dev.off()