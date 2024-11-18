#Auswertung qPCR
# LZ & HB 02.01.2024
setwd("/path/to/validation")


stars_gen <- function(x){
  stars <- c( "***", "**", "*","")
  vec <- c(0, 0.001, 0.01, 0.05,1.01)
  i <- findInterval(x, vec)
  stars[i]
}

library(readxl)
library(tidyverse)
library(readr)
library(ggpubr)

ZEB1 <- read_excel("Experiments_qPCR/2023-12-20 092926_ZEB1.xls", 
                   sheet = "Results", skip = 48, n_max = 89)

PDE10A <- read_excel("Experiments_qPCR/2023-12-20 122604_PDE10A.xls", 
                   sheet = "Results", skip = 48, n_max = 89)

KCNIP4 <- read_excel("Experiments_qPCR/2023-12-20 140254_KCNIP4.xls", 
                   sheet = "Results", skip = 48, n_max = 89)

CACNA1C <- read_excel("Experiments_qPCR/2023-12-20 153903_CACNA1C.xls", 
                   sheet = "Results", skip = 48, n_max = 89)

condition <- read_excel("Experiments_qPCR/condition.xlsx")


#identify outliers
table(CACNA1C$`Ct SD`>0.5) #no outlier
table(KCNIP4$`Ct SD`>0.5)  #no outlier
table(PDE10A$`Ct SD`>0.5)  #no outlier
table(ZEB1$`Ct SD`>0.5)    #no outlier


genes <- c("ZEB1", "CACNA1C", "KCNIP4","PDE10A")
ngenes <- 4

results <- data.frame(gene = genes,
                      W_3vs0 = rep(0, times = ngenes),
                      p_3vs0 = rep(0, times = ngenes),
                      W_3vsCtrl = rep(0, times = ngenes),
                      p_3vsCtrl = rep(0, times = ngenes),
                      W_0vsCtrl = rep(0, times = ngenes),
                      p_0vsCtrl = rep(0, times = ngenes))
rownames(results) <- genes


for (i in genes) {
  j <- get(i)
  j <- j[,c("Sample Name", "Target Name","Ct Mean")]
  j <- j[!is.na(j$`Ct Mean`),]
  j <- j[!duplicated(j),]
  j <- spread(j, key = "Target Name", value = "Ct Mean")
  
  j$delta_ct <- NA
  if(i=="CACNA1C"){j[,4] = j[,2]-j[,3]}
  if(i!="CACNA1C"){j[,4] = j[,3]-j[,2]}
  

  avgCT_ctrl <- mean(j$delta_ct[j$`Sample Name`%in%c("51","52","53","54","55","56")])
  j$ddCT <- j$delta_ct-avgCT_ctrl
  j$ddCT2 <- 2^(-j$ddCT)
  
  j <- merge(j, condition, by="Sample Name")
  
  #Wilcoxon-Tests
  test_3vs0 <- wilcox.test(j$ddCT2[j$condition=="3crit"],j$ddCT2[j$condition=="0crit"])
  test_3vsCtrl <- wilcox.test(j$ddCT2[j$condition=="3crit"],j$ddCT2[j$condition=="ctrl"])
  test_0vsCtrl <- wilcox.test(j$ddCT2[j$condition=="0crit"],j$ddCT2[j$condition=="ctrl"])

  results[i,2] <- test_3vs0$statistic
  results[i,3] <- test_3vs0$p.value
  results[i,4] <- test_3vsCtrl$statistic
  results[i,5] <- test_3vsCtrl$p.value  
  results[i,6] <- test_0vsCtrl$statistic
  results[i,7] <- test_0vsCtrl$p.value

  #Plot
  datplot <- data.frame(condition=c("control","0crit","3crit"),
                       mean=c(mean(j$ddCT2[j$condition=="ctrl"]),mean(j$ddCT2[j$condition=="0crit"]),mean(j$ddCT2[j$condition=="3crit"])))

  p <- ggbarplot(
    j, x = "condition", y = "ddCT2", 
    add = c("mean_se", "jitter"),
    fill= "condition") + ylim(0,2) +
    theme_classic() + ggtitle(i) + scale_fill_manual(values = c("gray60","gray30","gray90")) +
    xlab("condition") + ylab("2^(-delta delta CT)")+
    theme(text = element_text(size=20),axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'))
  
  ggsave(p,file=paste0("results/plots/",i,".pdf"))
  }

write_csv(results,"results_wilcoxtest.csv")

