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
                      t_3vs0 = rep(0, times = ngenes),
                      p_3vs0 = rep(0, times = ngenes),
                      t_3vsCtrl = rep(0, times = ngenes),
                      p_3vsCtrl = rep(0, times = ngenes),
                      t_0vsCtrl = rep(0, times = ngenes),
                      p_0vsCtrl = rep(0, times = ngenes))
rownames(results) <- genes

group_means <- data.frame(gene = genes,
                         m3_3vs0 = rep(0, times = ngenes),
                         m0_3vs0 = rep(0, times = ngenes),
                         m3_3vsCtrl = rep(0, times = ngenes),
                         mC_3vsCtrl = rep(0, times = ngenes),
                         m0_0vsCtrl = rep(0, times = ngenes),
                         mC_0vsCtrl = rep(0, times = ngenes))

rownames(group_means) <- genes

for (i in genes) {
  j <- get(i)
  j <- j[,c("Sample Name", "Target Name","Ct Mean")]
  j <- j[!is.na(j$`Ct Mean`),]
  j <- j[!duplicated(j),]
  j <- spread(j, key = "Target Name", value = "Ct Mean")
  
  j$delta_ct <- NA
  j[,4] = j[,2]-j[,3]
  
  j <- merge(j, condition, by="Sample Name")
  
  #t-Tests
  test_3vs0 <- t.test(j$delta_ct[j$condition=="3crit"],j$delta_ct[j$condition=="0crit"])
  test_3vsCtrl <- t.test(j$delta_ct[j$condition=="3crit"],j$delta_ct[j$condition=="ctrl"])
  test_0vsCtrl <- t.test(j$delta_ct[j$condition=="0crit"],j$delta_ct[j$condition=="ctrl"])
  
  results[i,2] <- test_3vs0$statistic
  results[i,3] <- test_3vs0$p.value
  results[i,4] <- test_3vsCtrl$statistic
  results[i,5] <- test_3vsCtrl$p.value  
  results[i,6] <- test_0vsCtrl$statistic
  results[i,7] <- test_0vsCtrl$p.value
  
  group_means[i,2] <- test_3vs0$estimate[1]
  group_means[i,3] <- test_3vs0$estimate[2]
  group_means[i,4] <- test_3vsCtrl$estimate[1]
  group_means[i,5] <- test_3vsCtrl$estimate[2]
  group_means[i,6] <- test_0vsCtrl$estimate[1]
  group_means[i,7] <- test_0vsCtrl$estimate[2]
  
  #Plot
  datplot <- data.frame(condition=c("control","0crit","3crit"),
                       mean=c(mean(j$delta_ct[j$condition=="ctrl"]),mean(j$delta_ct[j$condition=="0crit"]),mean(j$delta_ct[j$condition=="3crit"])))

  
  p <- ggbarplot(
    j, x = "condition", y = "delta_ct", 
    add = c("mean_se", "jitter"),
    fill= "condition") + ylim(0,7) +
    theme_classic() + ggtitle(i) + scale_fill_manual(values = c("gray60","gray30","gray90")) +
    xlab("condition") + ylab("delta CT")+
    theme(text = element_text(size=20),axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  
  
  if(i=="CACNA1C") {
    p <- ggbarplot(
      j, x = "condition", y = "delta_ct", 
      add = c("mean_se", "jitter"),
      fill= "condition") + ylim(-4,1) +
      theme_classic() + ggtitle(i) + scale_fill_manual(values = c("gray60","gray30","gray90")) +
      xlab("condition") + ylab("delta CT")+
      theme(text = element_text(size=20),axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
    
     }
  
  ggsave(p,file=paste0("results/plots/",i,".pdf"))
  }

write_csv(results,"results_ttest.csv")

group_means$delta_3vs0 <- log2(2^-(group_means$m3_3vs0 - group_means$m0_3vs0))
group_means$delta_3vsCtrl <- log2(2^-(group_means$m3_3vsCtrl - group_means$mC_3vsCtrl))
group_means$delta_0vsCtrl <- log2(2^-(group_means$m0_0vsCtrl - group_means$mC_0vsCtrl))

