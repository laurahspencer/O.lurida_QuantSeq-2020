# read in data 
histology <- read.csv("raw-data/histology_data.csv", header=T, stringsAsFactors = T, na.strings = "NA") %>%
    mutate_at(c("Female.Stage", "Male.Stage", "Sex.redo", "Dom.Stage.redo"), funs(factor(.))) %>% droplevels()

# Relevel a couple columns 
histology$Sex.redo <- factor(histology$Sex.redo, levels=c("I", "M", "HPM", "H", "HPF", "F"))
histology$PCO2 <- factor(histology$PCO2, levels = c("Pre","High","Amb")) 

# make an empty dataframe for test statistic results 

# ------- Compare pH treatments - effect of pH? 
#histology <- histology %>% filter(TEMPERATURE==6)

# See if there are any differences between previous temperature treatments within pCO2 treatments. Can I merge temp treat? 
histology.treats <- histology %>% filter(!TREATMENT %in% c("6", "10")) %>% droplevels()
(CT.sex.treat <- table(histology.treats$TREATMENT, histology.treats$Sex.redo, histology.treats$POPULATION))
fisher.test(CT.sex.treat[c(1,3),,"HL"], simulate.p.value = T, B = 10000) 
fisher.test(CT.sex.treat[c(2,4),,"HL"], simulate.p.value = T, B = 10000) 
fisher.test(CT.sex.treat[c(1,3),,"NF"], simulate.p.value = T, B = 10000) 
fisher.test(CT.sex.treat[c(2,4),,"NF"], simulate.p.value = T, B = 10000) 
fisher.test(CT.sex.treat[c(1,3),,"SN"], simulate.p.value = T, B = 10000) 
fisher.test(CT.sex.treat[c(2,4),,"SN"], simulate.p.value = T, B = 10000) 

(CT.malestage.treat <- table(histology.treats$TREATMENT, histology.treats$Male.Stage, histology.treats$POPULATION))
fisher.test(CT.malestage.treat[c(1,3),,"HL"], simulate.p.value = T, B = 10000) 
fisher.test(CT.malestage.treat[c(2,4),,"HL"], simulate.p.value = T, B = 10000) 
fisher.test(CT.malestage.treat[c(1,3),,"NF"], simulate.p.value = T, B = 10000) 
fisher.test(CT.malestage.treat[c(2,4),,"NF"], simulate.p.value = T, B = 10000) 
fisher.test(CT.malestage.treat[c(1,3),,"SN"], simulate.p.value = T, B = 10000) 
fisher.test(CT.malestage.treat[c(2,4),,"SN"], simulate.p.value = T, B = 10000) 

(CT.femstage.treat <- table(histology.treats$TREATMENT, histology.treats$Female.Stage, histology.treats$POPULATION))
fisher.test(CT.femstage.treat[c(1,3),,"HL"], simulate.p.value = T, B = 10000) 
fisher.test(CT.femstage.treat[c(2,4),,"HL"], simulate.p.value = T, B = 10000) 
fisher.test(CT.femstage.treat[c(1,3),,"NF"], simulate.p.value = T, B = 10000) 
fisher.test(CT.femstage.treat[c(2,4),,"NF"], simulate.p.value = T, B = 10000) 
fisher.test(CT.femstage.treat[c(1,3),,"SN"], simulate.p.value = T, B = 10000) 
fisher.test(CT.femstage.treat[c(2,4),,"SN"], simulate.p.value = T, B = 10000) 

# CONCLUSION: NO DIFFERENCES BETWEEN 6C OR 10C WITHIN PCO2 TREATMENTS & POPULATIONS. MERGE!

# Prepare contingency tables 
CT.sex.pop <- table(histology$PCO2, histology$Sex.redo, histology$POPULATION)
CT.sex.pop.plots <- table(histology$PCO2, histology$Sex.redo, histology$POPULATION)
CT.domsex.stage.pop <- table(histology$PCO2, histology$Dom.Stage.redo, histology$POPULATION)
CT.malestage.pop <- table(histology$PCO2, histology$Male.Stage, histology$POPULATION)
CT.femstage.pop <- table(histology$PCO2, histology$Female.Stage, histology$POPULATION)

# How many samples?
CT.sex.pop[,,"HL"] %>% rowSums()
CT.sex.pop[,,"NF"] %>% rowSums()
CT.sex.pop[,,"SN"] %>% rowSums()


# Extraneous tables across all popualtions 
#CT.sex <- table(histology$PCO2, histology$Sex.redo)
#CT.domsex.stage <- table(histology$PCO2, histology$Dom.Stage.redo)
#CT.malestage <- table(histology$PCO2, histology$Male.Stage)
#CT.femstage <- table(histology$PCO2, histology$Female.Stage)

# ======== STATS FOR PAPER ============== 
# FOR EACH SET OF GONAD DATA (SEX RATIO, SPERM, EGGS) - ALPHA= 0.05/3=0.0167 - THREE COMPARISONS ARE BEING MADE FOR EACH POPULATION. 

# SEX RATIOS

# DABOB BAY
chisq.test(CT.sex.pop[-1,,"HL"], simulate.p.value = T, B = 10000) # ambient vs high? NO diff <--- 
chisq.test(CT.sex.pop[-2,,"HL"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? NO diff 
chisq.test(CT.sex.pop[-3,,"HL"], simulate.p.value = T, B = 10000)  #Pre vs. high?  NO diff

# FIDALGO BAY
chisq.test(CT.sex.pop[-1,,"NF"], simulate.p.value = T, B = 10000) # ambient vs high? NO diff <--- 
chisq.test(CT.sex.pop[-2,,"NF"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? YES diff 
chisq.test(CT.sex.pop[-3,,"NF"], simulate.p.value = T, B = 10000)  #Pre vs. high?  NO to diff

# OYSTER BAY
chisq.test(CT.sex.pop[-1,,"SN"], simulate.p.value = T, B = 10000) # ambient vs high? YES diff <--- 
chisq.test(CT.sex.pop[-2,,"SN"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? YES diff 
chisq.test(CT.sex.pop[-3,,"SN"], simulate.p.value = T, B = 10000)  #Pre vs. high?  NO diff

#CONCLUSION: OYSTER BAY SEX RATIO DIFFERENT BETWEEN PCO2 EXPOSURE. MORE FEMALES. FIDALGO BAY AND OYSTER BAY SEX RATIO BOTH CHANGED IN AMBIENT PCO2, BUT DABOB BAY DID NOT CHANGE.  

## SPERM DEVELOPMENT 

# Stage of sperm differ by pCO2, by each population 

# DABOB BAY
fisher.test(CT.malestage.pop[-1,-1,"HL"], simulate.p.value = T, B = 10000) # ambient vs high? NO diff <---  
fisher.test(CT.malestage.pop[-2,-1,"HL"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? NO diff 
fisher.test(CT.malestage.pop[-3,-1,"HL"], simulate.p.value = T, B = 10000)  #Pre vs. high? NO diff

# FIDALGO BAY
chisq.test(CT.malestage.pop[-1,-1,"NF"], simulate.p.value = T, B = 10000) # ambient vs high? NO diff <---  
chisq.test(CT.malestage.pop[-2,-1,"NF"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? NO diff 
chisq.test(CT.malestage.pop[-3,-1,"NF"], simulate.p.value = T, B = 10000)  #Pre vs. high?  NO diff

# OYSTER BAY
chisq.test(CT.malestage.pop[-1,-1,"SN"], simulate.p.value = T, B = 10000) # ambient vs high? NO diff <---  
chisq.test(CT.malestage.pop[-2,-1,"SN"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? NO diff 
chisq.test(CT.malestage.pop[-3,-1,"SN"], simulate.p.value = T, B = 10000)  #Pre vs. high?  NO diff

# CONCLUSION: no sign. differences in developmental status of sperm (when present). NOTE: I removed "stage 0" b/c that reflects NO sperm. 

## EGG DEVELOPMENT 

# DABOB BAY
fisher.test(CT.femstage.pop[-1,-1,"HL"], simulate.p.value = T, B = 10000) # ambient vs high? NO diff <--- 
fisher.test(CT.femstage.pop[-2,-1,"HL"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? NO diff 
fisher.test(CT.femstage.pop[-3,-1,"HL"], simulate.p.value = T, B = 10000)  #Pre vs. high? NO diff

# FIDALGO BAY
fisher.test(CT.femstage.pop[-1,-1,"NF"], simulate.p.value = T, B = 10000) # ambient vs high? NO diff <--- 
fisher.test(CT.femstage.pop[-2,-1,"NF"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? YES diff
fisher.test(CT.femstage.pop[-3,-1,"NF"], simulate.p.value = T, B = 10000)  #Pre vs. high?  NO diff

CT.femstage.pop[,-1,"NF"]

# OYSTER BAY
fisher.test(CT.femstage.pop[-1,-1,"SN"], simulate.p.value = T, B = 10000) # ambient vs high? NO diff <--- 
fisher.test(CT.femstage.pop[-2,-1,"SN"], simulate.p.value = T, B = 10000)  #Pre vs. ambient? NO diff 
fisher.test(CT.femstage.pop[-3,-1,"SN"], simulate.p.value = T, B = 10000)  #Pre vs. high?  NO diff

# CONCLUSION: no sign. differences in developmental status of eggs (when present) between acid. and control. NOTE: I removed "stage 0" b/c that reflects NO eggs.
# BUT IN Fidalgo stage of eggs after control treatment differed compared to pre-treatment. But how? 
CT.femstage.pop[,,"NF"] # <-- seems like Fidalgo oocytes developed in ambient, but didn't change much in acidification. 

### PLOTS 

# Rename columns 
colnames(CT.sex.pop.plots) <- c("Indeterminate", "Male", "Male dominant", "Hermaphroditic", "Female dominant", "Female")
#colnames(CT.sex.plots) <- c("Indeterminate", "Male", "Male dominant", "Hermaphroditic", "Female dominant", "Female")
#colnames(CT.domsex.stage) <- c("Empty/No Follicles (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned/Regressing (4)")
#colnames(CT.domsex.stage.pop) <- c("Empty/No Follicles (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned/Regressing (4)")

plot.cohort <- c("HL","NF","SN")
plot.names <- c("Dabob Bay","Fidalgo Bay", "Oyster Bay C1")

# ----  Gonad sex for each cohort 

#pdf(file="Figures/gonad-sex-by-cohort", height = 5.75, width = 7.6)
par(mfrow=c(1,3), oma=c(5,3,1,2), mar=c(0,3,5,0), mgp=c(2.6,0.6,0))
for (i in 1:3) {
    barplot(t(prop.table(CT.sex.pop.plots[,,plot.cohort[i]], 1)), xlab=F, las=1, col=c("#f7f7f7", "#67a9cf", "#d1e5f0","gray85", "#fddbc7","#ef8a62" ), cex.lab=1.4, cex.axis = 1.2, col.axis = "gray30", col.lab = "gray30", legend.text = F)
    title(plot.names[i], line = 2, cex.main=1.5, col.main = "gray30", font.main = 1)
}
mtext(side=1,text=(expression(paste(pCO[2], " treatment"))), outer=T,line=3.5, col="gray30", font=1, cex=1.1)
mtext(side=2,text="Proportion Sampled", outer=T,line=0, col="gray30", font=1, cex=1, at=0.5)
mtext(side=3,outer=T,line=-1.5, col="gray30", font=3, cex=1.2, text=expression(paste("Gonad sex ratio by cohort & ", pCO[2])))
#dev.off()

# ---- Sperm stage by each cohort, treatment

#pdf(file="Figures/male-gonad-stage-by-cohort", height = 5.75, width = 7.6)
par(mfrow=c(1,3), oma=c(5,3,0,2), mar=c(0,3,5,0), mgp=c(2.6,0.6,0))
for (i in 1:3) {
    barplot(t(prop.table(CT.malestage.pop[,,plot.cohort[i]], 1)), xlab=F, las=1, col=c("#f7f7f7", "#cccccc", "#636363", "#252525",  "#969696"), cex.lab=1.4, cex.axis = 1.2, col.axis = "gray30", col.lab = "gray30", legend.text = F)
    title(plot.names[i], line = 1, cex.main=1.5, col.main = "gray30", font.main = 1)

}
mtext(side=1,text=(expression(paste(pCO[2], " treatment"))), outer=T,line=3.5, col="gray30", font=1, cex=1.1)
mtext(side=2,text="Proportion Sampled", outer=T,line=0, col="gray30", font=1, cex=1, at=0.5)
mtext(side=3,outer=T,line=-2, col="gray30", font=3, cex=1.2, text=expression(paste("Sperm stage by cohort & ", pCO[2])))
#dev.off()

# Egg stage by each cohort, treatment

#pdf(file="Figures/female-gonad-stage-by-cohort", height = 5.75, width = 7.6)
par(mfrow=c(1,3), oma=c(5,3,0,2), mar=c(0,3,5,0), mgp=c(2.6,0.6,0))
for (i in 1:3) {
    barplot(t(prop.table(CT.femstage.pop[,,plot.cohort[i]], 1)), xlab=F, las=1, col=c("#f7f7f7", "#cccccc", "#636363", "#252525",  "#969696"), cex.lab=1.4, cex.axis = 1.2, col.axis = "gray30", col.lab = "gray30", legend.text = F)
    title(plot.names[i], line = 1, cex.main=1.5, col.main = "gray30", font.main = 1)
    
}
mtext(side=1,text=(expression(paste(pCO[2], " treatment"))), outer=T,line=3.5, col="gray30", font=1, cex=1.1)
mtext(side=2,text="Proportion Sampled", outer=T,line=0, col="gray30", font=1, cex=1, at=0.5)
mtext(side=3,outer=T,line=-2, col="gray30", font=3, cex=1.2, text=expression(paste("Egg stage by cohort ", pCO[2])))
#dev.off()

# Anything interesting in these scatter plots?  hard to say. 
ggplot(histology, aes(x=Male.Stage, y=Female.Stage, shape=Sex.redo, col=PCO2)) + geom_jitter(size=3) + facet_wrap(~POPULATION) + theme_minimal()

# Percentage HPF or F 
#Dabob Bay: 
(5+2)/(0+7+3+1+5+2) #pre = 39%
(1+3)/(2+5+0+1+1+3) #High = 33%
(2+1)/(1+4+1+3+2+1) #Amb = 25%

#Fidalgo Bay: 
(5+2)/(0+7+3+1+5+2) #pre = 39%
(1+3)/(2+5+0+1+1+3) #High = 33%
(2+1)/(1+4+1+3+2+1) #Amb = 25%


# ALL POPS COMBINED COMPARISONS AMONG PH 
CT.sex.ph <- table(histology$PCO2, histology$Sex.redo)
fisher.test(CT.sex.ph[-1,], simulate.p.value = T, B = 10000) # all pops combined - sex ratio? NO 

CT.malestage.ph <- table(histology$PCO2, histology$Male.Stage)
fisher.test(CT.malestage.pop[-1,], simulate.p.value = T, B = 10000) # all pops combined - male stage? YES 

CT.femstage.ph <- table(histology$PCO2, histology$Female.Stage)
fisher.test(CT.femstage.pop[-1,], simulate.p.value = T, B = 10000) # all pops combined - female stage? NO 


# FIGURES FOR FINAL EXAM 
# PREVALENCE OF FEMALES 
as.data.frame(rbind(
    prop.table(CT.sex.pop[,,"HL"], margin = 1)[,5:6] %>% rowSums() %>% 
        t() %>% as.data.frame()  %>% mutate(Population="Dabob Bay"),
    prop.table(CT.sex.pop[,,"NF"], margin = 1)[,5:6] %>% rowSums() %>% 
        t() %>% as.data.frame()  %>% mutate(Population="Fidalgo Bay"),
    prop.table(CT.sex.pop[,,"SN"], margin = 1)[,5:6] %>% rowSums() %>% 
        t() %>% as.data.frame()  %>% mutate(Population="Oyster Bay"))) %>%
    mutate_at(vars(Population), factor) %>% 
    pivot_longer(names_to = "Treatment", cols = c("Pre", "High", "Amb")) %>%
    mutate(time=str_replace(Treatment, "High|Amb", "Post-pH")) %>%
    mutate(Treatment=str_replace(Treatment, "Pre", "Amb")) %>% 
    rbind(., subset(., time=="Pre") %>% 
              mutate(Treatment=str_replace(Treatment, "Amb", "High"))) %>%
    
    ggplot(aes(x=fct_rev(time), y=value, fill=Treatment, group=Treatment)) + 
    geom_line(aes(col=Treatment)) +
    geom_point(color="black", size=3, pch=22, alpha=0.75) + 
    facet_wrap(~Population) +
    scale_fill_manual(values=c("blue","red"), name="Treatment", 
                      labels=c("Amb"="Control", "High"="Acidification", "Pre"="Pre-Treatment")) +
    scale_color_manual(values=c("blue","red"), guide="none") +
    theme_bw(base_size = 12) + 
    theme(plot.title = element_text(size = 14, hjust = 0, colour = "gray30"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.title.x = element_blank(), legend.position="right", 
          axis.text.x = element_text(size=8)) +  
    labs(title=(expression(paste("The proportion of females by ", pCO[2], 
                                 " exposure & population"))), 
         y=expression("Proportion females")) +
    ylim(c(0,1))

# All populations combined - female prevalence 
as.data.frame(rbind(
    prop.table(CT.sex.ph, margin = 1)[,5:6] %>% rowSums() %>% 
        t() %>% as.data.frame())) %>%
    pivot_longer(names_to = "Treatment", cols = c("Pre", "High", "Amb")) %>%
    mutate(time=str_replace(Treatment, "High|Amb", "Post-pH")) %>%
    mutate(Treatment=str_replace(Treatment, "Pre", "Amb")) %>% 
    rbind(., subset(., time=="Pre") %>% 
              mutate(Treatment=str_replace(Treatment, "Amb", "High"))) %>%
    
    ggplot(aes(x=fct_rev(time), y=value, fill=Treatment, group=Treatment)) + 
    geom_line(aes(col=Treatment)) +
    geom_point(color="black", size=3.5, pch=21, alpha=0.75) + 
    scale_fill_manual(values=c("#67a9cf","#ef8a62"), name="Treatment", 
                      labels=c("Amb"="Control", "High"="Acidification", "Pre"="Pre-Treatment")) +
    scale_color_manual(values=c("#67a9cf","#ef8a62"), guide="none") +
    theme_bw(base_size = 12) + 
    theme(plot.title = element_text(size = 14, hjust = 0, colour = "gray30"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.title.x = element_blank(), legend.position="right", 
          axis.text.x = element_text(size=8)) +  
    labs(title=(expression(paste("The proportion of females by ", pCO[2], 
                                 " exposure & population"))), 
         y=expression("Proportion females")) +
    ylim(c(0,1))

# All populations combined - sperm stage 


h <- histology %>% filter(TEMPERATURE==6) %>% 
    group_by(PCO2) %>% mutate(Male.Stage=as.numeric(Male.Stage)) %>% 
    summarise(mean=mean(Male.Stage), sd=sd(Male.Stage)) %>% ungroup() %>%
    mutate(SAMPLING=str_replace(PCO2, "High|Amb", "APRIL")) %>%
    mutate(SAMPLING=str_replace(SAMPLING, "Pre", "FEBRUARY")) %>% 
    rbind(., subset(.,SAMPLING=="FEBRUARY") %>% mutate(PCO2=str_replace(PCO2, "Amb", "High"))) %>%
    mutate(PCO2=as.factor(PCO2), SAMPLING=as.factor(SAMPLING))


#  try to do a sperm stage figure 
histology  %>% mutate(Male.Stage=as.numeric(Male.Stage)) %>%
    ggplot(aes(x=PCO2, y=Male.Stage, fill=PCO2)) + 
    geom_violin(alpha=0.5) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=0.3, binwidth=0.05, col="black") + #binwidth=1, 
    #facet_wrap(~Population, labeller = labeller(Population=c("HL"="Dabob Bay","NF"="Fidalgo Bay","SN"="Oyster Bay"))) +
    scale_fill_manual(values=c("gray75","#ef8a62","#67a9cf"), name="Treatment",
                      labels=c("Ambient"="Control", "Low"="Acidification", "Pre-pH"="Pre-Treatment")) +
    # scale_fill_manual(values=c("#67a9cf","#ef8a62"), name="Treatment",
    #                    labels=c("Ambient"="Control", "Low"="Acidification", "Pre-pH"="Pre-Treatment")) +
    #scale_color_manual(values=c("#67a9cf","#ef8a62"), guide="none") +
    theme_bw(base_size = 12) + 
    theme(plot.title = element_text(size = 14, hjust = 0, colour = "gray30"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.title.x = element_blank(), legend.position="right", axis.text.x = element_text(size=8)) +  
    labs(title=(expression(paste("Sperm stage by ", pCO[2], " exposure"))), 
         y=expression("Sperm Developmental Stage")) #+
#    stat_summary(fun.y=mean, geom="point", shape=22, size=5, color="black", fill="white")



# Try showing male & female 

as.data.frame(rbind(
    prop.table(CT.sex.pop[,,"HL"], margin = 1)[,5:6] %>% rowSums() %>% 
        t() %>% as.data.frame()  %>% mutate(Population="Dabob Bay", Sex="Female"),
    prop.table(CT.sex.pop[,,"HL"], margin = 1)[,2:3] %>% rowSums() %>% 
        t() %>% as.data.frame()  %>% mutate(Population="Dabob Bay", Sex="Male"),
        prop.table(CT.sex.pop[,,"NF"], margin = 1)[,5:6] %>% rowSums() %>% 
        t() %>% as.data.frame()  %>% mutate(Population="Fidalgo Bay", Sex="Female"),
    prop.table(CT.sex.pop[,,"NF"], margin = 1)[,2:3] %>% rowSums() %>% 
        t() %>% as.data.frame()  %>% mutate(Population="Fidalgo Bay", Sex="Male"),
    prop.table(CT.sex.pop[,,"SN"], margin = 1)[,5:6] %>% rowSums() %>% 
        t() %>% as.data.frame()  %>% mutate(Population="Oyster Bay", Sex="Female"),
    prop.table(CT.sex.pop[,,"SN"], margin = 1)[,3:4] %>% rowSums() %>% 
    t() %>% as.data.frame()  %>% mutate(Population="Oyster Bay", Sex="Male"))) %>%
    mutate_at(vars(Population), factor) %>% 
    pivot_longer(names_to = "Treatment", cols = c("Pre", "High", "Amb")) %>%
    #dplyr::rename(Sex=Var1) %>% dplyr::rename(Treatment=Var2) %>% 
    mutate(time=str_replace(Treatment, "High|Amb", "Post-pH")) %>%
    mutate(Treatment=str_replace(Treatment, "Pre", "Amb")) %>% 
    rbind(., subset(., time=="Pre") %>% 
              mutate(Treatment=str_replace(Treatment, "Amb", "High"))) %>%
    ggplot(aes(x=fct_rev(time), y=value, fill=Treatment, group=Treatment)) + 
    geom_line(aes(col=Treatment)) +
    geom_point(color="black", size=4.5, pch=21, alpha=0.75) + 
    facet_wrap(~Sex+Population) +
    scale_fill_manual(values=c("blue","red"), name="Treatment",
                      labels=c("Amb"="Control", "High"="Acidification", "Pre-pH"="Pre-Treatment")) +
    scale_color_manual(values=c("blue","red"), guide="none") +
    theme_bw(base_size = 12) + 
    theme(plot.title = element_text(size = 14, hjust = 0, colour = "gray30"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.title.x = element_blank(), legend.position="right", 
          axis.text.x = element_text(size=8)) +  
    labs(title=(expression(paste("The proportion of females by ", pCO[2], 
                                 " exposure & population"))), 
         y=expression("Proportion females")) +
    ylim(c(0,1))

# Sex, all populations 

colnames(CT.sex.ph) <- c("Indeterminate", "Male", "Male dominant", "Hermaphroditic", "Female dominant", "Female")

#par(oma=c(5,3,1,2), mar=c(0,3,5,0), mgp=c(2.6,0.6,0))
barplot(t(prop.table(CT.sex.ph, 1)), xlab=NA, las=1, col=c("#f7f7f7", "#67a9cf", "#d1e5f0","gray85", "#fddbc7","#ef8a62" ), 
        cex.lab=1.4, cex.axis = 1.2, col.axis = "gray30", col.lab = "gray30", legend.text = F)
#title("Gonad Sex", line = 2, cex.main=1.5, col.main = "gray30", font.main = 1)
# mtext(side=1,text=(expression(paste(pCO[2], " treatment"))), outer=T,line=3.5, col="gray30", font=1, cex=1.1)
# mtext(side=2,text="Proportion Sampled", outer=T,line=0, col="gray30", font=1, cex=1, at=0.5)
# mtext(side=3,outer=T,line=-1.5, col="gray30", font=3, cex=1.2, text=expression(paste("Gonad sex ratio by cohort & ", pCO[2])))