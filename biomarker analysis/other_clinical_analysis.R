library(dplyr)
cancertype <- c("Blood", "Bone", "Brain", "Breast", "Cervix", "Colorectal", "Esophagus", "Head&Neck", "Kidney", "Liver", "Lung", "Nervous_system", "Ovary", "Pancreas", "Prostate", "Skin", "Stomach", "Uterus")
sample_data <- read.table(file = "ICGC_data_2017/Bladder/Bladder_sample.tsv", sep = '\t', header = TRUE)
sample_data <- sample_data %>% select(1, 4, 6)
sample_data$cancertype = "Bladder"
donor_data <- read.table(file = "ICGC_data_2017/Bladder/Bladder_donor.tsv", sep = '\t', header = TRUE)
spec_data <- read.table(file = "ICGC_data_2017/Bladder/Bladder_specimen.tsv", sep = '\t', header = TRUE)
processed_data <- merge(x = sample_data, y = donor_data, by = "icgc_donor_id", all.x = TRUE)
processed_data <- merge(x = processed_data, y = spec_data, by = "icgc_specimen_id", all.x = TRUE)
processed_data <- processed_data %>% select("icgc_sample_id", "cancertype", "donor_age_at_diagnosis", "donor_tumour_stage_at_diagnosis", "tumour_grade", "specimen_donor_treatment_type")

for (i in cancertype){
  sample_data <- read.table(file = paste0("ICGC_data_2017/", i, "/", i, "_sample.tsv"), sep = '\t', header = TRUE)
  sample_data <- sample_data %>% select(1, 4, 6)
  sample_data$cancertype = i
  donor_data <- read.table(file = paste0("ICGC_data_2017/", i, "/", i, "_donor.tsv"), sep = '\t', header = TRUE)
  spec_data <- read.table(file = paste0("ICGC_data_2017/", i, "/", i, "_specimen.tsv"), sep = '\t', header = TRUE)
  data <- merge(x = sample_data, y = donor_data, by = "icgc_donor_id", all.x = TRUE)
  data <- merge(x = data, y = spec_data, by = "icgc_specimen_id", all.x = TRUE)
  data <- data %>% select("icgc_sample_id", "cancertype", "donor_age_at_diagnosis", "donor_tumour_stage_at_diagnosis", "tumour_grade", "specimen_donor_treatment_type")
  processed_data <- rbind(processed_data, data)
}
View(processed_data)

samples <- read.table(file = "samples.csv", sep = ',', header = TRUE)
processed_data <- merge(x= samples, y= processed_data, by="icgc_sample_id", all.x=TRUE)
View(processed_data)
write.csv(processed_data, file="processed_data.csv")


processed_data <- read.csv("processed_data.csv")

View(unique(processed_data$specimen_donor_treatment_type))
treatment <- processed_data[!is.na(processed_data$specimen_donor_treatment_type) & processed_data$specimen_donor_treatment_type != "", ]
num <- treatment %>% group_by(cancertype, subtype, specimen_donor_treatment_type) %>% summarise(c = n())
View(num)
num <- num %>% mutate(treatment = specimen_donor_treatment_type)

for(i in unique(num$cancertype)){
  g <- num[num$cancertype == i, ]
  dir.create(i)
  for (j in unique(g$subtype)){
    f <- g[g$subtype == j, ]
    f <- f %>% ungroup() %>% mutate(treatment = factor(treatment, levels = c("surgery", "radiation therapy" , "other therapy", "no treatment" , "immunotherapy", "combined chemo+radiation therapy", "chemotherapy")),
                                    cumulative = cumsum(c),
                                    midpoint = cumulative - c / 2,
                                    label = paste0(treatment, " ", c ," ", sum(c)))

    ggplot(f, aes(x = 1, weight = c, fill = treatment)) +
      geom_bar(width = 1, position = "stack") +
      coord_polar(theta = "y") +
      geom_text(aes(x = 1.3, y = midpoint, label = label))+
      theme_nothing()

    ggsave(paste0(i,"/", j, ".jpg"))
  }
}
dir.create("all treatment")
for (i in num$cancertype){
g <- num[num$cancertype == i, ]
for (j in g$subtype){
dir.create(paste0("all treatment/", i, "-", j))
f <- g[g$subtype == j, ]
f <- f %>% ungroup() %>% mutate(treatment = factor(treatment, levels = c("surgery", "radiation therapy" , "other therapy", "no treatment" , "immunotherapy", "combined chemo+radiation therapy", "chemotherapy")),
                                cumulative = cumsum(c),
                                midpoint = cumulative - c / 2,
                                label = c)

ggplot(f, aes(x = 1, weight = c, fill = treatment)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") +
  geom_text(aes(x = 1.3, y = midpoint, label = label))+
  theme_void()
ggsave(paste0("all treatment/", i,"-", j,  "/", "in ", j, ".jpg"))

f1 <- g[g$subtype != j, ]
f1 <- f1 %>% group_by(treatment) %>% summarise(c2 = sum(c))
f1 <- f1 %>% ungroup() %>% mutate(treatment = factor(treatment, levels = c("surgery", "radiation therapy" , "other therapy", "no treatment" , "immunotherapy", "combined chemo+radiation therapy", "chemotherapy")),
                                cumulative = cumsum(c2),
                                midpoint = cumulative - c2 / 2,
                                label = c2)

ggplot(f1, aes(x = 1, weight = c2 , fill = treatment)) +
  geom_bar(width = 1,  position = "stack") +
  coord_polar("y") +
  geom_text(aes(x = 1.3, y = midpoint, label = label))+
  theme_void()
  #theme_nothing()
ggsave(paste0("all treatment/", i,"-", j,  "/", "not in ", j, ".jpg"))
  }
}
"Blood" "Esophagus" "kidney"
############ age ####################
View(processed_data)
for (i in unique(processed_data$cancertype)){
  dir.create(i)
  temp <- processed_data[processed_data$cancertype == i, ]
  for (j in unique(temp$subtype)){
  tmp <- temp[temp$subtype == j, ]
  ggplot(tmp, aes(x=donor_age_at_diagnosis)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666")
  ggsave(paste0(i, "/", j, ".jpg"))
  }
}
a = "Blood"
b = "C11"
age <- processed_data[processed_data$cancertype == a, ]
theage <- age[age$subtype == b, ]
the_age <- age[age$subtype != b, ]
ggplot(theage, aes(x=donor_age_at_diagnosis)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank())+
  labs(title = paste0(a, " samples in ", b))
ggsave(paste0(a,b,".jpg"))

ggplot(the_age, aes(x=donor_age_at_diagnosis)) + geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                 panel.background = element_blank())+
  labs(title = paste0(a, " samples not in ", b))
ggsave(paste0(a, "not", b, ".jpg"))



"Blood C7 age"
"Blood C10 age"
"Blood C11 age"
"Brain C1 age"
"Brain C5 age"










######################### Grade #############################
#tmp <- processed_data[processed_data$cancer_type == "Bladder", ]
tumor <- processed_data[!is.na(processed_data$tumour_grade) &(processed_data$tumour_grade == "1" |  processed_data$tumour_grade == "2"  |  processed_data$tumour_grade == "3" |  processed_data$tumour_grade == "4") , ]
write.csv(tumor, file="tumor_grade.csv")
tumor <- read.csv("tumor_grade.csv")
#tmp <- tumor[tumor$cancertype!=" ", ]
tmp <- tumor
library(ggplot2)
View(tmp)

temp <- tmp

num <- temp %>% group_by(tumour_grade, subtype, cancertype) %>% summarise(c = n())
num$tumour_grade <- sub("^", "grade", num$tumour_grade )
View(num)
library(waffle)
for (i in unique(num$subtype)){
  f <- num[num$subtype == i, ]
  f <- setNames(as.list(f$c), f$tumour_grade)
  f <- unlist(f)
  waffle(f, rows=5, size=0.6,
         colors=c("#44D2AC", "#E48B8B", "#B67093",
                  "#3A9ABD"),
         title="Age Groups bifurcation")
  ggsave(paste0(i, ".jpg"))
}
library(ggmap) # for theme_nothing
for(i in unique(num$cancertype)){
  g <- num[num$cancertype == i, ]
  dir.create(i)
  for (j in unique(g$subtype)){
    f <- g[g$subtype == j, ]
    #r <- setNames(as.list(f$c), f$tumour_grade)
    #r <- unlist(r)
    #  waffle(r, rows=5, size=0.6,
    #         colors=colScale,
    #         title="Age Groups bifurcation")
    #  ggsave(paste0(i,"/", j, ".jpg"))

    f <- f %>% ungroup() %>% mutate(tumour_grade = factor(tumour_grade, levels = c("grade4", "grade3", "grade2", "grade1")),
           cumulative = cumsum(c),
           midpoint = cumulative - c / 2,
           label = paste0(tumour_grade, " ", c ," ", sum(c)))

    ggplot(f, aes(x = 1, weight = c, fill = tumour_grade)) +
      geom_bar(width = 1, position = "stack") +
      coord_polar(theta = "y") +
      geom_text(aes(x = 1.3, y = midpoint, label = label))+
      theme_void()

    ggsave(paste0(i,"/", j, ".jpg"))
  }
}
dir.create("alL grade")
for (i in unique(num$cancertype)){
g <- num[num$cancertype == i, ]
for (j in unique(g$subtype)){
  dir.create(paste0("alL grade/", i, "-", j))
f <- g[g$subtype == j, ]
f <- f %>% ungroup() %>% mutate(tumour_grade = factor(tumour_grade, levels = c("grade4", "grade3", "grade2", "grade1")),
                                cumulative = cumsum(c),
                                midpoint = cumulative - c / 2,
                                label = c )

ggplot(f, aes(x = 1, weight = c, fill = tumour_grade)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") +
  geom_text(aes(x = 1.3, y = midpoint, label = label))+
  theme_void()
ggsave(paste0("alL grade/", i,"-", j,  "/", "in ", j, ".jpg", sep=""))

f <- g[g$subtype != j, ]
View(f)
f <- f %>% group_by(tumour_grade) %>% summarise(c2 = sum(c))
f <- f %>% ungroup() %>% mutate(tumour_grade = factor(tumour_grade, levels = c("grade4", "grade3", "grade2", "grade1")),
                                cumulative = cumsum(c2),
                                midpoint = cumulative - c2 / 2,
                                label =  c2 )

ggplot(f, aes(x = 1, weight = c2, fill = tumour_grade)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") +
  geom_text(aes(x = 1.3, y = midpoint, label = label))+
  theme_void()
ggsave(paste0("alL grade/", i,"-", j,  "/", "not in ", j, ".jpg"))
  }
}

bp<- ggplot(num, aes(x="", y=c, fill=tumour_grade))+ geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_brewer(palette="Dark2")+ facet_wrap(subtype ~.)
ggsave("grade.jpg")



#Breast C7



#brain grade









