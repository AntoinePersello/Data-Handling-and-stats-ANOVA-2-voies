# Library -----------------------------------------------------------------
library(tidyverse)
library(Rmisc)
library(multcomp)
library(boot)
library(car)
library(rstatix)
library(ggpubr)
library(emmeans)
library(agricolae)
library(DescTools)
library(FSA)
# Prep --------------------------------------------------------------------
PathWD = "C:\\Users\\Antoine\\OneDrive\\Papierss\\Traitements\\Scripts"
setwd(PathWD)
PathForDataImport = "C:\\Users\\Antoine\\OneDrive\\Papierss\\Traitements"
Data.Tab = read.csv(
  file = paste(PathForDataImport,
               "Table for R analysis .csv",
               sep = "/"),
  header = T,
  sep = ","
)
str(Data.Tab)# verify we have only numerical data

## First data filtering ----------------------------------------------------
Data.filter = Data.Tab |>
  dplyr::filter(c(Variable.type == "BGC" |
                    Variable.type == "Vitals"| 
                    Variable.type == "Bioch")) |>
  dplyr::filter(!(Variable ==  "BUN/Crea")) |>
  tidyr::pivot_longer(cols = dplyr::contains("_", ignore.case = T))
## Selection of variable to filter and filter ------------------------------
stopped_at <- NULL
VarList = unique(Data.filter$Variable)
for (i in 1:length(VarList)) {
  ToFilter = VarList[i]
  Data.Final = Data.filter |>
    dplyr::filter(Variable == ToFilter) |> # Second variable filter
    separate(name, into = c("Group", "number"), sep = "_")
  #---------Descriptive stats-----
  #aggregate( value ~ Group + Timing,data = Data.Final,FUN = mean) #mean/group only
  summary_data <-
    summarySE(Data.Final,
              measurevar = "value",
              groupvars = c("Timing", "Group"), na.rm = T) #mean se sd ci and group size by timing
  summary_data <- summary_data %>%
    mutate(lower_ci = value - ci,
           upper_ci = value + ci)
  if (length(unique(summary_data$Timing)) == 1) {
    dt = DunnettTest(
      x = Data.Final$value,
      g = Data.Final$Group,
      control = c("CEC", "Sham")
    )
    dtc = as.data.frame(dt$CEC)
    dts = as.data.frame(dt$Sham)
    dt.df=bind_rows(dtc, dts)
    write.csv(dt.df, paste(ToFilter, "KW-DT.csv", sep = "_"), row.names = F)
    rm(dtc, dts)
  } else {
    # Building linear model ---------------------------------------------------
    mod1 <-
      lm(
        value ~ Timing * Group,
        contrasts = list(Timing = contr.sum, Group = contr.sum),
        data = Data.Final
      )
    # 2 way ANOVA -------------------------------------------------------------
    #Note: aov() funciton is designed for balanced design (equal sample size) Anova function works fine with unbalanced design (unequal numbers of subkect per subgroup)
    #Note2: type II Anova is for not significant interaction and typeIII Anova is for significant interaction
    ANOVARES = Anova(mod1, type = 3) #type 3 because siginificant interaction (visible on plot draw before)
    #verification de : independance des residus, normalite et Homogeneites des variances
    #Data.Final$grp <- interaction(Data.Final$Timing, Data.Final$Group, sep="_")
    #Method2 ANOVA:
    mod = aov(value ~ Group * Timing, data = Data.Final) #the + sign is used to include independent variables without an interaction, the * sign is used to include independent variables with an interaction
    #-------Comparison 2 by 2------
    #Method 1:
    Timing_Group <- emmeans(mod1,  ~ Timing | Group)
    T_G = pairs(Timing_Group, adjust = "tukey")# 1) Between all Timing within each Group
    Group_Timing <-
      emmeans(mod1, ~ Group |
                Timing)# 2) Between all Groups for each Timing
    G_T = pairs(Group_Timing, adjust = "tukey")
    VarNumber = ifelse(length(unique(Data.Final$Timing)) > 3, yes = "Perop", no = "PreECC")
    TimingvsPerop = contrast(Timing_Group, method = "trt.vs.ctrl", ref = VarNumber)# 3) All timing compare to Start within each Group
    GroupVsCEC =  contrast(Group_Timing, method = "trt.vs.ctrl", ref = "CEC")# 4) Other groups vs CEC for each Timing
    # Save data ---------------------------------------------------------------
    write.csv(ANOVARES,
              paste(ToFilter, "resANOVA.csv", sep = "_"),
              row.names = F)
    write.csv(summary_data,
              paste(ToFilter, "ResTest.csv", sep = "_"),
              row.names = F)
    write.csv(as.data.frame(T_G),
              paste(ToFilter, "T_G.csv", sep = "_"),
              row.names = F)
    write.csv(as.data.frame(G_T),
              paste(ToFilter, "G_T.csv", sep = "_"),
              row.names = F)
    write.csv(
      as.data.frame(TimingvsPerop),
      paste(ToFilter, "Timingvs.csv", sep = "_"),
      row.names = F
    )
    write.csv(
      as.data.frame(GroupVsCEC),
      paste(ToFilter, "GroupVsCEC.csv", sep = "_"),
      row.names = F
    )
  }
print(ToFilter)
}

# Test code -----------------------------------------------------------

