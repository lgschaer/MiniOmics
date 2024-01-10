
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")

library("KEGGREST")
library(tidyverse)
# Pathway IDs from: https://www.genome.jp/pathway/ec00362

test <- keggGet(c("ec00624",   #Polycyclic aromatic hydrocarbon degradation
                  "ec00362",   #Benzoate degradation
                  "ec00627",   #Aminobenzoate degradation
                  "ec00020"))  #TCA cycle

e1 <- as.data.frame(test[[1]]$ENZYME) %>% 
  mutate(EC_number = test[[1]]$ENZYME, 
         General_Pathway = test[[1]]$NAME)

e2 <- as.data.frame(test[[2]]$ENZYME) %>% 
  mutate(EC_number = test[[2]]$ENZYME, 
         General_Pathway = test[[2]]$NAME)

e3 <- as.data.frame(test[[3]]$ENZYME) %>% 
  mutate(EC_number = test[[3]]$ENZYME, 
         General_Pathway = test[[3]]$NAME)

e4 <- as.data.frame(test[[4]]$ENZYME) %>% 
  mutate(EC_number = test[[4]]$ENZYME, 
         General_Pathway = test[[4]]$NAME)

head(e1)

interesting_enzymes <- e1 %>%
  full_join(e2, by = c("EC_number", "General_Pathway")) %>%
  full_join(e3, by = c("EC_number", "General_Pathway")) %>%
  full_join(e4, by = c("EC_number", "General_Pathway")) %>%
  dplyr::select(c(EC_number, General_Pathway))
head(interesting_enzymes)