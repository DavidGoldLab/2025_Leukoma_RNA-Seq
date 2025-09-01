# Load libraries
library(tidyr)
library(ggplot2)
library(dplyr)
library(AICcmodavg)

# Data wrangling ----------------------------
orig_pH <- read.csv("pH_11_16_2023.csv")
head(orig_pH)

orig_pH <- orig_pH[-grep("sump", orig_pH$Sample),]

# to mess with the names of the data...make new df
splitDF <- as.data.frame(orig_pH$Sample)

colnames(splitDF) <- "Name"
splitDF$Name_full <- ifelse(grepl("pore", splitDF$Name), 
                            paste0(splitDF$Name), paste0(splitDF$Name, "_over"))
splitDF$Name_toSplit <- splitDF$Name_full 
splitDF <- separate(splitDF, col = Name_toSplit, into = c('Date', 'Bucket', 'Position'), sep = '_')
splitDF <- separate(splitDF, col = Position, into = c('Position', 'Replicate'), sep='\\.')
splitDF$Replicate[is.na(splitDF$Replicate)] <- "3"
splitDF$Replicate <- as.numeric(splitDF$Replicate)

splitDF$index <- paste0(rownames(splitDF))
orig_pH$index <- paste0(rownames(splitDF))

# new data frame organized how we want:
DF_complete_pH <- left_join(splitDF, orig_pH, by = 'index')
View(DF_complete_pH)

## now make a new column with pooled buckets (C1+C2 = C)
DF_complete_pH$Bucket_pool <- gsub('[[:digit:]]+', '', DF_complete_pH$Bucket)

DF_complete_pH$Date <- gsub('D1','D01', DF_complete_pH$Date)
DF_complete_pH$Date <- gsub('D4','D04', DF_complete_pH$Date)
DF_complete_pH$Date <- gsub('D011','D11', DF_complete_pH$Date)
DF_complete_pH$Date <- gsub('D018','D18', DF_complete_pH$Date)
DF_complete_pH$Date <- gsub('D046','D46', DF_complete_pH$Date)

# to get the time variable for plotting values over time...
DF_complete_pH$Time <- DF_complete_pH$Date
DF_complete_pH$Time <- gsub('D', "", DF_complete_pH$Time)
DF_complete_pH$Time <- as.numeric(DF_complete_pH$Time)

### Rename to EA and EASH (pool + individual buckets)
DF_complete_pH <- DF_complete_pH %>%
  mutate(Bucket_pool = recode(Bucket_pool, "OA" = "EA", "OASH" = "EASH"))
DF_complete_pH <- DF_complete_pH %>%
  mutate(Bucket = recode(Bucket, "OA1" = "EA1", "OA2" = "EA2", "OASH1" = "EASH1", "OASH2" = "EASH2"))

# Data visualization for pH --------

# COLLAPSED REPLICATES: points are colored by pooled treatment (Bucket_pool), not by replicate Bucket
date_box_pH <- ggplot(
  transform(DF_complete_pH, Bucket_pool = factor(Bucket_pool, levels = c("C", "SH", "EA", "EASH"))),
  aes(x = Position, y = pH.corrected.TRIS)
) +
  facet_grid(~ Bucket_pool, scales = "fixed", space = "free") +
  theme(text = element_text(family = "Times New Roman", size = 15)) +
  geom_boxplot(color = '#5e5e5e', outlier.shape = NA) +
  xlab("Water Position") + ylab(bquote(pH[TOT])) +
  # NOTE: use pooled-treatment mapping here to avoid separating replicates
  geom_point(
    position = position_jitterdodge(jitter.width = 2, dodge.width = 0),
    pch = 21, aes(fill = Bucket_pool), size = 2, show.legend = TRUE
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  guides(fill = guide_legend(title = "Treatment"))

date_box_pH

# Timeline (already pooled by Bucket_pool via its aes and uses Position == "pore")
timeline_plot_pH <- ggplot(subset(DF_complete_pH, Position %in% "pore"),
                           aes(x = Time, y = pH.corrected.TRIS)) +
  geom_smooth(method = lm, aes(fill = Bucket_pool), se = FALSE) +
  ylab("pore water pH") +
  geom_point(aes(fill = Bucket_pool), size = 2, colour = "#5e5e5e", shape = 21) +
  guides(size = "none", fill = guide_legend(title = "Bucket")) +
  theme(text = element_text(family = "Times New Roman", size = 15))
timeline_plot_pH

# Save plots individually as PDFs (PDF-safe font override) ----------------

date_box_pH_safe      <- date_box_pH      + theme(text = element_text(family = "sans"))
timeline_plot_pH_safe <- timeline_plot_pH + theme(text = element_text(family = "sans"))

pdf("date_box_pH.pdf", width = 11, height = 8.5, family = "sans")
print(date_box_pH_safe)
dev.off()

pdf("timeline_plot_pH.pdf", width = 11, height = 8.5, family = "sans")
print(timeline_plot_pH_safe)
dev.off()

##### ANOVAS
# Using anova csv from desktop editing pH...---------
pH_anova <- read.csv("pH_09_06_2023.csv")
View(pH_anova)
str(pH_anova)

# MAKE FACTORS
pH_anova$sed_type  <- as.factor(pH_anova$sed_type)
pH_anova$over_pH   <- as.factor(pH_anova$over_pH)
pH_anova$bucket_ID <- as.factor(pH_anova$bucket_ID)

str(pH_anova)

## ANOVA MODEL 1 -------
two.way <- aov(log(pH_corrected_TRIS) ~ sed_type + over_pH, data = pH_anova) 
summary(two.way)

## ANOVA MODEL 2 -------
interaction <- aov(log(pH_corrected_TRIS) ~ sed_type * over_pH, data = pH_anova)
summary(interaction)

## ANOVA MODEL 3 -------
buckets <- aov(log(pH_corrected_TRIS) ~ sed_type * over_pH + bucket_ID, data = pH_anova)
summary(buckets)

# AIC calculation --------------------- 
model.set   <- list(two.way, interaction, buckets)
model.names <- c("two.way", "interaction", "buckets")
aictab(model.set, modnames = model.names)

### After AIC, need to do a Tukey HSD test
TukeyHSD(interaction)
