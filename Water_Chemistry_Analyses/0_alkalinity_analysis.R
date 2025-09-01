# Library ----------------------
library(tidyr)
library(ggplot2)
library(dplyr)
library(AICcmodavg)

# Data wrangling -----------------------------
orig_alkalinity <- read.csv("alkalinity_11_16_2023.csv")

# Delete sump samples
orig_alkalinity <- orig_alkalinity[-grep("sump", orig_alkalinity$Short.name),]

# Build name parts
splitDF <- as.data.frame(orig_alkalinity$Short.name)
colnames(splitDF) <- "Name"
splitDF$Name_full <- ifelse(grepl("pore", splitDF$Name), splitDF$Name, paste0(splitDF$Name, "_over"))
splitDF$Name_full <- make.names(splitDF$Name_full, unique = TRUE)

splitDF$Name_toSplit <- splitDF$Name_full
splitDF <- separate(splitDF, col = Name_toSplit, into = c('Date', 'Bucket', 'Position'), sep = '_')
splitDF <- separate(splitDF, col = Position, into = c('Position', 'Replicate'), sep='\\.')
splitDF$Replicate[is.na(splitDF$Replicate)] <- "3"
splitDF$Replicate <- as.numeric(splitDF$Replicate)

splitDF$index <- paste0(rownames(splitDF))
orig_alkalinity$index <- paste0(rownames(splitDF))

DF_complete <- left_join(splitDF, orig_alkalinity, by = 'index')
View(DF_complete)

# Pooled buckets (C1 + C2 -> C)
DF_complete$Bucket_pool <- gsub('[[:digit:]]+', '', DF_complete$Bucket)

# Normalize dates D1->D01 etc.
DF_complete$Date <- gsub('D1','D01', DF_complete$Date)
DF_complete$Date <- gsub('D4','D04', DF_complete$Date)
DF_complete$Date <- gsub('D011','D11', DF_complete$Date)
DF_complete$Date <- gsub('D018','D18', DF_complete$Date)
DF_complete$Date <- gsub('D046','D46', DF_complete$Date)

# Time for plotting
DF_complete$Time <- as.numeric(gsub('D', "", DF_complete$Date))

# Rename OA→EA and OASH→EASH
DF_complete <- DF_complete %>%
  mutate(Bucket = stringr::str_replace_all(Bucket, c("OA1"="EA1","OA2"="EA2","OASH1"="EASH1","OASH2"="EASH2")),
         Bucket_pool = stringr::str_replace_all(Bucket_pool, c("OA"="EA","OASH"="EASH")))

# (Optional) lock order of Position
DF_complete$Position <- factor(DF_complete$Position, levels = c("over", "pore"))

# Plots --------------------

# Faceted boxplot with pooled replicates (points filled by Bucket_pool)
date_box_alk <- ggplot(
  transform(DF_complete, Bucket_pool = factor(Bucket_pool, levels = c("C","SH","EA","EASH"))),
  aes(x = Position, y = Final.corrected.TA)
) +
  facet_grid(~ Bucket_pool, scales = "fixed", space = "free") +
  geom_boxplot(color = '#5e5e5e', outlier.shape = NA) +
  xlab("Water Position") + ylab("TA (umol/kg)") +
  theme(text = element_text(family = "Times New Roman", size = 15)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 2, dodge.width = 0),
    pch = 21, aes(fill = Bucket_pool), size = 2, show.legend = TRUE
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  guides(fill = guide_legend(title = "Treatment"))
date_box_alk

# Timeline (pore only), already pooled by Bucket_pool
timeline_plot_alk <- ggplot(subset(DF_complete, Position %in% "pore"),
                            aes(x = Time, y = Final.corrected.TA)) +
  geom_smooth(method = lm, aes(color = Bucket_pool), se = FALSE, linewidth = 0.7) +
  geom_point(aes(fill = Bucket_pool), size = 2, colour = "#5e5e5e", shape = 21) +
  ylab("TA (umol/kg)") +
  guides(color = "none", fill = guide_legend(title = "Treatment"), size = "none") +
  theme(text = element_text(family = "Times New Roman", size = 15))
timeline_plot_alk

# Save key plots as PDFs (PDF-safe font override) -------------

date_box_alk_safe      <- date_box_alk      + theme(text = element_text(family = "sans"))
timeline_plot_alk_safe <- timeline_plot_alk + theme(text = element_text(family = "sans"))

pdf("date_box_alk.pdf", width = 11, height = 8.5, family = "sans")
print(date_box_alk_safe)
dev.off()

pdf("timeline_plot_alk.pdf", width = 11, height = 8.5, family = "sans")
print(timeline_plot_alk_safe)
dev.off()

##### ANOVAS
alk_anova <- read.csv("alk_anova_09_07_2023.csv")
View(alk_anova)
str(alk_anova)

# Factors
alk_anova$sed_type  <- as.factor(alk_anova$sed_type)
alk_anova$over_pH   <- as.factor(alk_anova$over_pH)
alk_anova$bucket_ID <- as.factor(alk_anova$bucket_ID)
str(alk_anova)

# Models
two.way     <- aov(log(pore_TA) ~ sed_type + over_pH, data = alk_anova)
interaction <- aov(log(pore_TA) ~ sed_type * over_pH, data = alk_anova)
buckets     <- aov(log(pore_TA) ~ sed_type * over_pH + bucket_ID, data = alk_anova)

summary(two.way)
summary(interaction)
summary(buckets)

# AIC
model.set   <- list(two.way, interaction, buckets)
model.names <- c("two.way", "interaction", "buckets")
aictab(model.set, modnames = model.names)

# Tukey (as in your original)
TukeyHSD(buckets)
