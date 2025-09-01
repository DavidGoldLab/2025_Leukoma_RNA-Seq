# saturation state results
library(tidyr)
library(ggplot2)
library(dplyr)
library(AICcmodavg)

# Data wrangling -----------------------------
orig_saturation <- read.csv("omega_11_16_2023.csv")
View(orig_saturation)

# Build name parts
splitDF <- as.data.frame(orig_saturation$Sample)
colnames(splitDF) <- "Name"
splitDF$Name_full <- ifelse(grepl("pore", splitDF$Name), splitDF$Name, paste0(splitDF$Name, "_over"))
splitDF$Name_toSplit <- splitDF$Name_full
splitDF <- separate(splitDF, col = Name_toSplit, into = c('Date','Bucket','Position'), sep = '_')
splitDF <- separate(splitDF, col = Position, into = c('Position','Replicate'), sep='\\.')
splitDF$Replicate[is.na(splitDF$Replicate)] <- "3"
splitDF$Replicate <- as.numeric(splitDF$Replicate)

splitDF$index <- paste0(rownames(splitDF))
orig_saturation$index <- paste0(rownames(splitDF))

DF_complete_saturation <- left_join(splitDF, orig_saturation, by = 'index')
View(DF_complete_saturation)

# Pooled buckets (C1 + C2 -> C)
DF_complete_saturation$Bucket_pool <- gsub('[[:digit:]]+', '', DF_complete_saturation$Bucket)

# Normalize dates
DF_complete_saturation$Date <- gsub('D1','D01', DF_complete_saturation$Date)
DF_complete_saturation$Date <- gsub('D4','D04', DF_complete_saturation$Date)
DF_complete_saturation$Date <- gsub('D011','D11', DF_complete_saturation$Date)
DF_complete_saturation$Date <- gsub('D018','D18', DF_complete_saturation$Date)
DF_complete_saturation$Date <- gsub('D046','D46', DF_complete_saturation$Date)

# Time for plotting
DF_complete_saturation$Time <- as.numeric(gsub('D', "", DF_complete_saturation$Date))

# Rename OA→EA and OASH→EASH (both pooled and individual buckets)
DF_complete_saturation <- DF_complete_saturation %>%
  mutate(Bucket = gsub("OA1","EA1", Bucket),
         Bucket = gsub("OA2","EA2", Bucket),
         Bucket = gsub("OASH1","EASH1", Bucket),
         Bucket = gsub("OASH2","EASH2", Bucket))
DF_complete_saturation <- DF_complete_saturation %>%
  mutate(Bucket_pool = gsub("OA","EA", gsub("OASH","EASH", Bucket_pool)))

# (Optional) lock order of Position
DF_complete_saturation$Position <- factor(DF_complete_saturation$Position, levels = c("over", "pore"))

# Plots --------------------

# Faceted boxplot with pooled replicates (points filled by Bucket_pool)
date_box_saturation <- ggplot(
  transform(DF_complete_saturation, Bucket_pool = factor(Bucket_pool, levels = c("C","SH","EA","EASH"))),
  aes(x = Position, y = OmegaAragonite)
) +
  facet_grid(~ Bucket_pool, scales = "fixed", space = "free") +
  theme(text = element_text(family = "Times New Roman", size = 15)) +
  geom_boxplot(color = '#5e5e5e', outlier.shape = NA) +
  xlab("Water Position") + ylab("Aragonite Saturation State (mg/L)") +
  geom_point(
    position = position_jitterdodge(jitter.width = 2, dodge.width = 0),
    pch = 21, aes(fill = Bucket_pool), size = 2, show.legend = TRUE
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  guides(fill = guide_legend(title = "Treatment"))
date_box_saturation

# Timeline (pore only), already pooled by Bucket_pool
timeline_plot_saturation <- ggplot(subset(DF_complete_saturation, Position %in% "pore"),
                                   aes(x = Time, y = OmegaAragonite)) +
  geom_smooth(method = lm, aes(color = Bucket_pool), se = FALSE, linewidth = 0.7) +
  geom_point(aes(fill = Bucket_pool), size = 2, colour = "#5e5e5e", shape = 21) +
  ylab("Aragonite Saturation State") +
  guides(color = "none", fill = guide_legend(title = "Treatment"), size = "none") +
  theme(text = element_text(family = "Times New Roman", size = 15))
timeline_plot_saturation

# Save key plots as PDFs (PDF-safe font override) -------------

date_box_saturation_safe      <- date_box_saturation      + theme(text = element_text(family = "sans"))
timeline_plot_saturation_safe <- timeline_plot_saturation + theme(text = element_text(family = "sans"))

pdf("date_box_omega.pdf", width = 11, height = 8.5, family = "sans")
print(date_box_saturation_safe)
dev.off()

pdf("timeline_plot_omega.pdf", width = 11, height = 8.5, family = "sans")
print(timeline_plot_saturation_safe)
dev.off()

# ANOVAs -------------------------------------

omega_anova <- read.csv("omega_09_07_2023.csv")
View(omega_anova)
str(omega_anova)

omega_anova$sed_type  <- as.factor(omega_anova$sed_type)
omega_anova$over_pH   <- as.factor(omega_anova$over_pH)
omega_anova$bucket_ID <- as.factor(omega_anova$bucket_ID)

two.way     <- aov(log(pore_sat) ~ sed_type + over_pH, data = omega_anova)
interaction <- aov(log(pore_sat) ~ sed_type * over_pH, data = omega_anova)
buckets     <- aov(log(pore_sat) ~ sed_type * over_pH + bucket_ID, data = omega_anova)

summary(two.way)
summary(interaction)
summary(buckets)

library(AICcmodavg)
model.set   <- list(two.way, interaction, buckets)
model.names <- c("two.way", "interaction", "buckets")
aictab(model.set, modnames = model.names)

# Tukey (as in your original)
TukeyHSD(interaction)
