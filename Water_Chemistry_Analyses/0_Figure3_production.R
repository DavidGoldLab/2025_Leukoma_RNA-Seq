# ----- Libraries -----
library(tidyr)
library(dplyr)
library(ggplot2)

# =========================
# Helper: tidy sample names
# =========================
tidy_names <- function(df, sample_col) {
  splitDF <- as.data.frame(df[[sample_col]])
  colnames(splitDF) <- "Name"
  splitDF$Name_full <- ifelse(grepl("pore", splitDF$Name),
                              splitDF$Name,
                              paste0(splitDF$Name, "_over"))
  splitDF$Name_toSplit <- splitDF$Name_full
  splitDF <- separate(splitDF, col = Name_toSplit,
                      into = c("Date","Bucket","Position"), sep = "_")
  splitDF <- separate(splitDF, col = "Position",
                      into = c("Position","Replicate"), sep = "\\.")
  splitDF$Replicate[is.na(splitDF$Replicate)] <- "3"
  splitDF$Replicate <- as.numeric(splitDF$Replicate)
  splitDF$index <- paste0(seq_len(nrow(splitDF)))
  df$index <- paste0(seq_len(nrow(splitDF)))
  left_join(splitDF, df, by = "index")
}

# =========================
# Read & wrangle: pH
# =========================
pH_raw <- read.csv("pH_11_16_2023.csv")
pH_raw <- pH_raw[!grepl("sump", pH_raw$Sample), ]
DF_pH  <- tidy_names(pH_raw, "Sample")

# Pool buckets (C1+C2->C etc.) and normalize dates
DF_pH$Bucket_pool <- gsub("[[:digit:]]+", "", DF_pH$Bucket)
DF_pH$Date <- gsub("D1","D01", DF_pH$Date)
DF_pH$Date <- gsub("D4","D04", DF_pH$Date)
DF_pH$Date <- gsub("D011","D11", DF_pH$Date)
DF_pH$Date <- gsub("D018","D18", DF_pH$Date)
DF_pH$Date <- gsub("D046","D46", DF_pH$Date)
DF_pH$Time <- as.numeric(gsub("D","", DF_pH$Date))

# Rename OA→EA, OASH→EASH (both pooled & individual)
DF_pH <- DF_pH %>%
  mutate(Bucket_pool = recode(Bucket_pool, "OA"="EA", "OASH"="EASH"),
         Bucket = recode(Bucket,
                         "OA1"="EA1","OA2"="EA2","OASH1"="EASH1","OASH2"="EASH2"))

# =========================
# Read & wrangle: Alkalinity
# =========================
alk_raw <- read.csv("alkalinity_11_16_2023.csv")
alk_raw <- alk_raw[!grepl("sump", alk_raw$Short.name), ]
DF_alk  <- tidy_names(alk_raw, "Short.name")

DF_alk$Bucket_pool <- gsub("[[:digit:]]+", "", DF_alk$Bucket)
DF_alk$Date <- gsub("D1","D01", DF_alk$Date)
DF_alk$Date <- gsub("D4","D04", DF_alk$Date)
DF_alk$Date <- gsub("D011","D11", DF_alk$Date)
DF_alk$Date <- gsub("D018","D18", DF_alk$Date)
DF_alk$Date <- gsub("D046","D46", DF_alk$Date)
DF_alk$Time <- as.numeric(gsub("D","", DF_alk$Date))

DF_alk <- DF_alk %>%
  mutate(Bucket      = stringr::str_replace_all(Bucket,      c("OA1"="EA1","OA2"="EA2","OASH1"="EASH1","OASH2"="EASH2")),
         Bucket_pool = stringr::str_replace_all(Bucket_pool, c("OA"="EA","OASH"="EASH")))

# =========================
# Read & wrangle: Omega
# =========================
omega_raw <- read.csv("omega_11_16_2023.csv")
DF_om     <- tidy_names(omega_raw, "Sample")

DF_om$Bucket_pool <- gsub("[[:digit:]]+", "", DF_om$Bucket)
DF_om$Date <- gsub("D1","D01", DF_om$Date)
DF_om$Date <- gsub("D4","D04", DF_om$Date)
DF_om$Date <- gsub("D011","D11", DF_om$Date)
DF_om$Date <- gsub("D018","D18", DF_om$Date)
DF_om$Date <- gsub("D046","D46", DF_om$Date)
DF_om$Time <- as.numeric(gsub("D","", DF_om$Date))

DF_om <- DF_om %>%
  mutate(Bucket = gsub("OA1","EA1", Bucket),
         Bucket = gsub("OA2","EA2", Bucket),
         Bucket = gsub("OASH1","EASH1", Bucket),
         Bucket = gsub("OASH2","EASH2", Bucket),
         Bucket_pool = gsub("OA","EA", gsub("OASH","EASH", Bucket_pool)))

# =========================
# Standardize factors/order
# =========================
bucket_order  <- c("C","SH","EA","EASH")
position_order<- c("over","pore")

DF_pH$Bucket_pool   <- factor(DF_pH$Bucket_pool,   levels = bucket_order)
DF_alk$Bucket_pool  <- factor(DF_alk$Bucket_pool,  levels = bucket_order)
DF_om$Bucket_pool   <- factor(DF_om$Bucket_pool,   levels = bucket_order)

DF_pH$Position  <- factor(DF_pH$Position,  levels = position_order)
DF_alk$Position <- factor(DF_alk$Position, levels = position_order)
DF_om$Position  <- factor(DF_om$Position,  levels = position_order)

# =========================
# Build a single long table
# =========================
pH_long <- DF_pH %>%
  transmute(Bucket_pool, Position, Value = pH.corrected.TRIS, Measure = "pH")

alk_long <- DF_alk %>%
  transmute(Bucket_pool, Position, Value = Final.corrected.TA, Measure = "Alkalinity")

om_long <- DF_om %>%
  transmute(Bucket_pool, Position, Value = OmegaAragonite, Measure = "Omega")

all_long <- bind_rows(pH_long, alk_long, om_long)
all_long$Measure <- factor(all_long$Measure, levels = c("pH","Alkalinity","Omega"))

# =========================
# 3x4 matrix of boxplots with custom colors
# =========================

matrix_3x4 <- ggplot(all_long, aes(x = Position, y = Value)) +
  geom_boxplot(color = "#5e5e5e", outlier.shape = NA) +
  # jitter points, colored by treatment
  geom_point(
    aes(fill = Bucket_pool),
    position = position_jitter(width = 0.2, height = 0),
    shape = 21, size = 2, alpha = 0.8
  ) +
  facet_grid(Measure ~ Bucket_pool, scales = "free_y", space = "free_x") +
  xlab("Water Position") + ylab(NULL) +
  theme_bw(base_size = 10) +   # set global font size to 10 pt
  theme(
    legend.position = "none",               
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10)
  ) +
  # custom muted colors: C=red, SH=purple, EA=green, EASH=blue
  scale_fill_manual(values = c(
    "C"    = "#e41a1c",  # muted red
    "SH"   = "#984ea3",  # muted purple
    "EA"   = "#4daf4a",  # muted green
    "EASH" = "#377eb8"   # muted blue
  ))

# Show plot
matrix_3x4

# Export at Word-friendly size (fits text width ~6.9 in)
ggsave("combined_boxplots_3x4.pdf", matrix_3x4, 
       width = 6.9, height = 5.5, units = "in", dpi = 300)

ggsave("combined_boxplots_3x4.png", matrix_3x4, 
       width = 6.9, height = 5.5, units = "in", dpi = 300)
