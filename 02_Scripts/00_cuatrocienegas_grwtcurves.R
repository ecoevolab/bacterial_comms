# Author: Mayra Beatriz Mendoza Velazquez
# Title: Cuatro ciénegas growth curves

library(readr)
library(ggplot2)
library(gridExtra)

# Load data.frame 
individual_growth_curves <- read_tsv("01_RawData/individual_strains_growth_curves_filtered.tsv")
individual_growth <- as.data.frame(individual_growth_curves)

# Subset the data.frame for each strain 
CH23_GVals <- subset(individual_growth, Cepa == "CH23")
CH29_GVals <- subset(individual_growth, Cepa == "CH29")
CH90_GVals <- subset(individual_growth, Cepa == "CH90")
CH99b_GVals <- subset(individual_growth, Cepa == "CH99b")
CH111_GVals <- subset(individual_growth, Cepa == "CH111")
CH149a_GVals <- subset(individual_growth, Cepa == "CH149a")
CH154a_GVals <- subset(individual_growth, Cepa == "CH154a")
CH161d_GVals <- subset(individual_growth, Cepa == "CH161d")
CH447_GVals <- subset(individual_growth, Cepa == "CH447")
CH450_GVals <- subset(individual_growth, Cepa == "CH450")

#### Growth curves ####

#### CH23 ####
# CH23 - 30° 
CH23_30 <- subset(CH23_GVals, temp == 30)
View(CH23_30)

CH23_30_plot <- ggplot(CH23_30, aes(x = hr,
                           y = OD_real,
                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH23 - 30°C)")


# CH23 - 37° 
CH23_37 <- subset(CH23_GVals, temp == 37)

CH23_37_plot <- ggplot(CH23_37, aes(x = hr,
                    y = OD_real,
                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH23 - 37°C)")

# CH23 - 42°C 
CH23_42 <- subset(CH23_GVals, temp == 42)

CH23_42_plot <- ggplot(CH23_42, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH23 - 42°C)")

grid.arrange(CH23_30_plot, 
             CH23_37_plot, 
             CH23_42_plot, 
             nrow = 1,
             top = "CH23 growth curves") 

#### CH29 ####
# CH29 - 30° 
CH29_30 <- subset(CH29_GVals, temp == 30)
View(CH29_30)

CH29_30_plot <- ggplot(CH29_30, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH29 - 30°C)")


# CH29 - 37° 
CH29_37 <- subset(CH29_GVals, temp == 37)

CH29_37_plot <- ggplot(CH29_37, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH29 - 37°C)")

# CH29 - 42°C 
CH29_42 <- subset(CH29_GVals, temp == 42)

CH29_42_plot <- ggplot(CH29_42, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH29 - 42°C)")


grid.arrange(CH29_30_plot, 
             CH29_37_plot, 
             CH29_42_plot, 
             nrow = 1,
             top = "CH29 growth curves") 

#### CH90 ####
# CH90 - 30° 
CH90_30 <- subset(CH90_GVals, temp == 30)
View(CH90_30)

CH90_30_plot <- ggplot(CH90_30, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH90 - 30°C)")


# CH90 - 37° 
CH90_37 <- subset(CH90_GVals, temp == 37)

CH90_37_plot <- ggplot(CH90_37, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH90 - 37°C)")

# CH90 - 42°C 
CH90_42 <- subset(CH90_GVals, temp == 42)

CH90_42_plot <- ggplot(CH90_42, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH90 - 42°C)")


grid.arrange(CH90_30_plot, 
             CH90_37_plot, 
             CH90_42_plot, 
             nrow = 1,
             top = "CH90 growth curves") 

#### CH23 ####
# CH23 - 30° 
CH23_30 <- subset(CH23_GVals, temp == 30)
View(CH23_30)

CH23_30_plot <- ggplot(CH23_30, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH23 - 30°C)")


# CH23 - 37° 
CH23_37 <- subset(CH23_GVals, temp == 37)

CH23_37_plot <- ggplot(CH23_37, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH23 - 37°C)")

# CH23 - 42°C 
CH23_42 <- subset(CH23_GVals, temp == 42)

CH23_42_plot <- ggplot(CH23_42, aes(x = hr,
                                    y = OD_real,
                                    group = rep, colour = rep)) +
  geom_line(na.rm = TRUE) +
  theme_classic() +
  xlab("Time (hr)") + ylab("Absorbance (OD600 nm)") +
  labs(title = "Growth curve (CH23 - 42°C)")


grid.arrange(CH23_30_plot, 
             CH23_37_plot, 
             CH23_42_plot, 
             nrow = 1,
             top = "CH23 growth curves") 