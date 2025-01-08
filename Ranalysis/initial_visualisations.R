# Set up of data and required packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
install.packages('ggrepel')
install.packages("WebGestaltR")
install.packages("ggplot2")
install.packages("tidyr")


library(dplyr)
library(ggplot2)
library(tidyr)
library(limma)
library(ggrepel)
library(WebGestaltR)

exp_data <- read.csv('data/X.csv', header=FALSE)
metadata <-read.csv('data/obs.csv')
var_data<- read.csv('data/var.csv')

# Set matrix as readable dataframe with legible col and row names
rownames(exp_data)<-metadata$CellID
colnames(exp_data)<-var_data$Gene


# Visualise number of cells per patient to see over/under-represented patients

cells_per_patient <- metadata %>%
  group_by(Patient) %>%
  tally()

colnames(cells_per_patient) <- c('Patient', 'Total_Cells')

ggplot(cells_per_patient, aes(x = reorder(Patient, -Total_Cells), y = Total_Cells)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Total Number of Cells per Patient", x = "Patient ID", y = "Total Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Visualise number of cells per stage to identify over/under-represented stages
cells_per_stage <- metadata %>%
  group_by(Disease_stage) %>%
  tally()

colnames(cells_per_stage) <- c('Disease_stage', 'Total_Cells')

ggplot(cells_per_stage, aes(x = Disease_stage, y = Total_Cells, fill = Disease_stage)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Number of Cells per Disease Stage", x = "Disease Stage", y = "Total Number of Cells") +
  theme_minimal()


# Visualise cell types and disease stage associations, to investigate associations between cells types and disease
cell_counts <- metadata %>%
  group_by(Celltype, Disease_stage) %>%
  tally()

colnames(cell_counts) <- c('Celltype', 'Disease_stage', 'Cell_count')

ggplot(cell_counts, aes(x = Disease_stage, y = Cell_count, fill = Celltype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cell Type Distribution by Disease Stage", x = "Disease Stage", y = "Cell Count") +
  theme_minimal()
