if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BioQC")
library(BioQC)
install.packages("BayesFactor")
library(BayesFactor)
install.packages("dplyr")
library(dplyr)
install.packages("HDInterval")
library(HDInterval)
install.packages("ggplot2")
library(ggplot2)


# data prep
exp_data <- read.csv('data/X.csv', header=FALSE)
metadata <-read.csv('data/obs.csv')
var_data<- read.csv('data/var.csv')

rownames(exp_data)<-metadata$CellID
colnames(exp_data)<-var_data$Gene

# subset to cell lines
metadata_SE1 <- metadata[metadata$Celltype %in% c('Secretory Epithelial-1'),]
exp_data_SE1 <- exp_data[rownames(exp_data) %in% metadata_SE1$CellID,]

metadata_SE2 <- metadata[metadata$Celltype %in% c('Secretory Epithelial-2'),]
exp_data_SE2 <- exp_data[rownames(exp_data) %in% metadata_SE2$CellID,]

metadata_STIC <- metadata[metadata$Celltype %in% c('STIC lesion'),]
exp_data_STIC <- exp_data[rownames(exp_data) %in% metadata_STIC$CellID,]

# transform
exp_data_SE1 <- t(exp_data_SE1)
exp_data_SE2 <- t(exp_data_SE2)
exp_data_STIC <- t(exp_data_STIC)


# create compilation dataframes by sampling the subset expression data by cell line
# each df will take 100 samples of 500 genes
compilation_SE1 <- data.frame()

for (x in 1:100) {
  n <- sample(1:nrow(exp_data_SE1), 50, replace=FALSE)
  test <- exp_data_SE1[n,]
  entropy <- BioQC::entropy(test)
  samples <- append(n,entropy)
  compilation_SE1 <- rbind(compilation_SE1, samples)
}

compilation_SE2 <- data.frame()

for (x in 1:100) {
  n <- sample(1:nrow(exp_data_SE2), 50, replace=FALSE)
  test <- exp_data_SE2[n,]
  entropy <- BioQC::entropy(test)
  samples <- append(n,entropy)
  compilation_SE2 <- rbind(compilation_SE2, samples)
}

compilation_STIC <- data.frame()

for (x in 1:100) {
  n <- sample(1:nrow(exp_data_STIC), 50, replace=FALSE)
  test <- exp_data_STIC[n,]
  entropy <- BioQC::entropy(test)
  samples <- append(n,entropy)
  compilation_STIC <- rbind(compilation_STIC, samples)
}


colnames(compilation_SE1) <- c(paste0("f", 1:50), "Entropy")
boxplot(compilation_SE1$Entropy,
        main="Secretory Epithelial-1",
        ylab="Shannon Entropy")

colnames(compilation_SE2) <- c(paste0("f", 1:50), "Entropy")
boxplot(compilation_SE2$Entropy,
        main="Secretory Epithelial-2",
        ylab="Shannon Entropy")

colnames(compilation_STIC) <- c(paste0("f", 1:50), "Entropy")
boxplot(compilation_STIC$Entropy,
        main="STIC lesion",
        ylab="Shannon Entropy")

concat_df <- data.frame(Entropy = c(compilation_SE1$Entropy, compilation_SE2$Entropy, compilation_STIC$Entropy),
                        Group = c(rep('Secretory Epithelial-1', times=100), rep('Secretory Epithelial-2', times=100), rep('STIC lesion', times=100)))

concat_entropy <- ggplot(concat_df, aes(x=Group, y=Entropy, fill=Group)) +
  geom_violin(trim=F) +
  scale_color_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill='white') +
  labs(
    title = "Violin plot illustrating entropy across samples of cell types",
    x="Cell Type",
    y="Shannon Entropy"
  )
  
concat_entropy

# bayes
bayes_SE1 <-  ttestBF(x= compilation_SE1$Entropy, posterior = TRUE ,iterations = 1000)

bayes_STIC <-  ttestBF(x= compilation_STIC$Entropy, posterior = TRUE ,iterations = 1000)

bayes_SE1_STIC <-  ttestBF(x= compilation_STIC$Entropy, y=compilation_SE1$Entropy,
                           posterior = TRUE ,iterations = 1000, paired = TRUE)

post_SE1<- data.frame(mu = as.numeric(bayes_SE1[,"mu"]), Posterior = 'Secretory Epithelial-1')
post_STIC<- data.frame(mu = as.numeric(bayes_STIC[,"mu"]), Posterior = 'STIC lesion')
post_SE1_STIC<- data.frame(mu = as.numeric(bayes_SE1_STIC[,"mu"]),
                           Posterior = 'Secretory Epithelial-1 & STIC lesion')


ci95<-hdi(post_SE1_STIC)


plot_post_SE1 <- ggplot(post_SE1, aes(x=mu))+
  geom_histogram(color="darkblue", fill="lightblue") +
  ggtitle("Secretory Epithelial-1")

plot_post_STIC <- ggplot(post_STIC, aes(x=mu))+
  geom_histogram(color="red", fill="pink") +
  ggtitle("STIC lesion")

plot_post_SE1_STIC <- ggplot(post_SE1_STIC, aes(x=mu))+
  geom_histogram(color="darkgreen", fill="green") +
  geom_vline(xintercept = as.numeric(ci95[1:2]), linetype="dashed",
             color = "black", size=1, ) +
  ggtitle("Secretory Epithelial-1 & STIC lesion")

plot_post_SE1
plot_post_STIC
plot_post_SE1_STIC

install.packages("cowplot")
library(cowplot)

bayes_plot <- plot_grid(
  plot_post_SE1, plot_post_STIC, plot_post_SE1_STIC,
  ncol = 1,
  labels = 
)


posterior_distribution_plot <- ggdraw() +
  draw_label("Posterior Distributions of Entropy",
             fontface = 'bold', size = 16, x = 0.5, y = 0.95, hjust = 0.5) +
  draw_plot(bayes_plot, x=0, y=0, width=1, height=1)

posterior_distribution_plot <- plot_grid(
  ggdraw() + draw_label("Posterior Distributions of Entropy: SE-1, STIC",
                        fontface = 'bold', size = 16, x = 0.5, y = 0.6, hjust = 0.5),
  bayes_plot,
  ncol=1,
  rel_heights = c(2,15)
)

posterior_distribution_plot



# rerun, different choices
bayes_SE1 <-  ttestBF(x= compilation_SE1$Entropy, posterior = TRUE ,iterations = 1000)

bayes_SE2 <-  ttestBF(x= compilation_SE2$Entropy, posterior = TRUE ,iterations = 1000)

bayes_SE1_SE2 <-  ttestBF(x= compilation_SE2$Entropy, y=compilation_SE1$Entropy,
                           posterior = TRUE ,iterations = 1000, paired = TRUE)

post_SE1<- data.frame(mu = as.numeric(bayes_SE1[,"mu"]), Posterior = 'Secretory Epithelial-1')
post_SE2<- data.frame(mu = as.numeric(bayes_SE2[,"mu"]), Posterior = 'Secretory Epithelial-2')
post_SE1_SE2<- data.frame(mu = as.numeric(bayes_SE1_SE2[,"mu"]),
                           Posterior = 'Secretory Epithelial-1 & Secretory Epithalial-2')


ci95<-hdi(post_SE1_SE2)


plot_post_SE1 <- ggplot(post_SE1, aes(x=mu))+
  geom_histogram(color="darkblue", fill="lightblue") +
  ggtitle("Secretory Epithelial-1")

plot_post_SE2 <- ggplot(post_SE2, aes(x=mu))+
  geom_histogram(color="red", fill="pink") +
  ggtitle("Secretory Epithelial-2")

plot_post_SE1_SE2 <- ggplot(post_SE1_SE2, aes(x=mu))+
  geom_histogram(color="darkgreen", fill="green") +
  geom_vline(xintercept = as.numeric(ci95[1:2]), linetype="dashed",
             color = "black", size=1, ) +
  ggtitle("Secretory Epithelial-1 & Secretory Epithelial-2")

plot_post_SE1
plot_post_SE2
plot_post_SE1_SE2


bayes_plot <- plot_grid(
  plot_post_SE1, plot_post_SE2, plot_post_SE1_SE2,
  ncol = 1,
  labels = 
)


posterior_distribution_plot <- ggdraw() +
  draw_label("Posterior Distributions of Entropy",
             fontface = 'bold', size = 16, x = 0.5, y = 0.95, hjust = 0.5) +
  draw_plot(bayes_plot, x=0, y=0, width=1, height=1)

posterior_distribution_plot <- plot_grid(
  ggdraw() + draw_label("Posterior Distributions of Entropy: SE-1, SE-2",
                        fontface = 'bold', size = 16, x = 0.5, y = 0.6, hjust = 0.5),
  bayes_plot,
  ncol=1,
  rel_heights = c(2,15)
)

posterior_distribution_plot




