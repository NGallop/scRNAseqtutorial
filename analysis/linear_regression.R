secretorymeta<- subset(metadata, metadata$Celltype %in%c('Secretory Epithelial-1','Secretory Epithelial-2'))
secretoryexp<- subset(exp_data, rownames(exp_data) %in% secretorymeta$CellID)

secretorymeta$JUNB<-exp_data[row.names(secretoryexp),"JUNB"]
secretorymeta$IER2<-exp_data[row.names(secretoryexp),"IER2"]

filtered_meta <- secretorymeta %>%
  filter(JUNB >= 0 & IER2 >= 0)

model <- lm(JUNB ~ IER2, data = filtered_meta)
summary(model)

geneexpression <- ggplot(filtered_meta, aes(x = JUNB, y = IER2)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_minimal() +
  labs(title = "JUNB Expression vs.IER2 Expression",
       x = "JUNB Expression",
       y = "IER2 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


r_squared <- summary(model)$r.squared

geneexpression<- geneexpression + annotate("text", x = 2, y = 2.5, label = paste("R² = ", round(r_squared,3)), size=6)

stagegeneexpression <- geneexpression + geom_point(aes(color=Disease_stage))
stagegeneexpression

tissuegeneexpression <- geneexpression + geom_point(aes(color=Tissue))
tissuegeneexpression


