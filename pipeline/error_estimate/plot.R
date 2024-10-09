library(ggplot2)
data <- read.csv("out.csv",  header = T)
data$concordanceQv <- as.numeric(data$concordanceQv)
data$concordance <- as.numeric(data$concordance)
data$is_corr <- as.factor(data$is_corr)
data$mismatchBp <- as.numeric(data$mismatchBp)
data$nonHpDeletionBp <- as.numeric(data$nonHpDeletionBp)
data$nonHpInsertionBp <- as.numeric(data$nonHpInsertionBp)
data$hpInsertionBp <- as.numeric(data$hpInsertionBp)
data$hpDeletionBp <- as.numeric(data$hpDeletionBp)

pdf("concordance.diff.pdf")
#concordance
ggplot(data,aes(concordance,fill=is_corr)) + 
	geom_density(alpha=0.5, position = "identity") + 
	scale_fill_manual(values = c("#a3cd5b","#0fb9b1")) + 
	xlim(0.8, 1.1)
dev.off()
# QV

pdf("concordanceQV.diff.pdf")
ggplot(data,aes(x=is_corr, y=concordanceQv,fill=is_corr)) + 
	geom_boxplot() +
	scale_fill_manual(values = c("#a3cd5b","#0fb9b1")) 
dev.off()

pdf("mismatchBp.diff.pdf")
ggplot(data,aes(x=is_corr, y=mismatchBp,fill=is_corr)) + 
	geom_boxplot() +
	scale_fill_manual(values = c("#a3cd5b","#0fb9b1")) 
dev.off()


pdf("nonHpInsertionBp.diff.pdf")
ggplot(data,aes(x=is_corr, y=nonHpInsertionBp,fill=is_corr)) + 
	geom_boxplot() +
	scale_fill_manual(values = c("#a3cd5b","#0fb9b1")) 
dev.off()

pdf("nonHpDeletionBp.diff.pdf")
ggplot(data,aes(x=is_corr, y=nonHpDeletionBp,fill=is_corr)) + 
	geom_boxplot() +
	scale_fill_manual(values = c("#a3cd5b","#0fb9b1")) 
dev.off()

pdf("hpInsertionBp.diff.pdf")
ggplot(data,aes(x=is_corr, y=hpInsertionBp,fill=is_corr)) + 
	geom_boxplot() +
	scale_fill_manual(values = c("#a3cd5b","#0fb9b1")) 
dev.off()


pdf("hpDeletionBp.diff.pdf")
ggplot(data,aes(x=is_corr, y=hpDeletionBp,fill=is_corr)) + 
	geom_boxplot() +
	scale_fill_manual(values = c("#a3cd5b","#0fb9b1")) 
dev.off()
