library(ggplot2)
library(dplyr)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
imsk_go_df <- as.data.frame(read.table("IMSK_GO.csv",header = T, sep=",",stringsAsFactors=FALSE))
imsk_go_df$log10.p.value = -1 * imsk_go_df$log10.p.value
seleted_df = imsk_go_df %>% select(description, log10.p.value)
colnames(seleted_df)<-c("description", "-log10(p-value)")
library(ggpubr)
ggplot(seleted_df, aes(x = reorder(seleted_df$description, seleted_df$`-log10(p-value)`), y = seleted_df$`-log10(p-value)`, fill=seleted_df$`-log10(p-value)`)) + 
  geom_col(aes(fill = seleted_df$`-log10(p-value)`)) + 
  scale_fill_gradient2(name ="-log10(p-value)", low = "blue", high = "red") + #, mid = "white"
  coord_flip() + 
  labs(y= "-log10(p-value)", x="Enriched GO term") + theme(
    axis.title.x = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=12, face="bold")
  ) + font("xy.text", size = 10, color="black",  face = "bold")
ggsave("IMSK.svg", width =12, height = 8)