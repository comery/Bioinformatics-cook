```shell
input=$1
# -j 是线程数
seqkit fx2tab -j 4 -l  -n -i -H $input | cut -f 4 > Length.txt
# 查看Length.txt
head  Length.txt
```

```R
library(tidyverse)

length <- read_tsv("Length.txt") %>% group_by(length) %>%
  summarise(Count = n())
  length$length <- as.character(length$length)
  sum <- sum(length$Count)
  ggplot(length) + geom_col(aes(length, Count), width = 0.8) + 
    geom_line(aes(length, Count), group = 1) + geom_point(aes(length, Count)) + 
	  scale_y_continuous(sec.axis = sec_axis(~.*100/sum, name = "% Relative Abundance")) + xlab("Length") +
	    theme_bw() + theme(panel.grid = element_blank(), 
		                     axis.title = element_text(size = 15))

							 ggsave("Length.png", height = 5, width = 8)
							 ggsave("Length.pdf", height = 5, width = 8)
```
