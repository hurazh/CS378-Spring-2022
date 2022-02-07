orgs <- read.csv("organisms_genes.tsv", sep='\t')
head(orgs)
table(orgs$type)

library(ggplot2)
ggplot(orgs, aes(bp, coding, color=type)) +
      geom_point(size=3)
orgs[which.max(orgs$bp),]
mean(orgs$bp)

orgs_filter <- orgs[orgs$name != "Wheat",]
ggplot(orgs_filter, aes(bp, coding, color=type)) +
        geom_point(size=3)


orgs$log_bp <- log10(orgs$bp)
orgs$log_coding <- log10(orgs$coding)
library(dplyr)
orgs <- orgs %>% mutate(log_bp = log10(bp),
                        log_coding = log10(coding))
head(orgs)
ggplot(orgs, aes(log_bp, log_coding, color=type)) +
  geom_point(size=3)


orgs <- orgs %>% mutate(genes = coding + noncoding)
orgs %>%
  group_by(type) %>%
  summarize(n=n(),
            cr.avg=mean(coding/genes),
            cr.sd=sd(coding/genes))


orgs_filter %>%
  group_by(type) %>%
  summarize(basepairs=median(bp),
            min=min(bp),
            max=max(bp)) %>%
  ggplot(aes(type, basepairs, ymin=min, ymax=max)) +
          geom_pointrange()

ggplot(orgs_filter, aes(bp, coding, col=type)) +
        geom_point() + geom_smooth(se=FALSE, method="lm")
