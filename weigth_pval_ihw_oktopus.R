# BiocManager::install("IHWpaper")
library(IHWpaper)

dff

ihwRes <- ihw(PValue ~ logFC,  data = dff, alpha = 0.1)

ihwResDf <- as.data.frame(ihwRes)

dff2$baseMeanGroup2 <- groups_by_filter(dff2$baseMeanGroup, 8)

ggplot(dff, aes(x=PValue)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ baseMeanGroup2, nrow = 2)


baseMeanGroup <- sapply( levels(dds$group), function(group) rowMeans( counts(dds,normalized=TRUE)[,dds$group == group, drop=F] ) )

baseMeanGroup %>%
  as_tibble(rownames = 'ID') %>%
  pivot_longer(cols = colnames(baseMeanGroup), names_to = 'group', values_to = 'baseMeanGroup') %>%
  filter(baseMeanGroup > 0) %>%
  left_join(mtd %>% distinct(group, .keep_all = T)) %>%
  right_join(dff, by = c("ID", "Tissue")) %>%
  distinct(baseMeanGroup, .keep_all = T) %>%
  drop_na() -> dff2

ihwRes <- ihw(PValue ~ baseMeanGroup,  data = dff2, alpha = 0.1)
