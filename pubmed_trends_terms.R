# get pubmet trends!!!!
# https://guangchuangyu.github.io/cn/2017/11/pubmed-trend/

devtools::install_github('GuangchuangYu/yyplot')

require(yyplot)

term <- c('"H7N9"', '"H5N1"', '"RSV"')
pm <- pubmed_trend(term, year=2001:2014)
plot(pm)

pubmed_trend("Yu Guangchuang[Full Author Name]", 2010:2016)
