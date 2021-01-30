# error in gglayer (needed for yyplot:
# object 'setup_data_continuous_color' not found
# Error : unable to load R code in package ‘gglayer’

# get pubmet trends!!!!
# https://guangchuangyu.github.io/cn/2017/11/pubmed-trend/
# https://mp.weixin.qq.com/s/ne2AEUxmD7EDcWzh7i-jAQ

devtools::install_github('GuangchuangYu/yyplot', force = T, upgrade = F)
devtools::install_github('GuangchuangYu/gglayer', force = T, upgrade = F)


.git_packages <- c('gglayer', 'yyplot')
.inst <- .git_packages %in% installed.packages()
if (any(!.inst)) {   
  if (!require('devtools')) {
    install.packages("devtools", dep=TRUE, repos='http://cran.us.r-project.org')
  } else 
    sapply(c(.git_packages), install_github, username = 'GuangchuangYu')
    sapply(c(.git_packages), require, character.only = TRUE)
}


devtools::install_github("PhDMeiwp/basicPackages@master", force = TRUE)
# https://github.com/PhDMeiwp/basicPackages

basicPackages::install.yyplot()

require(yyplot)

term <- c('"H7N9"', '"H5N1"', '"RSV"')
pm <- pubmed_trend(term, year=2001:2014)
plot(pm)

pubmed_trend("Yu Guangchuang[Full Author Name]", 2010:2016)
