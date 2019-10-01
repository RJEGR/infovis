---
title: "Exploratory data Analysis (EDA)"
output:
  html_document:
    df_print: paged
  html_notebook: default
  pdf_document: default
---

Sometext with hyperlinks [here](http://rmarkdown.rstudio.com) Notebook.

# Data Visualization

- Plots of data easily communicate information that is difficult to extract from tables of raw values.

- Data visualization is a key component of exploratory data analysis (EDA), in which the properties of data are explored through visualization and summarization techniques.

- Data visualization can help discover biases, systematic errors, mistakes and other unexpected problems in data before those data are incorporated into potentially flawed analysis.

- This course covers the basics of data visualization and EDA in R using the ggplot2 package and motivating examples from world health, economics and infectious disease.


```{r}
install.packages('dslabs')
library(dslabs)
data(murders)
head(murders)
```

# Introduction to distribution
- The most basic statistical summary of a list of objects is its distribution.

- In some cases, data can be summarized by a two-number summary: the average and standard deviation. We will learn to use data visualization to determine when that is appropriate.

# Data types


- Categorical data are variables that are defined by a small number of groups.

- Ordinal categorical data have an inherent order to the categories (mild/medium/hot, for example).

- Non-ordinal categorical data have no order to the categories.

- Numerical data take a variety of numeric values.

- Continuous variables can take any value.

- Discrete variables are limited to sets of specific values.

# Describing numerical data to someone
- A distribution is a function or description that shows the possible values of a variable and how often those values occur.

- For categorical variables, the distribution describes the proportions of each category.

- A frequency table is the simplest way to show a categorical distribution. Use prop.table to convert a table of counts to a frequency table. Barplots display the distribution of categorical variables and are a way to visualize the information in frequency tables.

```{r}
# load the dataset
library(dslabs)
data(heights)

# make a table of category proportions
freq.sex <- prop.table(table(heights$sex))
# and visualize
barplot(freq.sex)
```
- For continuous numerical data, reporting the frequency of each unique entry is not an effective summary as many or most values are unique. Instead, a distribution function is required.

- The cumulative distribution function (CDF) is a function that reports the proportion of data below a value  𝑎  for all values of  𝑎 :  𝐹(𝑎)=Pr(𝑥≤𝑎) .

- The proportion of observations between any two values  𝑎  and  𝑏  can be computed from the CDF as  𝐹(𝑏)−𝐹(𝑎) .

- A histogram divides data into non-overlapping bins of the same size and plots the counts of number of values that fall in that interval.


```{r}
x <- heights$height
cummulative_distribution <- ecdf(x)
plot(cummulative_distribution)
```



