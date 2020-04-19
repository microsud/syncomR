
# syncomR: Synthetic Microbial Community Data Exploration, Analysis and Visualization using R     


<font color="red">This package is in experimental stage and should be used with caution!</font>

## About

The syncomR package is developed to analyse and visualize temporal microbial community profiles. The development was done for and using the human gut synthetic microbial community designed to investigate ecophysiological properties of core human gut microbes by [Shetty SA et al., 2020. Minimalist Approach for Deciphering the Ecophysiology of Human Gut Microbes](https://github.com/microsud/syncomR/tree/master). The package requires `phyloseq-class` object of microbial community profile for most of the analysis aimed at the ecological aspects. In addition, package data includes metatranscriptomics data, metabolites data, and gut micorbiota module expression data. The focus is to analyse temporal dynamics of  microbial community, test changes in abundances, overall community variability, resistance and resilience properties of the microbial communities. The ecological stability properties are calculated using the modified R code from [Liu et al., 2018. Ecological Stability Properties of Microbial Communities Assessed by Flow Cytometry](http://msphere.asm.org/content/3/1/e00564-17), while community level temporal variability and succession is calculated with modified R code from [Guittar J et al., 2019. Trait-based community assembly and succession of the infant gut microbiome](https://www.nature.com/articles/s41467-019-08377-w).
Additional functions are provided for data handling and formatting. Visualization tools are provided to explore the microbial community dynamics, including area, bar and line plot for composition, ordinaion plots for visualizing movement of community in and out of reference phase.


## Installing syncomR   

Before intalling `syncomR`, you will need to install R (3.5.1 or later) and RStudio (1.1.463 or later). Windows users also have to install [RTools](https://cran.r-project.org/bin/windows/Rtools/). Then open RStudio as administrator and in the console panel run the commands below to install syncomR:  

``` 
install.packages("devtools")

devtools::install_github("microsud/syncomR")
```

## Author
Sudarshan A. Shetty (sudarshanshetty9[@]gmail[dot]com)

