# Decrypting the molecular basis of cellular drug phenotypes by dose-resolved expression proteomics


A very interesting article that demonstrates the Omics (with a capital O) character of Proteomics is
[Decrypting the molecular basis of cellular drug phenotypes by dose-resolved expression proteomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC11919725/){.uri}, published in Nature Biotechnology in March 2024. Here, the authors measure the protein expression of more than 8000 protein groups after exposing cell lines to 144 different drugs for a period of 18 hours. They use this data to calculate dose response curves for all of the proteins.

Please familiarize yourself with the article, we'll try to reproduce some plots from it.

# Software needed

## File download software
For downloading the data you will need an FTP (file transfer protocol) client. For windows, you can obtain Filezilla client (take care not to install extra software that is sometimes offered in the installer). Macintosh users have a built-in FTP client under connect to server, prefix the URL with ftp://.

## R
Since you are reading this in Rmd, it should be obvious you need R and probably RStudio.
The following packages need to be installed:
drc for dose response curve calculation
tidyverse packages: readr, readxl, tidyr, dplyr, purrr, ggplot2


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

library(readr)
library(readxl)

library(tidyr)
library(dplyr)
library(purrr)

library(ggplot2)

library(drc)
```

# Data

We'd like to re-analyze some data from this article, to try to recreate some dose-response curves.
The data is stored in ProteomeXChange, but unfortunately the R package does not support the repository (MassIVE), so we'll need to do some manual work.

In the article, a reference to the dataset is available. Obtain the proteinGroups.txt file from the time dependent dataset.

A 500 Megabytes file with proteinGroups at an FDR of 0.01 can be found in the FTP server of MassIVE in this folder: "/v06/MSV000093659/search/Search results - Long-term performance test Jurkat"

Read the file into R. Rename the headers of the Gene names table to be more easily used in the rest of the script.

```{r read data}
doseresponse <- read_delim("proteinGroups_fdr0.01.txt") %>%
  rename(Genes = `Gene names`)
head(doseresponse)
```

The drug names are given as numbers, so we will need to obtain the names from the metadata data found in the MassIVE repository.

```{r read drug names}
drugs <- read_xlsx("MassIVE_MappingSheet.xlsx", sheet=2, skip=1, col_names = c("Drug", "Number")) %>% pull(Drug)
```

If you have a look at the dimensions of the dataset, it's clear that there are a massive amount of columns. We don't need all of them, use the tidyverse `select` function to obtain the colums that contain protein names, descriptions, and LFQ values.

```{r filter}
dmso <- doseresponse %>%
  select(Genes, starts_with("LFQ") & contains("DMSO"))

doseresponse <- doseresponse %>%
  select(c(Genes), starts_with("LFQ")) %>%
  select(-(contains("DMSO")))
```

### Tidying and normalization

Use dplyr to make a long-form dataset with the data points coming from the LFQ columns. We'll need to split the column names into their separate parts, this can be done with the `separate` command.

According to the methods, the protein values are divided by the mean DMSO expression value.

```{r tidy}

dmso_tidy <- dmso %>% pivot_longer(-Genes) %>% summarize(.by=Genes, mean_dmso = mean(value, na.rm=TRUE)) %>% select(Genes, mean_dmso)

tidy_doseresponse <- doseresponse %>%
  pivot_longer(starts_with("LFQ")) %>%
  separate(name, sep=" ", into=c(NA, NA, "compound", "dose")) %>%
  mutate(compound = drugs[as.numeric(compound)], dose = as.numeric(dose))

tidy_doseresponse %<>% left_join(dmso_tidy, by="Genes")

doseresponse_norm <- tidy_doseresponse %>%
  mutate(.by=c(Genes, compound), value = value / mean_dmso) %>%
  select(-mean_dmso)

```

We also need to add a zero-concentration value for the curves, so we have a common starting point.

```{r example plot}
zero_values <- doseresponse_norm %>% summarize(.by = c(Genes, compound), dose=0.1, value=1)

doseresponse_norm %<>% bind_rows(zero_values)
```

#### Example plot

Here is an example plot reproduction of the article.

```{r example plot}
doseresponse_norm %>%
  filter(Genes %in% c("BAG3", "DNAJB1", "DNAJB4", "HSPA1A;HSPA1B"), compound=="Geldanamycin") %>%
  ggplot(aes(x=dose, y=value, col=Genes)) +
    geom_smooth(method=drm, method.args=list(fct=L.4()), se=FALSE) +
    geom_point() +
    scale_x_log10()
```

## Assignment

Given this dataset, try to do the following:
- Reproduce some of the dose-response curves, or try to find some that have proteins that have relatively high or low values in the higher drug concentration regimes
- Create a heatmap for the proteins for one or more selected drugs, across the relative concentration range
- Perform hierarchical clustering separately for the clusters you observe, have a look over here: https://www.statology.org/hierarchical-clustering-in-r/
- use a GO annotation tool to find out if you observe some interesting enrichment in your clusters, and try to link the action of the drug(s) to the observed terms.
