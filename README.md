# Lecture and exercises on simulations in population genetics at University of Copenhagen 2023

*These lecture materials were originally part of the ["Summer course in analysis of high throughput data for population genetics 2023"](http://popgen.dk/popgen23/).*

------------------------------------------------------------------------

### You can find the slides [here](https://bodkan.quarto.pub/ku-popgen2023/) (lecture, exercises, bonus content).

### [Here](https://bodkan.quarto.pub/ku-popgen2023-onepage/) is a one-page version (useful for reference while doing the exercises).

### [Here](https://github.com/bodkan/ku-popgen2023/tree/main/solutions) are solutions to Exercises #1-#3 (#4 is a homework 'bonus').

------------------------------------------------------------------------

This README summarizes steps needed to set up your machine for the lecture and exercises. After you're done installing everything, make sure to [run a small testing simulation](#testing-the-setup) to know that everything works as needed.

**For course participants:** If you don't want to work on your own laptop, you can also use the [online RStudio server](http://emily.popgen.dk:3838/) provided by the course organizers.

------------------------------------------------------------------------

# Installation instructions

**In case you will be using RStudio (highly recommended), you should do [this](#workaround-for-an-rstudio-bug) first.** This fixes a questionable default setting of RStudio which not only complicates using Python from R, but can break reproducibility of your analyses.

## If you want to use the RStudio server provided for the course

**If you want to use** the [online RStudio server](http://emily.popgen.dk:3838/) provided by the course organizers, you should run this bit of code in the R console after you log in:

```
library(slendr)
setup_env(agree = TRUE)
```

This will automatically install and set up necessary Python modules. **The server is quite slow, so the process can easily take five or more minutes!**

## If you want to use your computer (the server can be slow)

You will need:

-   R version 4.x&mdash;installators for macOS and Windows are provided [here](https://cloud.r-project.org), Linux users will manage.
-   [RStudio](https://www.rstudio.com/products/rstudio/download/) (not crucial but highly recommended)

Getting *slendr* to work is critical. The whole lecture is dedicated to this package.

First, run this in your R console:

```         
install.packages("slendr")
```

Then load *slendr* itself and set it up by running this bit of code in your R console:

```         
library(slendr)
setup_env(agree = TRUE)
```

Finally, make sure you get a positive confirmation from the following check:

```         
check_env()
```

## Other R package dependencies

I will use some tidyverse packages for analysis and plotting. You can install them with:

```         
install.packages(c("dplyr", "ggplot2"))"
```

# Testing the setup

Copy the following script to your R session **after (!) you successfully installed _slendr_ and ran `setup_env()` as described above**.

```         
library(slendr)
init_env()

o <- population("outgroup", time = 1, N = 100)
b <- population("b", parent = o, time = 500, N = 100)
c <- population("c", parent = b, time = 1000, N = 100)
x1 <- population("x1", parent = c, time = 2000, N = 100)
x2 <- population("x2", parent = c, time = 2000, N = 100)
a <- population("a", parent = b, time = 1500, N = 100)

gf <- gene_flow(from = b, to = x1, start = 2100, end = 2150, rate = 0.1)

model <- compile_model(
  populations = list(a, b, x1, x2, c, o), gene_flow = gf,
  generation_time = 1, simulation_length = 2200
)

ts <- msprime(model, sequence_length = 10e6, recombination_rate = 1e-8)

ts
```

If this runs without error and you get a small summary table from the `ts` object, you're all set!

# Workaround for an RStudio bug

RStudio sometimes interferes with Python which is needed for simulations.

Go to `Tools` → `Global Options` in your RStudio and set the following options like this:

![](images/rstudio_setting.png)
