# ResultsHelper: Help Printing Results for APA Publications

Integrate statistical results into Markdown by referencing model objects; no need to precompute, just point to the model and parameter for the test.

# Installation

```{r}
install.packages("devtools") # if not installed, get devtools
devtools::install_github("https://github.com/Cmell/ResultsHelper")
```

# Example Use

```{r}
library(ResultsHelper)
library(lmSupport)

# a model to run

my_model = lm(
  mpg ~ wt * cyl, data = mtcars
)

summary(my_model)
```


## Results Section

Both weight, `r printStats('wt', my_model)`, and the number of cylinders, `r printStats('cyl', my_model)`, were significantly and inversely associated with fuel efficiency, controlling for the other. Heavier cars showed less fuel efficiency, as did cars with larger engines. In addition, an interaction emerged between the two, `r printStats('wt:cyl', my_model)`. 

## Customizing the Results String

See the details of `?printStats` for specifics. 

Change the format for a string:

Weight emerged as a significant predictor of fuel efficiency, `r printStats('wt', my_model, template='b_f_ci_p')`.
