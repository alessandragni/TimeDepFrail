# TimeDepFrail: Time-Dependent Shared Frailty Cox Models in R

TimeDepFrail is the ultimate R package for fitting and analyzing Time-Dependent Shared Frailty Cox Models. These models extend the traditional Shared (Gamma) Frailty Cox Models by incorporating a time-dependent frailty component, making it a robust tool for studying how unexplained heterogeneity in data evolves over time.

This package implements the methods discussed in "Centre-Effect on Survival After Bone Marrow Transplantation: Application of Time-Dependent Frailty Models" by C.M. Wintrebert et al. (2013).


## Installation
You can install the development version of the package from `GitHub`:

```{r, eval=FALSE}
devtools::install_github("alessandragni/TimeDepFrail")
```


## Dataset data_dropout
The `data_dropout` dataset is used to exemplify the package. 
It tracks the academic progress of students enrolled in 2012 over three academic years (six semesters). This dataset aims to explore the factors leading to student dropout.

The dataset is composed of four variables:
- `Gender`: Categorical covariate indicating gender (Male or Female).
- `CFUP`: Numeric covariate representing the standardized number of CFUs (Credito Formativo Universitario) passed by the student in the first semester.
- `time_to_event`: The time (in semesters) when a student decides to drop out. A value greater than 6.0 means the student did not drop out during the follow-up period.
- `group`: Categorical variable representing the student's course of study, with 16 levels from CosA to CosP.

Students are followed for a maximum of 6 semesters (3 academic years), from the start of lectures until they drop out or the follow-up ends.


## Model execution
To fit a Time-Dependent Shared Frailty model, the following elements are required:
- dataset: e.g. `data_dropout`
- time-domain vector: it can coincide with the follow-up or can be contained in it. For our dataset, no events are registered in the first semester and it would be useless to fit the model there. Thus, the time-domain starts at the end of the first semester (t=1) and it ends at the end of the third academic year (t=6). The temporal unit of measure is semester.
- categories_range_min and categories_range_max vectors: to maximize the log-likelihood function, we apply a constraint optimization method in multi-dimension and we need to provide the minimum and maximum existence range of each parameter (even though some parameters are not linked to any constraint). Since the latter may have the same meaning, we decide to group them in categories and to provide both the minimum and maximum range to each category.
According to the model, different categories are individuated.
- formula object: it specifies the relationship between the time-to-event, the dataset covariates and the group variable. The names reported in the formula must be contained also in the dataset, otherwise an error is thrown and the model stops the execution. 
Concerning the clustering variable (i.e. 'group'), it must be the arguments of a sort of cluster function: e.g. cluster(group).

Then we call the general method, specified which is the model we want to apply: 'AdPaikModel'. Code also include implementation of 'PowParModel' and 'StocTimeDepModel', but their efficiency could be improved and thus related functions are secondary.
The execution returns an S3 class object with all the provided information and computed variables.

Look at the 'Examples/ModelsApplication.R' script for an example on how to proceed.


## How to analyze the results
To analyze the results, we can take advantage of already implemented methods that plot:
- baseline hazard step-function
- frailty standard deviation (or variance)
- posterior frailty estimates
and we can call the summary method that prints the most relevant results of the model call.

Look at the 'Examples/ModelsApplication.R' script for an example on how to proceed.


## Which are the exportable functions?
Among all the functions saved in the workspace, some of them can be directly called by the user, while others cannot because they are internally used by other methods.
The callable ones are reported in the following list, without arguments for convenience:
- AdPaikModel(), (also available PowParModel(), StocTimeDepModel())
- AdPaik_1D(), (also available PowPar_1D(), StocTimeDep_1D())
- TimeDepFrailty()
- summary(), summary.AdPaik(), (also available summary.PowPar(), summary.StocTimeDep)
- frailty_sd.AdPaik()
- plot_bas_hazard(), plot_frailty_sd(), plot_post_frailty_est().


## To be aware of
- The AdPaikModel model is quite fast to optimize the log-likelihood function and to produce an output. However, considering the current dataset, it takes longer to estimate the 'Male' regression coefficient rather than the 'Female' one (reference). If reference category is switched, the estimated coefficient changes, but not the optimal log-likelihood value.
- The 'Centre-Specific Frailty Model with Power Parameter' is not as fast as the AdPaikModel, but it produces coherent and expected results.
- The 'Stochastic Time-Dependent' Centre-Specific Frailty Model' is really slow and it may not reach convergence.


## Authors and maintainers of the code
Alessandra Ragni (alessandra.ragni@polimi.it),
Giulia Romani (giulia.romani@mail.polimi.it),
Chiara Masci (chiara.masci@polimi.it).

