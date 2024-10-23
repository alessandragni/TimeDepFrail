# TimeDepFrail: Time-Dependent Shared Frailty Cox Models in R

TimeDepFrail is the ultimate R package for fitting and analyzing Time-Dependent Shared Frailty Cox Models. These models extend the traditional Shared (Gamma) Frailty Cox Models by incorporating a time-dependent frailty component, offering powerful tools for studying how unexplained heterogeneity in data evolves over time.

This package implements the models discussed in the influential paper "Centre-Effect on Survival After Bone Marrow Transplantation: Application of Time-Dependent Frailty Models" by C.M. Wintrebert et al. (2013).


## Installation
Get started by installing the development version from `GitHub`:

```{r, eval=FALSE}
devtools::install_github("alessandragni/TimeDepFrail")
```


## Dataset data_dropout
To exemplify the package, the model are applied on a dataset called 'data_dropout'.

This dataset comes from a university administrative database and tracks students enrolled in 2012 over three academic years (or 6 semesters). We are interested in understanding what factors lead to students dropping out.
Dropout students with a time-instant in the first semester have been removed, for internal reasons (the university cannot take preventive action to reduce or avoid their withdrawal).

, extracted from an administrative database provided by a non-specified university. 
The students are followed for at most 3 academic years or, equivalently, 6 semesters (follow-up periods), from the first day of lecture up to the time-instant of withdrawal (i.e. survival event) or the end of the academic year. We need to specify that the 

The dataset is composed of four variables:
- Gender: categorical covariate (Male or Female).
- CFUP: standardized numerical covariate indicating the number of CFU (Credito Formativo Universitario) passed by the students by the end of the first semester. 
- time_to_event: the time (in semesters) when a student decides to leave the university. A value greater than 6.0 indicates the student did not drop out during the follow-up (e.g. 6.1 semesters).
-group: categorical variable indicating the student's course of study, with 16 different levels from CosA, CosB, ... , CosP.


## How to execute a model
To execute a model several variables need to be specified:
- dataset: e.g. 'data_dropout'
- time-domain vector: it can coincide with the follow-up or can be contained in it. For our dataset, no events are registered in the first semester and it would be useless to fit the model there. Thus, the time-domain starts at the end of the first semester (t=1) and it ends at the end of the third academic year (t=6). The temporal unit of measure is semester.
- categories_range_min and categories_range_max vectors: to maximize the log-likelihood function, we apply a constraint optimization method in multi-dimension and we need to provide the minimum and maximum existence range of each parameter (even though some parameters are not linked to any constraint). Since the latter may have the same meaning, we decide to group them in categories and to provide both the minimum and maximum range to each category.
According to the model, different categories are individuated.
- formula object: it specifies the relationship between the time-to-event, the dataset covariates and the group variable. The names reported in the formula must be contained also in the dataset, otherwise an error is thrown and the model stops the execution. 
Concerning the clustering variable (i.e. 'group'), it must be the arguments of a sort of cluster function: e.g. cluster(group).

Then we call the general method, specified which is the model we want to apply: 'AdPaikModel', 'PowParModel' or 'StocTimeDepModel'.
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
- AdPaikModel(), PowParModel(), StocTimeDepModel()
- AdPaik_1D(),PowPar_1D(), StocTimeDep_1D()
- TimeDepFrailty()
- summary(), summary.AdPaik(), summary.PowPar(), summary.StocTimeDep(...coming)
- frailty_sd.AdPaik()
- plot_bas_hazard(), plot_frailty_sd(), plot_post_frailty_est().


## To be aware of
- The first model is quite fast to optimize the log-likelihood function and to produce an output. However, considering the current dataset, it takes longer to estimated the 'Male' regression coefficient rather than the 'Female' one (reference). Therefore, I suggest to switch them and to set the 'Male' to the reference level and to estimate the 'Female' one. The estimated coefficient changes, but not the optimal log-likelihood value.
- The 'Centre-Specific Frailty Model with Power Parameter' is not as fast as the first one, but it produces coherent and expected results.
- The 'Stochastic Time-Dependent' Centre-Specific Frailty Model' is really slow and it may not reach convergence.


## Authors and maintener of the codes
Alessandra Ragni (alessandra.ragni@polimi.it),
Giulia Romani (giulia.romani@mail.polimi.it),
Chiara Masci (chiara.masci@polimi.it).

