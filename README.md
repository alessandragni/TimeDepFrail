# TimeDepFrail: Time-Dependent Shared Frailty Cox Models in R

TimeDepFrail is the ultimate R package for fitting and analyzing Time-Dependent Shared Frailty Cox Models. 
These models extend the traditional Shared (Gamma) Frailty Cox Models by incorporating a time-dependent frailty component, 
making it a robust tool for studying how unexplained heterogeneity in data evolves over time.

This package implements the methods discussed in "Centre-Effect on Survival After Bone Marrow Transplantation: Application of Time-Dependent Frailty Models" by C.M. Wintrebert et al. (2004).


## Installation
You can install the development version of the package from `GitHub`:

```
devtools::install_github("alessandragni/TimeDepFrail")
```


## Dataset `data_dropout`
The `data_dropout` dataset is used to exemplify the package. 
It tracks the academic progress of students enrolled in 2012 over three academic years (six semesters). This dataset aims to explore the factors leading to student dropout.

The dataset is composed of four variables:

- `Gender`: Categorical covariate indicating gender (Male or Female).
- `CFUP`: Numeric covariate representing the standardized number of credits or CFUs (Credito Formativo Universitario) passed by the student in the first semester.
- `time_to_event`: The time (in semesters) when a student decides to drop out. A value greater than 6.0 means the student did not drop out during the follow-up period.
- `group`: Categorical variable representing the student's course of study, with 16 levels from CosA to CosP.

Students are followed for a maximum of 6 semesters (3 academic years), from the end of first semester until they drop out or the follow-up ends.


## Model execution
To fit a Time-Dependent Shared Frailty model, the following elements are required:

- data as `data.frame`, e.g. `data_dropout`
- `time_axis` vector: The time intervals for which the model is applied. For example, in the `data_dropout` dataset, no events occur in the first semester, so the `time_axis` starts at the end of the first semester (t = 1) and ends at the end of the third year (t = 6).
- `categories_range_min` and `categories_range_max` vectors: Provide minimum (`categories_range_min`) and maximum (`categories_range_max`) bounds for each parameter category to constrain the optimization.
- `formula` object: Specify the relationship between time-to-event, covariates, and group. For the clustering variable (`group`), it must be provided as `cluster(group)` in the formula.

Once these elements are prepared, you can call the desired model using the `AdPaikModel()` function. 

For full examples, refer to the `Examples/ReplicationCode.R` script.

Additionally, for guidance on selecting model parameters such as `time_axis`, `categories_range_min` and `categories_range_max`, 
we recommend base these choices on insights gained after fitting a Time-Unvarying Shared Frailty model. 
You can find a relevant example in the `Examples/ReplicationCode.R`, in the Appendix D Section.

## Analyzing results
Several built-in methods are available to analyze the results of the fitted model:

- Model Summaries: `summary()`, `plot()`, `print()`, `coef()`, `confint()`, `coefseAdPaik()`, `extractAIC()`, `logLik()`
- Baseline Hazard Step-Function: `bas_hazard()` and `plot_bas_hazard()`
- Frailty Standard Deviation/Variance: `frailty_sd()` and `plot_frailty_sd()`
- Posterior Frailty Estimates: `post_frailty_est()`, `plot_post_frailty_est()`, `post_frailty_var()`, `plot_post_frailty_var()`, `post_frailty_confint()`
- Conditional Survival Function: `survivalAdPaik()` and `plot_survivalAdPaik()`

These methods provide insightful visualizations and summaries to help you interpret your model results effectively.

Furthermore, also a support function suitable for the choice of the range of parameters and analysis of the 1D log-likelihood is available, `AdPaik_1D()`.


## To be aware of
- The `AdPaikModel` model is optimized for fast computation although estimating certain coefficients (e.g., `Male` versus `Female`) may vary slightly in computational time. 
  Note that changing the reference category (e.g., using `Male` as the baseline) alters the coefficient estimates but not the overall log-likelihood or model fit. 
  Users should choose reference categories based on interpretability rather than performance.


## Authors and maintainers of the code
Alessandra Ragni (alessandra.ragni@polimi.it),
Giulia Romani (giulia.romani@mail.polimi.it),
Chiara Masci (chiara.masci@polimi.it).

