# Time-Dependent Shared Frailty Cox Models-R
This Github repository is based on the paper "Centre-Effect on Survival After Bone Marrow Transplantation: Application of Time-Dependent Frailty Models", written by C.M. Wintrebert, H. Putter, A.H. Zwinderman and J.C. van Houwelingen.

Briefly, this paper describes three different "Time-Dependent Shared Frailty Cox Models", that are elaborated and complex Cox regression models, at which a time-dependent random term (called frailty) is added for studying the temporal behaviour of a portion of data heterogeneity that cannot be solely explained by the regression coefficients. More details are provided in my thesis. 

## Dataset data_dropout
So far, these are the only available codes implementing the models and they are applied on a dataset called 'data_dropout', extracted from an administrative database provided by a non-specified university. 
The students are followed for at most 3 academic years or, equivalently, 6 semesters (follow-up periods), from the first day of lecture up to the time-instant of withdrawal (i.e. survival event) or the end of the academic year. We need to specify that the dropout students with a time-instant in the first semester have been removed, for internal reasons (the university cannot take preventive action to reduce or avoid their withdrawal).


The dataset contains data of students enrolled at university in 2012 and it is composed of four variables:
- Gender: categorical covariate indicating the gender of an individual (Male or Female).
- CFUP: covariate indicating the number of CFU (Credito Formativo Universitario) passed by the students by the end of the first semester. This variable is standardized.
- time_to_event: variable indicating the time-instant during the follow-up in which a student decides to abandon university. If a student does not abandon university in the follow-up or does not abandon it at all, the variable assumes a default value greater to the end of the follow-up (e.g. 6.1 semesters).
-centre: categorical variable indicating the course of study each student belongs to. In this case, there are 16 different levels composed of 4 letters: CosA, CosB, ... , CosP.

In general, numerical covariates must be standardized and categorical covariates should not be transformed into dummy variables (each model does this operation).


## How to execute a model
To execute a model several variables need to be specified:
- dataset 'data_dropout'
- time-domain vector: it can coincide with the follow-up or can be contained in it. Indeed, in the case related to our dataset, no events are registered in the first semester and it would be useless to fit the model there. Thus, the time-domain starts at the end of the first semester (t=1) and it ends at the end of the third academic year (t=6). The temporal unit of measure is semester.
- categories_range_min and categories_range_max vectors: to maximize the log-likelihood function, we apply a constraint optimization method in multi-dimension and we need to provide the minimum and maximum existence range of each parameter. Since the latter may have the same meaning, we decide to group them in category and to provide both the minimum and maximum range to each category.
According to the model, different categories are individuated.
- formula object: it specifies the relationship between the time-to-event, the dataset covariates and the cluster variable. The names reported in the formula must be contained also in the dataset, otherwise an error is thrown and the model stops the execution. 
Concerning the clustering variable (i.e. 'centre'), it must be the arguments of a sort of cluster function: e.g. cluster(centre).

Then we call the general method, specified which is the model we want to apply: 'AdPaikModel', 'PowParModel' or 'StocTimeDepModel'.
The execution returns an S3 class object with all the provided information and computed variables.

Look at the 'ModelsApplication.R' script to observe how to proceed.


## How to analyze the results
To analyze the results, we can take advantage of already implemented methods that plot:
- baseline hazard step-function
- frailty standard deviation (or variance)
- posterior frailty estimates
and we can call the summary method that prints the most relevant results of the model call.

Look at the 'ModelsApplication.R' script to observe how to proceed.


## Which are the exportable functions?
Among all the functions saved in the workspace, some of them can be directly called by the user, while others cannot because they are internally used by other methods.
The callable ones are reported in the following list, without arguments for convenience:
- AdPaikModel(), PowParModel(), StocTimeDepModel()
- AdPaik_1D(),PowPar_1D(), StocTimeDep_1D()
- TimeDepFrailty()
- summary(), summary.AdPaik(), summary.PowPar(), summary.StocTimeDep(...coming)
- frailty.sd.AdPaik()
- plot.bas_hazard(), plot.frailty_sd(), plot.post_frailty_est().


## To be aware of
- The first model is quite fast to optimize the log-likelihood function and to produce an output. However, considering the current dataset, it takes longer to estimated the 'Male' regression coefficient rather than the 'Female' one (reference). Therefore, I suggest to switch them and to set the 'Male' to the reference level and to estimate the 'Female' one. How to do that? Simply add a special character in front of 'Male' so that it becomes the first one in alphabetical order (e.g. '_Male' and 'Female'). The estimated coefficient changes, but not the optimal log-likelihood value.
- The 'Centre-Specific Frailty Model with Power Parameter' is not as fast as the first one, but it produces coherent and expected results.
- The 'Stochastic Time-Dependent' Centre-Specific Frailty Model' is really slow and it may not reach convergence.


## Authors and maintener of the codes
Alessandra Ragni (alessandra.ragni@polimi.it)
Giulia Romani (giulia.romani@mail.polimi.it),
Chiara Masci (chiara.masci@polimi.it).

