# Posterior Frailty Confidence Intervals

Function for computing the posterior frailty confidence intervals of the
time-dependent shared frailty Cox model.

Recalling the structure of the frailty \\Z\_{jk} = \alpha_j +
\epsilon\_{jk}, \forall j,k\\ with \\k=1,\dots,L\\ and \\j=1,\dots,N\\
as being composed by the sum of two independent gamma distributions:

- \\\alpha_j \sim gamma(\mu_1/\nu, 1/\nu), \forall j\\

- \\\epsilon\_{jk} \sim gamma(\mu_2/\gamma_k, 1/\gamma_k), \forall j,k\\

The posterior frailty estimate is \\\hat{Z}\_{jk} =
\hat{\alpha}\_{j}/\hat{\alpha}\_{max} +
\hat{\epsilon}\_{jk}/\hat{\epsilon}\_{max}\\. This function allows to
get the confidence intervals of either the entire posterior frailty
estimates \\\hat{Z}\_{jk}\\ or its time-independent
\\\frac{\hat{\alpha}\_{j}}{\hat{\alpha}\_{\text{max}}}\\ or
time-dependent
\\\frac{\hat{\epsilon}\_{jk}}{\hat{\epsilon}\_{\text{max}}}\\
components. The user can control which components to display using the
flag_eps and flag_alpha parameters. Only one of these flags can be set
to TRUE at a time.

## Usage

``` r
post_frailty_confint(
  object,
  level = 0.95,
  flag_eps = FALSE,
  flag_alpha = FALSE
)
```

## Arguments

- object:

  S3 object of class 'AdPaik' returned by the main model output, that
  contains all the information for the computation of the frailty
  standard deviation.

- level:

  A numeric value representing the confidence level for the posterior
  frailty confidence intervals. Default is 0.95 for 95% confidence.

- flag_eps:

  Logical flag indicating whether to extract only the time-dependent
  posterior frailty estimates. Default is FALSE.

- flag_alpha:

  Logical flag indicating whether to extract only the time-independent
  posterior frailty estimates. Default is FALSE.

## Value

A list for posterior frailty confidence intervals, depending on the
flag_eps and flag_alpha values. Specifically:

- A list of length equal to the N containing posterior frailty
  confidence intervals for \\\alpha_j, \forall j\\. In this case the
  flag_eps must be FALSE and the flag_alpha must be TRUE.

- A list of length equal to the NxL containing posterior frailty
  confidence intervals for \\\epsilon\_{jk}, \forall j,k\\. In this case
  the flag_eps must be TRUE and the flag_alpha must be FALSE.

- A list of length equal to the NxL containing posterior frailty
  confidence intervals for \\Z\_{jk} \forall j,k\\. In this case the
  flag_eps must be FALSE and the flag_alpha must be FALSE.

## Examples

``` r
# Consider the 'Academic Dropout dataset'
data(data_dropout)

# Define the variables needed for the model execution
formula <- time_to_event ~ Gender + CFUP + cluster(group)
time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
eps <- 1e-10
categories_range_min <- c(-8, -2, eps, eps, eps)
categories_range_max <- c(-eps, 0, 1 - eps, 1, 10)

# \donttest{
# Call the main model
result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max)

post_frailty_confint(result)
#>                      2.5 %    97.5 %
#> CosA_Interval_1  0.9167325 1.2928835
#> CosB_Interval_1  0.6567874 1.0190160
#> CosC_Interval_1  0.7949313 1.5302039
#> CosD_Interval_1  1.1758024 1.6785402
#> CosE_Interval_1  0.7948309 1.2354165
#> CosF_Interval_1  0.9036411 1.4136082
#> CosG_Interval_1  1.0817985 1.9315772
#> CosH_Interval_1  0.8990500 1.4806323
#> CosI_Interval_1  0.8474523 1.5584086
#> CosJ_Interval_1  1.0559897 1.4271198
#> CosK_Interval_1  1.0190923 1.7410384
#> CosL_Interval_1  1.1942270 1.5658521
#> CosM_Interval_1  1.2707017 1.7122075
#> CosN_Interval_1  1.0814305 1.4260787
#> CosO_Interval_1  0.8972905 1.5261456
#> CosP_Interval_1  0.8873230 1.3955676
#> CosA_Interval_2  0.8351864 1.3823668
#> CosB_Interval_2  0.5196663 1.0481828
#> CosC_Interval_2  0.7097843 1.5556009
#> CosD_Interval_2  1.2266343 1.9104840
#> CosE_Interval_2  0.7082403 1.3133426
#> CosF_Interval_2  0.8125194 1.4623067
#> CosG_Interval_2  1.0445438 2.0024782
#> CosH_Interval_2  0.8088437 1.5217847
#> CosI_Interval_2  0.7880660 1.6177908
#> CosJ_Interval_2  1.0055593 1.5532949
#> CosK_Interval_2  0.9474188 1.7805086
#> CosL_Interval_2  1.1760371 1.7180985
#> CosM_Interval_2  1.1625246 1.7375353
#> CosN_Interval_2  0.9613203 1.4560053
#> CosO_Interval_2  0.8200018 1.5750614
#> CosP_Interval_2  0.8094129 1.4659200
#> CosA_Interval_3  0.8976099 1.3973100
#> CosB_Interval_3  0.4239758 0.8735753
#> CosC_Interval_3  0.6933422 1.5395084
#> CosD_Interval_3  1.2782355 1.9230717
#> CosE_Interval_3  0.6244366 1.1659660
#> CosF_Interval_3  0.8895560 1.5370812
#> CosG_Interval_3  1.0180628 1.9704374
#> CosH_Interval_3  0.7937533 1.4878784
#> CosI_Interval_3  0.7614071 1.5842778
#> CosJ_Interval_3  1.0262074 1.5019917
#> CosK_Interval_3  1.1021239 1.9758415
#> CosL_Interval_3  1.3215293 1.8137036
#> CosM_Interval_3  1.3456640 1.9037451
#> CosN_Interval_3  1.0086518 1.4352600
#> CosO_Interval_3  0.8165988 1.5619362
#> CosP_Interval_3  0.8090769 1.4391427
#> CosA_Interval_4  0.9180675 1.2938300
#> CosB_Interval_4  0.6639971 1.0258038
#> CosC_Interval_4  0.7949638 1.5296607
#> CosD_Interval_4  1.1725445 1.6746987
#> CosE_Interval_4  0.7928139 1.2325331
#> CosF_Interval_4  0.9163478 1.4262561
#> CosG_Interval_4  1.0789614 1.9281685
#> CosH_Interval_4  0.8899802 1.4705773
#> CosI_Interval_4  0.8481401 1.5585455
#> CosJ_Interval_4  1.0490884 1.4195207
#> CosK_Interval_4  1.0157052 1.7370400
#> CosL_Interval_4  1.2071205 1.5796493
#> CosM_Interval_4  1.2714213 1.7126931
#> CosN_Interval_4  1.0771341 1.4213907
#> CosO_Interval_4  0.8962495 1.5244639
#> CosP_Interval_4  0.8823747 1.3897368
#> CosA_Interval_5  0.8936260 1.3111611
#> CosB_Interval_5  0.6408342 1.0481439
#> CosC_Interval_5  0.7857410 1.5452513
#> CosD_Interval_5  1.1926882 1.7317784
#> CosE_Interval_5  0.7691959 1.2463506
#> CosF_Interval_5  0.9126772 1.4578471
#> CosG_Interval_5  1.0521840 1.9216008
#> CosH_Interval_5  0.8630941 1.4721552
#> CosI_Interval_5  0.8505229 1.5871209
#> CosJ_Interval_5  1.0109549 1.4203228
#> CosK_Interval_5  0.9942036 1.7391671
#> CosL_Interval_5  1.2014325 1.6164856
#> CosM_Interval_5  1.2426890 1.7186006
#> CosN_Interval_5  1.0518228 1.4384600
#> CosO_Interval_5  0.8765721 1.5321583
#> CosP_Interval_5  0.8732685 1.4148920
#> CosA_Interval_6  0.7367814 1.2441551
#> CosB_Interval_6  0.5252786 1.0792313
#> CosC_Interval_6  0.7138588 1.5996846
#> CosD_Interval_6  1.0834631 1.7339452
#> CosE_Interval_6  0.6683107 1.2755931
#> CosF_Interval_6  0.7310733 1.3696335
#> CosG_Interval_6  0.9645411 1.9308133
#> CosH_Interval_6  0.7902921 1.5260837
#> CosI_Interval_6  0.7618983 1.6196770
#> CosJ_Interval_6  0.9678294 1.4924210
#> CosK_Interval_6  0.9480190 1.8153680
#> CosL_Interval_6  1.3317374 1.9160571
#> CosM_Interval_6  1.2767537 1.8911833
#> CosN_Interval_6  1.1473118 1.6896530
#> CosO_Interval_6  0.8385430 1.6346456
#> CosP_Interval_6  0.8100919 1.4945953
#> CosA_Interval_7  0.6404724 1.4083561
#> CosB_Interval_7  0.3875845 1.1901538
#> CosC_Interval_7  0.6037430 1.6686329
#> CosD_Interval_7  0.9333248 1.7927660
#> CosE_Interval_7  0.5388991 1.3791834
#> CosF_Interval_7  0.6734710 1.5584523
#> CosG_Interval_7  0.8995554 2.0419026
#> CosH_Interval_7  0.7285191 1.7249530
#> CosI_Interval_7  0.6577245 1.7012948
#> CosJ_Interval_7  0.9586584 1.8858162
#> CosK_Interval_7  0.8253953 1.8706013
#> CosL_Interval_7  1.0126725 1.8345027
#> CosM_Interval_7  1.1527069 2.0849237
#> CosN_Interval_7  0.8749452 1.6731142
#> CosO_Interval_7  0.7432704 1.7750578
#> CosP_Interval_7  0.6487338 1.5375769
#> CosA_Interval_8  0.7199901 1.3231513
#> CosB_Interval_8  0.4084059 1.0325998
#> CosC_Interval_8  0.5637580 1.5161362
#> CosD_Interval_8  0.9651260 1.6664852
#> CosE_Interval_8  0.6678072 1.4281561
#> CosF_Interval_8  0.7306409 1.4984956
#> CosG_Interval_8  1.0090820 2.1449377
#> CosH_Interval_8  0.8983950 1.8826959
#> CosI_Interval_8  0.7643005 1.7943786
#> CosJ_Interval_8  0.8816486 1.4701200
#> CosK_Interval_8  0.9215460 1.9362668
#> CosL_Interval_8  1.0163122 1.5707044
#> CosM_Interval_8  1.5478614 2.4400334
#> CosN_Interval_8  0.8996862 1.4351760
#> CosO_Interval_8  0.8231881 1.7882959
#> CosP_Interval_8  0.6317781 1.3489353
#> CosA_Interval_9  0.9166334 1.2923978
#> CosB_Interval_9  0.6610701 1.0226819
#> CosC_Interval_9  0.7930939 1.5276722
#> CosD_Interval_9  1.1760018 1.6783422
#> CosE_Interval_9  0.7924958 1.2321801
#> CosF_Interval_9  0.9129983 1.4227321
#> CosG_Interval_9  1.0714834 1.9204424
#> CosH_Interval_9  0.8919920 1.4726328
#> CosI_Interval_9  0.8494692 1.5598652
#> CosJ_Interval_9  1.0517121 1.4224801
#> CosK_Interval_9  1.0200715 1.7414984
#> CosL_Interval_9  1.2052159 1.5778709
#> CosM_Interval_9  1.2682185 1.7094584
#> CosN_Interval_9  1.0838396 1.4288219
#> CosO_Interval_9  0.8950359 1.5231596
#> CosP_Interval_9  0.8847265 1.3921705
#> CosA_Interval_10 0.9132515 1.3099891
#> CosB_Interval_10 0.6583012 1.0419383
#> CosC_Interval_10 0.7841369 1.5293338
#> CosD_Interval_10 1.1806040 1.6993150
#> CosE_Interval_10 0.7868259 1.2443833
#> CosF_Interval_10 0.9105993 1.4360949
#> CosG_Interval_10 1.0808231 1.9397373
#> CosH_Interval_10 0.8890526 1.4834501
#> CosI_Interval_10 0.8529553 1.5749487
#> CosJ_Interval_10 1.0495297 1.4417032
#> CosK_Interval_10 1.0157300 1.7482060
#> CosL_Interval_10 1.2100878 1.6045445
#> CosM_Interval_10 1.2713018 1.7309957
#> CosN_Interval_10 1.1056320 1.4762193
#> CosO_Interval_10 0.8877955 1.5283644
#> CosP_Interval_10 0.8900068 1.4137797
# }
```
