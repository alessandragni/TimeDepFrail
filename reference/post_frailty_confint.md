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
#> CosA_Interval_1  0.9167341 1.2928857
#> CosB_Interval_1  0.6567888 1.0190181
#> CosC_Interval_1  0.7949324 1.5302050
#> CosD_Interval_1  1.1758042 1.6785425
#> CosE_Interval_1  0.7948325 1.2354187
#> CosF_Interval_1  0.9036419 1.4136089
#> CosG_Interval_1  1.0817996 1.9315785
#> CosH_Interval_1  0.8990520 1.4806354
#> CosI_Interval_1  0.8474539 1.5584110
#> CosJ_Interval_1  1.0559915 1.4271221
#> CosK_Interval_1  1.0190943 1.7410414
#> CosL_Interval_1  1.1942288 1.5658542
#> CosM_Interval_1  1.2707037 1.7122101
#> CosN_Interval_1  1.0814321 1.4260807
#> CosO_Interval_1  0.8972920 1.5261477
#> CosP_Interval_1  0.8873245 1.3955697
#> CosA_Interval_2  0.8351875 1.3823688
#> CosB_Interval_2  0.5196672 1.0481847
#> CosC_Interval_2  0.7097850 1.5556020
#> CosD_Interval_2  1.2266357 1.9104864
#> CosE_Interval_2  0.7082414 1.3133448
#> CosF_Interval_2  0.8125197 1.4623073
#> CosG_Interval_2  1.0445447 2.0024796
#> CosH_Interval_2  0.8088453 1.5217877
#> CosI_Interval_2  0.7880672 1.6177933
#> CosJ_Interval_2  1.0055605 1.5532970
#> CosK_Interval_2  0.9474204 1.7805116
#> CosL_Interval_2  1.1760383 1.7181004
#> CosM_Interval_2  1.1625260 1.7375376
#> CosN_Interval_2  0.9613212 1.4560068
#> CosO_Interval_2  0.8200028 1.5750635
#> CosP_Interval_2  0.8094140 1.4659221
#> CosA_Interval_3  0.8976105 1.3973107
#> CosB_Interval_3  0.4239764 0.8735763
#> CosC_Interval_3  0.6933425 1.5395087
#> CosD_Interval_3  1.2782360 1.9230724
#> CosE_Interval_3  0.6244372 1.1659672
#> CosF_Interval_3  0.8895556 1.5370804
#> CosG_Interval_3  1.0180631 1.9704379
#> CosH_Interval_3  0.7937545 1.4878807
#> CosI_Interval_3  0.7614080 1.5842796
#> CosJ_Interval_3  1.0262079 1.5019924
#> CosK_Interval_3  1.1021252 1.9758439
#> CosL_Interval_3  1.3215296 1.8137039
#> CosM_Interval_3  1.3456648 1.9037461
#> CosN_Interval_3  1.0086521 1.4352603
#> CosO_Interval_3  0.8165995 1.5619374
#> CosP_Interval_3  0.8090775 1.4391437
#> CosA_Interval_4  0.9180692 1.2938322
#> CosB_Interval_4  0.6639986 1.0258059
#> CosC_Interval_4  0.7949649 1.5296619
#> CosD_Interval_4  1.1725463 1.6747010
#> CosE_Interval_4  0.7928154 1.2325353
#> CosF_Interval_4  0.9163487 1.4262569
#> CosG_Interval_4  1.0789626 1.9281698
#> CosH_Interval_4  0.8899822 1.4705804
#> CosI_Interval_4  0.8481418 1.5585479
#> CosJ_Interval_4  1.0490902 1.4195231
#> CosK_Interval_4  1.0157072 1.7370430
#> CosL_Interval_4  1.2071223 1.5796517
#> CosM_Interval_4  1.2714234 1.7126958
#> CosN_Interval_4  1.0771358 1.4213928
#> CosO_Interval_4  0.8962510 1.5244659
#> CosP_Interval_4  0.8823763 1.3897390
#> CosA_Interval_5  0.8936288 1.3111622
#> CosB_Interval_5  0.6408367 1.0481446
#> CosC_Interval_5  0.7857425 1.5452514
#> CosD_Interval_5  1.1926897 1.7317785
#> CosE_Interval_5  0.7691986 1.2463520
#> CosF_Interval_5  0.9126783 1.4578462
#> CosG_Interval_5  1.0521860 1.9216018
#> CosH_Interval_5  0.8630973 1.4721579
#> CosI_Interval_5  0.8505247 1.5871222
#> CosJ_Interval_5  1.0109583 1.4203247
#> CosK_Interval_5  0.9942066 1.7391699
#> CosL_Interval_5  1.2014348 1.6164862
#> CosM_Interval_5  1.2426922 1.7186026
#> CosN_Interval_5  1.0518256 1.4384611
#> CosO_Interval_5  0.8765745 1.5321597
#> CosP_Interval_5  0.8732708 1.4148930
#> CosA_Interval_6  0.7367821 1.2441563
#> CosB_Interval_6  0.5252792 1.0792328
#> CosC_Interval_6  0.7138592 1.5996855
#> CosD_Interval_6  1.0834639 1.7339467
#> CosE_Interval_6  0.6683115 1.2755947
#> CosF_Interval_6  0.7310731 1.3696334
#> CosG_Interval_6  0.9645415 1.9308141
#> CosH_Interval_6  0.7902934 1.5260864
#> CosI_Interval_6  0.7618993 1.6196793
#> CosJ_Interval_6  0.9678303 1.4924223
#> CosK_Interval_6  0.9480203 1.8153709
#> CosL_Interval_6  1.3317383 1.9160584
#> CosM_Interval_6  1.2767549 1.8911852
#> CosN_Interval_6  1.1473124 1.6896540
#> CosO_Interval_6  0.8385438 1.6346474
#> CosP_Interval_6  0.8100927 1.4945970
#> CosA_Interval_7  0.6404733 1.4083583
#> CosB_Interval_7  0.3875852 1.1901561
#> CosC_Interval_7  0.6037434 1.6686344
#> CosD_Interval_7  0.9333257 1.7927683
#> CosE_Interval_7  0.5388999 1.3791858
#> CosF_Interval_7  0.6734710 1.5584533
#> CosG_Interval_7  0.8995558 2.0419042
#> CosH_Interval_7  0.7285205 1.7249566
#> CosI_Interval_7  0.6577255 1.7012976
#> CosJ_Interval_7  0.9586597 1.8858190
#> CosK_Interval_7  0.8253967 1.8706046
#> CosL_Interval_7  1.0126737 1.8345050
#> CosM_Interval_7  1.1527084 2.0849267
#> CosN_Interval_7  0.8749462 1.6731162
#> CosO_Interval_7  0.7432713 1.7750604
#> CosP_Interval_7  0.6487346 1.5375793
#> CosA_Interval_8  0.7199908 1.3231523
#> CosB_Interval_8  0.4084064 1.0326010
#> CosC_Interval_8  0.5637581 1.5161367
#> CosD_Interval_8  0.9651268 1.6664864
#> CosE_Interval_8  0.6678080 1.4281576
#> CosF_Interval_8  0.7306406 1.4984951
#> CosG_Interval_8  1.0090824 2.1449386
#> CosH_Interval_8  0.8983965 1.8826988
#> CosI_Interval_8  0.7643014 1.7943809
#> CosJ_Interval_8  0.8816493 1.4701210
#> CosK_Interval_8  0.9215473 1.9362695
#> CosL_Interval_8  1.0163129 1.5707053
#> CosM_Interval_8  1.5478624 2.4400347
#> CosN_Interval_8  0.8996868 1.4351766
#> CosO_Interval_8  0.8231888 1.7882975
#> CosP_Interval_8  0.6317787 1.3489365
#> CosA_Interval_9  0.9166351 1.2924001
#> CosB_Interval_9  0.6610716 1.0226840
#> CosC_Interval_9  0.7930950 1.5276733
#> CosD_Interval_9  1.1760036 1.6783445
#> CosE_Interval_9  0.7924974 1.2321823
#> CosF_Interval_9  0.9129993 1.4227329
#> CosG_Interval_9  1.0714846 1.9204437
#> CosH_Interval_9  0.8919939 1.4726358
#> CosI_Interval_9  0.8494708 1.5598677
#> CosJ_Interval_9  1.0517140 1.4224825
#> CosK_Interval_9  1.0200735 1.7415014
#> CosL_Interval_9  1.2052178 1.5778732
#> CosM_Interval_9  1.2682206 1.7094611
#> CosN_Interval_9  1.0838414 1.4288240
#> CosO_Interval_9  0.8950374 1.5231617
#> CosP_Interval_9  0.8847281 1.3921727
#> CosA_Interval_10 0.9132540 1.3099900
#> CosB_Interval_10 0.6583035 1.0419390
#> CosC_Interval_10 0.7841387 1.5293344
#> CosD_Interval_10 1.1806061 1.6993159
#> CosE_Interval_10 0.7868283 1.2443845
#> CosF_Interval_10 0.9106007 1.4360946
#> CosG_Interval_10 1.0808242 1.9397375
#> CosH_Interval_10 0.8890552 1.4834524
#> CosI_Interval_10 0.8529572 1.5749502
#> CosJ_Interval_10 1.0495322 1.4417041
#> CosK_Interval_10 1.0157324 1.7482083
#> CosL_Interval_10 1.2100902 1.6045452
#> CosM_Interval_10 1.2713044 1.7309970
#> CosN_Interval_10 1.1056333 1.4762186
#> CosO_Interval_10 0.8877977 1.5283658
#> CosP_Interval_10 0.8900086 1.4137805
# }
```
