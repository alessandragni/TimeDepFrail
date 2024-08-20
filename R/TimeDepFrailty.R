# Define a unique function that receives all the parameters and the name of the
# frailty model we want to apply.
TimeDepFrailty <- function (formula, data, time_axis,
                            categories_range_min, categories_range_max,
                            C_mult = 0, flag_fullsd = TRUE,
                            time_domain = 0, flag_time_domain = FALSE,
                            n_extrarun = 60, tol_ll = 1e-6, tol_optimize = 1e-6, h_dd = 1e-3,
                            print_previous_ll_values = c(TRUE, 3),
                            model_type = c('AdPaikModel','PowParModel', 'StocTimeDepModel')){
  call <- switch (model_type,
                   'AdPaikModel' = AdPaikModel(formula, data, time_axis,
                                               categories_range_min, categories_range_max,
                                               flag_fullsd = TRUE,
                                               n_extrarun = 60, tol_ll = 1e-6, tol_optimize = 1e-6, h_dd = 1e-3,
                                               print_previous_ll_values = c(TRUE, 3)),
                   'PowParModel' = PowParModel(formula, data, time_axis,
                                               categories_range_min, categories_range_max,
                                               C_mult,
                                               n_extrarun = 60, tol_ll = 1e-6, tol_optimize = 1e-6, h_dd = 1e-3,
                                               print_previous_ll_values = c(TRUE, 3)),
                   'StocTimeDepModel' = StocTimeDepModel(formula, data, time_axis,
                                                         categories_range_min, categories_range_max,
                                                         C_mult,
                                                         time_domain = 0, flag_time_domain = FALSE,
                                                         n_extrarun = 40, tol_ll = 1e-6, tol_optimize = 1e-6, h_dd = 1e-3,
                                                         print_previous_ll_values = c(TRUE, 3)))
  return (call)
}
