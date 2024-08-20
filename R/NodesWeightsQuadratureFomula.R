#' Nodes and weights for the Gauss_hermite quadrature formula 
#' for the 'Centre-Specific Frailty Model with Power Parameter'.
#' The nodes and weights have been extracted from the 'Handbook of Mathematical functions'
#' pag 940.

n_nodes <- 9
nodes9_ghqm <- matrix(c(0.000000000000000,
                       -0.723551018752838,
                       0.723551018752838,
                       -1.468553289216668,
                       1.468553289216668,
                       -2.266580584531843,
                       2.266580584531843,
                       -3.190993201781528,
                       3.190993201781528), nrow = n_nodes, ncol = 1)
weights9_ghqm <- matrix(c(7.202352156061e-1,
                         4.3265155900261e-1,
                         4.326515590026e-1,
                         8.847452739438e-2,
                         8.847452739438e-2,
                         4.943624275537e-3,
                         4.943624275537e-3,
                         3.960697726326e-5,
                         3.960697726326e-5), nrow = n_nodes, ncol = 1)


#' Nodes and weights for the Gauss-Hermite quadrature formula,
#' for the 'Stochastic Time-Dependent Centre-Specific Frailty Model'.
#' For the G function, the chosen nodes should not contain the zero (node)
#' since it appears at the denominator of a fraction.
#' Also in this case, the nodes and weights have been extracted from the 
#' 'Handbook of Mathematical functions', pag 940.

n_nodesG <- 10
nodes_ghqm <- matrix(c(0.342901327223705, 
                        -0.342901327223705,
                        1.036610829789514,
                        -1.036610829789514,
                        1.756683649299882,
                        - 1.756683649299882,
                        2.532731674232790,
                        -2.532731674232790,
                        3.436159118837738,
                        -3.436159118837738), nrow = n_nodesG, ncol = 1)
weightsG_ghqm <- matrix(c(6.108626337353e-1,
                          6.108626337353e-1,
                          2.401386110823e-1,
                          2.401386110823e-1,
                          3.387439445548e-2,
                          3.387439445548e-2,
                          1.343645746781e-3,
                          1.343645746781e-3,
                          7.640432855233e-6,
                          7.640432855233e-6), nrow = n_nodesG, ncol = 1)