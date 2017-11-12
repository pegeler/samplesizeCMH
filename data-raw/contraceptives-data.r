marginal_n <- matrix(
  c(273,716,2641,7260), 2, 2,
  dimnames = list(
    "OC Usage" = c("Exposed","Not Exposed"),
    "Disease Status" = c("Case","Control")
  )
)

(contraceptives_marginal <- as.table(marginal_n))

partial_n <- array(
  c(
    31,40,363,327,
    57,107,521,1146,
    69,196,760,1914,
    86,241,725,2503,
    30,132,272,1370
  ),
  dim = c(2,2,5),
  dimnames = list(
    "OC Usage" = c("Exposed","Not Exposed"),
    "Disease Status" = c("Case","Control"),
    "Age Group" = c("<= 34","35 - 39","40 - 44","45 - 49","50 - 55")
  )
)

(contraceptives <- as.table(partial_n))
