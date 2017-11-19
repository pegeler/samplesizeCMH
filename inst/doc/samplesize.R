## ----unstratified--------------------------------------------------------
library(samplesizeCMH)

sample_size_unstratified <- power.cmh.test(
  p2 = c(0.75,0.70,0.65,0.60),
  theta = 3,
  power = 0.9,
  t = c(0.10,0.40,0.35,0.15),
  alternative = "one",
  method = "unstratified"
)

print(sample_size_unstratified, detail = FALSE)

## ----binomial------------------------------------------------------------
sample_size_binomial <- power.cmh.test(
  p2 = c(0.75,0.70,0.65,0.60),
  theta = 3,
  power = 0.9,
  t = c(0.10,0.40,0.35,0.15),
  alternative = "one",
  method = "bin"
)

sample_size_binomial

## ----cc------------------------------------------------------------------
# Continuity corrected sample size estimate added by Nam
sample_size_corrected <- power.cmh.test(
  p2 = c(0.75,0.70,0.65,0.60),
  theta = 3,
  power = 0.9,
  t = c(0.10,0.40,0.35,0.15),
  alternative = "one",
  method = "cc.bin"
)

sample_size_corrected

## ----cc2-----------------------------------------------------------------
sample_size_corrected$N

## ----s-------------------------------------------------------------------
power.cmh.test(
  p2 = c(0.75,0.70,0.65,0.60),
  s = 1/3,
  theta = 3,
  power = 0.9,
  t = c(0.10,0.40,0.35,0.15),
  alternative = "one",
  method = "cc.bin"
)

## ----effect-size---------------------------------------------------------
effect.size(0.75,3)

## ----prop2odds-----------------------------------------------------------
# Control
prop2odds(0.75)

# Case
prop2odds(0.9)

## ----props2theta---------------------------------------------------------
props2theta(0.90,0.75)

