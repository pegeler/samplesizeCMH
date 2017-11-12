## Stata power cmh manual
## Should make our own examples
#  Example 1
power.cmh.test(
 p1 = NULL,
 p2 = c(.426,.444,.364),
 theta = 2.5,
 sig.level = 0.05,
 # power = .8,
 N = 171,
 alternative = "one",
 method = "cc"
) %>% unclass
#  Example 2
power.cmh.test(
 p1 = NULL,
 p2 = c(.426,.444,.364),
 theta = 2.5,
 sig.level = 0.05,
 power = .8,
 alternative = "two",
 t = c(4,1,4) / sum(c(4,1,4))
)
#  Example 3a
power.cmh.test(
 p1 = NULL,
 p2 = c(.426,.444,.364),
 theta = 2.5,
 sig.level = 0.05,
 power = .8,
 alternative = "two",
 t = c(4,1,4) / sum(c(4,1,4)),
 s = 1 - c(.47,.57,.51)
)
#  Example 3b
power.cmh.test(
 p1 = NULL,
 p2 = c(.426,.444,.364),
 theta = 2.5,
 sig.level = 0.05,
 power = .8,
 alternative = "two",
 t = c(4,1,4) / sum(c(4,1,4)),
 s = 1 - c(.8,.7,.3)
)

# From "Sample size determination for case-control studies and the comparison
# of stratified and unstratified analyses", Nam (1992). See references.

# Uncorrected sample size estimate first introduced
# by Woolson and others in 1986
sample_size_uncorrected <- power.cmh.test(
  p2 = c(0.75,0.70,0.65,0.60),
  theta = 3,
  power = 0.9,
  t = c(0.10,0.40,0.35,0.15),
  alternative = "one",
  method = "bin"
)

sample_size_uncorrected

# We see that the N.exact is 171, the same as calculated by Nam
sample_size_uncorrected$N.exact


# Continuity corrected sample size estimate added by Nam
sample_size_corrected <- samplesize.cmh.test(
  p2 = c(0.75,0.70,0.65,0.60),
  theta = 3,
  power = 0.9,
  t = c(0.10,0.40,0.35,0.15),
  alternative = "one",
  method = "cc.bin"
)

sample_size_corrected

# We see that the N.exact is indeed equal to that which is reported in the paper
sample_size_corrected$N.exact


# Hypergeometric
samplesize.cmh.test(
  p1 = 0.5,
  theta = 1/rep(2.5,5),
  t = c(.24,.23,.20,.18,.15),
  r = c(.021,.042,.034,.06,.08),
  s = c(.0007,.0014,.0037,.009,.0204),
  power = .8,
  alternative = "one",
  method = "h"
)
