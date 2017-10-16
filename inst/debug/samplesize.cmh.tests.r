## Stata power cmh manual
## Should make our own examples
#  Example 1
samplesize.cmh(
 p1 = NULL,
 p2 = c(.426,.444,.364),
 theta = 2.5,
 sig.level = 0.05,
 power = .8,
 alternative = "two",
 method = "c"
)
#  Example 2
samplesize.cmh(
 p1 = NULL,
 p2 = c(.426,.444,.364),
 theta = 2.5,
 sig.level = 0.05,
 power = .8,
 alternative = "two",
 t = c(4,1,4) / sum(c(4,1,4))
)
#  Example 3a
samplesize.cmh(
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
samplesize.cmh(
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
sample_size_uncorrected <- samplesize.cmh(
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
sample_size_corrected <- samplesize.cmh(
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

