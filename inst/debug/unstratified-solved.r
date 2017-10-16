z_a <- stats::qnorm(.05 / 2)
z_b <- stats::qnorm(1 - .8)
s = 0.5
t = 1 / 3
p2 = c(.426,.444,.364)
theta = 2.5
p1 = effect.size(p2, theta)

J_vec = rep(1, 3)

sum_ts <- sum(t * s * J_vec) # Need to ensure each element is getting replicated to J

sum_t1ms <- sum(t * (1 - s) * J_vec)

p2p <- sum(t * (1 - s) * p2 / sum_t1ms)

p1p <- sum(t * s * p1 / sum_ts) # Note to confirm

ppp <- p1p * sum_ts + p2p * sum_t1ms

(z_a * sqrt(sum_ts * sum_t1ms * ppp * (1 - ppp)) +
    z_b * sqrt(sum_ts^2 * sum_t1ms^2 *
                 (p1p * (1 - p1p) / sum_ts + p2p * (1 - p2p) / sum_t1ms)))^2 /
  ((p1p - p2p) * sum_ts * sum_t1ms)^2
