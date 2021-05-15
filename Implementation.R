### Data simulation ###

rm(list = ls())

# First we simulate ask prices:
# for a given stock with expected rate of return mu and volatility sigma, and 
# initial ask price ask0, and a time horizon T = 1, simulate n = 1000 trajectories 
# of the price Pt from time  t = 0 up until  t = T  through sampling_points many time periods, 
# each of length  delta t = T/sampling_points = 1, assuming the geometric Brownian motion model
library(sde)
set.seed(123)
# Set parameters
mu =  - 0.01; sigma = 0.01; ask0 = 40; T = 1 # a day
n = 1000; sampling_points = 1
# Generate n trajectories
dt = T/sampling_points; t = seq(0, T, by = dt)
X = matrix(rep(0, length(t) * n), nrow = n)
for (i in 1:n) {X[i,]= GBM(x = ask0, r = mu, sigma = sigma, T = T, N = sampling_points)}
# Plot evolution of ask price from t = 0 to t + 1
ymax = max(X); ymin = min(X) # bounds for simulated prices
plot(t, X[1,], t='l', ylim = c(ymin, ymax), col=1,
        ylab = "Ask price P(t)",xlab = "time t")
for(i in 2:n){lines(t,X[i,], t='l', ylim = c(ymin, ymax), col=i)}
rm(dt, i, mu, sampling_points, sigma, t, T, ymax, ymin)

ask0 <- X[,1]
ask1 <- X[,2]
# Let's now simulate bid-ask spread at time t = 0
# we have ask price which is 40
# the bid ask spread as a fraction of ask price  will be uniformly distributed
# with mean = 0.01368875 and sd = 0.006468388
set.seed(123)
spread0 <- runif(n, 0.1, 1)
spread0 <- spread0 / ask0   # bid-ask spread as a fraction of ask price
mean(spread0)
sd(spread0)

# we recover bid prices at time t = 0 from the spread
bid0 <- ask0 * (1 - spread0)

sum(bid0 <= ask0) # 1000

# Then we simulate discounts as a fraction of ask0
# which are just random
# but this decision should be influence by the bid ask spread at time t = 0
# but it can not be just function of bid ask spread at time t = 0,
# then it would not make sense to have both in the probit regression
set.seed(123)
gamma <- rnorm(n, 0.6, 0.2)
sum(gamma >= 0) # 1000
gamma <- gamma / ask0
mean(gamma)
sd(gamma)

orderPrice <- ask0 - gamma * ask0

sum(orderPrice <= ask0) # 1000
sum(orderPrice >= bid0)

summary(orderPrice)

# To be able to calculate returns, we need to know bid-ask spread at time t + 1
set.seed(123)
spread1 <- runif(n, 0.1, 1) # this is not expressed as a fraction of anything
bid1 <- ask1 - spread1
sum(bid1 <= ask1)

# now check the execution of the orders
# the order is executed at time t + 1 if orderPrice >= ask1
sum(orderPrice >= ask0) # no order is market order or marketable limit order
sum(orderPrice >= ask1) # 353 orders are executed
fill <- numeric(n)
fill[orderPrice >= ask1] <- 1
valuationPrice <- (ask1 + bid1)/2
r <- numeric(n)
r[orderPrice >= ask1] <- (valuationPrice[orderPrice >= ask1] - orderPrice[orderPrice >= ask1])/ask0[orderPrice >= ask1]


# Is the investor forced to trade at the end of the trading period if the limit order is not executed?
p <- 0.5
set.seed(123)
z <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(p, 1 - p))
# conditional simulation from the mixture components
r[(orderPrice < ask1) & z] <- (valuationPrice[(orderPrice < ask1) & z] - ask1[(orderPrice < ask1) & z])/ask0[(orderPrice < ask1) & z]

hist(r)

# we need to divide our sample into training and test data set
set.seed(123)
smp_size <- floor(0.75 * n)
train_ind <- sample(seq_len(n), size = smp_size)

###  Probit regression ###
# estimation of beta parameters
X <- as.data.frame(X)
names(X)[names(X) == "V1"] <- "ask0"
names(X)[names(X) == "V2"] <- "ask1"
X$bid0 <- bid0
X$bid1 <- bid1
X$spread0 <- spread0
X$gamma <- gamma
X$gamma2 <- gamma^2
X$gamma3 <- gamma^3
X$r <- r
X$r2 <- r^2
X$r3 <- r^3
X$fill <- fill
X$forced <- z

train <- X[train_ind, ]
test <- X[-train_ind, ]

myprobit <- glm(fill ~ spread0 + gamma + gamma2 + gamma3 + r + r2 + r3, 
                family = binomial(link = "probit"), data = train)

b0 <- myprobit$coefficients[1]
b1 <- myprobit$coefficients[2]
b2 <- myprobit$coefficients[3]
b3 <- myprobit$coefficients[4]
b4 <- myprobit$coefficients[5]
b5 <- myprobit$coefficients[6]
b6 <- myprobit$coefficients[7]
b7 <- myprobit$coefficients[8]

summary(myprobit)

# How good is our probit model?
prediction <- predict(myprobit, test, type = "response")
true <- test$fill
dif <- prediction - true
plot(test$gamma, dif)
rm(true, prediction, dif)
# all fills from test are predicted correctly

# Summary statistics
mean(fill)
sd(fill)
mean(spread0)
sd(spread0)
mean(gamma)
sd(gamma)
mean(r)
sd(r)

# Correlations
cor(spread0, fill)
cor(gamma, fill)
cor(gamma, spread0)
cor(r, fill)
cor(r, spread0)
cor(r, gamma)


library(writexl)
write_xlsx(X, "SimulatedData.xlsx")

# Plot bid prices
Y <- data.frame(X$bid0, X$bid1)
t <- 0:1
ymax = max(Y); ymin = min(Y) # bounds for simulated prices
plot(t, Y[1,], t='l', col=1,
     ylab = "Bid price P(t)",xlab = "time t")
for(i in 2:n){lines(t, Y[i,], t='l', col=i)}
rm(t,Y)


# we need to estimate the distribution of return as a normal distribution
# and get parameters mu and sigma
library(fitdistrplus)
FIT <- fitdist(r[train_ind], "norm")    ## note: it is "norm" not "normal"
class(FIT)
mu <- FIT$estimate[1] # mean(r)
sigma <- FIT$estimate[2] # sd(r)
plot(FIT)    ## use method `plot.fitdist`

I1 <- function(gamma){
  integrand <- function(r) {
    # here we need to have 1D spread, so I've chosen one
    dnorm(r, mu, sigma) * pnorm(b0 + b1 * spread0[1] + b2 * gamma + b3 * gamma^ 2 + b4 * gamma^3 + 
          b5 * r + b6 * r^2  + b7 * r^3)
  }
  integrate(integrand, lower = - Inf, upper = Inf)
}
# normal density x Phi is a function of return r and returns the function of gamma

# bid and ask prices are just inputs
# r_a is a variable I need to integrate over, but
# the result is function of gamma

