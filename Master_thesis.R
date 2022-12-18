### Load Libraries
library(tidyquant)
library(rmgarch)
library(tseries)
library(FinTS)
library(ggplot2)
require(gridExtra)
library(grid)
library(magrittr)
library(dplyr) 
library(tibble)
library(ks)
library(WeightedPortTest)

### Load data, return series
ret = read.csv("~/Desktop/Master/returns_R.csv")

### Univariate normal GARCH(1,1) for each return series
garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                          variance.model = list(garchOrder = c(1,1), 
                                                model = "sGARCH"), 
                          distribution.model = "norm")

### DCC estimation
dcc.garch11.spec = dccspec(uspec = multispec( replicate(3, garch11.spec) ), 
                           dccOrder = c(1,1), # dcc specification - DCC-GARCH(1,1) for conditional correlations
                           distribution = "mvnorm")
dcc.fit = dccfit(dcc.garch11.spec, data = ret, out.sample = 0)

### GO-GARCH estimation
go.garch11.spec <- gogarchspec( mean.model = list(armaOrder = c(0,0)),
                                variance.model     = list(model = "sGARCH", garchOrder = c(1,1)), 
                                distribution.model = "mvnorm")

go.fit <- gogarchfit(go.garch11.spec, ret, solver = "solnp", out.sample = 0)

### Model comparison, use information criteria and llh
## DCC-GARCH
dcc.fit@mfit$llh # Maximum loglikelihood
infocriteria(dcc.fit) # AIC and BIC

# Misspecification test of the residuals
skewness(dcc.fit@mfit$stdresid)
kurtosis(dcc.fit@mfit$stdresid)

dcc.fit@mfit$stdresid[,1] %>% jarque.bera.test()
dcc.fit@mfit$stdresid[,2] %>% jarque.bera.test()
dcc.fit@mfit$stdresid[,3] %>% jarque.bera.test()

Weighted.LM.test(dcc.fit@model$residuals[,1], dcc.fit@model$sigma[,1]^2, lag = 2,type = c("correlation", "partial"),fitdf = 1, weighted = FALSE)
Weighted.LM.test(dcc.fit@model$residuals[,2], dcc.fit@model$sigma[,2]^2, lag = 2,type = c("correlation", "partial"),fitdf = 1, weighted = FALSE)
Weighted.LM.test(dcc.fit@model$residuals[,3], dcc.fit@model$sigma[,3]^2, lag = 2,type = c("correlation", "partial"),fitdf = 1, weighted = FALSE)

## GO-GARCH
go.fit@mfit$llh # Maximum loglikelihood
# Information criteria
object = go.fit
npvar = 0
m = dim(object@model$modeldata$data)[2]
T = object@model$modeldata$T
itest = rugarch:::.information.test(object@mfit$llh, nObs = T, nPars =  npvar + (m^2 - m)/2 + length(object@mfit$matcoef[,1]))
itest$AIC
itest$BIC

# Misspecification test of the residuals
skewness(go.fit@mfit$residuals)
kurtosis(go.fit@mfit$residuals)

(go.fit@mfit$residuals[,1]/sqrt(rcov(go.fit)[1,1,])) %>% jarque.bera.test()
(go.fit@mfit$residuals[,2]/sqrt(rcov(go.fit)[2,2,])) %>% jarque.bera.test()
(go.fit@mfit$residuals[,3]/sqrt(rcov(go.fit)[3,3,])) %>% jarque.bera.test()

Weighted.LM.test(go.fit@mfit$residuals[,1], rcov(go.fit)[1,1,], lag = 2,type = c("correlation", "partial"),fitdf = 1, weighted = FALSE)
Weighted.LM.test(go.fit@mfit$residuals[,2], rcov(go.fit)[2,2,], lag = 2,type = c("correlation", "partial"),fitdf = 1, weighted = FALSE)
Weighted.LM.test(go.fit@mfit$residuals[,3], rcov(go.fit)[3,3,], lag = 2,type = c("correlation", "partial"),fitdf = 1, weighted = FALSE)

### Results from initial estimation
# Save results from DCC-GARCH
T = length(t(ret))/3
cov_JPM_MSFT_dcc = list()
for (t in 1:T){
  cov_JPM_MSFT_dcc = append(cov_JPM_MSFT_dcc, dcc.fit@mfit$H[1,2,t])}

cov_JPM_MSFT_dcc = (as.matrix(cov_JPM_MSFT_dcc ))%>% unlist

cov_JPM_KO_dcc = list()
for (t in 1:T){
  cov_JPM_KO_dcc = append(cov_JPM_KO_dcc, dcc.fit@mfit$H[1,3,t])}

cov_JPM_KO_dcc = (as.matrix(cov_JPM_KO_dcc ))%>% unlist

cov_MSFT_KO_dcc = list()
for (t in 1:T){
  cov_MSFT_KO_dcc = append(cov_MSFT_KO_dcc, dcc.fit@mfit$H[2,3,t])}

cov_MSFT_KO_dcc = (as.matrix(cov_MSFT_KO_dcc ))%>% unlist

# Save results from GO-GARCH
cov_go = rcov(go.fit)

cov_JPM_MSFT_go = list()
for (t in 1:T){
  cov_JPM_MSFT_go = append(cov_JPM_MSFT_go, cov_go[1,2,t])}

cov_JPM_MSFT_go = (as.matrix(cov_JPM_MSFT_go))%>% unlist

cov_JPM_KO_go = list()
for (t in 1:T){
  cov_JPM_KO_go = append(cov_JPM_KO_go, cov_go[1,3,t])}

cov_JPM_KO_go = (as.matrix(cov_JPM_KO_go ))%>% unlist

cov_MSFT_KO_go = list()
for (t in 1:T){
  cov_MSFT_KO_go = append(cov_MSFT_KO_go, cov_go[2,3,t])}

cov_MSFT_KO_go = (as.matrix(cov_MSFT_KO_go ))%>% unlist

# Import results from eigen-GARCH
cov_JPM_KO_eigen = read.csv("~/Desktop/Master/cov_JPM_KO_eigen.csv")
cov_JPM_MSFT_eigen = read.csv("~/Desktop/Master/cov_JPM_MSFT_eigen.csv")
cov_KO_MSFT_eigen = read.csv("~/Desktop/Master/cov_KO_MSFT_eigen.csv")

# Get dates
# Load data, time series closing prices
symbol.vec = c("JPM")
getSymbols(symbol.vec, from = "2001-01-02", to = "2009-12-31")

# Extract adjusted closing prices
JPM = `JPM`[, "JPM.Adjusted", drop=F][-2,]
length(JPM)

# Covariance plots
plot_4 <- ggplot()+
  geom_line(as.data.frame(cov_MSFT_KO_go), mapping = aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_go), colour="GO-GARCH"), size = 0.2) +
  geom_line(as.data.frame(cov_MSFT_KO_dcc), mapping=aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_dcc), colour="DCC-GARCH" ), size = 0.2) + 
  geom_line(as.data.frame(cov_KO_MSFT_eigen), mapping=aes(x = zoo::index(JPM), y = as.numeric(X0), colour="λ-GARCH"), size = 0.2) + 
  ggtitle("Microsoft & Coca-Cola") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  scale_color_manual( values = c("GO-GARCH" = "darkblue", "DCC-GARCH" = "red","λ-GARCH" = "green")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=10), legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_3 <- ggplot()+
  geom_line(as.data.frame(cov_MSFT_KO_go), mapping = aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_go)), size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_MSFT_KO_dcc), mapping=aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_dcc)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_KO_MSFT_eigen), mapping=aes(x = zoo::index(JPM), y = as.numeric(X0)), size = 0.2, colour = "green") + 
  ggtitle("(c) Microsoft & Coca-Cola") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15)) +
  labs(x = "", y ="")

plot_2 <- ggplot()+
  geom_line(as.data.frame(cov_JPM_KO_go), mapping = aes(x = zoo::index(JPM), y = as.numeric(cov_JPM_KO_go)), size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_JPM_KO_dcc), mapping=aes(x = zoo::index(JPM), y = as.numeric(cov_JPM_KO_dcc)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_JPM_KO_eigen), mapping=aes(x = zoo::index(JPM), y = as.numeric(X0)), size = 0.2, colour = "green") + 
  ggtitle("(b) JPMorgan & Coca-Cola") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15)) +
  labs(x = "", y ="")

plot_1 <- ggplot()+
  geom_line(as.data.frame(cov_JPM_MSFT_go), mapping = aes(x = zoo::index(JPM), y = as.numeric(cov_JPM_MSFT_go)), size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_JPM_MSFT_dcc), mapping=aes(x = zoo::index(JPM), y = as.numeric(cov_JPM_MSFT_dcc)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_JPM_MSFT_eigen), mapping=aes(x = zoo::index(JPM), y = as.numeric(X0)), size = 0.2, colour = "green") + 
  ggtitle("(a) JPMorgan & Microsoft") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15)) +
  labs(x = "", y ="")

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(plot_4)

grid.arrange(arrangeGrob(plot_1, plot_2, plot_3, ncol = 1),
             shared_legend, nrow = 2, heights = c(10, 1))

### Forecast
## DCC-GARCH
dccrolln <- dccroll(dcc.garch11.spec, 
                    data = ret, 
                    n.ahead = 1,
                    forecast.length = 755, 
                    refit.every = 1,
                    refit.window = "moving", 
                    fit.control=list(scale=TRUE),
                    solver.control=list(trace=1))

cov_dcc_f = rcov(dccrolln)

cov_JPM_MSFT_dcc_f = list()
for (t in 1:755){
  cov_JPM_MSFT_dcc_f = append(cov_JPM_MSFT_dcc_f, cov_dcc_f[1,2,t])}

cov_JPM_MSFT_dcc_f = (as.matrix(cov_JPM_MSFT_dcc_f))

cov_JPM_KO_dcc_f = list()
for (t in 1:755){
  cov_JPM_KO_dcc_f = append(cov_JPM_KO_dcc_f, cov_dcc_f[1,3,t])}

cov_JPM_KO_dcc_f = (as.matrix(cov_JPM_KO_dcc_f))

cov_MSFT_KO_dcc_f = list()
for (t in 1:755){
  cov_MSFT_KO_dcc_f = append(cov_MSFT_KO_dcc_f, cov_dcc_f[2,3,t])}

cov_MSFT_KO_dcc_f = (as.matrix(cov_MSFT_KO_dcc_f))

## GO-GARCH
gorolln <- gogarchroll(go.garch11.spec, 
                       data = ret, 
                       n.ahead = 1,
                       forecast.length = 755, 
                       refit.every = 1,
                       refit.window = "moving", 
                       fit.control=list(scale=TRUE),
                       solver.control=list(trace=1))

cov_go_f = rcov(gorolln)

cov_JPM_MSFT_go_f = list()
for (t in 1:755){
  cov_JPM_MSFT_go_f = append(cov_JPM_MSFT_go_f, cov_go_f[1,2,t])}

cov_JPM_MSFT_go_f = (as.matrix(cov_JPM_MSFT_go_f))

cov_JPM_KO_go_f = list()
for (t in 1:755){
  cov_JPM_KO_go_f = append(cov_JPM_KO_go_f, cov_go_f[1,3,t])}

cov_JPM_KO_go_f = (as.matrix(cov_JPM_KO_go_f))

cov_MSFT_KO_go_f = list()
for (t in 1:755){
  cov_MSFT_KO_go_f = append(cov_MSFT_KO_go_f, cov_go_f[2,3,t])}

cov_MSFT_KO_go_f = (as.matrix(cov_MSFT_KO_go_f))

## Eigen-GARCH
cov_JPM_KO_eigen_f = read.csv("~/Desktop/Master/cov_JPM_KO_forecast_eigen.csv")
cov_JPM_MSFT_eigen_f = read.csv("~/Desktop/Master/cov_JPM_MSFT_forcast_eigen.csv")
cov_KO_MSFT_eigen_f = read.csv("~/Desktop/Master/cov_KO_MSFT_forcast_eigen.csv")

## Results
plot_4 <- ggplot()+
  geom_line(as.data.frame(cov_MSFT_KO_go_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_MSFT_KO_go_f),colour="GO-GARCH" ),size = 0.2) +
  geom_line(as.data.frame(cov_MSFT_KO_dcc_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_MSFT_KO_dcc_f), colour="DCC-GARCH"), size = 0.2) + 
  geom_line(as.data.frame(cov_KO_MSFT_eigen_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(X0), colour="λ-GARCH"), size = 0.2) + 
  ggtitle("Microsoft & Coca Cola") +
  scale_color_manual(values = c("GO-GARCH" = "darkblue", "DCC-GARCH" = "red","λ-GARCH" = "green")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=10), legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_3 <- ggplot()+
  geom_line(as.data.frame(cov_MSFT_KO_go_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_MSFT_KO_go_f)),size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_MSFT_KO_dcc_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_MSFT_KO_dcc_f)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_KO_MSFT_eigen_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(X0)), size = 0.2, colour = "green") + 
  ggtitle("(c) Microsoft & Coca-Cola") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15), legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_2 <- ggplot()+
  geom_line(as.data.frame(cov_JPM_KO_go_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_JPM_KO_go_f)),size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_JPM_KO_dcc_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_JPM_KO_dcc_f)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_JPM_KO_eigen_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(X0), colour="Eigen-GARCH"), size = 0.2, colour = "green") + 
  ggtitle("(b) JPMorgan & Coca-Cola") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(x = "", y ="")

plot_1 <- ggplot()+
  geom_line(as.data.frame(cov_JPM_MSFT_go_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_JPM_MSFT_go_f)), size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_JPM_MSFT_dcc_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_JPM_MSFT_dcc_f)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_JPM_MSFT_eigen_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(X0)), size = 0.2, colour = "green") + 
  ggtitle("(a) JPMorgan & Microsoft") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(x = "", y ="")

# Apply user-defined function to extract legend
shared_legend <- extract_legend(plot_4)

grid.arrange(arrangeGrob(plot_1, plot_2, plot_3, ncol = 1),
             shared_legend, nrow = 2, heights = c(10, 1))


### Global Minimum Variance Portfolio
## DCC-GARCH

## Retun data
symbol.vec = c("JPM", "MSFT", "KO")
getSymbols(symbol.vec, from = "2001-01-02", to = "2009-12-31")

# Extract adjusted closing prices
JPM = `JPM`[, "JPM.Adjusted", drop=F]
MSFT = `MSFT`[, "MSFT.Adjusted", drop=F]
KO = `KO`[, "KO.Adjusted", drop=F]

# Calculate log-returns
JPM.ret = CalculateReturns(JPM)*100
MSFT.ret = CalculateReturns(MSFT)*100
KO.ret = CalculateReturns(KO)*100

# Remove first NA observation
JPM.ret = JPM.ret[-1,]
MSFT.ret = MSFT.ret[-1,]
KO.ret = KO.ret[-1,]

ret = na.omit(merge(merge(JPM.ret,MSFT.ret),KO.ret))
write.csv(ret, file="~/Desktop/Master/returns_R_simple.csv", row.names = FALSE)

pred_min_dcc <- apply(rcov(dccrolln), 3, function(x){ 
  w <- solve(x %>% as.matrix()) %*% rep(1, 3)
  return(w/sum(w))}) %>% 
  t()

apply(pred_min_dcc * ret %>% tail(755), 1, sum) %>% 
  enframe() %>% 
  summarise(mean = mean(value)*sqrt(250),
            sd = sqrt(250) * sd(value),
            sharpe = mean / sd)

ret_dcc = apply(pred_min_dcc * (ret)%>% tail(755), 1, sum)

# GO_GARCH
pred_min_go <- apply(rcov(gorolln), 3, function(x){ 
  w <- solve(x %>% as.matrix()) %*% rep(1, 3)
  return(w/sum(w))}) %>% 
  t()

apply(pred_min_go * ret %>% tail(755), 1, sum) %>% 
  enframe() %>% 
  summarise(mean = mean(value)*sqrt(250),
            sd = sqrt(250) * sd(value),
            sharpe = mean / sd)

ret_go = apply(pred_min_go * (ret)%>% tail(755), 1, sum)

# EW - portfolio
apply(rep(rep(1/3, 3),755) * ret %>% tail(755), 1, sum) %>% 
  enframe() %>% 
  summarise(mean = mean(value)*sqrt(250),
            sd = sqrt(250) * sd(value),
            sharpe = mean / sd)

EW = apply(rep(rep(1/3, 3),755) * (ret)%>% tail(755), 1, sum)

# Results
MVP_eigen = read.csv("~/Desktop/Master/MVP_eigen.csv")

plot_1 <- ggplot()+
  geom_line(as.data.frame(EW), mapping=aes(x = tail(zoo::index(ret),755), y = cumsum(EW), colour="EW"), size = 0.2) + 
  geom_line(as.data.frame(ret_dcc), mapping=aes(x = tail(zoo::index(ret),755), y = cumsum(ret_dcc), colour="DCC-GARCH"), size = 0.2) + 
  geom_line(as.data.frame(ret_go), mapping=aes(x = tail(zoo::index(ret),755), y = cumsum(ret_go),colour="GO-GARCH" ),size = 0.2) +
  geom_line(as.data.frame(MVP_eigen), mapping=aes(x = tail(zoo::index(ret),756), y = cumsum(X0), colour="λ-GARCH"), size = 0.2) + 
  scale_color_manual(name = "MGARCH-Models", values = c(EW ="lightblue", "GO-GARCH" = "darkblue", "DCC-GARCH" = "red","λ-GARCH" = "green")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_1

### Test for equal Sharpe Ratios
load("/Users/silkesommerfugl/Downloads/SharpeR/Sharpe.RData")
MVP_eigen_ret = read.csv("~/Desktop/Master/MVP_eigen.csv")

## Equal Sr Test: MGARCH vs EW
boot.time.inference(cbind(ret_dcc, EW), M = 755, b=4)
boot.time.inference(cbind(ret_go, EW), M = 755, b=4)

EW = apply(rep(rep(1/3, 3),756) * (ret)%>% tail(756), 1, sum)

boot.time.inference(cbind(MVP_eigen_ret[2], EW), M = 755, b=4)

## Equal Sr Test: Eigen vs MGARCH
boot.time.inference(cbind(tail(MVP_eigen_ret[2],755), ret_dcc), M = 755, b=4)
boot.time.inference(cbind(tail(MVP_eigen_ret[2],755), ret_go), M = 755, b=4)

### WEEKLY
## DCC estimation
garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                          variance.model = list(garchOrder = c(1,1), 
                                                model = "sGARCH"), 
                          distribution.model = "norm")
dcc.garch11.spec = dccspec(uspec = multispec( replicate(3, garch11.spec) ), 
                           dccOrder = c(1,1), 
                           distribution = "mvnorm")

ret_DCC_weekly <- list()
for (i in seq(from=0, to=754, by = 5)){
  dcc.fit = dccfit(dcc.garch11.spec, data = ret[(0+i):(1506+i),], out.sample = 0)
  dcc.focast = dccforecast(dcc.fit, n.ahead = 5, n.roll=0)
  x = rcov(dcc.focast)[[1]][,,1]+rcov(dcc.focast)[[1]][,,2]+rcov(dcc.focast)[[1]][,,3]+rcov(dcc.focast)[[1]][,,4]+rcov(dcc.focast)[[1]][,,5]
  w = solve(x %>% as.matrix()) %*% rep(1, 3)
  cov_DCC_weekly= w/sum(w)
  ret_DCC_weekly[i+1] = ret[(1506+i+1),1]*cov_DCC_weekly[1,1] + ret[(1506+i+1),2]*cov_DCC_weekly[2,1] + ret[(1506+i+1),3]*cov_DCC_weekly[3,1]
  ret_DCC_weekly[i+2] = ret[(1506+i+2),1]*cov_DCC_weekly[1,1] + ret[(1506+i+2),2]*cov_DCC_weekly[2,1] + ret[(1506+i+2),3]*cov_DCC_weekly[3,1]
  ret_DCC_weekly[i+3] = ret[(1506+i+3),1]*cov_DCC_weekly[1,1] + ret[(1506+i+3),2]*cov_DCC_weekly[2,1] + ret[(1506+i+3),3]*cov_DCC_weekly[3,1]
  ret_DCC_weekly[i+4] = ret[(1506+i+4),1]*cov_DCC_weekly[1,1] + ret[(1506+i+4),2]*cov_DCC_weekly[2,1] + ret[(1506+i+4),3]*cov_DCC_weekly[3,1]
  ret_DCC_weekly[i+5] = ret[(1506+i+5),1]*cov_DCC_weekly[1,1] + ret[(1506+i+5),2]*cov_DCC_weekly[2,1] + ret[(1506+i+5),3]*cov_DCC_weekly[3,1]
  }

mean(unlist(ret_DCC_weekly))*sqrt(250)
sqrt(250)*sd(unlist(ret_DCC_weekly))
mean(unlist(ret_DCC_weekly))/sd(unlist(ret_DCC_weekly))

## GO-GARCH estimation
go.garch11.spec <- gogarchspec( mean.model = list(armaOrder = c(0,0)),
                                variance.model     = list(model = "sGARCH", garchOrder = c(1,1)), 
                                distribution.model = "mvnorm")

ret_go_weekly <- list()
for (i in seq(from=0, to=754, by = 5)){
  go.fit <- gogarchfit(go.garch11.spec, ret[(0+i):(1506+i),], solver = "solnp", out.sample = 0)
  go.focast = gogarchforecast(go.fit, n.ahead = 5, n.roll = 0)

  x = rcov(go.focast)[[1]][,,1]+rcov(go.focast)[[1]][,,2]+rcov(go.focast)[[1]][,,3]+rcov(go.focast)[[1]][,,4]+rcov(go.focast)[[1]][,,5]
  w = solve(x %>% as.matrix()) %*% rep(1, 3)
  cov_go_weekly= w/sum(w)
  ret_go_weekly[i+1] = ret[(1506+i+1),1]*cov_go_weekly[1,1] + ret[(1506+i+1),2]*cov_go_weekly[2,1] + ret[(1506+i+1),3]*cov_go_weekly[3,1]
  ret_go_weekly[i+2] = ret[(1506+i+2),1]*cov_go_weekly[1,1] + ret[(1506+i+2),2]*cov_go_weekly[2,1] + ret[(1506+i+2),3]*cov_go_weekly[3,1]
  ret_go_weekly[i+3] = ret[(1506+i+3),1]*cov_go_weekly[1,1] + ret[(1506+i+3),2]*cov_go_weekly[2,1] + ret[(1506+i+3),3]*cov_go_weekly[3,1]
  ret_go_weekly[i+4] = ret[(1506+i+4),1]*cov_go_weekly[1,1] + ret[(1506+i+4),2]*cov_go_weekly[2,1] + ret[(1506+i+4),3]*cov_go_weekly[3,1]
  ret_go_weekly[i+5] = ret[(1506+i+5),1]*cov_go_weekly[1,1] + ret[(1506+i+5),2]*cov_go_weekly[2,1] + ret[(1506+i+5),3]*cov_go_weekly[3,1]
}

mean(unlist(ret_go_weekly))*sqrt(250)
sqrt(250)*sd(unlist(ret_go_weekly))
mean(unlist(ret_go_weekly))/sd(unlist(ret_go_weekly))

# Results
MVP_eigen_weekly = read.csv("~/Desktop/Master/MVP_eigen_week.csv")

plot_1 <- ggplot()+
  geom_line(as.data.frame(EW), mapping=aes(x = tail(zoo::index(ret),755), y = cumsum(EW), colour="EW"), size = 0.2) + 
  geom_line(as.data.frame(unlist(ret_DCC_weekly)), mapping=aes(x = tail(zoo::index(ret),755), y = cumsum(unlist(ret_DCC_weekly)), colour="DCC-GARCH"), size = 0.2) + 
  geom_line(as.data.frame(unlist(ret_go_weekly)), mapping=aes(x = tail(zoo::index(ret),755), y = cumsum(unlist(ret_go_weekly)),colour="GO-GARCH" ),size = 0.2) +
  geom_line(as.data.frame(MVP_eigen_weekly[2]), mapping=aes(x = tail(zoo::index(ret),756), y = cumsum(t(X0)), colour="λ-GARCH"), size = 0.2) + 
  scale_color_manual(name = "MGARCH-Models", values = c(EW ="lightblue", "GO-GARCH" = "darkblue", "DCC-GARCH" = "red","λ-GARCH" = "green")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_1

length(t(MVP_eigen_weekly[1]))
plot(cumsum(t(MVP_eigen_weekly[2])))

dev.off()

## Equal Sr Test: MGARCH vs EW

EW = apply(rep(rep(1/3, 3),755) * (ret)%>% tail(755), 1, sum)
boot.time.inference(cbind(as.numeric(ret_DCC_weekly), EW), M = 755, b=4)
boot.time.inference(cbind(as.numeric(ret_go_weekly), EW), M = 755, b=4)
boot.time.inference(cbind(as.numeric(unlist(MVP_eigen_weekly[2])), EW), M = 755, b=4)

## Equal Sr Test: Eigen vs MGARCH
boot.time.inference(cbind(as.numeric(tail(unlist(MVP_eigen_weekly[2]),755)), as.numeric(ret_DCC_weekly)), M = 755, b=4)
boot.time.inference(cbind(as.numeric(tail(unlist(MVP_eigen_weekly[2]),755)), as.numeric(ret_go_weekly)), M = 755, b=4)

## Discussion
skewness(ret_dcc)
skewness(ret_go)
skewness(EW)
skewness(MVP_eigen)

kurtosis(ret_dcc)
kurtosis(ret_go)
kurtosis(EW)
kurtosis(MVP_eigen)

# Plot weights
wieght_eigen_1 = read.csv("~/Desktop/Master/weights_eigen_0.csv")
wieght_eigen_2 = read.csv("~/Desktop/Master/weights_eigen_1.csv")
wieght_eigen_3 = read.csv("~/Desktop/Master/weights_eigen_2.csv")
wieght_eigen_1

plot_4 = ggplot()+
  geom_line(as.data.frame(rep(1/3, 755)), mapping=aes(x = tail(zoo::index(ret),755), y = rep(1/3, 755), colour="EW"), size = 0.2) + 
  geom_line(as.data.frame(pred_min_dcc), mapping=aes(x = tail(zoo::index(ret),755), y = V1, colour="DCC-GARCH"), size = 0.2) + 
  geom_line(as.data.frame(pred_min_go), mapping=aes(x = tail(zoo::index(ret),755), y = V1,colour="GO-GARCH" ),size = 0.2) +
  geom_line(tail(as.data.frame(wieght_eigen_1),755), mapping=aes(x = tail(zoo::index(ret),755), y = X0, colour="λ-GARCH"), size = 0.2) + 
  ggtitle("Weights: JPM") +
  scale_color_manual(name = "MGARCH-Models", values = c(EW ="lightblue", "GO-GARCH" = "darkblue", "DCC-GARCH" = "red","λ-GARCH" = "green")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_3 = ggplot()+
  geom_line(as.data.frame(rep(1/3, 755)), mapping=aes(x = tail(zoo::index(ret),755), y = rep(1/3, 755)),colour="lightblue", size = 0.2) + 
  geom_line(as.data.frame(pred_min_dcc), mapping=aes(x = tail(zoo::index(ret),755), y = V1), colour="red", size = 0.2) + 
  geom_line(as.data.frame(pred_min_go), mapping=aes(x = tail(zoo::index(ret),755), y = V1), colour="darkblue", size = 0.2) +
  geom_line(tail(as.data.frame(wieght_eigen_1),755), mapping=aes(x = tail(zoo::index(ret),755), y = X0),colour="green", size = 0.2) + 
  ggtitle("(c) Weights: Coca-Cola") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_2 = ggplot()+
  geom_line(as.data.frame(rep(1/3, 755)), mapping=aes(x = tail(zoo::index(ret),755), y = rep(1/3, 755)),colour="lightblue", size = 0.2) + 
  geom_line(as.data.frame(pred_min_dcc), mapping=aes(x = tail(zoo::index(ret),755), y = V2), colour="red", size = 0.2) + 
  geom_line(as.data.frame(pred_min_go), mapping=aes(x = tail(zoo::index(ret),755), y = V2), colour="darkblue" , size = 0.2) +
  geom_line(tail(as.data.frame(wieght_eigen_2),755), mapping=aes(x = tail(zoo::index(ret),755), y = X0), colour="green", size = 0.2) + 
  ggtitle("(b) Weights: Microsoft") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_1 = ggplot()+
  geom_line(as.data.frame(rep(1/3, 755)), mapping=aes(x = tail(zoo::index(ret),755), y = rep(1/3, 755)),  colour="lightblue", size = 0.2) + 
  geom_line(as.data.frame(pred_min_dcc), mapping=aes(x = tail(zoo::index(ret),755), y = V3), colour="red", size = 0.2) + 
  geom_line(as.data.frame(pred_min_go), mapping=aes(x = tail(zoo::index(ret),755), y = V3 ), colour="darkblue" ,size = 0.2) +
  geom_line(tail(as.data.frame(wieght_eigen_3),755), mapping=aes(x = tail(zoo::index(ret),755), y = X0, colour="λ-GARCH"), colour="green", size = 0.2) + 
  ggtitle("(a) Weights: JPM") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

# Apply user-defined function to extract legend
shared_legend <- extract_legend(plot_4)

grid.arrange(arrangeGrob(plot_1, plot_2, plot_3, ncol = 1),
             shared_legend, nrow = 2, heights = c(10, 1))

### Excluding the financial crisis
## DCC estimation
ew_fc = EW[0:251]
mean(ew_fc)*sqrt(250)
sqrt(250)*sd(ew_fc)
mean(ew_fc)/sd(ew_fc)

dcc_fc = ret_dcc[0:251]
mean(dcc_fc)*sqrt(250)
sqrt(250)*sd(dcc_fc)
mean(dcc_fc )/sd(dcc_fc )

go_fc = ret_go[0:251]
mean(go_fc)*sqrt(250)
sqrt(250)*sd(go_fc)
mean(go_fc)/sd(go_fc)

eigen_fc = MVP_eigen[0:251,]["X0"]
mean(unlist(eigen_fc))*sqrt(250)
sqrt(250)*sd(unlist(eigen_fc))
mean(unlist(eigen_fc))/sd(unlist(eigen_fc))

## Equal Sr Test: MGARCH vs EW

EW = apply(rep(rep(1/3, 3),755) * (ret)%>% tail(755), 1, sum)
boot.time.inference(cbind(as.numeric(dcc_fc), EW), M = 755, b=4)
boot.time.inference(cbind(as.numeric(go_fc ), EW), M = 755, b=4)
boot.time.inference(cbind(as.numeric(unlist(eigen_fc)), EW), M = 755, b=4)

## Equal Sr Test: Eigen vs MGARCH
boot.time.inference(cbind(as.numeric(tail(unlist(eigen_fc),755)), as.numeric(dcc_fc)), M = 755, b=4)
boot.time.inference(cbind(as.numeric(tail(unlist(eigen_fc),755)), as.numeric(go_fc)), M = 755, b=4)




