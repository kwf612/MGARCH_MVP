### Load Libraries
library(tidyquant)
library(rmgarch)
library(tseries)
library(FinTS)
library(ggplot2)
require(gridExtra)
library(magrittr)
library(dplyr) 
library(tibble)

### Load data, return series
ret = read.csv("~/Desktop/Master/returns_R.csv")
length(t(ret))/3
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

Box.test(dcc.fit@mfit$stdresid[,1]^2, lag = 10, type = "Ljung")
Box.test(dcc.fit@mfit$stdresid[,2]^2, lag = 10, type = "Ljung")
Box.test(dcc.fit@mfit$stdresid[,3]^2, lag = 10, type = "Ljung")


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
plot(go.fit@mfit$residuals[,3])

go.fit@mfit$residuals[,1] %>% jarque.bera.test()
go.fit@mfit$residuals[,2] %>% jarque.bera.test()
go.fit@mfit$residuals[,3] %>% jarque.bera.test()

Box.test(go.fit@mfit$residuals[,1]^2, lag = 10, type = "Ljung")
Box.test(go.fit@mfit$residuals[,2]^2, lag = 10, type = "Ljung")
Box.test(go.fit@mfit$residuals[,3]^2, lag = 10, type = "Ljung")

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
plot_3 <- ggplot()+
  geom_line(as.data.frame(cov_MSFT_KO_go), mapping = aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_go), colour="GO-GARCH"), size = 0.2) +
  geom_line(as.data.frame(cov_MSFT_KO_dcc), mapping=aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_dcc), colour="DCC-GARCH" ), size = 0.2) + 
  geom_line(as.data.frame(cov_KO_MSFT_eigen), mapping=aes(x = zoo::index(JPM), y = as.numeric(X0), colour="λ-GARCH"), size = 0.2) + 
  ggtitle("Microsoft & Coca Cola") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  scale_color_manual( values = c("GO-GARCH" = "darkblue", "DCC-GARCH" = "red","λ-GARCH" = "green")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=10), legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_2 <- ggplot()+
  geom_line(as.data.frame(cov_JPM_KO_go), mapping = aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_go)), size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_JPM_KO_dcc), mapping=aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_dcc)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_JPM_KO_eigen), mapping=aes(x = zoo::index(JPM), y = as.numeric(X0)), size = 0.2, colour = "green") + 
  ggtitle("JPMorgan & Coca Cola") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=10)) +
  labs(x = "", y ="")

plot_1 <- ggplot()+
  geom_line(as.data.frame(cov_JPM_MSFT_go), mapping = aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_go)), size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_JPM_MSFT_dcc), mapping=aes(x = zoo::index(JPM), y = as.numeric(cov_MSFT_KO_dcc)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_JPM_MSFT_eigen), mapping=aes(x = zoo::index(JPM), y = as.numeric(X0)), size = 0.2, colour = "green") + 
  ggtitle("JPMorgan & Microsoft") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=10)) +
  labs(x = "", y ="")

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(plot_3)

grid.arrange(arrangeGrob(plot_1, plot_2, plot_1, ncol = 1),
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
plot_3 <- ggplot()+
  geom_line(as.data.frame(cov_MSFT_KO_go_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_MSFT_KO_go_f),colour="GO-GARCH" ),size = 0.2) +
  geom_line(as.data.frame(cov_MSFT_KO_dcc_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_MSFT_KO_dcc_f), colour="DCC-GARCH"), size = 0.2) + 
  geom_line(as.data.frame(cov_KO_MSFT_eigen_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(X0), colour="λ-GARCH"), size = 0.2) + 
  ggtitle("Microsoft & Coca Cola") +
  scale_color_manual(values = c("GO-GARCH" = "darkblue", "DCC-GARCH" = "red","λ-GARCH" = "green")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=10), legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_2 <- ggplot()+
  geom_line(as.data.frame(cov_JPM_KO_go_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_JPM_KO_go_f)),size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_JPM_KO_dcc_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_JPM_KO_dcc_f)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_JPM_KO_eigen_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(X0), colour="Eigen-GARCH"), size = 0.2, colour = "green") + 
  ggtitle("JPMorgan & Coca Cola") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(x = "", y ="")

plot_1 <- ggplot()+
  geom_line(as.data.frame(cov_JPM_MSFT_go_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_JPM_MSFT_go_f)), size = 0.2, colour = "darkblue") +
  geom_line(as.data.frame(cov_JPM_MSFT_dcc_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(cov_JPM_MSFT_dcc_f)), size = 0.2, colour = "red") + 
  geom_line(as.data.frame(cov_JPM_MSFT_eigen_f), mapping=aes(x = tail(zoo::index(JPM),755), y = as.numeric(X0)), size = 0.2, colour = "green") + 
  ggtitle("JPMorgan & Microsoft") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(x = "", y ="")

# Apply user-defined function to extract legend
shared_legend <- extract_legend(plot_3)

grid.arrange(arrangeGrob(plot_1, plot_2, plot_1, ncol = 1),
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
  ggtitle("GMVP: Cummulated returns") +
  scale_color_manual(name = "MGARCH-Models", values = c(EW ="lightblue", "GO-GARCH" = "darkblue", "DCC-GARCH" = "red","λ-GARCH" = "green")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),legend.position = "bottom",legend.title = element_blank()) +
  labs(x = "", y ="")

plot_1

### Test for equal Sharpe Ratios
load("/Users/silkesommerfugl/Downloads/SharpeR/Sharpe.RData")
MVP_eigen_ret = read.csv("~/Desktop/Master/MVP_eigen.csv")
## Block size: MGARCH vs EW
#block.size.calibrate(cbind(ret_dcc, EW))
#block.size.calibrate(cbind(ret_go, EW))
#block.size.calibrate(cbind(MVP_eigen, EW))

## Block size: Eigen vs MGARCH
#block.size.calibrate(cbind(MVP_eigen, ret_dcc))
#block.size.calibrate(cbind(MVP_eigen, ret_go))

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
  geom_line(as.data.frame(EW), mapping=aes(x = tail(zoo::index(ret),756), y = cumsum(EW), colour="EW"), size = 0.2) + 
  geom_line(as.data.frame(unlist(ret_DCC_weekly)), mapping=aes(x = tail(zoo::index(ret),755), y = cumsum(unlist(ret_DCC_weekly)), colour="DCC-GARCH"), size = 0.2) + 
  geom_line(as.data.frame(unlist(ret_go_weekly)), mapping=aes(x = tail(zoo::index(ret),755), y = cumsum(unlist(ret_go_weekly)),colour="GO-GARCH" ),size = 0.2) +
  geom_line(as.data.frame(MVP_eigen_weekly[2]), mapping=aes(x = tail(zoo::index(ret),756), y = cumsum(t(X0)), colour="λ-GARCH"), size = 0.2) + 
  ggtitle("GMVP: Cummulated returns") +
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

