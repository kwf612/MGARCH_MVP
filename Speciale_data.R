# Load Libraries
library(tidyquant)
library(ggplot2)
library(PerformanceAnalytics)
library(tseries)
require(gridExtra)
library(grid)
library(FinTS)

# Load data, time series closing prices
symbol.vec = c("JPM", "MSFT", "KO")
getSymbols(symbol.vec, from = "2001-01-02", to = "2009-12-31")

# Extract adjusted closing prices
JPM = `JPM`[, "JPM.Adjusted", drop=F]
MSFT = `MSFT`[, "MSFT.Adjusted", drop=F]
KO = `KO`[, "KO.Adjusted", drop=F]

# Plot time series prices 
plot_1 <- ggplot(JPM, aes(x = zoo::index(JPM), y = JPM.Adjusted)) +
  geom_line(size = 0.2) + 
  ggtitle("JPM") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label= "")+
  labs(x = "", y ="Adjusted closing price")

plot_2 <-ggplot(MSFT, aes(x = zoo::index(JPM), y = `MSFT.Adjusted`)) +
  geom_line(size = 0.2) + 
  ggtitle("MSFT") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label= "")+
  labs(x = "", y ="Adjusted closing price") 

plot_3 <-ggplot(KO, aes(x = zoo::index(JPM), y = `KO.Adjusted`)) +
  geom_line(size = 0.2) + 
  ggtitle("KO") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label= "")+
  labs(x = "", y ="Adjusted closing price") 

# Calculate returns
JPM.ret = CalculateReturns(JPM, method = "discrete")*100
MSFT.ret = CalculateReturns(MSFT, method = "discrete")*100
KO.ret = CalculateReturns(KO, method = "discrete")*100

# Remove first NA observation
JPM.ret = JPM.ret[-1,]
MSFT.ret = MSFT.ret[-1,]
KO.ret = KO.ret[-1,]

# Plot returns
plot_4 <- ggplot(JPM.ret, aes(x = zoo::index(JPM.ret), y = JPM.Adjusted)) +
  geom_line(size = 0.2) + 
  ggtitle("") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y ="Returns (%)") 

plot_5 <-ggplot(MSFT.ret, aes(x = zoo::index(JPM.ret), y = `MSFT.Adjusted`)) +
  geom_line(size = 0.2) + 
  ggtitle("") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y ="Returns (%)") 

plot_6 <-ggplot(KO.ret, aes(x = zoo::index(JPM.ret), y = `KO.Adjusted`)) +
  geom_line(size = 0.2) + 
  ggtitle("") +
  scale_x_date(breaks = "2 year", date_labels ="%Y" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y ="Returns (%)") 

JPM_plot = grid.arrange(plot_1, plot_4, ncol=2, top = textGrob("(a) JPMorgan Chase & Co. (JPM)",gp=gpar(fontsize=15,font=3)))
MSFT_plot = grid.arrange(plot_2, plot_5, ncol=2, top = textGrob("(b) Microsoft Corporation (MSFT)",gp=gpar(fontsize=15,font=3)))
KO_plot = grid.arrange(plot_3, plot_6, ncol=2, top = textGrob("(c) The Coca-Cola Company (KO)",gp=gpar(fontsize=15,font=3)))

grid.arrange(JPM_plot, MSFT_plot, KO_plot, ncol=1)

# Create combined data series
ret = na.omit(merge(merge(JPM.ret,MSFT.ret),KO.ret))

write.csv(ret, file="~/Desktop/Master/returns_R.csv", row.names = FALSE)

# Summary stats
summary(ret)
statsDaily = rbind(apply(ret, 2, mean),
                   apply(ret, 2, var),
                   apply(ret, 2, sd),
                   apply(ret, 2, skewness),
                   apply(ret, 2, kurtosis))
rownames(statsDaily) = c("Mean", "Variance", "Std Dev",
                         "Skewness",  "Excess Kurtosis")
round(statsDaily, digits=4)

# Normality test
ret$JPM.Adjusted %>% jarque.bera.test()
ret$MSFT.Adjusted %>% jarque.bera.test()
ret$KO.Adjusted %>% jarque.bera.test()

# ADF test
ret$JPM.Adjusted %>% adf.test()
ret$MSFT.Adjusted  %>% adf.test()
ret$KO.Adjusted  %>% adf.test()

# ARCH test
na.omit(ret$JPM.Adjusted) %>% ArchTest(lags = 10)
na.omit(ret$MSFT.Adjusted)  %>% ArchTest(lags = 10)
na.omit(ret$KO.Adjusted)  %>% ArchTest(lags = 10)