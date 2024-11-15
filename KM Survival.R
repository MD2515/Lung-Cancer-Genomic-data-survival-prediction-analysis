dataset1 <- read.csv("D:/BHI/Spring 2024/BHI 699/lung/lasf.csv")
# Kaplan Meier survival estimate
Y <-Surv(ut$time, ut$status)
summary(survfit(Y~1, conf.type="plain"))
plot(survfit(Y~1, conf.type="plain"), col="blue", lwd=2)
survfit(Y~1)

