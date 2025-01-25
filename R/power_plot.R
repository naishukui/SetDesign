# supplemental figure 1
# compare calculated power and simulated power for true and misspecified models in four settings,
# binary/continuous outcomes and uncorrelated/correlated genotypes

library(ggplot2)
beta1list1<-c(-0.9,-0.75,-0.5,-0.25,-0.1,
              -0.9,-0.75,-0.5,-0.25,-0.1,
              -0.9,-0.75,-0.5,-0.25,-0.1,
              -0.9,-0.75,-0.5,-0.25,-0.1,
              -0.9,-0.75,-0.5,-0.25,-0.1)
beta2list1<-c( rep(0.9,5),rep(0.75,5),rep(0.5,5),rep(0.25,5),rep(0.1,5))

compute_power_data <- function(result_path, derive_function, derive_args) {
  result <- get(load(result_path))

  powerE2 <- sapply(result, function(x) x$powerE2)
  powerCombined <- sapply(result, function(x) x$powerCombined)
  powerSeparated <- sapply(result, function(x) x$powerSeparated)

  powerDerived <- do.call(derive_function, derive_args)

  data <- data.frame(
    Derived = powerDerived,
    Power = c(powerE2, powerCombined, powerSeparated),
    Group = rep(c("C_misspecified", "S_misspecified", "S_true"), each = length(powerDerived))
  )

  return(data)
}

generate_power_plot <- function(data, subtitle) {
  ggplot(data, aes(x = Derived, y = Power, color = Group)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "violet", linetype = "dashed", linewidth = 1) +
    labs(
      title = subtitle,
      x = "C_true",
      y = "Other Power"
    ) +
    mytheme
}


data1 <- compute_power_data(
  result_path = "//Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/SKATsimulation/output/snpC_1.Rdata",
  derive_function = powerC_derive,
  derive_args = list(k = 50, n = 2000, alpha = 0.05, p = 0.01, list1 = beta1list1, list2 = beta2list1)
)

data2 <- compute_power_data(
  result_path = "//Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/SKATsimulation/output/snpCcor_1.Rdata",
  derive_function = powerC_deriveCor,
  derive_args = list(k = 50, n = 2000, alpha = 0.05, p = 0.05, rho = 0.15, list1 = beta1list1, list2 = beta2list1)
)

data3 <- compute_power_data(
  result_path = "//Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/SKATsimulation/output/new/snpD_1000_1.Rdata",
  derive_function = powerD_derive,
  derive_args = list(k = 50, n = 2000, alpha = 0.05, p = 0.01, list1 = beta1list1, list2 = beta2list1)
)


data4 <- compute_power_data(
  result_path = "//Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/SKATsimulation/output/new/snpDcor_1.Rdata",
  derive_function = powerD_deriveCor,
  derive_args = list(k = 50, n = 2000, alpha = 0.05, p = 0.05, rho = 0.15, list1 = beta1list1, list2 = beta2list1)
)


p1 <- generate_power_plot(data1,subtitle="Binary without correlation")
p2 <- generate_power_plot(data2,subtitle="Binary with correlation")
p3 <- generate_power_plot(data3,subtitle="Logistic without correlation")
p4 <- generate_power_plot(data4,subtitle="Logistic with correlation")


ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


