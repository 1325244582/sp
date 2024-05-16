#' ODE solution
#'
#' This function solves an ODE system based on the given parameter list.
#' @param user created parameters
#' @return solution of an ODE function system
#' @export

#showtext_auto()
# 表示之后用上同样的字体
#font_add("kaishu", "simkai.ttf")
# 如果系统已经安装了该字体，可以通过这种方式调用

model <- function(t, y, param) {
  # 从参数列表中获得5种状态SEIR
  S <- y[1]
  E <- y[2]
  I1 <- y[3]
  I2 <- y[4]
  R <- y[5]
  N <- param["N"]
  # 参数列表中提取模型参数
  beta <- param['beta']
  mu <- param['mu']
  gamma1 <- param['gamma1']
  gamma2 <- param['gamma2']
  lamda <- param['lamda']
  # 基于模型参数和状态信息创建传染并传播教学
  dSt <- mu * (N - S) - beta * S * (I1 + I2) / N
  dEt <- beta * S * (I1 + I2) / N - mu * E - lamda * E
  dI1t <- lamda * E - (mu + gamma1) * I1
  dI2t <- gamma1 * I1 - (mu + gamma2) * I2
  dRt <- gamma2 * I2 - mu * R

  # 汇总模型结果
  outcome <- c(dSt, dEt, dI1t, dI2t, dRt)

  list(outcome)
}

# 设置仿真参数
times <- seq(0, 156, by = 1/7)
param <- c(mu = 0.0000, lamda = 0.04, beta = 5, gamma1 = 0.2, gamma2 = 0.06, N = 1)
init <- c(S = 0.9999, E = 0.00000, I1 = 0.00003, I2 = 0.00000, R = 0)

# 微分方程求解函数code
result <- deSolve::ode(y = init, times = times, func = model, parms = param)
result <- as.data.frame(result)

# 尾部数据变化情况
tail(round(result, 6), 10)

# 绘制图形
#' @export
seirplot <- ggplot2::ggplot(data = result) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = S, col = 'S'), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = E, col = 'E'), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = I1, col = 'I1'), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = I2, col = 'I2'), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = R, col = 'R'), lwd = 2) +
  ggplot2::labs(x = "时间", y = "比率") +
  ggplot2::scale_color_manual(name = "传染病模型仿真",
                     values = c("S" = "orange", "E" = "purple", "I1" = "pink", "I2" = "green", "R" = "skyblue"))

# 绘制图形对象绘制图形结果
seirplot

# 保存为矢量图
ggplot2::ggsave(seirplot, file = "lwx.pdf", width = 7, height = 6)
ggplot2::ggsave(seirplot, file = "lwx.svg", width = 7, height = 6)


