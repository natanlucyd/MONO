# =============================================================
# Monografia II - Estimacao definitiva
#   1. Series mensais e semanais (memoria + proxy de demanda de IA)
#   2. Estacionariedade (ADF e KPSS)
#   3. VAR mensal bivariado (DDR4, DDR5) com dummy de crise
#   4. VAR trivariado (NVDA, DDR5, DDR4) mensal e semanal
#   5. IRF, IRF acumulada e decomposicao de variancia (FEVD)
#   6. Robustez: ordenacoes de Cholesky, janelas amostrais,
#      especificacoes alternativas (defasagens, media, segmentos,
#      Toda-Yamamoto, Johansen), janelas longas da NVIDIA e
#      intervalos de previsao por bootstrap (Liu, 2007)
#
# Bases (pasta raiz do projeto):
#   Ultimate_Memory_Shortage_Crisis_Dataset_10k.csv
#   nvidia_stock_data_1999_2026.csv
#   ram_pricing_intelligence_2026.csv
#
# Rodar a partir da pasta raiz do projeto:
#   Rscript monografia_II/estimacao_definitiva.R
# Saidas: monografia_II/resultados_definitivos.txt,
#         monografia_II/figuras/*.png, monografia_II/tabelas/*.csv
# =============================================================

suppressMessages({
  library(vars)        # VAR, VARselect, irf, fevd, causality
  library(tseries)     # adf.test, kpss.test
  library(urca)        # ca.jo (Johansen)
  library(strucchange) # Fstats, breakpoints
})

invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))  # acentos nos graficos
set.seed(2026)
dir_out <- "monografia_II"
dir.create(file.path(dir_out, "figuras"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_out, "tabelas"), showWarnings = FALSE, recursive = TRUE)
sink(file.path(dir_out, "resultados_definitivos.txt"), split = TRUE)

# Paleta (azul, laranja, verde-agua) e estilo dos graficos
cor1 <- "#2a78d6"; cor2 <- "#eb6834"; cor3 <- "#1baf7a"
tinta <- "#3a3a38"; grade <- "#e4e3df"
png_fig <- function(nome, w = 8, h = 4.6) {
  png(file.path(dir_out, "figuras", nome), width = w, height = h,
      units = "in", res = 200)
  par(mar = c(4.2, 4.6, 2.6, 1), col.axis = tinta, col.lab = tinta,
      col.main = tinta, fg = tinta, las = 1, cex.axis = 0.85, cex.main = 1)
}
grid_h <- function() abline(h = axTicks(2), col = grade, lwd = 0.8)

titulo <- function(x) cat("\n\n", strrep("=", 70), "\n", x, "\n",
                          strrep("=", 70), "\n", sep = "")

# Data de corte: ultimo mes completo antes da redacao (out/2026).
# Registros com data posterior sao tratados como projecoes da base
# e ficam fora da amostra principal (usados apenas na robustez).
corte <- as.Date("2026-09-30")

# =============================================================
# 1. DADOS
# =============================================================
titulo("1. CONSTRUCAO DAS SERIES")

crise_all <- read.csv("Ultimate_Memory_Shortage_Crisis_Dataset_10k.csv",
                      stringsAsFactors = FALSE, na.strings = c("", "NA"))
crise_all$data <- as.Date(crise_all$timestamp)
crise_all <- subset(crise_all, generation %in% c("DDR4", "DDR5") &
                      price_status != "Listing-Error" &
                      data >= as.Date("2024-01-01"))
crise_all$mes <- as.Date(format(crise_all$data, "%Y-%m-01"))
crise_all$sem <- as.Date(cut(crise_all$data, breaks = "week"))
crise_all$seg <- ifelse(crise_all$market_segment == "Enterprise/AI", "ent", "con")

crise <- subset(crise_all, data <= corte)
cat("Registros validos (jan/2024-set/2026):", nrow(crise), "\n")
print(table(crise$generation))
cat("Registros datados apos o corte (excluidos da amostra principal):",
    sum(crise_all$data > corte), "\n")

# Agrega por periodo (mediana ou media) e gera dummy de crise endogena
agrega <- function(df, per, fun = median) {
  a <- aggregate(df$price_per_gb, list(per = df[[per]], g = df$generation), fun)
  a <- reshape(a, idvar = "per", timevar = "g", direction = "wide")
  names(a) <- c("per", "ddr4", "ddr5")
  sh <- aggregate(df$price_status == "Shortage-Inflated", list(per = df[[per]]), mean)
  names(sh)[2] <- "share"
  a <- merge(a, sh)
  a$crise <- as.numeric(a$share > 0.5)
  na.omit(a[order(a$per), ])
}

mens <- agrega(crise, "mes")
sema <- agrega(crise, "sem")
cat("Meses:", nrow(mens), "| Semanas:", nrow(sema), "\n")
cat("Primeiro mes de crise:", format(min(mens$per[mens$crise == 1])), "\n")

# Series DDR5 por segmento de mercado (Enterprise/AI x Consumer)
seg <- aggregate(price_per_gb ~ mes + seg, subset(crise, generation == "DDR5"), median)
seg <- reshape(seg, idvar = "mes", timevar = "seg", direction = "wide")
names(seg) <- c("per", "ddr5_con", "ddr5_ent")
mens <- merge(mens, seg)

# Indicadores de oferta da propria base (media mensal)
of <- aggregate(cbind(fab_utilization_rate, global_inventory_weeks,
                      gpu_hbm_trend_gb) ~ mes, crise, mean)
names(of) <- c("per", "fab_util", "estoque_sem", "hbm_gb")
mens <- merge(mens, of)

# NVIDIA: retorno log diario -> mensal e semanal
nvda <- read.csv("nvidia_stock_data_1999_2026.csv", stringsAsFactors = FALSE,
                 na.strings = c("", "NA"))
nvda$date <- as.Date(nvda$date, format = "%d/%m/%Y")
nvda <- nvda[order(nvda$date), ]
nvda$ret <- c(NA, diff(log(nvda$close)))
nvda$mes <- as.Date(format(nvda$date, "%Y-%m-01"))
nvda$sem <- as.Date(cut(nvda$date, breaks = "week"))
cat("NVIDIA: de", format(min(nvda$date)), "a", format(max(nvda$date)), "\n")

nv_m <- aggregate(cbind(ret, close) ~ mes, nvda,
                  function(x) c(sum = sum(x), last = tail(x, 1)))
nv_m <- data.frame(per = nv_m$mes, nvda = nv_m$ret[, "sum"],
                   nvda_close = nv_m$close[, "last"])
nv_m$nvda_n <- as.vector(table(nvda$mes)[as.character(nv_m$per)])
nv_s <- aggregate(ret ~ sem, nvda, sum); names(nv_s) <- c("per", "nvda")
nv_s$n <- as.vector(table(nvda$sem)[as.character(nv_s$per)])

# Apenas periodos completos de negociacao (mar/2026 so tem 8 pregoes)
nv_m <- subset(nv_m, nvda_n >= 15)
nv_s <- subset(nv_s, n >= 4 | per < max(per))

base_m <- merge(mens, nv_m)          # jan/2024 - fev/2026
base_s <- merge(sema, nv_s)          # jan/2024 - mar/2026
cat("Base trivariada mensal:", format(min(base_m$per)), "a",
    format(max(base_m$per)), "(", nrow(base_m), "meses )\n")
cat("Base trivariada semanal:", format(min(base_s$per)), "a",
    format(max(base_s$per)), "(", nrow(base_s), "semanas )\n")

write.csv(mens, file.path(dir_out, "tabelas", "serie_mensal_memoria.csv"), row.names = FALSE)
write.csv(base_m, file.path(dir_out, "tabelas", "base_trivariada_mensal.csv"), row.names = FALSE)
write.csv(base_s, file.path(dir_out, "tabelas", "base_trivariada_semanal.csv"), row.names = FALSE)

# --- Validacao cruzada com a base de marketplace -------------
titulo("1b. VALIDACAO COM A BASE DE MARKETPLACE (mar-abr/2026)")
mkt <- read.csv("ram_pricing_intelligence_2026.csv", stringsAsFactors = FALSE,
                na.strings = c("", "NA"))
mkt$data <- as.Date(substr(mkt$lastUpdated, 1, 10))
mkt$ppg <- mkt$price / mkt$capacity_gb
mkt <- subset(mkt, ram_generation %in% c("DDR4", "DDR5") & is.finite(ppg) &
                ppg > 0 & data >= as.Date("2026-03-01"))
val <- data.frame(
  geracao = c("DDR4", "DDR5"),
  n_mkt = as.vector(table(mkt$ram_generation)[c("DDR4", "DDR5")]),
  mediana_mkt = tapply(mkt$ppg, mkt$ram_generation, median)[c("DDR4", "DDR5")],
  mediana_mkt_novo = tapply(mkt$ppg[mkt$condition == "New"],
                            mkt$ram_generation[mkt$condition == "New"], median)[c("DDR4", "DDR5")],
  mediana_crise = tapply(crise$price_per_gb[crise$mes %in% as.Date(c("2026-03-01", "2026-04-01"))],
                         crise$generation[crise$mes %in% as.Date(c("2026-03-01", "2026-04-01"))],
                         median)[c("DDR4", "DDR5")])
print(val, row.names = FALSE)
cat("Razao DDR5/DDR4 (marketplace):", round(val$mediana_mkt[2] / val$mediana_mkt[1], 2),
    "| (base de crise):", round(val$mediana_crise[2] / val$mediana_crise[1], 2), "\n")
write.csv(val, file.path(dir_out, "tabelas", "validacao_marketplace.csv"), row.names = FALSE)

# --- Descritivas mensais -------------------------------------
titulo("1c. DESCRITIVAS MENSAIS")
desc <- aggregate(cbind(ddr4, ddr5, ddr5_ent, ddr5_con) ~ crise, mens, mean)
print(desc)
cat("Variacao % entre regimes:\n")
print(round(100 * (desc[2, -1] / desc[1, -1] - 1), 1))
var_a <- function(x, d1, d0) 100 * (x[mens$per == d1] / x[mens$per == d0] - 1)
cat("Variacao % set/2026 vs dez/2024 -> DDR4:",
    round(var_a(mens$ddr4, as.Date("2026-09-01"), as.Date("2024-12-01")), 1),
    "| DDR5:", round(var_a(mens$ddr5, as.Date("2026-09-01"), as.Date("2024-12-01")), 1), "\n")
cat("Salto dez/2024 -> jan/2025 (DDR5, %):",
    round(var_a(mens$ddr5, as.Date("2025-01-01"), as.Date("2024-12-01")), 1), "\n")
cat("Correlacoes (niveis) entre preco DDR5 e indicadores de oferta:\n")
print(round(cor(mens[, c("ddr5", "fab_util", "estoque_sem", "hbm_gb")]), 3))
cat("Desvio-padrao da 1a diferenca dos indicadores de oferta:\n")
print(round(sapply(mens[, c("fab_util", "estoque_sem", "hbm_gb")], function(x) sd(diff(x))), 4))
d_of <- sapply(mens[, c("fab_util", "estoque_sem", "hbm_gb")], diff)
cat("R2 de cada indicador regredido em tendencia linear:\n")
print(round(sapply(mens[, c("fab_util", "estoque_sem", "hbm_gb")],
                   function(x) summary(lm(x ~ seq_along(x)))$r.squared), 4))

# Figura 8: series mensais
png_fig("fig08_series_mensais.png", h = 4.4)
plot(mens$per, mens$ddr5, type = "n", ylim = c(0, max(mens$ddr5_ent) * 1.05),
     xlab = "Mês", ylab = "Preço mediano por GB (US$)",
     main = "Preço mensal mediano por GB: DDR4 e DDR5 (jan/2024 - set/2026)",
     frame.plot = FALSE)
grid_h()
rect(min(mens$per[mens$crise == 1]) - 15, -10, max(mens$per) + 20, 1e3,
     col = adjustcolor("#9a9890", 0.10), border = NA)
lines(mens$per, mens$ddr5_ent, col = cor1, lwd = 1.4, lty = 3)
lines(mens$per, mens$ddr5_con, col = cor1, lwd = 1.4, lty = 2)
lines(mens$per, mens$ddr5, col = cor1, lwd = 2.2)
lines(mens$per, mens$ddr4, col = cor2, lwd = 2.2)
text(min(mens$per[mens$crise == 1]), max(mens$ddr5_ent) * 1.03,
     " regime de crise", adj = 0, cex = 0.8, col = tinta)
legend("topleft", bty = "n", cex = 0.8, lwd = c(2.2, 1.4, 1.4, 2.2),
       lty = c(1, 3, 2, 1), col = c(cor1, cor1, cor1, cor2),
       legend = c("DDR5 (todas)", "DDR5 Enterprise/IA", "DDR5 Consumidor", "DDR4"))
dev.off()

# Figura 9: NVIDIA (indice) x DDR5 x DDR4 (indice jan/2024 = 100)
idx <- function(x) 100 * x / x[1]
png_fig("fig09_indices_nvda_memoria.png", h = 4.4)
plot(base_m$per, idx(base_m$ddr5), type = "n",
     ylim = range(c(idx(base_m$ddr5), idx(base_m$ddr4), idx(base_m$nvda_close))),
     xlab = "Mês", ylab = "Índice (jan/2024 = 100)",
     main = "Proxy de demanda de IA e preços de memória (índices)",
     frame.plot = FALSE)
grid_h()
lines(base_m$per, idx(base_m$nvda_close), col = cor3, lwd = 2.2)
lines(base_m$per, idx(base_m$ddr5), col = cor1, lwd = 2.2)
lines(base_m$per, idx(base_m$ddr4), col = cor2, lwd = 2.2)
legend("topleft", bty = "n", cex = 0.8, lwd = 2.2, col = c(cor3, cor1, cor2),
       legend = c("NVIDIA (fechamento)", "DDR5 (US$/GB)", "DDR4 (US$/GB)"))
dev.off()

# =============================================================
# 2. ESTACIONARIEDADE
# =============================================================
titulo("2. ESTACIONARIEDADE (series mensais)")
testa <- function(x, nome) {
  a <- suppressWarnings(adf.test(x, k = 1))
  k <- suppressWarnings(kpss.test(x, null = "Level"))
  data.frame(serie = nome, n = length(x), adf_est = round(a$statistic, 3),
             adf_p = round(a$p.value, 3), kpss_est = round(k$statistic, 3),
             kpss_p = round(k$p.value, 3))
}
est <- rbind(
  testa(log(mens$ddr4), "log DDR4 (nivel)"),
  testa(log(mens$ddr5), "log DDR5 (nivel)"),
  testa(log(base_m$nvda_close), "log NVDA (nivel)"),
  testa(diff(log(mens$ddr4)), "dlog DDR4"),
  testa(diff(log(mens$ddr5)), "dlog DDR5"),
  testa(base_m$nvda, "retorno NVDA"))
print(est, row.names = FALSE)
cat("Obs.: ADF com 1 defasagem (amostra curta); H0 ADF = raiz unitaria;",
    "H0 KPSS = estacionaria em nivel.\n")
cat("dlog DDR5 com quebra (ADF na serie descontada da dummy de crise):\n")
r_d5 <- residuals(lm(diff(log(mens$ddr5)) ~ mens$crise[-1]))
print(suppressWarnings(adf.test(r_d5, k = 1)))
nv_long <- aggregate(ret ~ mes, nvda, sum)
nv_long <- subset(nv_long, mes > min(mes) & mes <= as.Date("2026-02-01"))
est <- rbind(est, testa(nv_long$ret, "retorno NVDA (1999-2026)"))
print(tail(est, 1), row.names = FALSE)
write.csv(est, file.path(dir_out, "tabelas", "estacionariedade_mensal.csv"), row.names = FALSE)

# =============================================================
# 3. VAR MENSAL BIVARIADO (DDR4, DDR5) - jan/2024 a set/2026
# =============================================================
titulo("3. VAR MENSAL BIVARIADO (d_ddr4, d_ddr5) + dummy de crise")
y2 <- data.frame(d_ddr4 = diff(log(mens$ddr4)), d_ddr5 = diff(log(mens$ddr5)))
ex2 <- matrix(mens$crise[-1], ncol = 1, dimnames = list(NULL, "crise"))
sel2 <- VARselect(y2, lag.max = 3, type = "const", exogen = ex2)$selection
print(sel2)
p2 <- as.integer(sel2["SC(n)"])   # todos os criterios indicam p = 2
var2 <- VAR(y2, p = p2, type = "const", exogen = ex2)
print(summary(var2))
print(serial.test(var2, lags.pt = 8, type = "PT.asymptotic"))
print(normality.test(var2)$jb.mul$JB)
print(arch.test(var2, lags.multi = 3)$arch.mul)
cat("Raizes:", round(roots(var2), 3), "\n")
print(causality(var2, cause = "d_ddr4"))
print(causality(var2, cause = "d_ddr5"))

# =============================================================
# 4. VAR TRIVARIADO - ordenacao NVDA -> DDR5 -> DDR4
# =============================================================
# Justificativa: a demanda de IA (NVDA) e a variavel mais exogena
# ao mercado de memoria no impacto; DDR5 esta na linha de frente
# de servidores/IA; DDR4 (geracao legada) responde por ultimo.
y3m <- data.frame(nvda = base_m$nvda[-1],
                  d_ddr5 = diff(log(base_m$ddr5)),
                  d_ddr4 = diff(log(base_m$ddr4)))
ex3m <- matrix(base_m$crise[-1], ncol = 1, dimnames = list(NULL, "crise"))
y3s <- data.frame(nvda = base_s$nvda[-1],
                  d_ddr5 = diff(log(base_s$ddr5)),
                  d_ddr4 = diff(log(base_s$ddr4)))
ex3s <- matrix(base_s$crise[-1], ncol = 1, dimnames = list(NULL, "crise"))

diagn <- function(m, lags) {
  s <- serial.test(m, lags.pt = lags, type = "PT.asymptotic")$serial
  j <- normality.test(m)$jb.mul$JB
  a <- arch.test(m, lags.multi = 2)$arch.mul
  data.frame(portmanteau_p = round(s$p.value, 3), jb_p = round(j$p.value, 3),
             arch_p = round(a$p.value, 3), max_raiz = round(max(roots(m)), 3))
}
granger_tab <- function(m) {
  vs <- colnames(m$y)
  do.call(rbind, lapply(vs, function(v) {
    g <- causality(m, cause = v)
    data.frame(causa = v, F = round(g$Granger$statistic, 3),
               p_granger = round(g$Granger$p.value, 3),
               chi2_inst = round(g$Instant$statistic, 3),
               p_inst = round(g$Instant$p.value, 3))
  }))
}
coef_tab <- function(m) {
  do.call(rbind, lapply(names(m$varresult), function(eq) {
    cf <- summary(m$varresult[[eq]])$coefficients
    data.frame(equacao = eq, termo = rownames(cf), coef = round(cf[, 1], 4),
               ep = round(cf[, 2], 4), p = round(cf[, 4], 3))
  }))
}

write.csv(coef_tab(var2), file.path(dir_out, "tabelas", "coef_var2_mensal.csv"), row.names = FALSE)
write.csv(granger_tab(var2), file.path(dir_out, "tabelas", "granger_bivariado_mensal.csv"), row.names = FALSE)

titulo("4a. VAR TRIVARIADO MENSAL (jan/2024 - fev/2026)")
print(VARselect(y3m, lag.max = 3, type = "const", exogen = ex3m))
var3m <- VAR(y3m, p = 1, type = "const", exogen = ex3m)
print(summary(var3m))
d3m <- diagn(var3m, 6); print(d3m)
g3m <- granger_tab(var3m); print(g3m, row.names = FALSE)
write.csv(coef_tab(var3m), file.path(dir_out, "tabelas", "coef_var3_mensal.csv"), row.names = FALSE)

titulo("4b. VAR TRIVARIADO SEMANAL (jan/2024 - mar/2026)")
sel_s <- VARselect(y3s, lag.max = 8, type = "const", exogen = ex3s)
print(sel_s)
p_s <- as.integer(sel_s$selection["SC(n)"])
var3s <- VAR(y3s, p = p_s, type = "const", exogen = ex3s)
print(summary(var3s))
d3s <- diagn(var3s, 12); print(d3s)
g3s <- granger_tab(var3s); print(g3s, row.names = FALSE)
write.csv(coef_tab(var3s), file.path(dir_out, "tabelas", "coef_var3_semanal.csv"), row.names = FALSE)

write.csv(rbind(cbind(freq = "mensal", g3m), cbind(freq = "semanal", g3s)),
          file.path(dir_out, "tabelas", "granger_trivariado.csv"), row.names = FALSE)
write.csv(rbind(cbind(freq = "mensal", p = 1, n = var3m$obs, d3m),
                cbind(freq = "semanal", p = p_s, n = var3s$obs, d3s),
                cbind(freq = "mensal_bivar", p = p2, n = var2$obs, diagn(var2, 8))),
          file.path(dir_out, "tabelas", "diagnosticos.csv"), row.names = FALSE)

# =============================================================
# 5. IRF E FEVD
# =============================================================
titulo("5. FUNCOES IMPULSO-RESPOSTA E DECOMPOSICAO DE VARIANCIA")
H_m <- 12; H_s <- 12
irf_m  <- irf(var3m, impulse = "nvda", n.ahead = H_m, boot = TRUE, runs = 1000, ci = 0.90)
irfc_m <- irf(var3m, impulse = "nvda", n.ahead = H_m, boot = TRUE, runs = 1000,
              ci = 0.90, cumulative = TRUE)
irf_s  <- irf(var3s, impulse = "nvda", n.ahead = H_s, boot = TRUE, runs = 1000, ci = 0.90)
irfc_s <- irf(var3s, impulse = "nvda", n.ahead = H_s, boot = TRUE, runs = 1000,
              ci = 0.90, cumulative = TRUE)
irf_own_m <- irf(var3m, impulse = "d_ddr5", response = c("d_ddr5", "d_ddr4"),
                 n.ahead = H_m, boot = TRUE, runs = 1000, ci = 0.90)

tab_irf <- function(ir, resp, hs) {
  data.frame(h = hs, resposta = resp,
             irf = round(ir$irf$nvda[hs + 1, resp], 4),
             inf = round(ir$Lower$nvda[hs + 1, resp], 4),
             sup = round(ir$Upper$nvda[hs + 1, resp], 4))
}
hs <- c(0, 1, 2, 3, 6, 12)
ti <- rbind(cbind(freq = "mensal", tipo = "pontual", rbind(tab_irf(irf_m, "d_ddr5", hs), tab_irf(irf_m, "d_ddr4", hs))),
            cbind(freq = "mensal", tipo = "acumulada", rbind(tab_irf(irfc_m, "d_ddr5", hs), tab_irf(irfc_m, "d_ddr4", hs))),
            cbind(freq = "semanal", tipo = "pontual", rbind(tab_irf(irf_s, "d_ddr5", hs), tab_irf(irf_s, "d_ddr4", hs))),
            cbind(freq = "semanal", tipo = "acumulada", rbind(tab_irf(irfc_s, "d_ddr5", hs), tab_irf(irfc_s, "d_ddr4", hs))))
print(ti, row.names = FALSE)
write.csv(ti, file.path(dir_out, "tabelas", "irf_choque_nvda.csv"), row.names = FALSE)
cat("Desvio-padrao do choque NVDA (Cholesky) - mensal:",
    round(sqrt(summary(var3m)$covres["nvda", "nvda"]), 4),
    "| semanal:", round(sqrt(summary(var3s)$covres["nvda", "nvda"]), 4), "\n")
cat("IRF propria DDR5 (mensal) h=0..3:\n")
print(round(cbind(irf_own_m$irf$d_ddr5[1:4, ], irf_own_m$Lower$d_ddr5[1:4, ],
                  irf_own_m$Upper$d_ddr5[1:4, ]), 4))

fevd_tab <- function(m, hs, freq) {
  f <- fevd(m, n.ahead = max(hs))
  do.call(rbind, lapply(names(f), function(v)
    data.frame(freq = freq, variavel = v, h = hs, round(100 * f[[v]][hs, , drop = FALSE], 1))))
}
fv <- rbind(fevd_tab(var3m, c(1, 3, 6, 12), "mensal"),
            fevd_tab(var3s, c(1, 4, 12), "semanal"))
print(fv, row.names = FALSE)
write.csv(fv, file.path(dir_out, "tabelas", "fevd_trivariado.csv"), row.names = FALSE)

# Figura 10: IRF (pontual e acumulada) - mensal
plot_irf <- function(ir, resp, titulo_g, cor, xlab) {
  h <- 0:(nrow(ir$irf$nvda) - 1)
  lo <- ir$Lower$nvda[, resp]; hi <- ir$Upper$nvda[, resp]; md <- ir$irf$nvda[, resp]
  plot(h, md, type = "n", ylim = range(c(lo, hi, 0)), xlab = xlab,
       ylab = "Resposta (log-pontos)", main = titulo_g, frame.plot = FALSE)
  grid_h()
  polygon(c(h, rev(h)), c(lo, rev(hi)), col = adjustcolor(cor, 0.18), border = NA)
  abline(h = 0, col = "#9a9890", lwd = 1)
  lines(h, md, col = cor, lwd = 2.2)
}
png_fig("fig10_irf_mensal.png", w = 8, h = 6.2)
par(mfrow = c(2, 2), mar = c(4, 4.6, 2.6, 1))
plot_irf(irf_m, "d_ddr5", "DDR5: resposta mensal", cor1, "Meses")
plot_irf(irf_m, "d_ddr4", "DDR4: resposta mensal", cor2, "Meses")
plot_irf(irfc_m, "d_ddr5", "DDR5: resposta acumulada", cor1, "Meses")
plot_irf(irfc_m, "d_ddr4", "DDR4: resposta acumulada", cor2, "Meses")
dev.off()
png_fig("fig11_irf_semanal.png", w = 8, h = 6.2)
par(mfrow = c(2, 2), mar = c(4, 4.6, 2.6, 1))
plot_irf(irf_s, "d_ddr5", "DDR5: resposta semanal", cor1, "Semanas")
plot_irf(irf_s, "d_ddr4", "DDR4: resposta semanal", cor2, "Semanas")
plot_irf(irfc_s, "d_ddr5", "DDR5: resposta acumulada", cor1, "Semanas")
plot_irf(irfc_s, "d_ddr4", "DDR4: resposta acumulada", cor2, "Semanas")
dev.off()

# Figura 12: FEVD (barras empilhadas) - mensal e semanal
png_fig("fig12_fevd.png", w = 8, h = 4.8)
par(mfrow = c(1, 2), mar = c(4.2, 4.6, 2.6, 0.5), oma = c(2, 0, 0, 0))
for (fr in list(list(m = var3m, nm = "Mensal", h = 12, lb = "Meses"),
                list(m = var3s, nm = "Semanal", h = 12, lb = "Semanas"))) {
  f <- fevd(fr$m, n.ahead = fr$h)$d_ddr5 * 100
  barplot(t(f), col = c(cor3, cor1, cor2), border = "white", space = 0.25,
          names.arg = 1:fr$h, xlab = fr$lb, ylab = "% da variância (DDR5)",
          main = paste0(fr$nm, ": FEVD de DDR5"), cex.names = 0.7)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom", horiz = TRUE, bty = "n", cex = 0.8, fill = c(cor3, cor1, cor2),
       border = "white", legend = c("Choque NVDA", "Choque DDR5", "Choque DDR4"))
dev.off()

# =============================================================
# 6. ROBUSTEZ
# =============================================================
titulo("6. ROBUSTEZ")

# --- 6.1 Ordenacoes de Cholesky ------------------------------
titulo("6.1 Ordenacoes alternativas de Cholesky (FEVD h=12, % devido a NVDA)")
perms <- list(c("nvda", "d_ddr5", "d_ddr4"), c("nvda", "d_ddr4", "d_ddr5"),
              c("d_ddr5", "nvda", "d_ddr4"), c("d_ddr4", "nvda", "d_ddr5"),
              c("d_ddr5", "d_ddr4", "nvda"), c("d_ddr4", "d_ddr5", "nvda"))
ord <- do.call(rbind, lapply(perms, function(o) {
  fm <- fevd(VAR(y3m[, o], p = 1, type = "const", exogen = ex3m), n.ahead = 12)
  fs <- fevd(VAR(y3s[, o], p = p_s, type = "const", exogen = ex3s), n.ahead = 12)
  im <- irf(VAR(y3m[, o], p = 1, type = "const", exogen = ex3m), impulse = "nvda",
            n.ahead = 12, boot = FALSE, cumulative = TRUE)
  data.frame(ordenacao = paste(o, collapse = " > "),
             m_ddr5 = round(100 * fm$d_ddr5[12, "nvda"], 1),
             m_ddr4 = round(100 * fm$d_ddr4[12, "nvda"], 1),
             s_ddr5 = round(100 * fs$d_ddr5[12, "nvda"], 1),
             s_ddr4 = round(100 * fs$d_ddr4[12, "nvda"], 1),
             irf_acum12_ddr5 = round(im$irf$nvda[13, "d_ddr5"], 4))
}))
print(ord, row.names = FALSE)
write.csv(ord, file.path(dir_out, "tabelas", "rob_ordenacoes.csv"), row.names = FALSE)

# --- 6.2 Janelas amostrais e especificacoes -------------------
titulo("6.2 Janelas amostrais e especificacoes alternativas")
resumo_spec <- function(nome, y, ex, p = 1) {
  m <- if (is.null(ex)) VAR(y, p = p, type = "const") else VAR(y, p = p, type = "const", exogen = ex)
  cf5 <- summary(m$varresult$d_ddr5)$coefficients
  cri <- if ("crise" %in% rownames(cf5)) cf5["crise", c(1, 4)] else c(NA, NA)
  fv <- fevd(m, n.ahead = 12)
  vs <- colnames(y)
  nv <- vs[1]
  g <- causality(m, cause = nv)$Granger
  data.frame(especificacao = nome, n = m$obs, p = p,
             granger_nvda_p = round(g$p.value, 3),
             crise_ddr5 = round(cri[1], 4), crise_ddr5_p = round(cri[2], 3),
             fevd12_nvda_em_ddr5 = round(100 * fv[[vs[2]]][12, nv], 1),
             fevd12_nvda_em_ddr4 = round(100 * fv[[vs[3]]][12, nv], 1),
             portmanteau_p = round(serial.test(m, lags.pt = ifelse(m$obs > 60, 12, 6))$serial$p.value, 3))
}
# (a) base
rb <- list(resumo_spec("Base mensal (jan24-fev26)", y3m, ex3m))
# (b) sem dummy de crise
rb[[length(rb) + 1]] <- resumo_spec("Mensal sem dummy de crise", y3m, NULL)
# (c) apenas regime de crise
ic <- which(base_m$crise[-1] == 1)
rb[[length(rb) + 1]] <- resumo_spec("Mensal so regime de crise (jan25-fev26)", y3m[ic, ], NULL)
cat("Regime de crise (mensal): coeficientes da equacao d_ddr5\n")
print(round(summary(VAR(y3m[ic, ], p = 1, type = "const")$varresult$d_ddr5)$coefficients, 4))
cat("Regime de crise (mensal): coeficientes da equacao d_ddr4\n")
print(round(summary(VAR(y3m[ic, ], p = 1, type = "const")$varresult$d_ddr4)$coefficients, 4))
# (d) duas defasagens
rb[[length(rb) + 1]] <- resumo_spec("Mensal com p = 2", y3m, ex3m, p = 2)
# (e) media em vez de mediana
mens_mu <- agrega(crise, "mes", mean); bm_mu <- merge(mens_mu, nv_m)
y_mu <- data.frame(nvda = bm_mu$nvda[-1], d_ddr5 = diff(log(bm_mu$ddr5)), d_ddr4 = diff(log(bm_mu$ddr4)))
rb[[length(rb) + 1]] <- resumo_spec("Mensal com media (nao mediana)", y_mu,
                                     matrix(bm_mu$crise[-1], dimnames = list(NULL, "crise")))
# (f) segmentos DDR5: Enterprise/IA e Consumidor no lugar de DDR5/DDR4
y_seg <- data.frame(nvda = base_m$nvda[-1], d_ddr5 = diff(log(base_m$ddr5_ent)),
                    d_ddr4 = diff(log(base_m$ddr5_con)))
rb[[length(rb) + 1]] <- resumo_spec("Segmentos: DDR5 Enterprise/IA e DDR5 Consumidor", y_seg, ex3m)
# (g) semanal base
rb[[length(rb) + 1]] <- resumo_spec(paste0("Base semanal (p = ", p_s, ")"), y3s, ex3s, p = p_s)
# (h) semanal so crise
ics <- which(base_s$crise[-1] == 1)
rb[[length(rb) + 1]] <- resumo_spec("Semanal so regime de crise", y3s[ics, ], NULL, p = p_s)
# (i) semanal com retorno da NVDA defasado 4 semanas (antecipacao de pedidos)
nv_lag <- c(rep(NA, 4), head(base_s$nvda, -4))
y_l4 <- na.omit(data.frame(nvda = nv_lag[-1], d_ddr5 = diff(log(base_s$ddr5)),
                           d_ddr4 = diff(log(base_s$ddr4)), crise = base_s$crise[-1]))
rb[[length(rb) + 1]] <- resumo_spec("Semanal, NVDA defasada 4 semanas", y_l4[, 1:3],
                                     matrix(y_l4$crise, dimnames = list(NULL, "crise")), p = p_s)
rob <- do.call(rbind, rb)
print(rob, row.names = FALSE)
write.csv(rob, file.path(dir_out, "tabelas", "rob_especificacoes.csv"), row.names = FALSE)
cat("Obs.: na linha de segmentos, as colunas *_ddr5 referem-se a DDR5 Enterprise/IA",
    "e *_ddr4 a DDR5 Consumidor.\n")

# Bivariado mensal na amostra completa da base (inclui registros apos set/2026)
mens_full <- agrega(crise_all[crise_all$data <= as.Date("2026-12-31"), ], "mes")
y2f <- data.frame(d_ddr4 = diff(log(mens_full$ddr4)), d_ddr5 = diff(log(mens_full$ddr5)))
v2f <- VAR(y2f, p = 1, type = "const",
           exogen = matrix(mens_full$crise[-1], dimnames = list(NULL, "crise")))
cat("\nBivariado mensal, amostra completa jan/2024-dez/2026 (n =", v2f$obs, "):\n")
print(round(summary(v2f$varresult$d_ddr5)$coefficients, 4))
cat("Granger DDR4 -> DDR5 p =", round(causality(v2f, cause = "d_ddr4")$Granger$p.value, 3),
    "| DDR5 -> DDR4 p =", round(causality(v2f, cause = "d_ddr5")$Granger$p.value, 3), "\n")

# --- 6.3 Causalidade em niveis: Toda-Yamamoto e Johansen ------
titulo("6.3 Niveis: Toda-Yamamoto e cointegracao de Johansen")
L <- data.frame(l_nvda = log(base_m$nvda_close), l_ddr5 = log(base_m$ddr5),
                l_ddr4 = log(base_m$ddr4))
print(VARselect(L, lag.max = 3, type = "both")$selection)
toda_yamamoto <- function(L, p = 1, dmax = 1) {
  m <- VAR(L, p = p + dmax, type = "both")
  out <- list()
  for (eq in names(L)) for (cs in setdiff(names(L), eq)) {
    cf <- coef(m$varresult[[eq]]); V <- vcov(m$varresult[[eq]])
    idx <- which(names(cf) %in% paste0(cs, ".l", 1:p))
    W <- as.numeric(t(cf[idx]) %*% solve(V[idx, idx]) %*% cf[idx])
    out[[length(out) + 1]] <- data.frame(causa = cs, efeito = eq, wald = round(W, 3),
                                         p = round(1 - pchisq(W, length(idx)), 3))
  }
  do.call(rbind, out)
}
ty <- toda_yamamoto(L)
print(ty, row.names = FALSE)
write.csv(ty, file.path(dir_out, "tabelas", "rob_toda_yamamoto.csv"), row.names = FALSE)
jo <- ca.jo(L, type = "trace", K = 2, ecdet = "trend")
print(summary(jo))
jo_tab <- data.frame(h0 = c("r <= 2", "r <= 1", "r = 0"), estat = round(jo@teststat, 2),
                     cv10 = jo@cval[, 1], cv5 = jo@cval[, 2], cv1 = jo@cval[, 3])
write.csv(jo_tab, file.path(dir_out, "tabelas", "rob_johansen.csv"), row.names = FALSE)

# --- 6.4 Janelas longas da proxy de demanda (NVIDIA 1999-2026) -
titulo("6.4 NVIDIA em janelas longas: pre-pandemia, pandemia e pos-2022")
nm_all <- aggregate(ret ~ mes, nvda, sum)
nm_all <- subset(nm_all, mes <= as.Date("2026-02-01"))
nm_all$janela <- cut(nm_all$mes, breaks = as.Date(c("1999-01-01", "2020-01-01", "2023-01-01", "2026-03-01")),
                     labels = c("Pre-pandemia (1999-2019)", "Pandemia (2020-2022)", "Pos-2022 / IA generativa (2023-2026)"),
                     right = FALSE)
jan <- do.call(rbind, lapply(split(nm_all, nm_all$janela), function(d) {
  ar <- arima(d$ret, order = c(1, 0, 0))
  data.frame(janela = d$janela[1], meses = nrow(d),
             ret_medio_mensal = round(100 * mean(d$ret), 2),
             vol_anualizada = round(100 * sd(d$ret) * sqrt(12), 1),
             ret_acumulado = round(100 * (exp(sum(d$ret)) - 1), 1),
             ar1 = round(coef(ar)[1], 3))
}))
print(jan, row.names = FALSE)
write.csv(jan, file.path(dir_out, "tabelas", "rob_janelas_nvda.csv"), row.names = FALSE)
# Quebras estruturais na media do retorno mensal (Bai-Perron)
rts <- ts(nm_all$ret, start = c(1999, 1), frequency = 12)
fs <- Fstats(rts ~ 1, from = 0.10)
print(sctest(fs, type = "supF"))
bp <- breakpoints(rts ~ 1, h = 0.10)
print(summary(bp)$RSS)
cat("Datas de quebra (BIC):", format(nm_all$mes[bp$breakpoints]), "\n")

# --- 6.5 Intervalos de previsao por bootstrap (Liu, 2007) ------
titulo("6.5 Intervalos de previsao por bootstrap (Liu, 2007)")
# Modelo AR(p) estimado apenas no regime pre-crise (dlog semanal);
# residuos reamostrados geram B trajetorias do log-preco no regime de crise.
boot_pi <- function(lp, n_pre, B = 2000, pmax = 4) {
  dy <- diff(lp[1:n_pre])
  p <- max(1, ar(dy, order.max = pmax, aic = TRUE)$order)
  fit <- arima(dy, order = c(p, 0, 0), method = "CSS-ML")
  phi <- coef(fit)[1:p]; mu <- coef(fit)["intercept"]
  e <- na.omit(residuals(fit)); e <- e - mean(e)
  H <- length(lp) - n_pre
  sims <- matrix(NA, B, H)
  for (b in 1:B) {
    # bootstrap de residuos condicionado aos parametros estimados
    hist <- tail(dy, p) - mu; lvl <- lp[n_pre]
    for (h in 1:H) {
      dnew <- sum(phi * rev(tail(hist, p))) + sample(e, 1)
      hist <- c(hist, dnew); lvl <- lvl + mu + dnew; sims[b, h] <- lvl
    }
  }
  q <- apply(sims, 2, quantile, c(0.025, 0.5, 0.975))
  obs <- lp[(n_pre + 1):length(lp)]
  list(p = p, q = q, obs = obs, fora = mean(obs > q[3, ] | obs < q[1, ]),
       acima = mean(obs > q[3, ]))
}
n_pre_s <- sum(sema$crise == 0 & sema$per < as.Date("2025-01-06"))
bp5 <- boot_pi(log(sema$ddr5), n_pre_s)
bp4 <- boot_pi(log(sema$ddr4), n_pre_s)
cat("Semanas pre-crise usadas na estimacao:", n_pre_s, "\n")
cat("DDR5: AR(", bp5$p, ") | % semanas de crise acima da banda de 95%:",
    round(100 * bp5$acima, 1), "| fora da banda:", round(100 * bp5$fora, 1), "\n")
cat("DDR4: AR(", bp4$p, ") | % semanas de crise acima da banda de 95%:",
    round(100 * bp4$acima, 1), "| fora da banda:", round(100 * bp4$fora, 1), "\n")
cat("Preco DDR5 observado na ultima semana vs mediana prevista (US$/GB):",
    round(exp(tail(bp5$obs, 1)), 2), "vs", round(exp(tail(bp5$q[2, ], 1)), 2),
    "| limite superior 95%:", round(exp(tail(bp5$q[3, ], 1)), 2), "\n")
cat("Preco DDR4 observado na ultima semana vs mediana prevista (US$/GB):",
    round(exp(tail(bp4$obs, 1)), 2), "vs", round(exp(tail(bp4$q[2, ], 1)), 2),
    "| limite superior 95%:", round(exp(tail(bp4$q[3, ], 1)), 2), "\n")
# Primeira semana em que o preco observado sai da banda
cat("Primeira semana fora da banda - DDR5:",
    format(sema$per[n_pre_s + which(bp5$obs > bp5$q[3, ])[1]]),
    "| DDR4:", format(sema$per[n_pre_s + which(bp4$obs > bp4$q[3, ])[1]]), "\n")

# NVIDIA: AR nos retornos mensais 2010-2022, previsao 2023-2026
nl <- aggregate(close ~ mes, nvda, function(x) tail(x, 1))
nl <- subset(nl, mes >= as.Date("2010-01-01") & mes <= as.Date("2026-02-01"))
n_pre_nv <- sum(nl$mes < as.Date("2023-01-01"))
bpn <- boot_pi(log(nl$close), n_pre_nv)
cat("NVIDIA: AR(", bpn$p, ") | % meses 2023-2026 acima da banda de 95%:",
    round(100 * bpn$acima, 1), "| fora:", round(100 * bpn$fora, 1), "\n")
write.csv(data.frame(serie = c("DDR5 semanal", "DDR4 semanal", "NVIDIA mensal"),
                     ar_p = c(bp5$p, bp4$p, bpn$p),
                     n_estimacao = c(n_pre_s, n_pre_s, n_pre_nv) - 1,
                     horizonte = c(length(bp5$obs), length(bp4$obs), length(bpn$obs)),
                     pct_acima = round(100 * c(bp5$acima, bp4$acima, bpn$acima), 1),
                     pct_fora = round(100 * c(bp5$fora, bp4$fora, bpn$fora), 1)),
          file.path(dir_out, "tabelas", "rob_bootstrap_liu.csv"), row.names = FALSE)

# Figura 13: bandas de previsao por bootstrap
plot_fan <- function(datas, lp, n_pre, bpx, tit, cor, xl) {
  y <- exp(lp); q <- exp(bpx$q); dfc <- datas[(n_pre + 1):length(datas)]
  plot(datas, y, type = "n", ylim = range(c(y, q)), log = "y", xlab = xl,
       ylab = "US$ (escala log)", main = tit, frame.plot = FALSE)
  grid_h()
  polygon(c(dfc, rev(dfc)), c(q[1, ], rev(q[3, ])), col = adjustcolor("#9a9890", 0.25), border = NA)
  lines(dfc, q[2, ], col = "#6f6e69", lwd = 1.4, lty = 2)
  lines(datas, y, col = cor, lwd = 2.2)
  abline(v = datas[n_pre], col = "#9a9890", lty = 3)
}
png_fig("fig13_bootstrap_liu.png", w = 8, h = 7.4)
par(mfrow = c(3, 1), mar = c(4, 4.6, 2.4, 1))
plot_fan(sema$per, log(sema$ddr5), n_pre_s, bp5,
         "DDR5: preço observado vs. banda de 95% prevista com dinâmica pré-crise", cor1, "Semana")
plot_fan(sema$per, log(sema$ddr4), n_pre_s, bp4,
         "DDR4: preço observado vs. banda de 95% prevista com dinâmica pré-crise", cor2, "Semana")
plot_fan(nl$mes, log(nl$close), n_pre_nv, bpn,
         "NVIDIA: fechamento observado vs. banda de 95% prevista com dinâmica 2010-2022", cor3, "Mês")
dev.off()

titulo("FIM")
sessionInfo()
sink()
