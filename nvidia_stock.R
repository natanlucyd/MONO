# =============================================================
# Monografia - NVIDIA: do GPU gamer a era da IA generativa
# Base: nvidia_stock_data_1999_2026.csv
#   Cotacao diaria (OHLCV) de 1999 a 2026, com market cap,
#   receita trimestral, medias moveis, era e eventos-chave.
#
# Objetivo: criar o data frame e medir o impacto da entrada
# no mercado de IA sobre a acao nos ultimos anos.
# =============================================================

# --- 1. Carregar os dados ------------------------------------
nvda <- read.csv("nvidia_stock_data_1999_2026.csv",
                 stringsAsFactors = FALSE,
                 na.strings = c("", "NA"))

# Data vem como dd/mm/aaaa
nvda$date <- as.Date(nvda$date, format = "%d/%m/%Y")

# Era como fator ORDENADO cronologicamente
nvda$era <- factor(nvda$era,
                   levels = c("Pre-GPU Era", "Gaming GPU Era",
                              "Deep Learning Era", "Data Center Growth",
                              "Generative AI Era"))

nvda <- nvda[order(nvda$date), ]

# --- 2. Inspecao ---------------------------------------------
str(nvda)
summary(nvda[, c("date", "close", "volume", "market_cap_usd_bn",
                 "quarterly_revenue_usd_bn")])
nrow(nvda)
table(nvda$era)

# Periodo coberto por cada era
tapply(nvda$date, nvda$era, range)

# --- 3. Retornos ---------------------------------------------
# Retorno log diario (some-os para obter retorno acumulado)
nvda$ret <- c(NA, diff(log(nvda$close)))

# Retorno anual (%): soma dos retornos log dentro de cada ano
nvda$ano <- as.numeric(format(nvda$date, "%Y"))
ret_anual <- tapply(nvda$ret, nvda$ano,
                    function(r) 100 * (exp(sum(r, na.rm = TRUE)) - 1))
round(ret_anual, 1)

# Retorno acumulado POR ERA (%)
ret_era <- tapply(nvda$ret, nvda$era,
                  function(r) 100 * (exp(sum(r, na.rm = TRUE)) - 1))
round(ret_era, 1)

# Volatilidade anualizada por era (%): desvio-padrao * sqrt(252)
vol_era <- tapply(nvda$ret, nvda$era,
                  function(r) 100 * sd(r, na.rm = TRUE) * sqrt(252))
round(vol_era, 1)

# --- 4. Historia completa (escala log) -----------------------
plot(nvda$date, nvda$close, type = "l", log = "y",
     xlab = "Ano", ylab = "Preco de fechamento (USD, escala log)",
     main = "NVIDIA 1999-2026: cada era em escala log")
# Linhas verticais no inicio de cada era
inicio_eras <- tapply(nvda$date, nvda$era, min)
abline(v = inicio_eras, lty = 3, col = "gray40")

# --- 5. Zoom na era da IA generativa -------------------------
ia <- subset(nvda, era == "Generative AI Era")
range(ia$date)                       # quando comeca a era da IA
head(ia[, c("date", "close", "market_cap_usd_bn", "key_event")])

# Valorizacao na era da IA
100 * (tail(ia$close, 1) / ia$close[1] - 1)            # preco (%)
range(ia$market_cap_usd_bn)                            # market cap (US$ bi)

plot(ia$date, ia$close, type = "l", col = "darkgreen",
     xlab = "Data", ylab = "Fechamento (USD)",
     main = "NVIDIA na era da IA generativa")

# Eventos-chave da era da IA (1a data de cada evento)
eventos <- aggregate(date ~ key_event, data = ia, FUN = min)
eventos <- eventos[order(eventos$date), ]
eventos
abline(v = eventos$date, lty = 3, col = "gray50")

# Maior queda diaria da era (ex.: choque DeepSeek)
ia[which.min(ia$ret), c("date", "close", "ret", "key_event")]

# --- 6. Market cap e receita: o salto da IA ------------------
plot(nvda$date, nvda$market_cap_usd_bn, type = "l", log = "y",
     xlab = "Ano", ylab = "Market cap (US$ bi, escala log)",
     main = "Valor de mercado da NVIDIA")
abline(v = min(ia$date), lty = 2, col = "red")
text(min(ia$date), max(nvda$market_cap_usd_bn),
     "  era da IA", adj = 0, cex = 0.8, col = "red")

# Receita trimestral media por era (US$ bi)
tapply(nvda$quarterly_revenue_usd_bn, nvda$era, mean, na.rm = TRUE)

# =============================================================
# PROXIMOS PASSOS (ligacao com a crise de RAM):
# A demanda de IA (NVIDIA) e candidata a CAUSA da escassez de
# memoria. Para testar, agregue o retorno da NVDA por semana e
# junte com a serie semanal de precos de RAM do VAR_crise_ram.R:
#
#   nvda$semana <- as.Date(cut(nvda$date, breaks = "week"))
#   nvda_sem <- aggregate(ret ~ semana, data = nvda,
#                         FUN = function(r) sum(r, na.rm = TRUE))
#   base <- merge(serie, nvda_sem, by = "semana")  # 'serie' do outro script
#   # -> VAR trivariado: d_ddr4, d_ddr5, ret_nvda
# =============================================================
