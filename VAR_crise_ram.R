# =============================================================
# Monografia - VAR: precos de DDR4 vs DDR5 na crise de memoria
# Base: Ultimate_Memory_Shortage_Crisis_Dataset_10k.csv
#   10.000 kits de RAM (2024-2027) com preco por GB, geracao,
#   status do preco (Normal / Shortage-Inflated / Listing-Error),
#   regiao, segmento etc.
#
# Objetivo: medir como os precos inflaram na crise de escassez
# (que nos dados comeca em jan/2025 e se intensifica a partir
# de jun/2025) e como DDR4 e DDR5 interagem dinamicamente.
# =============================================================

library(vars)     # VAR, VARselect, irf, fevd
library(tseries)  # adf.test

# --- 1. Carregar e filtrar -----------------------------------
crise <- read.csv("Ultimate_Memory_Shortage_Crisis_Dataset_10k.csv",
                  stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))

crise$data <- as.Date(crise$timestamp, format = "%Y-%m-%d")

# Range pedido: 2024-2026 | so DDR4 e DDR5 | descarta erros de anuncio
crise <- subset(crise,
                data >= as.Date("2024-01-01") &
                data <= as.Date("2026-12-31") &
                generation %in% c("DDR4", "DDR5") &
                price_status != "Listing-Error")

crise$generation <- factor(crise$generation)

str(crise[, c("data", "generation", "price_per_gb",
              "price_usd", "price_status", "market_segment")])
nrow(crise)
table(crise$generation)
table(crise$price_status)

# --- 2. Agregacao SEMANAL ------------------------------------
# Mediana semanal do preco por GB (mediana e robusta a outliers,
# importante numa base com anuncios especulativos)
crise$semana <- as.Date(cut(crise$data, breaks = "week"))

semanal <- aggregate(price_per_gb ~ semana + generation,
                     data = crise, FUN = median)

serie <- reshape(semanal,
                 idvar     = "semana",
                 timevar   = "generation",
                 direction = "wide")
names(serie) <- c("semana", "ddr4", "ddr5")
serie <- serie[order(serie$semana), ]
serie <- na.omit(serie)
nrow(serie)    # numero de semanas (T do VAR)

# Dummy de crise DERIVADA DOS DADOS: semana em que mais da metade
# dos anuncios esta marcada como 'Shortage-Inflated'
share_inflado <- aggregate(price_status ~ semana, data = crise,
                           FUN = function(x) mean(x == "Shortage-Inflated"))
names(share_inflado) <- c("semana", "share_inflado")
serie <- merge(serie, share_inflado, by = "semana")
serie$crise <- as.numeric(serie$share_inflado > 0.5)

inicio_crise <- min(serie$semana[serie$crise == 1])
inicio_crise   # primeira semana de crise

# --- 3. Visualizacao: a inflacao dos precos ------------------
plot(serie$semana, serie$ddr5, type = "l", col = "red",
     ylim = range(c(serie$ddr4, serie$ddr5)),
     xlab = "Semana", ylab = "Preco mediano por GB (USD)",
     main = "Crise de memoria: preco semanal por GB (2024-2026)")
lines(serie$semana, serie$ddr4, col = "blue")
abline(v = inicio_crise, lty = 2)
text(inicio_crise, max(serie$ddr5), "  inicio da crise", adj = 0, cex = 0.8)
legend("topleft", legend = c("DDR5", "DDR4"),
       col = c("red", "blue"), lty = 1)

# --- 4. Quanto os precos inflaram? ---------------------------
# Comparacao pre-crise vs crise (medias semanais)
impacto <- aggregate(cbind(ddr4, ddr5) ~ crise, data = serie, FUN = mean)
impacto
# Variacao percentual entre os dois regimes
100 * (impacto[2, c("ddr4", "ddr5")] / impacto[1, c("ddr4", "ddr5")] - 1)

# --- 5. Estacionariedade -------------------------------------
adf.test(serie$ddr4)   # em nivel: tendencia forte -> nao estacionaria
adf.test(serie$ddr5)

# Log-diferenca semanal (~ variacao percentual do preco)
y <- data.frame(d_ddr4 = diff(log(serie$ddr4)),
                d_ddr5 = diff(log(serie$ddr5)))
adf.test(y$d_ddr4)
adf.test(y$d_ddr5)

# --- 6. VAR com a crise como variavel exogena ----------------
# A dummy de crise entra como exogena: captura o deslocamento
# de regime sem contaminar a dinamica DDR4 <-> DDR5
exo <- matrix(serie$crise[-1], ncol = 1,
              dimnames = list(NULL, "crise"))

VARselect(y, lag.max = 8, type = "const", exogen = exo)

modelo <- VAR(y, p = 1, type = "const", exogen = exo)
summary(modelo)
# O coeficiente de 'crise' em cada equacao mede o quanto a
# variacao semanal media dos precos sobe durante a crise

# Diagnostico dos residuos
serial.test(modelo, lags.pt = 12, type = "PT.asymptotic")

# Causalidade de Granger: quem lidera o repasse de precos?
causality(modelo, cause = "d_ddr4")
causality(modelo, cause = "d_ddr5")

# Impulso-resposta e decomposicao da variancia
plot(irf(modelo, n.ahead = 12, boot = TRUE))
plot(fevd(modelo, n.ahead = 12))

# =============================================================
# LEITURA PARA A MONOGRAFIA:
# - Secao 4 da o tamanho da inflacao de precos (pre vs crise)
# - O coef. da dummy 'crise' da o efeito medio semanal do regime
# - Granger/IRF mostram se o choque se propaga de uma geracao
#   para a outra (ex.: escassez de DDR5 empurrando demanda e
#   preco de DDR4, geracao antiga)
# =============================================================
