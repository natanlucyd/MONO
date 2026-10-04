# =============================================================
# Monografia - Analise de precos de memoria RAM
# Base: ram_pricing_intelligence_2026.csv
#   ~2.900 anuncios de RAM (DDR3/DDR4/DDR5) em marketplace,
#   com preco, capacidade, condicao, vendedor e data do anuncio.
#
# ATENCAO (metodologia): a base e CROSS-SECTION (anuncios),
# nao serie temporal. Para rodar VAR e preciso AGREGAR os
# anuncios no tempo (ex.: preco medio diario por geracao).
# A secao 4 faz exatamente isso.
# =============================================================

# --- Pacotes -------------------------------------------------
# install.packages(c("vars", "tseries"))   # rodar 1x se nao tiver
library(vars)     # VAR, VARselect, irf, fevd
library(tseries)  # adf.test (raiz unitaria)

# --- 1. Carregar os dados ------------------------------------
ram <- read.csv("ram_pricing_intelligence_2026.csv",
                stringsAsFactors = FALSE,
                na.strings = c("", "NA"))

# --- 2. Limpeza / tipos --------------------------------------
ram$ram_generation <- factor(ram$ram_generation)
ram$condition      <- factor(ram$condition)
ram$brand          <- factor(ram$brand)
ram$lastUpdated    <- as.POSIXct(ram$lastUpdated,
                                 format = "%Y-%m-%d %H:%M:%S")
ram$data           <- as.Date(ram$lastUpdated)

# Preco por GB (util para comparar pentes de tamanhos diferentes)
ram$preco_por_gb <- ram$price / ram$capacity_gb

# --- 3. Primeira inspecao ------------------------------------
str(ram)          # estrutura do data frame (tipos das colunas)
head(ram)         # primeiras 6 linhas
summary(ram)      # estatisticas descritivas
nrow(ram)         # numero de observacoes

table(ram$ram_generation)              # quantos anuncios por geracao
tapply(ram$price, ram$ram_generation, summary)        # preco por geracao
tapply(ram$preco_por_gb, ram$ram_generation, median, na.rm = TRUE)

# Boxplot do preco por GB, por geracao (escala log)
boxplot(preco_por_gb ~ ram_generation, data = ram,
        log  = "y",
        xlab = "Geracao",
        ylab = "Preco por GB (USD, escala log)",
        main = "Preco por GB por geracao de RAM")

# --- 4. Agregacao em serie temporal --------------------------
# Preco medio diario por GB, separado em DDR4 e DDR5
# (DDR3 tem so 2 anuncios -> descartado)
diario <- aggregate(preco_por_gb ~ data + ram_generation,
                    data   = subset(ram, ram_generation %in% c("DDR4", "DDR5")),
                    FUN    = mean)

# Passa para formato largo: 1 linha por dia, 1 coluna por geracao
serie <- reshape(diario,
                 idvar     = "data",
                 timevar   = "ram_generation",
                 direction = "wide")
names(serie) <- c("data", "ddr4", "ddr5")
serie <- serie[order(serie$data), ]

# Mantem apenas dias com observacao nas DUAS series
serie <- na.omit(serie)
nrow(serie)   # quantos dias sobraram? VAR precisa de um T razoavel
head(serie)

plot(serie$data, serie$ddr4, type = "l", col = "blue",
     ylim = range(c(serie$ddr4, serie$ddr5)),
     xlab = "Data", ylab = "Preco medio por GB (USD)",
     main = "Preco medio diario por GB: DDR4 vs DDR5")
lines(serie$data, serie$ddr5, col = "red")
legend("topleft", legend = c("DDR4", "DDR5"),
       col = c("blue", "red"), lty = 1)

# --- 5. Estacionariedade -------------------------------------
# VAR padrao assume series estacionarias. Teste ADF:
# H0 = raiz unitaria (nao estacionaria)
adf.test(serie$ddr4)
adf.test(serie$ddr5)

# Se nao estacionarias, usar log-diferenca (~ variacao percentual)
y <- data.frame(d_ddr4 = diff(log(serie$ddr4)),
                d_ddr5 = diff(log(serie$ddr5)))
adf.test(y$d_ddr4)
adf.test(y$d_ddr5)

# --- 6. VAR --------------------------------------------------
# Escolha do numero de defasagens por criterio de informacao
VARselect(y, lag.max = 4, type = "const")

# Estima o VAR (ajuste p conforme o VARselect acima)
modelo <- VAR(y, p = 1, type = "const")
summary(modelo)

# Diagnostico: autocorrelacao serial dos residuos
serial.test(modelo, lags.pt = 8, type = "PT.asymptotic")

# Causalidade de Granger
causality(modelo, cause = "d_ddr4")
causality(modelo, cause = "d_ddr5")

# Funcao impulso-resposta (IRF)
plot(irf(modelo, n.ahead = 8, boot = TRUE))

# Decomposicao da variancia do erro de previsao (FEVD)
plot(fevd(modelo, n.ahead = 8))

# =============================================================
# PARA BASES FUTURAS:
#   nova <- read.csv("outra_base.csv", stringsAsFactors = FALSE,
#                    na.strings = c("", "NA"))
#   - repetir secoes 2-3 (tipos + inspecao)
#   - se ja for serie temporal, pular a agregacao (secao 4)
#     e ir direto para estacionariedade + VAR (secoes 5-6)
# =============================================================
