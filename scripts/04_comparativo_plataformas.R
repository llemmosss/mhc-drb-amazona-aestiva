# PIPELINE MHC-DRB (Amazona aestiva) - SCRIPT Comparativo entre plataformas - Autor Lemos

library(Biostrings)
library(tidyverse)

# 1.Carregar Arquivos Brutos
path_usp <- list.files("data/raw", pattern = "raw_amazona_aestiva_usp", full.names = TRUE)[1]
path_fiocruz <- list.files("data/raw", pattern = "raw_amazona_aestiva_fiocruz", full.names = TRUE)[1]

fasta_usp <- readDNAStringSet(path_usp)
fasta_fiocruz <- readDNAStringSet(path_fiocruz)

# 2.Métricas
calcular_metricas <- function(fasta_obj, nome_plataforma) {
  seqs <- as.character(fasta_obj)
  tibble(
    plataforma = nome_plataforma,
    cabecalho = names(fasta_obj),
    tamanho_pb = nchar(seqs),
    qtd_iupac = str_count(seqs, "[RYSWKM]"),
    porcentagem_iupac = round((qtd_iupac / tamanho_pb) * 100, 2)
  )
}

df_comparativo <- bind_rows(
  calcular_metricas(fasta_usp, "USP"),
  calcular_metricas(fasta_fiocruz, "Fiocruz")
)

# 3.Resumo Média
df_resumo <- df_comparativo %>%
  group_by(plataforma) %>%
  summarise(media_iupac = mean(porcentagem_iupac), .groups = "drop")

# 4.Gráfico (Média + Pontos Individuais + Rótulo)
grafico_iupac <- ggplot() +
  geom_col(data = df_resumo, aes(x = plataforma, y = media_iupac, fill = plataforma), 
           width = 0.45, alpha = 0.75, color = "black", linewidth = 0.3) +
  geom_jitter(data = df_comparativo, aes(x = plataforma, y = porcentagem_iupac, color = plataforma),
              width = 0.12, size = 3, alpha = 0.85) +
  geom_text(data = df_resumo, aes(x = plataforma, y = media_iupac + 0.7, label = sprintf("%.1f%%", media_iupac)),
            fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("USP" = "#2b5c8f", "Fiocruz" = "#e66101")) +
  scale_color_manual(values = c("USP" = "#112438", "Fiocruz" = "#803300")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Taxa Média de Ambiguidades (IUPAC) por Plataforma",
    subtitle = "Comparação da frequência de picos duplos/heterozigose (USP vs Fiocruz)",
    x = "Plataforma de Sequenciamento",
    y = "Média de Sítios Ambíguos (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold")
  )

# 5. Salvar
write_csv(df_comparativo, "results/comparativo_qualidade_plataformas.csv")
ggsave("results/grafico_comparativo_iupac.png", plot = grafico_iupac, width = 6.5, height = 4.5, dpi = 300)