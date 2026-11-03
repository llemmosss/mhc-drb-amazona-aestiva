# PIPELINE MHC-DRB (Amazona aestiva) - SCRIPT 03 - Autor: Lemos
# Alinhamento Múltiplo Multiespécie (Amazona aestiva + Forpus + Painel Externo de Aves)

# 1.Pacotes

library(Biostrings)
library(DECIPHER)
library(tidyverse)

# Estrutura de Diretórios
dir_processed <- "data/processed"
dir_ref       <- "data/references"
dir_results   <- "results"

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

# 2. Carregar Matriz do Script 02 
path_amae_forpus <- file.path(dir_processed, "02_alinhamento_amostras_e_refs.fasta")
if (!file.exists(path_amae_forpus)) {
  stop("ERRO: Arquivo '02_alinhamento_amostras_e_refs.fasta' não encontrado em data/processed/")
}

seqs_amae_forpus <- readDNAStringSet(path_amae_forpus)
message("-> Sequências do Script 02 (A. aestiva + Forpus) carregadas: ", length(seqs_amae_forpus))

# 3.Carregar Painel Extensível de Aves
path_refs_aves <- file.path(dir_ref, "ref_todas_aves.fasta")
if (!file.exists(path_refs_aves)) {
  stop("ERRO: Arquivo 'ref_todas_aves.fasta' não encontrado em: ", path_refs_aves)
}

seqs_outras_aves <- readDNAStringSet(path_refs_aves)

# Remove duplicata do Forpus no arquivo externo, se houver
idx_dup_forpus <- grep("EF710734", names(seqs_outras_aves), ignore.case = TRUE)
if (length(idx_dup_forpus) > 0) {
  seqs_outras_aves <- seqs_outras_aves[-idx_dup_forpus]
}

message("-> Sequências do painel externo de aves carregadas: ", length(seqs_outras_aves))

# 4.Unificação do Banco Multiespécie e Remoção de Gaps Prévios
banco_multiespecie_raw <- c(seqs_amae_forpus, seqs_outras_aves)

# Remove hífens (-) para que o AlignSeqs possa re-alinhar sem erros
banco_multiespecie <- RemoveGaps(banco_multiespecie_raw)

message("-> Total de sequências unificadas (sem gaps) para o alinhamento: ", length(banco_multiespecie))

# 5.Execução do Alinhamento Múltiplo Multiespécie (DECIPHER)
message("\n-> Executando Alinhamento Múltiplo de Perfis (Multiespécie)...")
alinhamento_bruto <- AlignSeqs(banco_multiespecie, verbose = FALSE)

# 6.Trimming Rígido Ancorado
matriz_dna <- as.matrix(alinhamento_bruto)
idx_forpus <- grep("EF710734", rownames(matriz_dna))[1]

# Filtra exclusivamente as colunas onde Forpus possui nucleotídeo real
cols_validas_forpus <- which(matriz_dna[idx_forpus, ] != "-")
matriz_ancorada     <- matriz_dna[, cols_validas_forpus, drop = FALSE]

# Converte a matriz recortada de volta para DNAStringSet
seqs_multiespecie_limpas <- DNAStringSet(apply(matriz_ancorada, 1, paste, collapse = ""))

# 7.Auditoria da Matriz Multiespécie e Cálculo de Similaridade
message("\n-> Gerando Relatório da Matriz Comparativa Multiespécie...")

str_forpus <- as.character(seqs_multiespecie_limpas[idx_forpus])
vetor_forpus_char <- strsplit(str_forpus, "")[[1]]

calcular_metricas_multiespecie <- function(nome_seq, seq_str) {
  seq_vetor <- strsplit(seq_str, "")[[1]]
  
  cat_seq <- if_else(str_detect(nome_seq, "EF710734"), "Gabarito_Forpus",
                     if_else(str_detect(nome_seq, "Amae_|REF-"), "Amazona_aestiva", "Referencia_Externa_Ave"))
  
  qtd_gaps <- sum(seq_vetor == "-")
  tamanho_efetivo_pb <- length(seq_vetor) - qtd_gaps
  
  # Percentual de Identidade contra o Forpus
  pos_comparaveis <- which(seq_vetor != "-" & vetor_forpus_char != "-")
  maches <- sum(seq_vetor[pos_comparaveis] == vetor_forpus_char[pos_comparaveis])
  pct_identidade <- round((maches / max(1, length(pos_comparaveis))) * 100, 2)
  
  tibble(
    identificador = nome_seq,
    grupo_taxonomico = cat_seq,
    comprimento_matriz_pb = nchar(seq_str),
    bases_cobertas_pb = tamanho_efetivo_pb,
    gaps_na_janela = qtd_gaps,
    identidade_forpus_pct = pct_identidade
  )
}

df_multiespecie <- map2_dfr(names(seqs_multiespecie_limpas), as.character(seqs_multiespecie_limpas), calcular_metricas_multiespecie)

# 8.Exibição do Resumo no Console
message(" RESUMO DO ALINHAMENTO MULTIESPÉCIE (ÉXON 2)")
resumo_console <- df_multiespecie %>%
  group_by(grupo_taxonomico) %>%
  summarise(
    total_sequencias = n(),
    media_cobertura_pb = round(mean(bases_cobertas_pb), 1),
    media_identidade_forpus = round(mean(identidade_forpus_pct), 1),
    .groups = "drop"
  )
print(resumo_console)

# 9.Salvar Arquivos de Saída
path_fasta_out <- file.path(dir_processed, "03_alinhamento_multiespecie_aves.fasta")
path_csv_out   <- file.path(dir_results, "relatorio_alinhamento_multiespecie.csv")

writeXStringSet(seqs_multiespecie_limpas, filepath = path_fasta_out)
write_csv(df_multiespecie, path_csv_out)