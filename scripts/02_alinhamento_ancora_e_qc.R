# PIPELINE MHC-DRB (Amazona aestiva) - SCRIPT 02 - Autor: Lemos
# Alinhamento Ancorado Exclusivo: Forpus passerinus (270 pb) + Amostras de Amazona aestiva
# Gerador de Relatório de Qualidade e Integridade do Alinhamento

# 1.Carregar Pacotes
library(Biostrings)
library(DECIPHER)
library(tidyverse)

# Estrutura de Diretórios
dir_processed <- "data/processed"
dir_results   <- "results"

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

# 2.Carregar Amostras(Script 01)
path_amostras <- file.path(dir_processed, "01_sequencias_orientadas.fasta")
if (!file.exists(path_amostras)) {
  stop("ERRO: Arquivo '01_sequencias_orientadas.fasta' não encontrado em: ", path_amostras)
}

seqs_amostras <- readDNAStringSet(path_amostras)
message("-> Amostras de Amazona aestiva carregadas: ", length(seqs_amostras))

# 3.Definir a Sequência Âncora (Forpus passerinus EF710734.1)
str_forpus <- "GGATATTCCAGGAGGCATTTAAGGCCGAGTGTCACTTCCTCAATGGCACCGAGCACGTCAGGTATCTGGAGAAGTATATCTACAACCGGGAGCAGTACGTGCACTTCGACAGCGATGTGGGGCACTATGTGGGTGATAACCCCGTGGGGGAGTTCCAAGCCAAGATCTGGAACAGCAAGGAGGAAATACTGGAGAATGCACGGGCTGAGGTGGACAGGTACTGCCGGCACAACTATGGGGTGTCGACCCCTTTCCTCGTAGATCGGAGGG"
seq_forpus <- DNAStringSet(setNames(str_forpus, "EF710734.1_Forpus_passerinus_Exon2_Anchor"))

# 4.Unificação: Âncora+ Amostras
banco_alinhamento <- c(seq_forpus, seqs_amostras)
message("-> Banco para alinhamento montado: 1 âncora (Forpus) + ", length(seqs_amostras), " amostras de A. aestiva")

# 5.Executar Alinhamento Múltiplo (DECIPHER)
message("\n-> Executando Alinhamento Múltiplo Ancorado...")
alinhamento_bruto <- AlignSeqs(banco_alinhamento, verbose = FALSE)

# 6.Trimming Rígido Ancorado
matriz_dna <- as.matrix(alinhamento_bruto)
idx_forpus <- grep("EF710734", rownames(matriz_dna))[1]

# Filtra colunas onde Forpus possui nucleotídeo real (sem gaps)
cols_validas_forpus <- which(matriz_dna[idx_forpus, ] != "-")
matriz_ancorada     <- matriz_dna[, cols_validas_forpus, drop = FALSE]

# Converte matriz filtrada de volta para DNAStringSet
seqs_alinhadas <- DNAStringSet(apply(matriz_ancorada, 1, paste, collapse = ""))

# 7.Auditoria de Qualidade e Integridade da Matriz
message("\n-> Calculando métricas de qualidade e integridade do alinhamento...")

vetor_forpus_char <- strsplit(str_forpus, "")[[1]]

calcular_metricas_qualidade <- function(nome_seq, seq_str) {
  seq_vetor <- strsplit(seq_str, "")[[1]]
  
  cat_seq  <- if_else(str_detect(nome_seq, "EF710734"), "Gabarito_Forpus", "Amostra_Amazona")
  qtd_gaps <- sum(seq_vetor == "-")
  tamanho_efetivo_pb <- length(seq_vetor) - qtd_gaps
  
  qtd_iupac <- sum(grepl("[RYSWKM]", seq_vetor, ignore.case = TRUE))
  pct_iupac <- round((qtd_iupac / max(1, tamanho_efetivo_pb)) * 100, 2)
  
  pos_comparaveis <- which(seq_vetor != "-" & vetor_forpus_char != "-")
  maches <- sum(seq_vetor[pos_comparaveis] == vetor_forpus_char[pos_comparaveis])
  pct_identidade <- round((maches / max(1, length(pos_comparaveis))) * 100, 2)
  
  #Verificação de Stop Codons
  seq_limpa <- gsub("-", "", seq_str)
  resto <- nchar(seq_limpa) %% 3
  if (resto != 0) seq_limpa <- paste0(seq_limpa, paste(rep("N", 3 - resto), collapse = ""))
  
  has_stop_codon <- FALSE
  if (nchar(seq_limpa) >= 3) {
    prot_str <- as.character(suppressWarnings(translate(DNAString(seq_limpa), if.fuzzy.codon = "solve")))
    prot_interna <- substr(prot_str, 1, max(1, nchar(prot_str) - 1))
    has_stop_codon <- str_detect(prot_interna, "\\*")
  }
  
  tibble(
    identificador = nome_seq,
    categoria = cat_seq,
    comprimento_matriz_pb = nchar(seq_str),
    bases_cobertas_pb = tamanho_efetivo_pb,
    quantidade_gaps = qtd_gaps,
    quantidade_iupac = qtd_iupac,
    porcentagem_iupac = pct_iupac,
    identidade_forpus_pct = pct_identidade,
    stop_codon_prematuro = has_stop_codon
  )
}

df_qualidade <- map2_dfr(names(seqs_alinhadas), as.character(seqs_alinhadas), calcular_metricas_qualidade)

# 8.Resumo
message(" RESUMO DA QUALIDADE DO ALINHAMENTO (FORPUS + AMAZE)")
resumo_console <- df_qualidade %>%
  group_by(categoria) %>%
  summarise(
    total_sequencias = n(),
    media_cobertura_pb = round(mean(bases_cobertas_pb), 1),
    media_identidade_forpus = round(mean(identidade_forpus_pct), 1),
    media_iupac_pct = round(mean(porcentagem_iupac), 2),
    com_stop_codon = sum(stop_codon_prematuro),
    .groups = "drop"
  )
print(resumo_console)

# 9.Salvar Arquivos de Saída
path_fasta_out <- file.path(dir_processed, "02_alinhamento_amostras_e_refs.fasta")
path_csv_out   <- file.path(dir_results, "relatorio_qualidade_alinhamento.csv")

writeXStringSet(seqs_alinhadas, filepath = path_fasta_out)
write_csv(df_qualidade, path_csv_out)