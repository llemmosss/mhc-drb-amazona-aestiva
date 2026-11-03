# PIPELINE MHC-DRB (Amazona aestiva) - SCRIPT DE TESTE IN SILICO (BENCHMARK) - Autor Lemos

library(Biostrings)
library(DECIPHER)
library(tidyverse)

# 1.Sequência (Forpus passerinus - EF710734.1 / 270 pb)
str_forpus <- "GGATATTCCAGGAGGCATTTAAGGCCGAGTGTCACTTCCTCAATGGCACCGAGCACGTCAGGTATCTGGAGAAGTATATCTACAACCGGGAGCAGTACGTGCACTTCGACAGCGATGTGGGGCACTATGTGGGTGATAACCCCGTGGGGGAGTTCCAAGCCAAGATCTGGAACAGCAAGGAGGAAATACTGGAGAATGCACGGGCTGAGGTGGACAGGTACTGCCGGCACAACTATGGGGTGTCGACCCCTTTCCTCGTAGATCGGAGGG"

# 2.Construção In Silico dos Controles
seq_pos_pura       <- str_forpus                                       
seq_com_flancos    <- paste0("ATCGATCGATCGATCGATCGATCGATCGATC", str_forpus, "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA") 
seq_frameshift     <- paste0(substr(str_forpus, 1, 100), substr(str_forpus, 102, 270)) 
seq_stop_prematuro <- paste0(substr(str_forpus, 1, 50), "TAA", substr(str_forpus, 54, 270)) 

# Banco de Benchmark
benchmark_set <- DNAStringSet(c(
  "CTRL_POSITIVO_PURO"      = seq_pos_pura,
  "TESTE_TRIMMING_FLANCOS"  = seq_com_flancos,
  "TESTE_MUT_FRAMESHIFT"    = seq_frameshift,
  "TESTE_STOP_CODON_PSEUDO" = seq_stop_prematuro
))

# 3.Processamento no Pipeline (Alinhamento + Trimming)
ancora <- DNAStringSet(setNames(str_forpus, "EF710734.1_Anchor"))
banco_alinhamento <- c(ancora, benchmark_set)

# Alinhamento e Trimming
aln <- AlignSeqs(banco_alinhamento, verbose = FALSE)
matriz <- as.matrix(aln)
idx_ancora <- grep("EF710734.1_Anchor", rownames(matriz))[1]
cols_validas <- which(matriz[idx_ancora, ] != "-")
matriz_recortada <- matriz[, cols_validas, drop = FALSE]

# Matriz Recortada Final
seqs_processadas <- DNAStringSet(apply(matriz_recortada, 1, paste, collapse = ""))

# 4.Avaliação de Desempenho
seqs_teste <- seqs_processadas[!str_detect(names(seqs_processadas), "Anchor")]
banco_teste <- banco_alinhamento[!str_detect(names(banco_alinhamento), "Anchor")]

relatorio_stress <- tibble(
  ID = names(seqs_teste),
  Tamanho_Original_pb = width(banco_teste),
  Tamanho_Recortado_pb = width(seqs_teste)
) %>%
  mutate(
    Trim_Exato_270pb = (Tamanho_Recortado_pb == 270),
    Tem_Stop_Codon = map_lgl(as.character(seqs_teste), function(s) {
      s_clean <- gsub("-", "", s)
      resto <- nchar(s_clean) %% 3
      if (resto != 0) s_clean <- paste0(s_clean, paste(rep("N", 3 - resto), collapse = ""))
      prot <- as.character(suppressWarnings(translate(DNAString(s_clean), if.fuzzy.codon = "solve")))
      prot_interna <- substr(prot, 1, max(1, nchar(prot) - 1))
      str_detect(prot_interna, "\\*")
    })
  )

# 5.Definição do Diretório de Saída e Exportação
dir_saida <- "C:/Users/llemm/Downloads/Pasta_Compartilhada/LEMOS_TCC/results"

if (!dir.exists(dir_saida)) {
  dir.create(dir_saida, recursive = TRUE)
}

caminho_arquivo <- file.path(dir_saida, "relatorio_teste_estresse_mhc.csv")
write_csv(relatorio_stress, caminho_arquivo)