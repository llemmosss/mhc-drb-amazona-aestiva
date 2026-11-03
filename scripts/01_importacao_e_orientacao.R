# PIPELINE MHC-DRB (Amazona aestiva) - SCRIPT 01 
# Autor original: Lemos 
# Descrição: Leitura, Padronização, Tag de Primer, Dicionário e Orientação (5' -> 3')

# 1.Pacotes

library(Biostrings)
library(tidyverse)

# 2.Definição de Estrutura de Diretórios Relativa
dir_raw       <- file.path("data", "raw")
dir_processed <- file.path("data", "processed")

if (!dir.exists(dir_processed)) {
  dir.create(dir_processed, recursive = TRUE)
}

# 3.Dicionário de Mapeamento (USP)
mapa_usp <- tibble::tribble(
  ~ref_code, ~animal_id,
  "REF-1",  "338-24",
  "REF-2",  "338-24",
  "REF-3",  "19-25",
  "REF-4",  "19-25",
  "REF-5",  "337-24",
  "REF-6",  "337-24",
  "REF-7",  "336-24",
  "REF-8",  "336-24",
  "REF-9",  "344-24",
  "REF-10", "344-24"
)


# 4.Processar Amostras USP

processar_usp <- function(caminho_arquivo, tabela_mapa) {
  if (is.na(caminho_arquivo) || !file.exists(caminho_arquivo)) {
    warning(" Arquivo USP não encontrado em: ", dir_raw)
    return(list(seqs = DNAStringSet(), metadata = tibble()))
  }
  
  fasta_raw <- readDNAStringSet(caminho_arquivo)
  headers_orig <- names(fasta_raw)
  
  seqs_list <- vector("list", length(fasta_raw))
  meta_list <- vector("list", length(fasta_raw))
  
  for (i in seq_along(fasta_raw)) {
    h_orig <- headers_orig[i]
    seq_dna <- fasta_raw[[i]]
    
    # Extrai padrão REF-X ou REF_X
    ref_match <- stringr::str_extract(h_orig, "REF[-_][0-9]+")
    
    if (!is.na(ref_match)) {
      ref_formatted <- gsub("_", "-", ref_match)
      num_ref <- as.numeric(gsub("REF-", "", ref_formatted))
      
      info <- tabela_mapa %>% filter(ref_code == ref_formatted)
      animal_id <- if (nrow(info) > 0) info$animal_id[1] else "UNKNOWN"
    } else {
      ref_formatted <- "REF-UNK"
      num_ref <- NA
      animal_id <- "UNKNOWN"
    }
    
    # Determina orientação por regra par/ímpar
    is_reverse <- !is.na(num_ref) && (num_ref %% 2 == 0)
    
    if (is_reverse) {
      seq_orientada <- reverseComplement(seq_dna)
      direcao <- "Reverse"
    } else {
      seq_orientada <- seq_dna
      direcao <- "Forward"
    }
    
    # Header limpo e conciso para o FASTA (evita truncamento em softwares)
    clean_id <- paste0("Amae|USP|", animal_id, "|", ref_formatted)
    
    seqs_list[[i]] <- seq_orientada
    
    # Tabela de Metadados detalhada
    meta_list[[i]] <- tibble(
      fasta_header_clean = clean_id,
      header_original    = h_orig,
      fonte              = "USP",
      animal_id          = animal_id,
      ref_code           = ref_formatted,
      direcao_inferida   = direcao,
      invertida_rc       = is_reverse,
      tamanho_pb         = length(seq_orientada)
    )
  }
  
  res_seqs <- DNAStringSet(seqs_list)
  names(res_seqs) <- map_chr(meta_list, "fasta_header_clean")
  
  return(list(seqs = res_seqs, metadata = bind_rows(meta_list)))
}

# 5.Processar Amostras FIOCRUZ

processar_fiocruz <- function(caminho_arquivo) {
  if (is.na(caminho_arquivo) || !file.exists(caminho_arquivo)) {
    warning(" Arquivo Fiocruz não encontrado em: ", dir_raw)
    return(list(seqs = DNAStringSet(), metadata = tibble()))
  }
  
  fasta_raw <- readDNAStringSet(caminho_arquivo)
  headers_orig <- names(fasta_raw)
  
  seqs_list <- vector("list", length(fasta_raw))
  meta_list <- vector("list", length(fasta_raw))
  
  for (i in seq_along(fasta_raw)) {
    h_orig <- headers_orig[i]
    h_limpo <- gsub(" \\(Amazona aestiva\\)", "", h_orig)
    seq_dna <- fasta_raw[[i]]
    
    # Identifica 'R' ou 'r' para Reverse
    is_reverse <- grepl("^[0-9]+[Rr]_", h_limpo)
    
    if (is_reverse) {
      seq_orientada <- reverseComplement(seq_dna)
      direcao <- "Reverse"
    } else {
      seq_orientada <- seq_dna
      direcao <- "Forward"
    }
    
    # Extração simples do ID do animal/amostra
    sample_id <- stringr::str_extract(h_limpo, "^[0-9]+[RrFf]?_[0-9]+") %>% 
      coalesce(h_limpo)
    
    clean_id <- paste0("Amae|FIOCRUZ|", sample_id)
    
    seqs_list[[i]] <- seq_orientada
    
    meta_list[[i]] <- tibble(
      fasta_header_clean = clean_id,
      header_original    = h_orig,
      fonte              = "FIOCRUZ",
      animal_id          = sample_id,
      ref_code           = "N/A",
      direcao_inferida   = direcao,
      invertida_rc       = is_reverse,
      tamanho_pb         = length(seq_orientada)
    )
  }
  
  res_seqs <- DNAStringSet(seqs_list)
  names(res_seqs) <- map_chr(meta_list, "fasta_header_clean")
  
  return(list(seqs = res_seqs, metadata = bind_rows(meta_list)))
}

# 6. Execução do Fluxo de Dados
# Localiza arquivos brutos
path_usp <- list.files(dir_raw, pattern = "raw_amazona_aestiva_usp", full.names = TRUE)[1]
path_fio <- list.files(dir_raw, pattern = "raw_amazona_aestiva_fiocruz", full.names = TRUE)[1]

# Processa ambos os conjuntos
res_usp     <- processar_usp(path_usp, mapa_usp)
res_fiocruz <- processar_fiocruz(path_fio)

# Une sequências e metadados
todas_seqs <- c(res_usp$seqs, res_fiocruz$seqs)
tabela_metadados <- bind_rows(res_usp$metadata, res_fiocruz$metadata)


# 7.Salvando os Resultados
out_fasta <- file.path(dir_processed, "01_sequencias_orientadas.fasta")
out_meta  <- file.path(dir_processed, "01_metadados_amostras.csv")

# Salva FASTA
writeXStringSet(todas_seqs, filepath = out_fasta)

# Salva Tabela de Metadados
write_csv(tabela_metadados, file = out_meta)