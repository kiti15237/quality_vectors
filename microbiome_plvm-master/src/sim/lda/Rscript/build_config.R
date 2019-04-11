library("jsonlite")
df <- expand.grid(
  D = c(20, 100),
  V = c(325, 650),
  K = 2,
  N = c(1625, 3250, 6500),
  alpha0 = 1,
  gamma0 = 1
)

MICROBIOME_PLVM_DIR = "C:/Users/ctata/Documents/Lab/quality_vectors/microbiome_plvm-master"
df$id <- seq_len(nrow(df))
lda_dir <- file.path(
  MICROBIOME_PLVM_DIR,
  "src",
  "sim",
  "lda"
)

cat(
  toJSON(df, auto_unbox = TRUE),
  file = file.path(lda_dir, "pipeline", "experiment.json")
)
