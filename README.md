# How to embed your data
## All code in this tutorial is meant as a guide/example. It's not a complete wrapper, and will likely not be plug-and-play for you and your system. A complete package is in the works, however, it is not available at this time.

## 1. Begin with a sample by ASV table as processed by Dada2, and a fasta file of the representative sequences.
Use this tutorial: https://benjjneb.github.io/dada2/tutorial.html to process your raw fastq reads. At the end of the tutorial, you are left with a variable called seqtab.nochim which is a sample by ASV table. Write this table to a file. Also write a fasta file with each full length sequence. See the end of scripts/embed_your_data/dada2.R or below:
To write a fasta file containing the sequences: 
 ```
 Rscript
    #################################################
    #### write sequence table #######################
    #################################################
  ```
    write.table(seqtab, paste(data_dir, "seqtab.txt", sep = ""),
                quote = F, sep = "\t", row.names = T, col.names = T)

    ##################################################
    ###  Write representative sequences  #############
    ##################################################

    fasta_file_name <- paste(data_dir, "queryseqs.fasta", sep = "")
    taxa_seqs <- colnames(seqtab)
    headers <- paste(">seq", seq(1, length(taxa_seqs)), sep = "")
    fasta <- paste(headers, taxa_seqs, sep = "\n", collapse = "\n")
    write(fasta, fasta_file_name)
  ```
 
  ```
## 2. Blastn against the embedding database
  ```
    #!/bin/bash
    query_path="" #CHANGE ME to the directory containing the fasta file of your query sequences, created in step 1.
    out_path="" #CHANGE ME to the directory you want to output your blast results to
    blast_path="" #CHANGE ME to the directory holding your downloaded blastn software
    db_path="" #CHANGE ME to the directory holding 
    $blast_dir/blastn -db data/blastdb/embedding_db -query "$query_path"/queryseqs.fasta -out "$out_path"/blast_hits.tsv  -outfmt "6 qseqid sseqid qseq sseq evalue bitscore"
    cat "$out_path"/blast_hits.tsv | sort -k1,1 -k5,5g -k6,6nr | sort -u -k1,1 --merge > "$out_path"/best_hits.tsv
  ```
## 3. Assign the appropriate embedding sequence id to each of your sequences. If there is no match, toss out that sequence. See scripts/embed_your_data/embed_asv_table.R 
```
   data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/"

   #########################
   ### Read ASV table ######
   #########################
   asv_file = paste(data_dir, "/halfvarson/seqtab.rds", sep = "")
   seqtab <- readRDS(asv_file)
   #seqtab <- read.table(asv_file, quote = F, sep = "\t", row.names = T, col.names = T)

   #### change column names to id
   query_seqs_file = paste(data_dir, "/halfvarson/queryseqs.fasta", sep ="")
   query_fasta <- read.table(query_seqs_file)
   query_headers <- query_fasta[seq(1, nrow(query_fasta), by = 2), 1]
   query_headers <- gsub(">", "", query_headers)
   query_seqs <- as.character(query_fasta[seq(2, nrow(query_fasta), by = 2), 1])

   query_headers <- query_headers[query_seqs %in% colnames(seqtab)]
   query_seqs <- query_seqs[query_seqs %in% colnames(seqtab)]

   colnames(seqtab)[match(colnames(seqtab), query_seqs)] == query_seqs #check
   colnames(seqtab)[match(colnames(seqtab), query_seqs)] <- query_headers

   ##############################
   ## Read quality vector  ######
   ## transformation table ######
   ##############################
   transform_mat_file <- paste(data_dir, "/embed/embed_.07_100dim.txt", sep = "") #CHANGE TO EMBED IN A DIFFERENT DIMENSIONAL SPACE
   qual_vecs <- read.table(transform_mat_file)

   #############################
   ## Read best hits table #####
   ## from blast ###############
   #############################

   blast_output_file = paste(data_dir, "halfvarson/embed/best_hits.tsv", sep = "")
   best_hits <- read.delim(blast_output_file, header=FALSE, row.names = 1)
   colnames(best_hits) <- c("hit_id", "query_seq", "hit_seq", "evalue", "bitscore")

   #Filter best_hits table to only include hits that pass the e-value threshold
   best_hits <- best_hits[best_hits$evalue < 1*10^(-29), ]

   #################################
   ### Assign nearest neighbor id ##
   #################################

   #Drop any ASVs from the table that don't have near enough hits in the transformation matrix
   seqtab <- seqtab[ , colnames(seqtab) %in% rownames(best_hits)] #17784 taxa left

   #Assign the id of each ASV's nearest hit in the embedding transformation table.
   colnames(seqtab) <- best_hits[colnames(seqtab), "hit_id"]


   ###############################
   ### Match orders ##############
   ###############################
   qual_vecs <- qual_vecs[colnames(seqtab), ]


   ###############################
   ### Take dot product ##########
   ###############################

   embedded <- as.matrix(seqtab) %*% as.matrix(qual_vecs)
   embedded_file = paste(data_dir, "halfvarson/embed/seqtab_embedded_.07_100dim", sep = "")
   saveRDS(embedded, paste(embedded_file, ".rds", sep = ""))
   write.table(embedded, paste(embedded_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)

```

You have successfully embedded your data! Each row represents that sample's microbiome signature.
