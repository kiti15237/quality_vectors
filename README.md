# How to embed your data
## All code in this tutorial is meant as a guide/example. It's not a complete wrapper, and will likely not be plug-and-play for you and your system. A complete R package is in the works, however, is not available at this time.

## 1. Begin with a sample by ASV table as processed by Dada2, and a fasta file of the representative sequences.
Use this tutorial: https://benjjneb.github.io/dada2/tutorial.html to process your raw fastq reads. At the end of the tutorial, you are left with a variable called seqtab.nochim which is a sample by ASV table.
To write a fasta file containing the sequences: 
 ```
 Rscript
    fasta_file_name <- "queryseqs.fasta" 
    taxa_seqs <- colnames(seqtab.nochim)
    headers <- paste(">seq", seq(1, length(taxa_seqs)), sep = "")
    fasta <- paste(headers, taxa_seqs, sep = "\n", collapse = "\n")
    write(fasta, fasta_file_name)
  ```
 
To make future processing clearer, we should also substitute the full length sequence for the sequence id in the column names of our ASV table
  ```
    colnames(seqtab.nochim) <- headers
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
## 3. Assign the appropriate embedding sequence id to each of your sequences. If there is no match, toss out that sequence. 
  ```
  Rscript
    blast_output_dir = "" #CHANGE ME
    best_hits <- read.delim(paste(blast_output_dir, "best_hits.tsv", sep = ""), header=FALSE, row.names = 1)
    colnames(best_hits) <- c("hit_id", "query_seq", "hit_seq", "evalue", "bitscore")

    #Filter best_hits table to only include hits that pass the e-value threshold
    best_hits <- best_hits[best_hits$evalue < 1*10^(-29), ]

    #Drop any ASVs from the table that don't have near enough hits in the transformation matrix
    seqtab.nochim <- seqtab.nochim[ , colnames(seqtab.nochim) %in% rownames(best_hits)]

    #Assign the id of each ASV's nearest hit in the embedding transformation table.
    colnames(seqtab.nochim) <- best_hits[colnames(seqtab.nochim), "hit_seq"]
  ```
## 4. Reorder your ASV table columns to match the embedding transformation matrix rows. 
  ```
    #Read in transformation matrix
    qual_vecs <- read.table("data/embed/embed_.07_100.txt", row.names=1, quote="\"", comment.char="")

    #Reorder the transformation matrix to match order of query table
    qual_vecs <- qual_vecs[colnames(seqtabe.nochim), :]
  ```
## 5. Take the dot product between the ASV table and the transformation matrix. This is your embedded data. Use this for machine learning.
  ```
    #Normalize the taxa counts
    seqtab.norm = asinh(as.matrix(seqtab.nochim))
    
    #Transform by dot product
    seqtab_embedded <- as.data.frame(seqtab.norm) %*% as.matrix(qual_vecs))
    rownames(seqtab_embedded) <- rownames(seqtab.nochim)
    colnames(seqtab_embedded) <- colnames(qual_vecs)
  ```
You have successfully embedded your data! Each row represents that sample's microbiome signature.
