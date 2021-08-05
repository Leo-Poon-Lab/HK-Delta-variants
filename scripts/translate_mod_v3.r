library(Biostrings)

getCodons <- function(myAln_new) {
    seqs <- as.character(myAln_new)
    len <- width(myAln_new)
    starts <- lapply(len, function(x){
        seq(from=1, to=x, by=3)
    })
    ends <- lapply(starts, function(x){
        x+2
    })   
    myViews <- lapply(seq_along(myAln_new), function(i) { 
        x <- myAln_new[[i]]
        starts_t <- starts[[i]]
        ends_t <- ends[[i]]
        Views(x, starts_t, ends_t)
    })
    myCodons <- lapply(myViews, function(x) {
        as.character(DNAStringSet(x))
    })
    myCodons
}
translateCodons <- function(myCodons, unknownCodonTranslatesTo="X") {
    ## make new genetic code
    gapCodon <- "-"
    # names(gapCodon) <- "---"
    names(gapCodon) <- "NNN"
    my_GENETIC_CODE <- c(GENETIC_CODE, gapCodon)
    
    ## translate the codons
    pep <- my_GENETIC_CODE[myCodons]
    
    ## check for codons that were not possible to translate, e.g. frameshift codons
    if (sum(is.na(pep))>0) {
        # cat("\nwarning - there were codons I could not translate. Using this character", unknownCodonTranslatesTo, "\n\n")
        pep[ which(is.na(pep)) ] <- unknownCodonTranslatesTo
    }
    
    ## prep for output
    pep <- paste(pep, collapse="")
    return(pep)
}

## wrap the getCodons and translateCodons functions together into one:
translateGappedAln <- function(myAln, unknownCodonTranslatesTo="X") {
    myAln_new <- myAln
    sapply(seq_along(myAln), function(i){
        # print(i)
        seq_tmp <- myAln[i]
        gap_codon <- "NNN"
        seq_tmp_new <- gsub(gap_codon, "", as.character(seq_tmp), fixed = T)
        myAln_new[i] <<- DNAStringSet(seq_tmp_new)
    })
    myCodons <- getCodons(myAln_new)
    myAAaln <- AAStringSet(unlist(lapply(myCodons, translateCodons, unknownCodonTranslatesTo=unknownCodonTranslatesTo)))
    names(myAAaln) <- names(myAln)

    file_in <- "./tmp.fasta"
    file_out <- "./tmp_aln.fasta"
    file_aa_ref <- "./tmp_ref_aa.fasta"
    writeXStringSet(myAAaln, file_in)
    system(paste0("mafft --6merpair --thread -10 --keeplength --addfragments ", file_in, " ", file_aa_ref, " > ", file_out)) #MAFFT
    AAaln <- readAAStringSet(file_out)
    file.remove(file_in)
    file.remove(file_out)
    file.remove(file_aa_ref)
    return(AAaln[-1])
}

get_nt <- function(seq_nt, gene = "spike"){
    stopifnot(all(width(seq_nt) == 29903))
    orf <- read.csv("./ORF_SCoV2.csv")
    idx <- which(tolower(orf$sequence) == tolower(gene))
    stopifnot(length(idx)>0)
    subseq(seq_nt, orf$start[idx], orf$stop[idx])
}

get_aa <- function(seq_nt, gene = "spike", unknownCodonTranslatesTo="X"){
    myAln <- get_nt(seq_nt, gene = gene)
    myAln <- chartr("-", "N", myAln)
    
    ref_seq <- readDNAStringSet("./reference.fasta")
    ref_seq_nt <- get_nt(ref_seq, gene)
    ref_seq_aa <- translate(ref_seq_nt)
    file_aa_ref <- "./tmp_ref_aa.fasta"
    writeXStringSet(ref_seq_aa, file_aa_ref)

    translateGappedAln(myAln, unknownCodonTranslatesTo = unknownCodonTranslatesTo)
}
