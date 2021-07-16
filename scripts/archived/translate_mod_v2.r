
getCodons <- function(myAln) {
    seqs <- as.character(myAln)
    len <- width(myAln)[1]
    starts <- seq(from=1, to=len, by=3)
    ends <- starts + 2
    myViews <- lapply(myAln, function(x) { 
        Views(x, starts, ends)
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
    return(myAAaln)
}
