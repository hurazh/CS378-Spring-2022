library(Biostrings)

# DNAString object
dna_str <- DNAString("AAACGCG")
dna_str
length(dna_str)
reverseComplement(dna_str)

# DNAStringSet object
dna_set <- DNAStringSet(c("AAAAA","CGCG","TCG"))
dna_set
reverseComplement(dna_set)
length(dna_set)
width(dna_set)

# occurence/frequency
letterFrequency(dna_str, 'C')
letterFrequency(dna_set, 'C')
letterFrequency(dna_set, "C", as.prob=TRUE)
letterFrequency(dna_str, "CG")
letterFrequency(dna_set, "CG")
# matching
# single pattern + multiple strings
dna <- DNAStringSet(c("AACTCTC","CTCTAAA","AAAGAG"))
matches <- vmatchPattern("CTC", dna)
matches
elementNROWS(matches)
# multiple pattern + single string
dna <- DNAString("AAACTCAAAGAGAAATTTAAA")
pd <- PDict(c("CTC","GAG","TTT","AAA"))
matches <- matchPDict(pd, dna)
matches
elementNROWS(matches)

# Full genome sequences for Homo sapiens
# BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
library(BSgenome.Hsapiens.NCBI.GRCh38)
Hsapiens
# seqence names
seqnames(Hsapiens)
seqinfo(Hsapiens)
chr1 <- Hsapiens[[1]]
length(chr1)
alphabetFrequency(chr1)
matchPattern("CTC", chr1)
n_pattern = "AAGCCTAAGCCTAAGCCTAA"
matchPattern(seq_pattern, chr1)
matches <- matchPattern(n_pattern, chr1, max.mismatch=2)
matches
mismatch(n_pattern, matches)
gr <- GRanges("1", IRanges(1e6 + c(1,101,201), width=100))
getSeq(Hsapiens, gr)

