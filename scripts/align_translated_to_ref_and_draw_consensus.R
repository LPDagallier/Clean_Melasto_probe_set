##############################################################################
# Description: ALIGN TRANSLATED NUCLEOTIDES TO REFERENCE AND DRAW CONSENSUS

# Author: Leo-Paul Dagallier
# Date created: 2022-09-27
##############################################################################

# Packages -------------------------------------------------------------------
# Package dependencies: DECIPHER, Biostrings
# Inputs/parameters ----------------------------------------------------------
# Needs to be run with 5 arguments:
#   - args[1]: input fasta file containing the hits to align ($ID"_hits_to_align.FNA")
#   - args[2]: output filename for the aligned hits to the reference ($ID"_hits_aligned_to_ref.FNA")
#   - args[3]: text file containing the name of the reference sequence to remove before drawing the consensus (IDs_to_remove_after_alignment.txt)
#   - args[4]: character string containing the name of the consensus sequence ($ID"_hits_consensus")
#   - args[5]: output filename for the consensus sequence ($ID"_hits_consensus.FNA")
# 
# Example: args = c("Tibouchina-AT5G67530_hits_to_align.FNA", "Tibouchina-AT5G67530_hits_aligned_to_ref.FNA", "IDs_to_remove_after_alignment.txt", "Tibouchina-AT5G67530_hits_consensus", "Tibouchina-AT5G67530_hits_consensus.FNA")

# Main code ------------------------------------------------------------------
# load the arguments
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
  stop("Five arguments must be supplied (see script file for details).", call.=FALSE)
} else {
  print("Script called with the following arguments:")
  print(args)
}

to_align <- Biostrings::readDNAStringSet(args[1], format="fasta")
# Align:
aligned <- DECIPHER::AlignTranslation(to_align, readingFrame = 1, gapOpening = -50)
# Export alignment:
Biostrings::writeXStringSet(aligned, filepath = args[2])
# Remove reference sequence:
aligned_no_ref <- aligned[- grep(pattern = paste(readLines(args[3]), collapse = "|"), x = names(aligned))]
# Draw consensus:
consensus <- DECIPHER::ConsensusSequence(aligned_no_ref,
                                         threshold = 0.05,
                                         ambiguity = TRUE,
                                         noConsensusChar = "+",
                                         ignoreNonBases = T,
                                         includeTerminalGaps = FALSE)
names(consensus) <- args[4]
# Export consensus:
Biostrings::writeXStringSet(consensus, filepath = args[5])

# End of script
##############################################################################