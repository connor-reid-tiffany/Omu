#Method for gathering metadata from KEGG API


KEGG_Gather <- function(x) UseMethod("KEGG_Gather", x)

KEGG_Gather.cpd <- function(x) "cpd"

KEGG_Gather.rxn <- function(x) "rxn"

KEGG_Gather.genes <- function(x) "KO"
