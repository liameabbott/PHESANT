
library(data.table)
source("summarise_phenotypes.r")

outcome_info <- read.table("./variable-info/outcome_info_final.tsv", sep='\t', quote="", comment.char="", header=TRUE)

for (idx in seq(1, 4)) {

    tsv_data_both <- read.table(paste0('../ukb31063.phesant.both_sexes.', idx, '.tsv'), header=TRUE, sep='\t')
    tsv_data_females <- read.table(paste0('../ukb31063.phesant.female.', idx, '.tsv'), header=TRUE, sep='\t')
    tsv_data_males <- read.table(paste0('../ukb31063.phesant.male.', idx, '.tsv'), header=TRUE, sep='\t')

    log_file_both <- paste0('../ukb31063.phesant.both_sexes.', idx, '.log')
    log_file_females <- paste0('../ukb31063.phesant.female.', idx, '.log')
    log_file_males <- paste0('../ukb31063.phesant.male.', idx, '.log')

    notes_both <- get_notes(tsv_data_both, log_file_both, outcome_info, codings_tables)
    write.table(notes_both, file=paste0('../ukb31063.both_sexes.tmp.', idx, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_females <- get_notes(tsv_data_females, log_file_females, outcome_info, codings_tables)
    write.table(notes_females, file=paste0('../ukb31063.female.tmp.', idx, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_males <- get_notes(tsv_data_males, log_file_males, outcome_info, codings_tables)
    write.table(notes_males, file=paste0('../ukb31063.male.tmp.', idx, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_codings_included_both <- include_PHESANT_reassignment_names(paste0('../ukb31063.both_sexes.tmp.', idx, '.tsv'), outcome_info)
    write.table(notes_codings_included_both, file=paste0('../ukb31063.both_sexes.phenosummary.', idx, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_codings_included_females <- include_PHESANT_reassignment_names(paste0('../ukb31063.female.tmp.', idx, '.tsv'), outcome_info)
    write.table(notes_codings_included_both, file=paste0('../ukb31063.female.phenosummary.', idx, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_codings_included_males <- include_PHESANT_reassignment_names(paste0('../ukb31063.male.tmp.', idx, '.tsv'), outcome_info)
    write.table(notes_codings_included_both, file=paste0('../ukb31063.male.phenosummary.', idx, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

}
