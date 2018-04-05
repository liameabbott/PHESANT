
source("summarise_phenotypes.r")

outcome_info <- read.table("./variable-info/outcome_info_final.tsv",
                       sep='\t', quote="", comment.char="", header=TRUE)

for (i in seq(1, 4)) {

    tsv_data_both <- read.table(paste0('../ukb31063.phesant.both_sexes.', i, '.log', header=TRUE, sep='\t')
    tsv_data_females <- read.table(paste0('../ukb31063.phesant.female.', i, '.log', header=TRUE, sep='\t')
    tsv_data_males <- read.table(paste0('../ukb31063.phesant.male.', i, '.log', header=TRUE, sep='\t')

    log_file_both <- read.table(paste0('../ukb31063.phesant.both_sexes.', i, '.log', header=TRUE, sep='\t')
    log_file_females <- read.table(paste0('../ukb31063.phesant.female.', i, '.log', header=TRUE, sep='\t')
    log_file_males <- read.table(paste0('../ukb31063.phesant.male.', i, '.log', header=TRUE, sep='\t')

    notes_both <- get_notes(tsv_data_both, log_file_both, outcome_info, codings_tables)
    write.table(notes_both, file=paste0('../ukb31063.both_sexes.tmp.', i, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_females <- get_notes(tsv_data_females, log_file_females, outcome_info, codings_tables)
    write.table(notes_females, file=paste0('../ukb31063.female.tmp.', i, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_males <- get_notes(tsv_data_males, log_file_males, outcome_info, codings_tables)
    write.table(notes_males, file=paste0('../ukb31063.male.tmp.', i, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_codings_included_both <- include_PHESANT_reassignment_names(paste0('../ukb31063.both_sexes.tmp.', i, '.tsv'), outcome_info)
    write.table(notes_codings_included_both, file=paste0('../ukb31063.both_sexes.phenosummary.', i, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_codings_included_females <- include_PHESANT_reassignment_names(paste0('../ukb31063.female.tmp.', i, '.tsv'), outcome_info)
    write.table(notes_codings_included_both, file=paste0('../ukb31063.female.phenosummary.', i, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    notes_codings_included_males <- include_PHESANT_reassignment_names(paste0('../ukb31063.male.tmp.', i, '.tsv'), outcome_info)
    write.table(notes_codings_included_both, file=paste0('../ukb31063.male.phenosummary.', i, '.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

}


#log_both <- read.table('../')
#filename <- "../ukb31063.output"
#tsv_filename <- paste(filename, ".tsv", sep="")
#log_file <- paste(filename, ".log", sep="")


#qc_data <- read.table("~/results/ukb1859_qc.tsv", header=TRUE)

#outcome_info <- read.table("~/Repositories/PHESANT/variable-info/outcome_info_final.tsv",
#					   sep='\t', quote="", comment.char="", header=TRUE)

#samples_for_removal <- as.character(read.table("~/results/w1859_20170726_participantwithdrawallist.csv")$V1)

#notes_for_manny_1859 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, qc_data, samples_for_removal)

#write.table(notes_for_manny_1859, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)

#pheno_summary_PHESANT <- "~/results/ukb1859_phenosummary_final.tsv"
#notes_for_manny_1859_PHESANT_codings_included <- include_PHESANT_reassignment_names(pheno_summary, outcome_info)
#write.table(notes_for_manny_1859_PHESANT_codings_included, file=pheno_summary_PHESANT, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
