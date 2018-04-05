trim <- function (x) {
	x <- gsub("^\\s+|\\s+$", "", x)
	x <- gsub("^\\\"|\"$", "", x)
	return(x)
}

remove_excess_whitespace <- function(x) x <- gsub("\\s+", " ", x)

# Read in all the data table codings
to_read <- paste("./WAS/codings/",
	dir("./WAS/codings"), sep="")
codings_tables <- list()

for(i in to_read) {
	name <- gsub("./WAS/codings/coding(.*).tsv", "\\1", i)
	codings_tables[[name]] <- read.table(i, header=TRUE, sep='\t', quote="", stringsAsFactors=FALSE)
}

get_subtype <- function(x) paste(x[2:length(x)],collapse=" ")

get_notes <- function(tsv_data, log_file, outcome_info, codings_tables, start_column=4)
{	

	if (ncol(tsv_data) > (start_column-1)) {
		# Create character matrix 'notes' that we will write to disk and pass to Manny.
		notes <- matrix(nrow=(ncol(tsv_data)-start_column+1),ncol=10)
		colnames(notes) <- c("FieldID", "Field", "N.non.missing", "N.missing", "N.cases", "N.controls",
			"Notes", "PHESANT.notes", "PHESANT.reassignments", "warning.for.case.control")
		# Rownames are the FieldIDs
		notes[,1] = colnames(tsv_data)[start_column:ncol(tsv_data)]
		samples <- nrow(tsv_data)
		k <- 1
		for(col in colnames(tsv_data)[start_column:ncol(tsv_data)]){
			type <- class(tsv_data[col][,1])
			if(type == "numeric") {
				var <- substr(col, 2, nchar(col))
				where <- which(outcome_info$FieldID == var)
				col_name <- paste(trim(outcome_info$Field[where]),"-",col)
				# The second column is the Field name.
				notes[k,2] <- trim(outcome_info$Field[where])
				# The seventh column is the 'Notes' field in variable-info file:
				notes[k,7] <- remove_excess_whitespace(as.character(trim(outcome_info$Notes[where])))
				
				if (log_file != FALSE){
					# The eighth column is information about how the data is parsed using PHESANT:
					matching_line <- trim(grep(paste('^', var, '_', sep=""), readLines(log_file), value=TRUE))
					notes[k,8] <- matching_line

					# The ninth column is usually empty, unless there are PHESANT reassignments, in which
					# case we detail those reassignments here:
					if (length(grep("reassignments", matching_line)) > 0) {
						notes[k,9] <- trim(gsub("^.*(reassignments: .*?)\\|\\|.*", "\\1", matching_line))
					}
				}

			} else if (type=="logical" | type=="integer") {
				var <- strsplit(col, split="_")[[1]][1]
				var <- substr(var, 2, nchar(var))
				subvar <- strsplit(col, split="_")[[1]][2]
				where <- which(outcome_info$FieldID == var)
				
				# Add the notes:
				# The seventh column is the 'notes' field in variable-info file:
				notes[k,7] <- remove_excess_whitespace(as.character(trim(outcome_info$Notes[where])))
				
				if (log_file != FALSE) {
					# The eight column is information about how the data is parsed using PHESANT:
					matching_line <- trim(grep(paste('^', var, '_', sep=""), readLines(log_file), value=TRUE))
					
					if (length(grep("^.*CAT-MULTIPLE", matching_line)) > 0) {
						new_matching_line <- paste(gsub(paste("^(.*?)\\|.*(CAT-MUL-BINARY-VAR", subvar, ".*?)CAT.*"), "\\1 \\|\\| \\2", matching_line))
						if(matching_line == new_matching_line){
							new_matching_line <- paste(gsub(paste("^(.*?)\\|.*(CAT-MUL-BINARY-VAR", subvar, ".*)"), "\\1 \\|\\| \\2", matching_line))
						}
						notes[k,8] <- trim(new_matching_line)
					} else if (length(grep("CAT-SINGLE-UNORDERED", matching_line)) > 0) {
						new_matching_line <- gsub(paste("^(.*?)\\|.*(Inc.*?: ", subvar, "\\([0-9]+\\)).*", sep=""),
							paste("\\1 \\|\\| CAT-SINGLE \\|\\| CAT-SINGLE-BINARY-VAR:", subvar, " \\|\\| \\2 \\|\\|"), matching_line)
						notes[k,8] <- trim(new_matching_line)
					} else {
						notes[k,8] <- trim(matching_line)
					}

					# The ninth column is usually empty, unless there are PHESANT reassignments, in which
					# case we detail those reassignments here:
					if (length(grep("reassignments", matching_line)) > 0) {
						notes[k,9] <- trim(gsub("^.*(reassignments: .*?)\\|\\|.*", "\\1", matching_line))
					}
				}

				col_name <- paste(trim(outcome_info$Field[where]), "-", col)

				col_subname <- as.character(outcome_info$Coding[where])
				if(nchar(col_subname)>0 && !is.na(subvar)) {
					to_parse <- strsplit(col_subname, split="\\|")[[1]]
					if(length(to_parse) > 1) {
						where_subvar <- which(sapply(strsplit(trim(to_parse)," "), "[[",1) == subvar) 
						col_subname <- ifelse(length(where_subvar) > 0,
							get_subtype(strsplit(trim(to_parse[where_subvar])," ")[[1]]),
							"PHESANT recoding")
						# The first column is the Field name - if we've entered this if statement,
						# we need to include the subfield (as the categorical has been split into 
						# loads of booleans).
						notes[k,2] <- paste(trim(outcome_info$Field[where]),": ",col_subname, sep="")
					} else {
						notes[k,2] <- trim(outcome_info$Field[where])
					}
				} else {
					# The first column is the Field name.
					notes[k,2] <- trim(outcome_info$Field[where])
				}

				if(length(grep("Too many", col_subname))>0) {
					# Then we need to look in the corresponding tsv_data-coding table to get the label.
					coding <- trim(gsub("^.*id=","",col_subname))
					where_coding <- which(codings_tables[coding][[1]]$coding == subvar)
					col_subname <- codings_tables[coding][[1]]$meaning[where_coding]
					# The first column is the Field name - if we've entered this if statement,
					# we need to include the subfield (as the categorical has been split into 
					# loads of booleans) - this is accessed from the associated coding table for the 
					# phenotype.
					notes[k,2] <- paste(trim(outcome_info$Field[where]),": ",col_subname, sep="")
				}
				y <- table(tsv_data[col][,1])

				# Finally, include the case and control numbers for the logical and integer phenotypes
				# with just two categories.
				if(type=="logical") {
					notes[k,5] <- y[names(y) == "TRUE"]
					notes[k,6] <- y[names(y) == "FALSE"]
				} else {
					if(length(y) == 2) {
						# I assume that 1 encodes a positive here:
						if(all(names(y) == c("0", "1"))) {
							notes[k,5] <- y[names(y) == "1"]
							notes[k,6] <- y[names(y) == "0"]
						} else {
							# All bets are off if the variable is not 0,1, 
							# but I report the numbers in each of the two categories.
							# This should be the larger number encoding the positive, according to the 
							# laws of the table function.
							notes[k,5] <- y[1]
							notes[k,6] <- y[2]
							notes[k,10] <- "YES"
						}
					}
				}

			} else {
				print("Error: not one of the specified types!")
				print(type)
				print(k)
			}
			# The third column is the number of missing data for this phenotype.
			n_miss <- sum(is.na(tsv_data[col][,1]))
			notes[k,4] <- n_miss
			# The second column is the number of non-missing data for this phenotype.
			n_non_miss <- sum(!is.na(tsv_data[col][,1]))
			notes[k,3] <- n_non_miss
			k <- k+1
		}
		# Get rid of the X that's appended at the start of each variable.
		notes[,1] = substr(notes[,1], 2, nchar(notes[,1]))
		return(notes)
	}
}

include_PHESANT_reassignment_names <- function(pheno_summary, outcome_info)
{	
	pheno_summary <- read.table(pheno_summary, sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE) 

	# Find the variables that have had a reassignment.
	where_reassignment <- which(!is.na(pheno_summary[,9]))
	reassignments <- as.character(pheno_summary[where_reassignment,9])
	
	if (length(reassignments) > 0) {
		reassignments_before_after_list <- list()
		reassignments <- gsub("reassignments: ", "", reassignments)
		reassignments_list <- strsplit(reassignments, split="\\|")
		for(i in 1:length(reassignments_list)) {
			reassignments_before_after_list[[i]] <- matrix(unlist(strsplit(reassignments_list[[i]], "=")),ncol=2, byrow=TRUE)
		}

		for(i in 1:length(where_reassignment)) {
			if(length(grep("_", pheno_summary[,1][where_reassignment[i]])) > 0) {
				for(j in 1:nrow(reassignments_before_after_list[[i]])) {
					if(length(grep(paste("_", reassignments_before_after_list[[i]][j,2], sep=""),
						pheno_summary[,1][where_reassignment[i]])) > 0) {
						# We've found this recoded variable in this row.
						# I want to now find what this codes for - so need to find the coding table for this variable
						# from the outcome-info file, then look at what the uncoded version of the phenotype is.
						var <- strsplit(pheno_summary[,1][where_reassignment[i]], split="_")[[1]][1]
						where_var <- which(outcome_info$FieldID == var)

						i_subname <- as.character(outcome_info$Coding[where_var])
						subvar <- reassignments_before_after_list[[i]][j,1]

						if(nchar(i_subname)>0 && !is.na(subvar)) {
							to_parse <- strsplit(i_subname, split="\\|")[[1]]
							if(length(to_parse) > 1) {
								where_subvar <- which(sapply(strsplit(trim(to_parse)," "), "[[",1) == subvar) 
								i_subname <- ifelse(length(where_subvar) > 0,
									get_subtype(strsplit(trim(to_parse[where_subvar])," ")[[1]]),
									"PHESANT recoding")
								
								pheno_summary[where_reassignment[i],2] <- gsub("PHESANT recoding", i_subname, pheno_summary[where_reassignment[i],2])
							}
						}
					
						if(length(grep("Too many", i_subname))>0) {
							# Then we need to look in the corresponding tsv_data-coding table to get the label.
							coding <- trim(gsub("^.*id=","",i_subname))
							where_coding <- which(codings_tables[coding][[1]]$coding == reassignments_before_after_list[[i]][j,1])
							i_subname <- codings_tables[coding][[1]]$meaning[where_coding]

							pheno_summary[where_reassignment[i],2] <- paste(pheno_summary[where_reassignment[i],2], i_subname, sep="")
						}
					}
				}
			}
		}
	}
	return(pheno_summary)
}
