
source("predict_functions.R")

#-------------#

#------------------#
# LOAD BETA MATRIX #
#------------------#

# read in example marsit beta matrix
beta <- read.delim(gzfile("marsit_beta.txt.gz"))

# Many beta matrices will use CpGs as rows and Samples as columns
# Transpose beta matrix to put samples as rows and CpGs as columns
beta <- as.data.frame(t(beta))
dim(beta)

#----------------------------#
# LOAD MODEL REFERENCE TABLE #
#----------------------------#

# read in the reference table for the desired model
ref <- read.delim("../data/model1_EPIC_reference_table.txt")

#---------#
# PREDICT #
#---------#

# calculate predictions on the new data
# returns psi, norm_psi, prob_smoker, prob_nonsmoker, predicted_class
res <- psi_predict(beta, ref)

res2 <- psi_predict(beta, ref, fill_missing_cpgs=TRUE)
