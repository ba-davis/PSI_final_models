


# normalize betas function
#   beta_sub: samples as rows, cpgs as columns
#             samples are rownames(beta_sub), no sample column
#   ref: model reference table with cols <cpg> <estimate> <mean> <weighted_mean> <sd> <y_int>
norm_betas <- function(beta_sub, ref) {
    # centering (subtract the mean from each cpg) and
    # scaling (divide by standard deviation)

    # subset ref table to contain only cpgs found in beta_sub
    ref_sub <- ref[ref$cpg %in% colnames(beta_sub), ]
    beta_sub <- beta_sub[, match(ref_sub$cpg, colnames(beta_sub))]

    # sanity check
    if (!(identical(colnames(beta_sub), ref_sub$cpg))) {
    print(paste0("WARNING: cpgs are not in proper order.",
    " Problem with PSI calculation."))
    }

    # for each cpg/column in the input beta matrix
    #   - subtract the mean for that cpg in ref from
    #     each row in the cpg column in user input beta matrix.
    #   - Then, divide this new value by the sd value for the cpg
    for (i in seq_len(ncol(beta_sub))) {
        # subtract the mean
        beta_sub[[colnames(beta_sub)[i]]] <- beta_sub[[colnames(beta_sub)[i]]] -
        ref_sub[ref_sub$cpg == colnames(beta_sub)[i], "mean"]
        # divide by the standard deviation
        beta_sub[[colnames(beta_sub)[i]]] <- beta_sub[[colnames(beta_sub)[i]]] /
        ref_sub[ref_sub$cpg == colnames(beta_sub)[i], "sd"]
    }

    return(beta_sub)
}

# multiply beta value by coefficient and store in matrix
#   beta_sub: samples as rows, cpgs as columns
#             samples are rownames(beta_sub), no sample column
#   ref: model reference table with cols <cpg> <estimate> <mean> <weighted_mean> <sd> <y_int>
calculate_psi <- function(beta_sub, ref) {
    # put the cpg rows of transposed beta_sub in the same order as
    # the cpg rows in ref table ref$cpg
    B <- t(beta_sub)
    # subset ref table to contain only cpgs found in beta_sub
    ref_sub <- ref[ref$cpg %in% rownames(B), ]
    B <- B[match(ref_sub$cpg, rownames(B)), ]
    # sanity check
    if (!(identical(rownames(B), ref_sub$cpg))) {
    print(paste0("WARNING: cpgs are not in proper order.",
    " Problem with PSI calculation."))
    }

    y <- NULL
    i <- NULL
    A <- ref_sub$estimate
    y <- matrix(nrow = nrow(B), ncol = ncol(B))
    # for each row (cpg) of B, multiply the entire row by
    #   the coefficient for that cpg
    for (i in 1:nrow(B)){
        y[i, ] <- rbind(as.numeric(B[i, ]) * A[i])
    }

    s <- colSums(y)
    # store psi score per sample in dataframe
    d2 <- data.frame(colnames(B), s)

    return(d2)
}

# convert logit probability value to probability value
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# take input beta matrix and predict to get output summary table
# input_beta: dataframe of beta values
#             samples as rows, cpgs as columns
#             samples are rownames(input_beta), there is no sample column
# ref: model reference table with cols:
#        <cpg> <estimate> <mean> <weighted_mean> <sd> <y_int>
# fill_missing_cpgs: TRUE or FALSE, whether to fill beta values of missing
#                    CpGs with weighted mean of training data, default FALSE
# prob_cutoff: cutoff value to use for classifying smoker or nonsmoker
#              default 0.5
#
# returns output df of results with columns:
#   <sample> <psi> <logit_prob> <probability_smoker> <probability_nonsmoker> <predicted_class>

psi_predict <- function(input_beta, ref, fill_missing_cpgs = FALSE, prob_cutoff = 0.5) {
    # check how many of the required cpgs for the model are present as columns
    # in user input beta matrix
    print("Checking input CpGs for those required by the model.")
    required_cpg_num <- length(ref$cpg)
    present_cpg_num <- length(ref$cpg[ref$cpg %in% colnames(input_beta)])
    print(paste0(present_cpg_num, " out of ", required_cpg_num,
        " required CpGs are present."))
    if (present_cpg_num != required_cpg_num) {
        print(paste0(required_cpg_num - present_cpg_num, " CpGs are missing."))
    }

    # subset input beta matrix to contain only present required cpgs
    beta_sub <- input_beta[, colnames(input_beta) %in% ref$cpg]

    # if you want to use weighted mean values for missing cpgs...
    # add the weighted mean beta values for the missing cpgs
    if (fill_missing_cpgs == TRUE) {
        print("Using weighted mean values for missing cpgs.")
        missing_cpgs <- ref$cpg[!(ref$cpg %in% colnames(beta_sub))]
        for (i in seq_along(missing_cpgs)) {
            beta_sub[[missing_cpgs[i]]] <- ref$weighted_mean[ref$cpg == missing_cpgs[i]]
        }
    }
    else {
        print("Using only present required CpGs for prediction.")
    }

    # Calculate PSI from unnormalized beta matrix
    #print("Calculating PSI.")
    #d2 <- calculate_psi(beta_sub, mytab)
    #colnames(d2) <- c("sample", "psi")

    # normalize the beta matrix according to training data
    #print("Normalizing betas.")
    beta_norm <- norm_betas(beta_sub, ref)

    # Calculate PSI score from the normalized betas
    d2 <- calculate_psi(beta_norm, ref)
    colnames(d2) <- c("Sample_Name", "PSI")

    # add y-intercept to the PSI
    y_int <- unique(ref$y_int)
    d2$logit_prob <- d2$PSI + y_int

    # add predicted probability columns
    d2$probability_smoker <- sapply(d2$logit_prob, logit2prob)
    d2$probability_nonsmoker <- 1 - d2$probability_smoker

    # add predicted class column
    print(paste0("Using ", prob_cutoff, " as probability threshold for class prediction."))
    d2$predicted_class <- ifelse(d2$probability_smoker >= prob_cutoff,
        "smoker", "nonsmoker")

    # remove logit_prob column
    d2$logit_prob <- NULL

    return(d2)
}
