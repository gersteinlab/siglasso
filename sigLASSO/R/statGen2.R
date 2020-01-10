statGen2 <-
function (N, mut_number, spectrum, sigs, sig_number, noise, priorPer, 
    Nsubsample) 
{
    if (missing(priorPer)) {
        priorPer = 0
    }
    if (missing(Nsubsample)) {
        Nsubsample = 0
    }
    spectrum.deconstruct <- data.frame(t(spectrum)/mut_number)
    colnames(spectrum.deconstruct) <- colnames(sigs.input)
    rownames(spectrum.deconstruct) <- seq(1, N)
    coef_deconstruct <- rep(0, dim(spectrum)[2])
    coef_all = rep(0, dim(spectrum)[2])
    for (i in c(1:dim(spectrum)[2])) {
        predictor <- rep(1, 30)
        if (priorPer > 0) {
            priorN = round(sig_number * priorPer)
            priorSet = sample(as.numeric(which(sigs[, i] > 0)), 
                priorN)
            predictor[priorSet] = 0
        }
        cv.glmnet <- cv.glmnet(data.matrix(sig[, -1]), spectrum[, 
            i]/mut_number, alpha = 1, intercept = F, lower.limit = 0, 
            lambda = seq(1e-05, 0.01, by = 1e-05), penalty.factor = predictor, 
            nfolds = 12)
        fit <- glmnet(data.matrix(sig[, -1]), spectrum[, i]/mut_number, 
            alpha = 1, intercept = F, lower.limit = 0, lambda = cv.glmnet$lambda.1se, 
            penalty.factor = predictor)
        if (Nsubsample > 0) {
            fit_sub_coef <- matrix(0, Nsubsample, 30)
            sum = NULL
            for (j in c(1:96)) {
                sum = c(sum, rep(j, spectrum[j, i]))
            }
            for (j in c(1:Nsubsample)) {
                set <- sample(sum, round(mut_number/2))
                set_data = matrix(0, nrow = 96, ncol = 1)
                for (k in c(1:length(set))) {
                  set_data[set[k], 1] = set_data[set[k], 1] + 
                    1
                }
                set_data <- set_data/length(set)
                cv.glmnet <- cv.glmnet(data.matrix(sig[, -1]), 
                  set_data, alpha = 1, intercept = F, lower.limit = 0, 
                  lambda = seq(1e-05, 0.01, by = 1e-05), penalty.factor = predictor, 
                  nfolds = 12)
                fit_set <- glmnet(data.matrix(sig[, -1]), set_data, 
                  alpha = 1, intercept = F, lower.limit = 0, 
                  lambda = cv.glmnet$lambda.1se, penalty.factor = predictor)
                fit_sub_coef[j, ] <- coef(fit_set)[-1, ]
            }
            fit_sub_coef[fit_sub_coef > 0] <- 1
            fit_sub_select <- colSums(fit_sub_coef)/Nsubsample
            fit_sub_select <- ifelse(fit_sub_select < 0.6, 0, 
                1)
        }
        coef <- coef(fit)[-1, ]
        coef[fit_sub_select == 0] <- 0
        coef_all = cbind(coef_all, data.matrix(coef))
        coef_deconstruct <- cbind(coef_deconstruct, t(whichSignatures(tumor.ref = spectrum.deconstruct, 
            signatures.ref = signatures.cosmic, sample.id = i, 
            contexts.needed = TRUE, tri.counts.method = "default")[1]$weights))
    }
    plot_coef <- coef_all[, -1]
    coef_deconstruct <- coef_deconstruct[, -1]
    F_neg = rep(NA, N)
    F_poz = rep(NA, N)
    T_poz = rep(NA, N)
    for (i in 1:N) {
        inferred_set <- as.numeric(which(plot_coef[, i] > 0))
        true_set <- as.numeric(which(sigs[, i] > 0))
        F_neg[i] <- length(which(!(true_set %in% inferred_set)))
        F_poz[i] <- length(which(!(inferred_set %in% true_set)))
        T_poz[i] <- length(which(inferred_set %in% true_set))
    }
    F_neg_d = rep(NA, N)
    F_poz_d = rep(NA, N)
    T_poz_d = rep(NA, N)
    for (i in 1:N) {
        inferred_set <- as.numeric(which(coef_deconstruct[, i] > 
            0))
        true_set <- as.numeric(which(sigs[, i] > 0))
        F_neg_d[i] <- length(which(!(true_set %in% inferred_set)))
        F_poz_d[i] <- length(which(!(inferred_set %in% true_set)))
        T_poz_d[i] <- length(which(inferred_set %in% true_set))
    }
    return(list((plot_coef - sigs)^2, cbind(T_poz, F_poz, F_neg), 
        (coef_deconstruct - sigs)^2, cbind(T_poz_d, F_poz_d, 
            F_neg_d), plot_coef, coef_deconstruct, sigs))
}
