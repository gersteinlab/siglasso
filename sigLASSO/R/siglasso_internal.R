siglasso_internal <-
function (spectrum, sig, prior, adaptive, elastic_net, gamma, 
    alpha_min, iter_max, sd_multiplier, conf) 
{
    p_star_hat <- data.matrix(spectrum/sum(spectrum))
    p_star_hat_now <- p_star_hat
    p_hat <- rep(0, nrow(sig))
    l_last <- -Inf
    iter_number = 0
    predictor <- data.matrix(sig)
    penalty <- rep(1, ncol(sig))
    las_pm <- estimate_lasso_parameters(data.matrix(sig), p_star_hat, 
        gamma, penalty, prior, adaptive = T, sd_multiplier, elastic_net, 
        conf)
    predictor <- las_pm$predictor
    lambda_min1se <- las_pm$lambda_min1se
    penalty <- las_pm$penalty
    lambda_seq <- las_pm$lambda_seq
    ols_fit <- las_pm$ols_fit
    alpha_min1se <- las_pm$alpha_min1se
    while (TRUE) {
        iter_number = iter_number + 1
        if (is.vector(predictor) || ncol(predictor) < 2) {
            coef_hat <- rep(0, ncol(sig))
            coef_hat[ols_fit$x > 0] <- coef(lm(predictor ~ p_star_hat - 
                1))
            lambda_min1se <- 0
            K = length(predictor)
        }
        else {
            K = nrow(predictor)
            fit <- glmnet(predictor, p_star_hat, alpha = alpha_min1se, 
                intercept = F, lower.limit = 0, lambda = lambda_seq, 
                standardize = T, penalty.factor = penalty)
            if (adaptive) {
                coef_hat <- rep(0, ncol(sig))
                coef_hat[ols_fit$x > 0] <- coef(fit, s = lambda_min1se, 
                  x = predictor, y = p_star_hat, penalty.factor = penalty)[-1, 
                  ]
            }
            else {
                coef_hat <- coef(fit, s = lambda_min1se, x = predictor, 
                  y = p_star_hat, penalty.factor = penalty)[-1, 
                  ]
            }
        }
        p_hat <- sig %*% coef_hat
        alpha = max(1/sum((p_hat - p_star_hat)^2) * (K - length(which(coef_hat > 
            0)))/conf, alpha_min)
        f <- function(lambda) {
            sum(p_hat) + (K/alpha) * lambda - 2 + sum(sqrt((p_hat + 
                lambda/alpha)^2 + 4 * spectrum/alpha))
        }
        low_lim = -sum(spectrum)
        uniroot_iter_number = 0
        while (f(low_lim) > 0) {
            uniroot_iter_number = uniroot_iter_number + 1
            if (uniroot_iter_number > 20) {
                stop("uniroot cannot find the root!")
            }
            low_lim = low_lim * 2
        }
        lambda_hat <- uniroot(f, c(low_lim, alpha/K))$root
        p_star_hat <- (p_hat + lambda_hat/alpha + sqrt((p_hat + 
            lambda_hat/alpha)^2 + 4 * spectrum/alpha))/2
        l = sum(spectrum * log(p_star_hat), na.rm = T) - alpha * 
            sum((p_star_hat - p_hat)^2) - lambda_min1se * sum(coef_hat) * 
            alpha
        if (sum((p_star_hat_now - p_star_hat)^2) < 1e-08 || iter_number > 
            iter_max) {
            print(iter_number)
            break
        }
        else if (sum((p_star_hat_now - p_star_hat)^2) > 0.001) {
            las_pm <- estimate_lasso_parameters(data.matrix(sig), 
                p_star_hat, gamma, penalty, prior, adaptive = T, 
                sd_multiplier, elastic_net, conf)
            predictor <- las_pm$predictor
            lambda_min1se <- las_pm$lambda_min1se
            alpha_min1se <- las_pm$alpha_min1se
            penalty <- las_pm$penalty
            lambda_seq <- las_pm$lambda_seq
            ols_fit <- las_pm$ols_fit
        }
        l_last = l
        p_star_hat_now = p_star_hat
    }
    names(coef_hat) <- colnames(sig)
    return(coef_hat)
}
