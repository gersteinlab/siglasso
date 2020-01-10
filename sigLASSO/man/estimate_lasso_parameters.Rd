\name{estimate_lasso_parameters}
\alias{estimate_lasso_parameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate_lasso_parameters(predictor, p_star_hat, gamma, penalty, prior, adaptive, sd_multiplier, elastic_net, conf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictor}{
%%     ~~Describe \code{predictor} here~~
}
  \item{p_star_hat}{
%%     ~~Describe \code{p_star_hat} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{penalty}{
%%     ~~Describe \code{penalty} here~~
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
}
  \item{adaptive}{
%%     ~~Describe \code{adaptive} here~~
}
  \item{sd_multiplier}{
%%     ~~Describe \code{sd_multiplier} here~~
}
  \item{elastic_net}{
%%     ~~Describe \code{elastic_net} here~~
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (predictor, p_star_hat, gamma, penalty, prior, adaptive, 
    sd_multiplier, elastic_net, conf) 
{
    lambda_min1se = 0
    alpha_min1se <- 1
    alpha_seq = seq(1, 0.1, -0.1)
    mysd <- function(y) sqrt(sum((y - mean(y))^2)/length(y))
    if (adaptive) {
        ols_fit <- nnls(predictor, p_star_hat)
        predictor <- data.matrix(predictor[, ols_fit$x > 0])
        penalty <- 1/(ols_fit$x[ols_fit$x > 0])^gamma * prior[ols_fit$x > 
            0]
    }
    if (!is.vector(predictor) && ncol(predictor) > 1) {
        Y = p_star_hat
        X = predictor
        n = ncol(predictor)
        if (elastic_net) {
            lambda_alpha = NULL
        }
        predict_sse <- NULL
        K = nrow(predictor)
        for (h in seq(1, 20)) {
            while (TRUE) {
                group <- NULL
                for (g in seq(1, 6)) {
                  group <- c(group, sample(rep(seq(1, 8), 2)))
                }
                flag = 0
                for (k in seq(1, 8)) {
                  if (var(p_star_hat[which(group != g)]) == 0) {
                    flag = 1
                    break
                  }
                  if (0 \%in\% apply(predictor[which(group != g), 
                    ], 2, var)) {
                    flag = 1
                    break
                  }
                }
                if (flag == 0) {
                  break
                }
            }
            for (g in seq(1, 8)) {
                train_set <- p_star_hat[which(group != g)]
                test_set <- p_star_hat[which(group == g)]
                train_predictor <- predictor[which(group != g), 
                  ]
                test_predictor <- predictor[which(group == g), 
                  ]
                if (elastic_net) {
                  lambda_alpha_search = NULL
                  predict_sse_alpha_search = NULL
                  for (a in alpha_seq) {
                    lambda_max = log10(max(glmnet(X, Y, alpha = a, 
                      intercept = F, lower.limit = 0, penalty.factor = penalty)$lambda))
                    lambda_seq <- 10^seq(lambda_max, lambda_max - 
                      8, -0.05)
                    fit <- glmnet(train_predictor, train_set, 
                      alpha = a, intercept = F, lower.limit = 0, 
                      lambda = lambda_seq, penalty.factor = penalty)
                    fit_predict <- predict(fit, test_predictor, 
                      alpha = a, intercept = F, lower.limit = 0, 
                      penalty.factor = penalty)
                    predict_sse_alpha_search <- c(predict_sse_alpha_search, 
                      apply(fit_predict, 2, function(x) sum(x - 
                        test_set)^2))
                    lambda_alpha_search <- c(lambda_alpha_search, 
                      lambda_seq)
                  }
                  lambda_alpha <- rbind(lambda_alpha, lambda_alpha_search)
                  predict_sse <- rbind(predict_sse, predict_sse_alpha_search)
                }
                else {
                  a = 1
                  lambda_max = log10(max(glmnet(X, Y, alpha = a, 
                    intercept = F, lower.limit = 0, penalty.factor = penalty)$lambda))
                  lambda_seq <- 10^seq(lambda_max, lambda_max - 
                    8, -0.05)
                  fit <- glmnet(train_predictor, train_set, alpha = 1, 
                    intercept = F, lower.limit = 0, lambda = lambda_seq, 
                    penalty.factor = penalty)
                  fit_predict <- predict(fit, test_predictor, 
                    alpha = 1, intercept = F, lower.limit = 0, 
                    penalty.factor = penalty)
                  predict_sse <- rbind(predict_sse, apply(fit_predict, 
                    2, function(x) sum(x - test_set)^2))
                }
            }
        }
        predict_mse <- apply(predict_sse, 2, mean)
        predict_mse_1sd <- apply(predict_sse, 2, sd)
        mse_max_allowed <- (predict_mse + sd_multiplier * predict_mse_1sd)[which.min(predict_mse)]
        if (elastic_net) {
            idx <- which(predict_mse <= mse_max_allowed)
            alpha_min1se <- max(sort(rep(alpha_seq, length(lambda_seq)), 
                decreasing = T)[idx])
            idx_col <- which(sort(rep(alpha_seq, length(lambda_seq)), 
                decreasing = T) == alpha_min1se)
            idx <- intersect(idx, idx_col)
            lambda_min1se <- max(lambda_alpha[1, idx])
            print("alaph min")
            print(alpha_min1se)
        }
        else {
            lambda_min1se <- max(lambda_seq[which(predict_mse <= 
                mse_max_allowed)])
        }
        return(list(lambda_min1se = lambda_min1se, penalty = penalty, 
            predictor = predictor, lambda_seq = lambda_seq, alpha_min1se = alpha_min1se, 
            ols_fit = ols_fit))
    }
    return(list(lambda_min1se = lambda_min1se, penalty = penalty, 
        predictor = predictor, lambda_seq = NULL, alpha_min1se = alpha_min1se, 
        ols_fit = ols_fit))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
