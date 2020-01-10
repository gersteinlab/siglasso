getmin <-
function (lambda, cvm, cvsd) 
{
    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    lambda.min = max(lambda[idmin], na.rm = TRUE)
    idmin = match(lambda.min, lambda)
    semin = (cvm + 3 * cvsd)[idmin]
    idmin = cvm <= semin
    lambda.1se = max(lambda[idmin], na.rm = TRUE)
    list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}
