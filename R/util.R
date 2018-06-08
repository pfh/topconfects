
#
# Ensure a numeric vector has a certain length. Repeat input if input length is 1.
#
broadcast <- function(vec, n) {
    assert_that(is.numeric(vec))
    names(vec) <- NULL
    if (length(vec) == 1) 
        vec <- rep(vec, n)
    assert_that(length(vec) == n)
    vec
}
