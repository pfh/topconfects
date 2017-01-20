
#
# Differential exon usage, 3' end points, 5' start points, etc
#
#



#' Group confects (differential exon usage, etc)
#'
#' Find differential exon usage, etc.
#'
#' If the order of the members of a group is important, ensure the order of rows is correct in the original matrix.
#'
#' Groups with less than two members will be ignored.
#'
#' To construct a group design matrix, the design from \code{fit} will be repeated in a block diagonal matrix.
#'
#' @param fit An edgeR DGEGLM object.
#'
#' @param group_id A factor of length \code{nrow(fit)}, assigning items to groups (eg genes).
#'
#' @param group_effect A \code{function(n)} returning an effect size object for a group with \code{n} members.
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
# edger_group_confects <- function(fit, group_id, group_effect, fdr=0.05, max=30.0, step=0.01) {
#     group_id <- factor(group_id)
#     group_effect <- memoise(group_effect)
#
#     assert_that(is(fit, "DGEGLM"))
#     assert_that(length(group_id) == nrow(fit))
#
#     members <- split(seq_len(nrow(fit)), group_id)
#     members <- members[ map_int(members,length) >= 2 ]
#     n <- length(members)
#
#     group_design <- function(m) {
#         design <- fit$design
#         zero_block <- matrix(0,nrow=nrow(design),ncol=ncol(design))
#         big_design <- matrix(0,nrow=0,ncol=n*ncol(design))
#         for(i in seq_len(m)) {
#             big_design <- rbind(big_design,
#                 do.call(cbind,c(
#                     rep(list(zero_block), i-1),
#                     list(design),
#                     rep(list(zero_block), m-i)
#                 ))
#             )
#         }
#
#         big_design
#     }
#     group_design <- memoise(group_design)
#
#     fit <- function(i, cons=NULL, initial=NULL) {
#         m <- length(members[[i]])
#         ...
#     }
# }








