#' Data frames row binding.
#'
#' BindFillRows
#' @keywords internal
#' @description Bind data frames by rows after filling missing columns with NA.
#' @param df_Lst <data.frames or list[data.frame]>: Data frames to bind or list of data.frames. If is a data.frame create a list with arguments `df_Lst` and `...`, else `...` are ignored.
#' @param ... <data.frames or list[data.frame]>: Data frames to bind or list of data.frames.
#' @return The binded data frame
#' @examples
#' df1 <- data.frame(a = seq_len(5), b = c(6:10))
#' df2 <- data.frame(a = c(11:15), b = c(16:20), c = LETTERS[seq_len(5)])
#' BindFillRows(df1, df2)
#' BindFillRows(list(df1, df2))
#'
BindFillRows <- function(
    df_Lst, ...
) {
    if (is.data.frame(df_Lst)) {
        df_Lst <- list(df_Lst, ...)
    }
    df_Lst <- lapply(
        seq_along(df_Lst),
        function(data.ndx) {
            data.df <- df_Lst[[data.ndx]]
            dataNames.chr <- lapply(df_Lst[-data.ndx], names) |>
                unlist() |>
                unique()
            data.df[setdiff(dataNames.chr, names(data.df))] <- NA
            return(data.df)
        }
    )
    return(do.call(rbind, df_Lst))
}
