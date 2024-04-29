#' Aggregation plot
#'
#' ggAPA
#' @description Create a ggplot object used for plot aggregation.
#' @param aggregatedMtx <matrix> : The matrix to plot.
#' (Default NULL)
#' @param title <character> : The title of plot.
#' (Default NULL)
#' @param trim <numeric> : A number between 0 and 100
#' that gives the percentage of trimming. (Default 0)
#' @param tails <character> : Which boundary must be trimmed?
#' If it's both, trim half of the percentage in inferior and superior.
#' see `QtlThreshold`. (Default "both")
#' @param colMin <numeric> : Minimal value of Heatmap,
#' force color range. If `NULL` automatically find. (Default NULL)
#' @param colMid <numeric> : Center value of Heatmap,
#' force color range. If `NULL` automatically find. (Default NULL)
#' @param colMax <numeric> : Maximal value of Heatmap,
#' force color range. If `NULL` automatically find. (Default NULL)
#' @param colBreaks <numeric> : Repartition of colors.
#'  If `NULL` automatically find. (Default NULL)
#' @param blurPass <numeric> : Number of blur pass. (Default 0)
#' @param boxKernel <numeric> : If `NULL` automatically
#' compute for 3 Sd. (Default NULL)
#' @param kernSize <numeric> : Size of box applied to blurr.
#' If `NULL` automatically compute for 3 Sd. (Default NULL)
#' @param stdev <numeric> : SD of gaussian smooth. (Default 0.5)
#' @param loTri <numeric> : The value that replace all value in
#' the lower triangle of matrice (Usefull when blur is apply).(Default NULL)
#' @param colors <character> : Heatmap color list.
#' If `NULL`, automatically compute. (Default NULL)
#' @param na.value <character> : Color of NA values.
#' (Default "#F2F2F2")
#' @param colorScale <character> : Shape of color scale on of
#' "linear" or "density" based. (Default "linear")
#' @param bias <numeric> : A positive number. Higher values give
#' more widely spaced colors at the high end. See `?grDevices::colorRamp`
#' for more details. (Default 1)
#' @param paletteLength <numeric> : The number of color in the
#' palette. (Default 51)
#' @param annotate <logical> : Should there be axis ticks?
#' (Default: TRUE)
#' @param anchor.name <character> Name of anchor for
#' annotation. (Default "Anchor")
#' @param bait.name <character> Name of bait for
#' annotation. (Default "Bait")
#' @param fixCoord <logical> Fix axes coordinates?
#' (Default TRUE)
#' @return A ggplot object.
#' @importFrom checkmate assertMatrix
#' @export
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#'
#' # Index Beaf32
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRangeList = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSizes = data.frame(seqnames = c("2L", "2R"),
#'         seqlengths = c(23513712, 25286936)),
#'     binSize = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] 
#' # subset 2000 first for exemple
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_PF.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_Ctrl.cmx_lst,
#'     referencePoint = "pf"
#' )
#'
#' # Aggregate matrices in one matrix
#' aggreg.mtx <- Aggregation(interactions_PF.mtx_lst)
#'
#' # Visualization
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx
#' )
#'
#' # Add Title
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA"
#' )
#'
#' # Trim values
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA 30% trimmed on upper tail of distribution",
#'     trim = 30,
#'     tails = "upper"
#' )
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA 30% trimmed on lower tail of distribution",
#'     trim = 30,
#'     tails = "lower"
#' )
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA 15% trimmed on each tail of distribution",
#'     trim = 30,
#'     tails = "both"
#' )
#'
#' # Change Minimal, Central and Maximal Colors scale value
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA [min 200, center 300, max 600]",
#'     colMin = 200,
#'     colMid = 300,
#'     colMax = 600
#' )
#'
#' # Change Color
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA",
#'     colors = viridis(6),
#'     na.value = "black"
#' )
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA",
#'     colors = c("black", "white"),
#' )
#'
#' # Change Color distribution
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA [100,150,200,250,300,350,600]",
#'     colBreaks = c(100, 150, 200, 250, 300, 350, 600) # Choosen Breaks
#' )
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA",
#'     colorScale = "density" # color distribution based on density
#' )
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA",
#'     bias = 2 # (>1 wait on extremums)
#' )
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA",
#'     bias = 0.5 # (<1 wait on center)
#' )
#'
#' # Apply a Blurr
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA",
#'     blurPass = 1,
#'     stdev = 0.5
#' )
#'
#' # ggplot2 object modifications
#' # Since the function returns a ggplot object, it is possible
#' # to modify it following the ggplot2 grammar.
#' ggAPA(
#'     aggregatedMtx = aggreg.mtx,
#'     title = "APA",
#' ) +
#'     ggplot2::labs(
#'         title = "New title",
#'         subtitle = "and subtitle"
#'     )
ggAPA <- function(
    aggregatedMtx = NULL, title = NULL,
    trim = 0, tails = "both",
    colMin = NULL, colMid = NULL,
    colMax = NULL, colBreaks = NULL,
    blurPass = 0, boxKernel = NULL,
    kernSize = NULL,
    stdev = 0.5, loTri = NULL,
    colors = NULL,
    na.value = "#F2F2F2",
    colorScale = "linear",
    bias = 1, paletteLength = 51,
    annotate = TRUE,
    anchor.name = "Anchor",
    bait.name = "Bait",
    fixCoord = TRUE
) {
    checkmate::assertMatrix(
        x = aggregatedMtx,
        null.ok = FALSE
    )
    # Trimming
    if (!is.null(colBreaks)) {
        colMin <- min(colBreaks)
        colMax <- max(colBreaks)
    }
    vec.num <- c(aggregatedMtx)
    if (is.null(trim)) {
        trim <- 0
    }
    if (trim != 0 ||
        !is.null(colMin) ||
        !is.null(colMax)
    ) {
        bounds.num_vec <- vec.num |>
            QtlThreshold(
                prctThr = trim,
                tails = tails
            ) |>
            stats::setNames(NULL)
        bounds.num_lst <- list(
            bounds.num_vec,
            list(
                colMin,
                colMax
            )
        )
        bounds.num_lst <- TransposeList(bounds.num_lst)
        bounds.num_vec <- c(
            max(
                unlist(bounds.num_lst[1]),
                na.rm = TRUE
            ),
            min(
                unlist(bounds.num_lst[2]),
                na.rm = TRUE
            )
        )
    } else {
        bounds.num_vec <- NULL
    }
    if (!is.null(bounds.num_vec)) {
        aggregatedMtx <- TrimOutliers(
            x = aggregatedMtx,
            thr = bounds.num_vec,
            clip = TRUE
        )
        vec.num <- c(aggregatedMtx)
    }
    # Smoothing
    if (blurPass) {
        for (i in seq_len(blurPass)) {
            aggregatedMtx <- BoxBlur(
                mtx = aggregatedMtx,
                stdev = stdev,
                boxKernel = boxKernel,
                kernSize = kernSize
            )|> `attributes<-`(attributes(aggregatedMtx))
        }
        if (!is.null(loTri)) {
            aggregatedMtx[lower.tri(
                aggregatedMtx,
                diag = FALSE
            )] <- loTri
        }
        vec.num <- c(aggregatedMtx)
    }
    # Breaks
    if (is.null(colBreaks)) {
        colBreaks <- BreakVector(
            x = vec.num,
            x_min = colMin,
            center = colMid,
            x_max = colMax,
            n = paletteLength,
            method = colorScale
        )
        colMin <- min(colBreaks)
        colMax <- max(colBreaks)
    }
    # Colors
    if (is.null(colors)) {
        colors <- dplyr::case_when(
            !is.null(colMid) && max(colBreaks) <= colMid ~
                rev(YlGnBu(
                    paletteLength = paletteLength,
                    bias = bias
                )),
            !is.null(colMid) && colMid <= min(colBreaks) ~
                YlOrRd(
                    paletteLength = paletteLength,
                    bias = bias
                ),
            TRUE ~
                c(
                    rev(YlGnBu(
                        paletteLength = floor((paletteLength -1)/2),
                        bias = bias
                    )),
                    "#FFFFD8",
                    YlOrRd(
                        paletteLength = ceiling((paletteLength -1)/2),
                        bias = bias
                    )
                )
        )
    }
    # Raster
    data.dtf <- MeltSpm(aggregatedMtx)
    plot.ggp <- ggplot2::ggplot(
        data.dtf, ggplot2::aes(
            .data$j,
            .data$i
        )
    ) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$x)) +
    ggplot2::scale_fill_gradientn(
        colours  = colors,
        values   = MinMaxScale(colBreaks),
        na.value = na.value,
        limits   = c(
            colMin,
            colMax
        )
    ) +
    ggplot2::labs(
        title = title,
        y = dimnames(aggregatedMtx)[[2]],
        x = dimnames(aggregatedMtx)[[2]]
    ) +
    ggplot2::theme_classic()
    if(annotate){
        if(attributes(aggregatedMtx)$referencePoint == "rf"){
            regionSize <- attributes(
                aggregatedMtx)$matriceDim/(
                    (attributes(aggregatedMtx)$shift*2)+1)
            breaks_y <- c(1,floor(
                (attributes(aggregatedMtx)$matriceDim-regionSize)/2)+1,
                    floor(
                (attributes(aggregatedMtx)$matriceDim+regionSize)/2),
                attributes(aggregatedMtx)$matriceDim)
            breaks_x <- breaks_y
            labels_y <- c(paste0("-",
                    attributes(aggregatedMtx)$shift, " x ROI\nwidth"),
                    anchor.name, bait.name, paste0("+",
                    attributes(aggregatedMtx)$shift, " x ROI\nwidth"))
            labels_x <- labels_y
        } else {
           if(is.null(colnames(aggregatedMtx))){
                extent <- (floor(
                    nrow(aggregatedMtx)/2) * attributes(
                        aggregatedMtx)$resolution)/1000
                breaks_y <- c(1,ceiling(nrow(aggregatedMtx)/2),
                    nrow(aggregatedMtx))
                breaks_x <- c(1,ceiling(ncol(aggregatedMtx)/2),
                    ncol(aggregatedMtx))
                labels_y <- c(paste0("-",
                    extent, "KB"), anchor.name, paste0("+",
                    extent, "KB"))
                labels_x <- c(paste0("-",
                    extent, "KB"), bait.name, paste0("+",
                    extent, "KB"))
            } else{
                breaks_y <- seq_along(rownames(aggregatedMtx))
                labels_y <- rownames(aggregatedMtx)
                breaks_x <- seq_along(colnames(aggregatedMtx))
                labels_x <- colnames(aggregatedMtx)
            }
        }
        plot.ggp <- plot.ggp+
            ggplot2::scale_y_reverse(
                breaks = breaks_y,
                labels = labels_y
            ) +
            ggplot2::scale_x_continuous(
                breaks = breaks_x,
                labels = labels_x
            ) +
            ggplot2::theme(
                axis.line.y  = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_line(colour="black"),
                axis.line.x  = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(angle=45,hjust=1),
                axis.ticks.x = ggplot2::element_line(colour="black"),
                legend.title = ggplot2::element_blank()
            )

    }else{
        breaks_y <- seq_along(rownames(aggregatedMtx))
        labels_y <- rownames(aggregatedMtx)
        breaks_x <- seq_along(colnames(aggregatedMtx))
        labels_x <- colnames(aggregatedMtx)
        plot.ggp <- plot.ggp+
            ggplot2::scale_y_reverse(
                breaks = breaks_y,
                labels = labels_y
            ) +
            ggplot2::scale_x_continuous(
                breaks = breaks_x,
                labels = labels_x
            ) +
            ggplot2::theme(
                axis.line.y  = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.line.x  = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                legend.title = ggplot2::element_blank()
            )
    }
    if(fixCoord) plot.ggp <- plot.ggp+ggplot2::coord_fixed(ratio = 1)
    return(plot.ggp)
}
