#' .GInteractionFormatting
#' Correcting the genomic features for which to extract
#' hic contacts. adding distance, anchor.bin, bait.bin, submatrix name
#' 
#' @keywords internal
#' @param genomicFeature GRanges, GInteractions object to format
#' @param hicResolution resolution
#'
#' @return returns A GInteraction object
#' @noRd
.GInteractionFormatting <- function(
        genomicFeature, hicResolution
    ) {
    if (inherits(genomicFeature, "GRanges")) {
        feature.gni <- InteractionSet::GInteractions(
            GenomicRanges::resize(genomicFeature, 1, "start"),
            GenomicRanges::resize(genomicFeature, 1, "end")
        )
    } else if (inherits(genomicFeature, "Pairs") &&
        inherits(genomicFeature@first, "GRanges") &&
        inherits(genomicFeature@second, "GRanges")) {
        feature.gni <- InteractionSet::GInteractions(
            genomicFeature@first,
            genomicFeature@second
        )
    } else if (inherits(genomicFeature, "GInteractions")) {
        feature.gni <- genomicFeature
    } else {
        stop("genomicFeature object is not a GRanges, Pairs
        or GInteractions object!")
    }
    S4Vectors::mcols(feature.gni) <- S4Vectors::mcols(genomicFeature)
    if (is.null(GenomeInfoDb::seqinfo(feature.gni))) {
        GenomeInfoDb::seqinfo(feature.gni) <-
            GenomeInfoDb::seqinfo(genomicFeature)
    }
    if (is.null(S4Vectors::mcols(feature.gni)$distance)) {
        S4Vectors::mcols(feature.gni)$distance <-
        InteractionSet::pairdist(feature.gni)
    }
    if (is.null(S4Vectors::mcols(feature.gni)$orientation)) {
        S4Vectors::mcols(feature.gni)$orientation <-
            (feature.gni == InteractionSet::swapAnchors(feature.gni))
    }
    if (is.null(S4Vectors::mcols(feature.gni)$anchor.bin)) {
        S4Vectors::mcols(feature.gni)$anchor.bin <- paste0(
            GenomeInfoDb::seqnames(
                InteractionSet::anchors(feature.gni)$first
            ), ":",
            ceiling(
                InteractionSet::anchors(feature.gni)$first@ranges@start/
                hicResolution
            )
        )
    }
    if (is.null(S4Vectors::mcols(feature.gni)$bait.bin)) {
        S4Vectors::mcols(feature.gni)$bait.bin <- paste0(
            GenomeInfoDb::seqnames(
                InteractionSet::anchors(feature.gni)$second
            ), ":",
            ceiling(
                InteractionSet::anchors(feature.gni)$second@ranges@start/
                hicResolution
            )
        )
    }
    if (is.null(S4Vectors::mcols(feature.gni)$name)) {
        S4Vectors::mcols(feature.gni)$name <- paste0(
            S4Vectors::mcols(feature.gni)$anchor.bin, "_",
            S4Vectors::mcols(feature.gni)$bait.bin
        )
    }
    if (is.null(S4Vectors::mcols(feature.gni)$submatrix.name)) {
        S4Vectors::mcols(feature.gni)$submatrix.name <- paste0(
            S4Vectors::mcols(feature.gni)$anchor.bin, "_",
            S4Vectors::mcols(feature.gni)$bait.bin
        )
        S4Vectors::mcols(feature.gni)$submatrix.name[
            !S4Vectors::mcols(feature.gni)$orientation
            ] <- paste0(
                S4Vectors::mcols(feature.gni)$bait.bin[
                    !S4Vectors::mcols(feature.gni)$orientation
                ], "_",
                S4Vectors::mcols(feature.gni)$anchor.bin[
                    !S4Vectors::mcols(feature.gni)$orientation
                ]
            )
    }
    if (!sum(
        S4Vectors::mcols(feature.gni)$anchor.bin !=
        S4Vectors::mcols(feature.gni)$bait.bin
    )) {
        S4Vectors::mcols(feature.gni)$bait.bin <- NULL
        S4Vectors::mcols(feature.gni)$bin <-
            S4Vectors::mcols(feature.gni)$anchor.bin
        S4Vectors::mcols(feature.gni)$anchor.bin <- NULL
    }
    return(feature.gni)
}

#' Submatrix extraction.
#'
#' ExtractSubmatrix
#' @description Extract matrices in the HiC maps list around genomic features.
#' @param genomicFeature <GRanges or Pairs[GRanges] or GInteractions>:
#'  The genomic coordinates on which compute the extraction of HiC submatrix.
#' @param hicLst <List[ContactMatrix][InteractionSet::ContactMatrix()]>:
#'  The HiC maps list.
#' @param referencePoint <character>: Type of
#' extracted submatrices. "rf" for "region feature" to extract triangle-shaped
#' matrices around regions or "pf" for "point feature" to extract
#' square-shaped matrices around points. (Default "rf")
#' @param hicResolution <numeric>: The resolution
#' used in hicLst. If `NULL`, automatically find in resolution attributes
#' of hicLst. (Default NULL)
#' @param matriceDim <numeric>: The size of matrices in output.
#' (Default 21).
#' @param shift <numeric>: Only when "referencePoint" is "rf".
#' Factor defining how much of the distance between anchor and bait is
#' extracted before and after the region (Default 1).
#' Ex: for shift=2, extracted matrices will be: 
#' `2*regionSize+regionSize+2*regionSize`.
#' @param remove_duplicates <Logical>: 
#' remove duplicated submatrices ?
#' This avoids duplicated submatrices when both anchor and bait bins are from
#'  the same feature. ex. BEAF32-BEAF32, same submatrix twice with opposite
#'  orientations(Default TRUE)
#' @param cores <integer> : An integer
#' to specify the number of cores. (Default 1)
#' @param verbose <logical>: If TRUE, 
#' show the progression in console. (Default FALSE)
#' @return A matrices list.
#' @importFrom checkmate assertLogical assertChoice assert checkNumeric
#' @importFrom S4Vectors mcols runValue runLength metadata
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
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#'
#' # Matrices extractions of regions defined between
#' # Beaf32 <-> Beaf32 interactions
#' interactions_RF.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_Ctrl.cmx_lst,
#'     referencePoint = "rf"
#' )
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 pointinteraction
#' interactions_PF.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_Ctrl.cmx_lst,
#'     referencePoint = "pf"
#' )
#'
#'


ExtractSubmatrix <- function(
    genomicFeature = NULL, hicLst = NULL, referencePoint = "pf",
    hicResolution = NULL, matriceDim = 21, shift = 1, 
    remove_duplicates = TRUE, cores = 1,verbose = FALSE
) {
    .validHicMatrices(matrices = hicLst)
    checkmate::assertLogical(
        x = remove_duplicates
    )
    checkmate::assertChoice(
        x = referencePoint,
        choices = c('rf','pf'),
        null.ok = FALSE)
    checkmate::assert(
        checkmate::checkNumeric(
            matriceDim, lower = 5,
            any.missing = FALSE),
        checkmate::checkNumeric(
            shift, lower = 1,
            any.missing = FALSE),
        checkmate::checkNumeric(
            cores, lower = 1,
            any.missing = FALSE),
        combine = "and")

    # Run Check Resolution
    if (!is.null(attr(hicLst, "resolution"))) {
        hicResolution <- attr(hicLst, "resolution")
    } else if (is.character(hicResolution)) {
        hicResolution <- GenomicSystem(hicResolution)
    }
    checkmate::assertNumeric(
        x = hicResolution,
        lower = 1, null.ok = FALSE,
        .var.name = "hicResolution"
    )
    # Formatting
    genomicFeature <- .GInteractionFormatting(
        genomicFeature = genomicFeature,
        hicResolution = hicResolution
    )
    if (is.null(GenomeInfoDb::seqinfo(genomicFeature)) ||
    any(is.na(GenomeInfoDb::seqlengths(
        GenomeInfoDb::seqinfo(genomicFeature))))) {
        GenomeInfoDb::seqlevels(genomicFeature) <-
            as.character(attributes(attributes(hicLst)$chromSize[[1]]))
        GenomeInfoDb::seqinfo(genomicFeature) <-
            GenomeInfoDb::Seqinfo(
                seqnames = attributes(attributes(hicLst)$chromSize[[1]]),
                seqlengths = attributes(attributes(hicLst)$chromSize[[2]]))
    }
    if (!sum(genomicFeature$anchor.bin != genomicFeature$bait.bin)) {
        referencePoint <- "pf"
    }
    # Resize Features according reference Points
    referencePoint <- tolower(referencePoint)
    if (referencePoint == "rf") {
        cis.lgk <- ReduceRun(
            GenomeInfoDb::seqnames(
                InteractionSet::anchors(genomicFeature)$first),
            GenomeInfoDb::seqnames(
                InteractionSet::anchors(genomicFeature)$second),
            reduceMethod = "paste", sep = "_"
            ) |>
            as.character() |>
            lapply(function(combinaison) {
                split <- unlist(strsplit(combinaison, "_"))
                return(split[1] == split[2])
            }) |>
            unlist()
        genomicFeature <- genomicFeature[cis.lgk]
        genomicFeature <- genomicFeature[
            which(genomicFeature$distance >= (3 * hicResolution))]
        ranges.lst_dtf <- lapply(
            c("first", "second"),
            function(anchorName.chr) {
                anchor.gnr <- InteractionSet::anchors(genomicFeature)[[
                    anchorName.chr
                ]]
                anchor.dtf <- data.frame(IRanges::ranges(anchor.gnr)) |>
                    dplyr::select("start", "end")
                colnames(anchor.dtf) <- paste0(
                    anchorName.chr, ".",
                    names(anchor.dtf)
                )
                return(anchor.dtf)
            }
        )
        ranges.dtf <- do.call(cbind, ranges.lst_dtf)
        ranges.dtf <- dplyr::mutate(
            ranges.dtf,
            start = pmin(ranges.dtf$first.start, ranges.dtf$second.start)
        )
        ranges.dtf <- dplyr::mutate(
            ranges.dtf,
            end = pmax(ranges.dtf$first.end, ranges.dtf$second.end)
        )
        feature.gnr <- GenomicRanges::GRanges(
            seqnames = GenomeInfoDb::seqnames(
                InteractionSet::anchors(genomicFeature)$first
            ),
            ranges = IRanges::IRanges(
                start = ranges.dtf$start,
                end = ranges.dtf$end
            ),
            seqlengths = GenomeInfoDb::seqlengths(genomicFeature),
            seqinfo = GenomeInfoDb::seqinfo(genomicFeature)
        )
        featureResize.gnr <- GenomicRanges::resize(
            feature.gnr,
            width = feature.gnr@ranges@width +
                genomicFeature$distance * shift * 2,
            fix = "center"
        )
        featureResize.gni <- InteractionSet::GInteractions(
            featureResize.gnr,
            featureResize.gnr
        )
        S4Vectors::mcols(featureResize.gni) <- S4Vectors::mcols(genomicFeature)
    } else if (referencePoint == "pf") {
        featureResize.gni <- suppressWarnings(
            GenomicRanges::resize(
                genomicFeature, width = hicResolution * (matriceDim - 1) + 1,
                fix = "center"
            )
        )
    }
    # Filt Out Of Bound
    featureFilt.gni <- featureResize.gni[which(
        1L <=
        data.frame(
            InteractionSet::anchors(featureResize.gni)$first@ranges
        )[,"start"] &
        data.frame(
            InteractionSet::anchors(featureResize.gni)$first@ranges
        )[,"end"] <=
        SeqEnds(InteractionSet::anchors(featureResize.gni)$first) &
        1L <=
        data.frame(
            InteractionSet::anchors(featureResize.gni)$second
        )[,"start"] &
        data.frame(
            InteractionSet::anchors(featureResize.gni)$second
        )[,"end"] <=
        SeqEnds(InteractionSet::anchors(featureResize.gni)$second)
    )]
    ## Filt Duplicated Submatrix before extraction
    featureNoDup.gni <-
        featureFilt.gni[!duplicated(featureFilt.gni$submatrix.name)]
    ## Order according Chromosomes combinaison
    ## ReduceRun takes longer than just paste
    chromosomesCombinaison.rle <- S4Vectors::Rle(paste(
        GenomeInfoDb::seqnames(
            InteractionSet::anchors(featureNoDup.gni)$first
        ),
        GenomeInfoDb::seqnames(
            InteractionSet::anchors(featureNoDup.gni)$second
        ),sep = "_"
    ))
    order.num <- unlist(
        lapply(
            unique(S4Vectors::runValue(chromosomesCombinaison.rle)),
            function(comb) {
                which(as.character(chromosomesCombinaison.rle) == comb)
            }
        )
    )
    featureNoDup.gni <- featureNoDup.gni[order.num]
    chromosomesCombinaison.rle <- chromosomesCombinaison.rle[order.num]
    matAnchors.gnr_lst <- lapply(hicLst, InteractionSet::anchors)
    # Keep combinaison that are present in hicLst
    present.ndx <- which(
        as.vector(chromosomesCombinaison.rle) %in% names(matAnchors.gnr_lst) |
            as.vector(chromosomesCombinaison.rle) %in% unlist(
                lapply(
                    names(matAnchors.gnr_lst),
                    function(name.chr) {
                        strsplit(name.chr, "_") |>
                        unlist() |>
                        rev() |>
                        paste(collapse = "_")
                    }
                )
            )
    )
    chromosomesCombinaison.rle <- chromosomesCombinaison.rle[present.ndx]
    featureNoDup.gni <- featureNoDup.gni[present.ndx]
    # Separate anchors
    anchors.gnr <- InteractionSet::anchors(featureNoDup.gni)$first
    baits.gnr <- InteractionSet::anchors(featureNoDup.gni)$second
    # Extraction
    multicoreParam <- MakeParallelParam(cores = cores,verbose = FALSE)
    submatrix.spm_lst <- BiocParallel::bplapply(
        BPPARAM = multicoreParam,
        seq_along(S4Vectors::runValue(chromosomesCombinaison.rle)),
        function(combinaison.ndx) {
            combinaisonName.chr <-
                S4Vectors::runValue(
                    chromosomesCombinaison.rle
                )[[combinaison.ndx]]
            combinaisonStart.ndx <- cumsum(c(1, S4Vectors::runLength(
                chromosomesCombinaison.rle)))[[combinaison.ndx]]
            combinaisonEnd.ndx <-cumsum(S4Vectors::runLength(
                chromosomesCombinaison.rle))[[combinaison.ndx]]
            if (combinaisonName.chr %in% names(matAnchors.gnr_lst)) {
                mat.ndx <- which(names(hicLst) == combinaisonName.chr)
                ovl_row <- data.frame(GenomicRanges::findOverlaps(
                    anchors.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],
                    matAnchors.gnr_lst[[mat.ndx]]$row
                ))
                ovl_row <- dplyr::group_by(
                    ovl_row,
                    queryHits = ovl_row$queryHits
                )
                ovl_row <- tidyr::nest(ovl_row)
                ovl_col <- data.frame(GenomicRanges::findOverlaps(
                    baits.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],
                    matAnchors.gnr_lst[[mat.ndx]]$col
                ))
                ovl_col <- dplyr::group_by(
                    ovl_col,
                    queryHits = ovl_col$queryHits
                )
                ovl_col <- tidyr::nest(ovl_col)
            } else if ({
                combinaisonName.chr |>
                strsplit("_") |>
                unlist() |>
                rev() |>
                paste(collapse = "_")} %in% names(matAnchors.gnr_lst)) {
                    mat.ndx <- which(
                        names(hicLst) == {
                            combinaisonName.chr |>
                            strsplit("_") |>
                            unlist() |>
                            rev() |>
                            paste(collapse = "_")
                        }
                    )
                    ovl_row <- data.frame(GenomicRanges::findOverlaps(
                        baits.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],
                        matAnchors.gnr_lst[[mat.ndx]]$row
                    ))
                    ovl_row <- dplyr::group_by(
                        ovl_row,
                        queryHits = ovl_row$queryHits
                    )
                    ovl_row <- tidyr::nest(ovl_row)
                    ovl_col <- data.frame(GenomicRanges::findOverlaps(
                        anchors.gnr[combinaisonStart.ndx:combinaisonEnd.ndx],
                        matAnchors.gnr_lst[[mat.ndx]]$col
                    ))
                    ovl_col <- dplyr::group_by(
                        ovl_col,
                        queryHits = ovl_col$queryHits
                    )
                    ovl_col <- tidyr::nest(ovl_col)
            }
            tempSubmatrix.spm_lst <- lapply(
                seq_along(combinaisonStart.ndx:combinaisonEnd.ndx),
                function(range.ndx) {
                    row.ndx <- unlist(
                        ovl_row[[range.ndx, "data"]],
                        use.names = FALSE
                    )
                    col.ndx <- unlist(
                        ovl_col[[range.ndx, "data"]],
                        use.names = FALSE
                    )
                    if (S4Vectors::metadata(hicLst[[mat.ndx]])$type ==
                        "cis"
                    ) {
                        gap.num <-
                            stats::median(col.ndx) - stats::median(row.ndx)
                    } else {
                        gap.num <- Inf
                    }
                    if (gap.num < 0) {
                        spMtx <-
                            hicLst[[mat.ndx]][col.ndx, row.ndx]@matrix
                    } else {
                        spMtx <-
                            hicLst[[mat.ndx]][row.ndx, col.ndx]@matrix
                    }
                    if (dim(spMtx)[1] != matriceDim) {
                        spMtx <- ResizeMatrix(
                            mtx = spMtx,
                            newDim = c(matriceDim, matriceDim)
                        )
                    }
                    if (abs(gap.num) < matriceDim) {
                        if (abs(gap.num) > 0) {
                            spMtx <- PadMtx(
                                mtx = spMtx, padSize = abs(gap.num),
                                val = 0, side = c("left", "bot")
                            )
                        }
                        spMtx[lower.tri(spMtx)] <- NA
                        if (abs(gap.num) > 0) {
                            spMtx <- spMtx[
                                seq_len(matriceDim),
                                (abs(gap.num)+1):(matriceDim+abs(gap.num))
                            ]
                        }
                    }
                    return(as.matrix(spMtx))
                }
            )
            if(verbose){
                message("extracted : ", length(tempSubmatrix.spm_lst), 
            " Matrices on ", combinaisonName.chr)
            }
            
            return(tempSubmatrix.spm_lst)
        }
    ) |>
        do.call(what = c) |>
        stats::setNames(featureNoDup.gni$submatrix.name)
    if(verbose){
        message(" FILTERED extracted submat: ", length(submatrix.spm_lst))
    }
    
    if(!remove_duplicates){
        submatrix.spm_lst <- submatrix.spm_lst[featureFilt.gni$submatrix.name]
        if(verbose){
            message("You have chosen to keep duplicated submatrices.")
        }
    }
    if(verbose){
        message(" TOTAL extracted submat: ", length(submatrix.spm_lst))
    }
    
    interactions.ndx <- seq_along(genomicFeature$name) |>
        stats::setNames(genomicFeature$name)
    if(remove_duplicates){
        interactions.ndx <- interactions.ndx[featureNoDup.gni$name]
    }else{
        interactions.ndx <- interactions.ndx[featureFilt.gni$name]
    }
    
    attributes(submatrix.spm_lst)$interactions <- genomicFeature[
        interactions.ndx]
    attributes(submatrix.spm_lst)$resolution <- hicResolution
    attributes(submatrix.spm_lst)$referencePoint <- referencePoint
    attributes(submatrix.spm_lst)$matriceDim <- matriceDim
    if (referencePoint == "rf") {
        attributes(submatrix.spm_lst)$shift <- shift
    }
    return(submatrix.spm_lst)
}
