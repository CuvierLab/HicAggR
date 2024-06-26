#' Import Hic data
#'
#' ImportHiC
#' @description Import ..hic, .cool, .mcool or .bedpe data
#' @param file <GRanges or Pairs[GRanges] or GInteractions>:
#'  The genomic feature on which compute the extraction of HiC submatrix.
#'  Extension should be .hic, .cool, .mcool, .h5, .hdf5, .HDF5 or .bedpe"
#'  assuming .h5 and .hdf5 are only for cool (not mcool).
#' @param hicResolution <numeric>: The HiC resolution.
#' @param chromSizes <data.frame>: A data.frame where first
#' colum correspond to the chromosomes names, and the second column correspond
#' to the chromosomes lengths in base pairs.
#' @param chrom_1 <numeric>: The seqnames of first chromosmes
#' (rows in matrix).
#' @param chrom_2 <numeric>: The seqnames of second chromosmes
#' (col in matrix).
#' If `NULL` variable will be assigned value of chrom_1 (Defalt NULL).
#' @param cores <numerical> : An integer to specify the number
#' of cores. (Default 1)
#' @param verbose <logical>: Show the progression in console?
#'  (Default FALSE)
#' @param hic_norm <character>: "norm" argument to supply to
#' [strawr::straw()].
#'  This argument is for .hic format data only. Available norms can be obtained
#'  through [strawr::readHicNormTypes()].
#'  (Default "NONE").
#' @param hic_matrix <character>: "matrix" argument to supply to
#' [strawr::straw()].
#'  This argument is for .hic format data only.
#'  Other options can be: "oe", "expected". (Default "observed").
#' @param cool_balanced <logical> Import already balanced matrix?
#'  (Default: FALSE)
#' @param cool_weight_name <character> Name of the correcter in
#' the cool file. (Default: weight).
#' [rhdf5::h5ls()] to see the available correctors.
#' @param cool_divisive_weights <logical> Does the correcter
#' vector contain divisive biases as in hicExplorer or multiplicative as in
#' cooltools? (Default: FALSE)
#' @param h5_fill_upper <logical> Do the matrix in h5 format
#' need to be transposed? (Default: TRUE)
#' 
#' @details If you request "expected" values when importing .hic format data,
#' you must do yourself the "oe" by importing manually the observed counts
#' as well.
#' 
#' Prior to v.0.9.0 cooltools had multiplicative weight only, so make sure
#' your correcters are divisive or multiplicative.
#' https://cooler.readthedocs.io/en/stable/releasenotes.html#v0-9-0
#'  
#' When loading hic matrix in h5 format make sure you have enough momory
#' to load the full matrix with all chromosomes regardless of values for
#' chrom_1 and chrom_2 arguments. The function first loads the whole matrix,
#' then extracts matrices per chromosome for the time being, it's easier ;).
#' 
#' @return A matrices list.
#' @export
#' @importFrom checkmate checkChoice assertCharacter assertFileExists
#' @importFrom withr local_options
#' @examples
#' \donttest{
#'
#' # Prepare Temp Directory
#' options(timeout = 3600)
#' temp.dir <- file.path(tempdir(), "HIC_DATA")
#' dir.create(temp.dir)
#'
#' # Download .hic file
#' Hic.url <- paste0(
#'     "https://4dn-open-data-public.s3.amazonaws.com/",
#'     "fourfront-webprod/wfoutput/",
#'     "7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic"
#' )
#' HicOutput.pth <- file.path(temp.dir, "Control_HIC.hic")
#' HicOutput.pth <- normalizePath(HicOutput.pth)
#' if(.Platform$OS.type == "windows"){
#'     download.file(Hic.url, HicOutput.pth, method = "auto",
#'     extra = "-k",mode="wb")
#' }else{
#'     download.file(Hic.url, HicOutput.pth, method = "auto", extra = "-k")
#' }
#'
#' # Import .hic file
#' HiC_Ctrl.cmx_lst <- ImportHiC(
#'     file = HicOutput.pth,
#'     hicResolution = 100000,
#'     chrom_1 = c("2L", "2L", "2R"),
#'     chrom_2 = c("2L", "2R", "2R")
#' )
#'
#' # Download .mcool file
#' Mcool.url <- paste0(
#'     "https://4dn-open-data-public.s3.amazonaws.com/",
#'     "fourfront-webprod/wfoutput/",
#'     "4f1479a2-4226-4163-ba99-837f2c8f4ac0/4DNFI8DRD739.mcool"
#' )
#' McoolOutput.pth <- file.path(temp.dir, "HeatShock_HIC.mcool")
#' HicOutput.pth <- normalizePath(McoolOutput.pth)
#' if(.Platform$OS.type == "windows"){
#'     download.file(Mcool.url, McoolOutput.pth, method = "auto",
#'     extra = "-k",mode="wb")
#' }else{
#'     download.file(Mcool.url, McoolOutput.pth, method = "auto",
#'     extra = "-k")
#' }
#'
#' # Import .mcool file
#' HiC_HS.cmx_lst <- ImportHiC(
#'     file = McoolOutput.pth,
#'     hicResolution = 100000,
#'     chrom_1 = c("2L", "2L", "2R"),
#'     chrom_2 = c("2L", "2R", "2R")
#' )
#' }
#' # Import .h5 file
#' h5_path <- system.file("extdata",
#'     "Control_HIC_10k_2L.h5",
#'     package = "HicAggR", mustWork = TRUE
#' )
#' binSize=10000
#' hicLst <- ImportHiC(
#'   file      = h5_path,
#'   hicResolution       = binSize,
#'   chromSizes = data.frame(seqnames = c("2L"), 
#'   seqlengths = c(23513712)),
#'   chrom_1   = c("2L")
#' )

ImportHiC <- function(
    file = NULL, hicResolution = NULL, chromSizes = NULL, chrom_1 = NULL,
    chrom_2 = NULL, verbose = FALSE, cores = 1,
    hic_norm="NONE", hic_matrix = "observed", cool_balanced = FALSE,
    cool_weight_name = "weight", cool_divisive_weights= FALSE,
    h5_fill_upper = TRUE
) {
    # Resolution Format
    withr::local_options(list(scipen = 999))

    checkmate::assertFileExists(
        x = file,
        access = "r",
        .var.name = "Hic_file")
    if (inherits(hicResolution, "character")) {
        hicResolution <- GenomicSystem(hicResolution)
    }
    # Chromosomes Format
    if (is.null(chrom_2)) {
        chrom_2 <- chrom_1
    } else if (length(chrom_1) != length(chrom_2)) {
        stop("chrom_1 and chrom_2 must have the same length")
    }
    
    if (GetFileExtension(file) == "hic") {
        if(is.null(chromSizes)){
            chromSizes <- strawr::readHicChroms(file)
        }
        if("index" %in% colnames(chromSizes)){
            chromSizes <- chromSizes |>
                dplyr::select(-"index")}
        colnames(chromSizes) <- c("name", "length")
    } else if (GetFileExtension(file) %in%
        c("cool", "mcool", "HDF5", "hdf5")
    ) {
        # Define HDF5groups
        chr.group <- ifelse(
            GetFileExtension(file) %in% c("cool", "HDF5", "hdf5"),
            yes = "/chroms",
            no = paste("resolutions", hicResolution, "chroms", sep = "/")
        )
        # Get SeqInfo
        chromSizes <- data.frame(rhdf5::h5read(file, name = chr.group))
    } else if ((GetFileExtension(file) == "h5")){
        # Get SeqInfo
        chromSizes <- data.frame(rhdf5::h5read(file, name = "intervals")) |>
            dplyr::group_by(.data$chr_list) |> 
            dplyr::summarise("length" = max(.data$end_list))|>
            dplyr::rename("name"="chr_list")
        i.group <- "/matrix/indices"
        p.group <- "/matrix/indptr"
        x.group <- "/matrix/data"
        ## got the structure from these
        ## https://github.com/deeptools/HiCMatrix/blob/master/
        ## hicmatrix/lib/h5.py lines 38:39
        ## https://docs.scipy.org/doc/scipy/reference/
        ## generated/scipy.sparse.csc_matrix.html
        ## https://github.com/rstudio/reticulate/blob/main/R/conversion.R
        ## lines 499-507
        ## indices in csr_matrix correspond to i, indptr to p and data to x
        hic_spm_full_h5 <- Matrix::sparseMatrix(
            i=as.vector(rhdf5::h5read(
            file,
            name = i.group
            )),
            p=as.vector(rhdf5::h5read(
            file,
            name = p.group
            )),
            x = as.vector(rhdf5::h5read(
            file,
            name = x.group
            )), dims = c(sum(ceiling(chromSizes$length/hicResolution)),
                sum(ceiling(chromSizes$length/hicResolution))),
            index1 = FALSE, repr = "C"
        )
        if(h5_fill_upper){
            hic_spm_full_h5 <- Matrix::t(hic_spm_full_h5)
        }
    } else if (GetFileExtension(file) == "bedpe" &&
        !is.null(chromSizes)
    ) {
        colnames(chromSizes) <- c("name", "length")
        hic.gnp <- rtracklayer::import(file, format = "bedpe")
        megaHic.dtf <- data.frame(
            chrom_1 = as.vector(hic.gnp@first@seqnames),
            i = ceiling(hic.gnp@first@ranges@start/hicResolution),
            chrom_2 = as.vector(hic.gnp@second@seqnames),
            j = ceiling(hic.gnp@second@ranges@start/hicResolution),
            counts = hic.gnp@elementMetadata$score
        )
    } else if (GetFileExtension(file) == "bedpe" &&
        is.null(chromSizes)
    ) {
        stop("To import HiC data in bedpe format, non null chromSizes
            argument is required!")
    } else {
        stop("file must be .hic, .cool, .mcool, .hdf5, .HDF5 or .bedpe")
    }
    # to avoid warning if chromSizes is tibble, cause setting row.names
    # for tibble is deprecated
    chromSizes <- as.data.frame(chromSizes)
    rownames(chromSizes) <- chromSizes$name
    
    if(is.null(chrom_1) && is.null(chrom_2)){
        if(verbose){
            message("chrom_1 and chrom_2 are NULL,
            so all chromosomes are chosen")
        }        
        chrom_1 <- chromSizes$name
        chrom_2 <- chromSizes$name
    }
    if ("ALL" %in% toupper(chrom_1)){
        chrom_1 <- chrom_1[-which(toupper(chrom_1) == "ALL")]
        if(verbose){
            message("ALL removed from chrom_1")
        }
    }
    if ("ALL" %in% toupper(chrom_2)){
        chrom_2 <- chrom_2[-which(toupper(chrom_2) == "ALL")]
        if(verbose){
            message("ALL removed from chrom_2")
        }        
    }
    chrom.chr <- c(chrom_1, chrom_2) |>
        unlist() |>
        unique()
    if (grepl(pattern = "chr", chrom.chr[1], fixed = TRUE)) {
        seqlevelsStyleHiC <- "UCSC"
    } else {
        seqlevelsStyleHiC <- "ensembl"
    }
    if (grepl("chr", rownames(chromSizes)[
        length(rownames(chromSizes))],fixed = TRUE) &&
        seqlevelsStyleHiC == "ensembl"
    ) {
        chromSizes$name <- unlist(lapply(
            strsplit(rownames(chromSizes),"chr"),`[[`, 2
        ))
        rownames(chromSizes) <- unlist(lapply(
            strsplit(rownames(chromSizes),"chr"),`[[`, 2
        ))
    } else if (!grepl("chr", rownames(chromSizes)[
        length(rownames(chromSizes))], fixed = TRUE) &&
        seqlevelsStyleHiC == "UCSC"
    ) {
        chromSizes$name <- paste0("chr", rownames(chromSizes))
        rownames(chromSizes) <- paste0("chr", rownames(chromSizes))
    }
    chromSizes <- dplyr::mutate(
        chromSizes,
        dimension = ceiling(chromSizes$length/hicResolution)
    )
    chrom.chr <- chrom.chr[chrom.chr %in% chromSizes$name]
    chrom_1 <- chrom_1[chrom_1 %in% chromSizes$name]
    chrom_2 <- chrom_2[chrom_2 %in% chromSizes$name]
    chromSizes <- chromSizes[which(chromSizes$name%in%chrom.chr),]
    # Create Genome as GRanges
    binnedGenome.grn <- chromSizes |>
        dplyr::pull("length") |>
        stats::setNames(chromSizes$name) |>
        GenomicRanges::tileGenome(
            tilewidth = hicResolution,
            cut.last.tile.in.chrom = TRUE
        )
    GenomeInfoDb::seqlengths(binnedGenome.grn) <- chromSizes$length |>
        stats::setNames(chromSizes$name)
    binnedGenome.grn <- GenomeInfoDb::sortSeqlevels(binnedGenome.grn)
    
    chromComb.lst <- paste(chrom_1, chrom_2, sep = "_")
    matrixSymmetric.bln <- strsplit(chromComb.lst, "_") |>
        lapply(function(name.chr) {name.chr[[1]] == name.chr[[2]]}) |>
        unlist()
    matrixType.str <- chromComb.lst
    matrixType.str[which(matrixSymmetric.bln)] <- "cis"
    matrixType.str[which(!matrixSymmetric.bln)] <- "trans"
    matrixKind <- chromComb.lst
    matrixKind[which(matrixSymmetric.bln)] <- "U"
    matrixKind[which(!matrixSymmetric.bln)] <- NA
    attributes.tbl <- dplyr::bind_cols(
        name = chromComb.lst, type = matrixType.str, kind = matrixKind,
        symmetric = matrixSymmetric.bln
    )
    # Dump file
    multicoreParam <- MakeParallelParam(
        cores = cores,
        verbose = verbose
    )
    hic.lst_cmx <- BiocParallel::bplapply(
        BPPARAM = multicoreParam, seq_along(chromComb.lst),
        function(ele.ndx) {
            # Chromosomes
            ele.lst <- unlist(strsplit(chromComb.lst[[ele.ndx]], "_"))
            chrom_1 <- ele.lst[[1]]
            chrom_2 <- ele.lst[[2]]
            # Dimension
            dims.num <- ele.lst |>
                lapply(
                    function(chrom) {
                        dplyr::filter(
                            chromSizes,
                            chromSizes$name == chrom) |>
                        dplyr::pull("dimension")
                    }
                ) |>
                unlist()
            if (GetFileExtension(file) == "hic") {

                if(!checkmate::checkChoice(
                    x = hic_norm,
                    choices = c("NONE", "VC", "VCSQRT", "KR", "ICE")
                )){
                    warning("hic_norm has unusual value,
                        please make sure norm type exists
                        in your .hic file")
                }
                if(!checkmate::checkChoice(
                    x = hic_matrix,
                    choices = c("observed", "oe", "expected")
                )){
                    warning("hic_matrix has unusual value,
                        please make sure matrix type exists
                        in your .hic file")
                }
                # Read .hic file
                hic.dtf <- strawr::straw(
                    hic_norm,
                    file,
                    chrom_1,
                    chrom_2,
                    "BP",
                    hicResolution,
                    hic_matrix
                )
                hic.dtf$j <- ceiling((hic.dtf$y + 1)/hicResolution)
                hic.dtf$i <- ceiling((hic.dtf$x + 1)/hicResolution)
            } else if (GetFileExtension(file) %in%
                c("cool", "mcool", "HDF5", "hdf5")
            ) {
                # Define HDF5groups
                indexes.group <- ifelse(
                    GetFileExtension(file) %in%
                        c("cool", "HDF5", "hdf5"),
                    yes = "/indexes",
                    no = paste("resolutions",hicResolution,"indexes",sep = "/")
                )
                pixels.group <- ifelse(
                    GetFileExtension(file) %in%
                        c("cool", "HDF5", "hdf5"),
                    yes = "/pixels",
                    no = paste("resolutions", hicResolution,"pixels",sep = "/")
                )
                # Define start and end of chromosomes
                ends.ndx <- chromSizes$dimension |>
                    cumsum() |>
                    stats::setNames(chromSizes$name)
                starts.ndx <- 1 + c(0, ends.ndx[-length(ends.ndx)]) |>
                    stats::setNames(chromSizes$name)
                # Read .mcool file
                bin1.ndx <- as.vector(rhdf5::h5read(
                    file,
                    name = paste(indexes.group, "bin1_offset", sep = "/"),
                    index = list(starts.ndx[chrom_1]:ends.ndx[chrom_1])
                ))
                slice.num <- sum(
                    bin1.ndx[-1] - bin1.ndx[-length(bin1.ndx)]
                    ) - 1
                chunk.num <- seq(bin1.ndx[1] + 1, bin1.ndx[1] + 1 + slice.num)
                hic.dtf <- data.frame(
                    i = as.vector(rhdf5::h5read(
                        file,
                        name = paste(pixels.group, "bin1_id", sep = "/"),
                        index = list(chunk.num)
                    )) + 1,
                    j = as.vector(rhdf5::h5read(
                        file,
                        name = paste(pixels.group, "bin2_id", sep = "/"),
                        index = list(chunk.num)
                    )) + 1,
                    counts = as.vector(rhdf5::h5read(
                        file,
                        name = paste(pixels.group, "count", sep = "/"),
                        index = list(chunk.num)
                    ))
                )
                if(cool_balanced){
                    checkmate::assertCharacter(
                        x = cool_weight_name, null.ok = FALSE)
                    bins.group <- ifelse(
                        GetFileExtension(file) %in%
                            c("cool", "HDF5", "hdf5"),
                        yes = "/bins",
                        no = paste("resolutions",
                            hicResolution,"bins",sep = "/")
                    )
                    bins <- data.frame(weight=
                        rhdf5::h5read(file, 
                        name = paste(bins.group,
                            get("cool_weight_name"),sep="/"),
                        bit64conversion = 'double',
                        index = list(
                            starts.ndx[chrom_1]:ends.ndx[chrom_1]))) |>
                        dplyr::mutate(
                            bin_id=c(starts.ndx[chrom_1]:ends.ndx[chrom_1])-1)
                    bins2 <- data.frame(weight=
                        rhdf5::h5read(file, 
                        name = paste(bins.group,
                            get("cool_weight_name"),sep="/"),
                        bit64conversion = 'double',
                        index = list(
                            starts.ndx[chrom_2]:ends.ndx[chrom_2]))) |>
                        dplyr::mutate(
                            bin_id=c(starts.ndx[chrom_2]:ends.ndx[chrom_2])-1)
                    hic.dtf <- dplyr::left_join(
                        (dplyr::mutate(hic.dtf, "bin1_id" = hic.dtf$i-1)),
                        (dplyr::rename(bins, "weight1" = "weight")),
                        dplyr::join_by("bin1_id"=="bin_id"))
                    hic.dtf <- dplyr::left_join(
                        (dplyr::mutate(hic.dtf, "bin2_id" = hic.dtf$j-1)),
                        (dplyr::rename(bins2, "weight2" = "weight")),
                        dplyr::join_by("bin2_id"=="bin_id"))
                    if(cool_divisive_weights){
                        hic.dtf <- dplyr::mutate(
                            hic.dtf,
                            "weight1"=1/hic.dtf$weight1
                        )
                        hic.dtf <- dplyr::mutate(
                            hic.dtf,
                            "weight2"=1/hic.dtf$weight2
                        )
                    }
                }
                filter.bin2 <- hic.dtf$j %in%
                    starts.ndx[chrom_2]:ends.ndx[chrom_2]
                hic.dtf <- hic.dtf[filter.bin2, ]
                # I have no idea what's the intent of the following lines
                # just leaving them as they were...
                hic.dtf <- dplyr::mutate(
                    hic.dtf,
                    i = hic.dtf$i - starts.ndx[chrom_1] + 1
                )
                hic.dtf <- dplyr::mutate(
                    hic.dtf,
                    j = hic.dtf$j - starts.ndx[chrom_2] + 1
                )
            }else if (GetFileExtension(file) == "h5") {
                # Define start and end of chromosomes
                ends.ndx <- chromSizes$dimension |>
                    cumsum() |>
                    stats::setNames(chromSizes$name)
                starts.ndx <- 1 + c(0, ends.ndx[-length(ends.ndx)]) |>
                    stats::setNames(chromSizes$name)
                hic.spm <- hic_spm_full_h5[
                    starts.ndx[chrom_1]:ends.ndx[chrom_1],
                    starts.ndx[chrom_2]:ends.ndx[chrom_2]]
            } else if (GetFileExtension(file) == "bedpe") {
                hic.dtf <- dplyr::filter(
                    megaHic.dtf,
                    # https://github.com/tidyverse/dplyr/issues/3139
                    # conflicting variable name
                    megaHic.dtf$chrom_1 == (!!chrom_1) &
                    megaHic.dtf$chrom_2 == (!!chrom_2)
                )
            }
            # Create Contact matrix
            if(cool_balanced && GetFileExtension(file)%in%c("cool","mcool")){
                hic.spm <- Matrix::sparseMatrix(
                    i = hic.dtf$i,
                    j = hic.dtf$j,
                    x = (hic.dtf$counts*hic.dtf$weight1*hic.dtf$weight2),
                    dims = dims.num
                )
            }else if(GetFileExtension(file) != "h5"){
                hic.spm <- Matrix::sparseMatrix(
                    i = hic.dtf$i,
                    j = hic.dtf$j,
                    x = hic.dtf$counts,
                    dims = dims.num
                )
            }
            
            row.regions <- binnedGenome.grn[which(
                as.vector(binnedGenome.grn@seqnames) == chrom_1
            )]
            col.regions <- binnedGenome.grn[which(
                as.vector(binnedGenome.grn@seqnames) == chrom_2
            )]
            hic <- InteractionSet::ContactMatrix(
                hic.spm,
                row.regions,
                col.regions
            )
            # Metadata
            hic@metadata <- dplyr::filter(
                attributes.tbl,
                attributes.tbl$name == paste(ele.lst, collapse = "_")
            ) |>
            tibble::add_column(resolution = hicResolution) |>
            as.list()
            if(hic_norm!="NONE" && GetFileExtension(file)=="hic"){
                hic@metadata <- append(hic@metadata,
                list(observed = hic.dtf$counts,
                normalizer = NULL,
                mtx = "norm"))
            }else if (cool_balanced && 
            GetFileExtension(file)%in%c("cool","mcool")) {
                hic@metadata <- append(hic@metadata,
                    list(observed = hic.dtf$counts,
                    normalizer = (hic.dtf$weight1 * hic.dtf$weight2),
                    mtx = "norm"))
            }
            if(hic_matrix!="obs" && GetFileExtension(file)=="hic"){
                hic@metadata <- append(hic@metadata,
                list(expected=hic_matrix))
                }
            
            return(hic)
        }
    )
    if(hic_matrix!="obs" && GetFileExtension(file)=="hic"){
        # Add attributes
        hic.lst_cmx <- hic.lst_cmx |>
            stats::setNames(chromComb.lst) |>
            AddAttr(
                attrs = list(
                    resolution = hicResolution,
                    chromSize = tibble::as_tibble(chromSizes),
                    matricesKind = attributes.tbl,
                    mtx = hic_matrix
                )
            )
    }else if (cool_balanced && 
            GetFileExtension(file)%in%c("cool","mcool")) {
        hic.lst_cmx <- hic.lst_cmx |>
            stats::setNames(chromComb.lst) |>
            AddAttr(
                attrs = list(
                    resolution = hicResolution,
                    chromSize = tibble::as_tibble(chromSizes),
                    matricesKind = attributes.tbl,
                    mtx = "norm"
                )
            )
    }else {
        # Add attributes
        hic.lst_cmx <- hic.lst_cmx |>
            stats::setNames(chromComb.lst) |>
            AddAttr(
                attrs = list(
                    resolution = hicResolution,
                    chromSize = tibble::as_tibble(chromSizes),
                    matricesKind = attributes.tbl,
                    mtx = "obs"
                )
            )
    }
    return(hic.lst_cmx)
}
