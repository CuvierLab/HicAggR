#' Draw aggregation plot.
#'
#' PlotAPA
#' @description Draw aggregation plot from aggregation matrices.
#' @param aggregatedMtx <matrix>: The aggregated matrix.
#' @param trim <numeric>: A number between 0 and 100 that gives
#' the percentage of triming in matrices.
#' @param colMin <matrix>: The minimal value in color scale.
#'  If Null automaticaly find.
#' @param colMid <matrix>: The middle value in color scale.
#'  If Null automaticaly find.
#' @param colMax <matrix>: The mximal value in color scale.
#'  If Null automaticaly find.
#' @param colMinCond <matrix>: Avalaible for plotting
#' differential aggregation. The minimal value in color scale in the
#' classsical aggregation plot. If `NULL` automaticaly find.
#' @param colMaxCond <matrix>: Avalaible for plotting
#' differantial aggregation. The maxiaml value in color scale in the
#' classsical aggregation plot. If `NULL` automaticaly find.
#' @param extra_info <logical> do you want to have a recall of
#' your arguments values? (Default FALSE)
#' @param ... additional arguments to [ggAPA()]
#' @return None
#' @export
#' @importFrom gridExtra grid.table
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
#' PlotAPA(
#'     aggregatedMtx = aggreg.mtx
#' )
#'
#' PlotAPA(
#'     aggregatedMtx= aggreg.mtx,
#'     trim= 20,
#'     colMin= -2,
#'     colMid= 0,
#'     colMax= 2,
#'     colMinCond = 0,
#'     colMaxCond = 2
#' )
#' 

PlotAPA <- function(aggregatedMtx = NULL, trim=0, colMin=NULL, colMid=NULL,
    colMax=NULL, colMinCond=NULL, colMaxCond=NULL, extra_info = FALSE,...){
    .ggDensity <- function(data.lst=NULL, colour.col=NULL, mean.bln=TRUE,
    title=NULL){
        data.lst_tbl <- lapply(seq_along(data.lst),function(element.ndx){
            return(tibble::tibble(
                value = data.lst[[element.ndx]],
                class = factor(names(data.lst)[[element.ndx]])
            ))
            }) 
        data.tbl <- dplyr::bind_rows(data.lst_tbl)
        if(is.null(colour.col)){
            colour.col <- Hue(length(data.lst)) |> 
                stats::setNames(names(data.lst))
        }
        plot.gp <- ggplot2::ggplot(data.tbl, 
            ggplot2::aes(x=data.tbl$value, fill=data.tbl$class,
                colour=data.tbl$class)) + 
            ggplot2::geom_density(alpha=0.1) +
            ggplot2::scale_color_manual(values = colour.col)+
            ggplot2::scale_fill_manual(values = colour.col) +
            ggplot2::labs(title=title) 
        if (mean.bln){
            data.tbl <- dplyr::group_by(data.tbl, class = data.tbl$class)
            mu.tbl <-  dplyr::summarise(data.tbl,
                grp.mean = mean(data.tbl$value))
            plot.gp <- plot.gp + 
                ggplot2::geom_vline(data = mu.tbl,
                    ggplot2::aes(xintercept = mu.tbl$grp.mean,
                        colour = mu.tbl$class),
                    linetype = "dashed")
        }
        return(plot.gp)
    }
    # Differential or not?
        differential.bln <- !is.null(attributes(aggregatedMtx)$matrices)
        if(differential.bln){
            colors <- NULL
        }else{
            colors <- viridis(255)
        }
    # Plot
    # Auto Scale
    plot.gp <- ggAPA(
        aggregatedMtx=aggregatedMtx, 
        colors=colors,
        title=ifelse(differential.bln,
            yes="Agregation of differential matrices",
    no="Agregation"),...
    ) + ggplot2::labs(subtitle="scale (auto), center()")
    plot(plot.gp)
# Auto Scale + Center
    if(!is.null(colMid)){
        plot.gp <- ggAPA(
            aggregatedMtx=aggregatedMtx, 
            colors=colors,
            colMid=colMid,
            title=ifelse(differential.bln,
                yes="Agregation of differential matrices",
                no="Agregation"),...
        ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",colMid,")"))
        plot(plot.gp)
    }
    # Trim Scale + Center
    if(!is.null(trim) && 0<trim){
        plot.gp <- ggAPA(
            aggregatedMtx=aggregatedMtx, 
            colors=colors,
            trim=trim,
            colMid=colMid,
            title=ifelse(differential.bln,
                yes="Agregation of differential matrices",
                no="Agregation"),...
        ) + ggplot2::labs(subtitle=
            paste0("scale (rm ",trim,"%), center(",colMid,")"))
        plot(plot.gp)
    }
    # MinMax Scale + Center
    if(!is.null(colMin) || !is.null(colMax)){
        plot.gp <- ggAPA(
            aggregatedMtx=aggregatedMtx, 
            colors=colors,
            colMin=colMin,
            colMid=colMid,
            colMax=colMax,
            title=ifelse(differential.bln,
                yes="Agregation of differential matrices",
                no="Agregation"),...
        ) + ggplot2::labs(subtitle=
            paste0("scale (",colMin,";",colMax,"), center(",colMid,")"))
        plot(plot.gp)
    }
    if (differential.bln){
        # Pval + Auto Scale
        if(!is.null(attributes(aggregatedMtx)$matrices$pVal) &&
        sum(!is.na(attributes(aggregatedMtx)$matrices$pVal))>=3){
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$pVal,
                colors=YlOrRd(9),
                title = "-log10(p.values)",...
            ) + ggplot2::labs(subtitle="scale (auto), center()")
        }else{
            plot.gp <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::annotate("text", x = 1, y = 1,
                    label = "Not enough pval computed to plot 
                    a pval matrix (<3) or nothing significant")
        }
        plot(plot.gp)
        # FiltPval + Auto Scale + Center
        if(!is.null(attributes(aggregatedMtx)$matrices$aggDiffPvalFilt) &&
        sum(!is.na(attributes(aggregatedMtx)$matrices$pVal))>=3){
            plot.gp <- ggAPA(
                aggregatedMtx=
                attributes(aggregatedMtx)$matrices$aggDiffPvalFilt,
                colors=colors,
                colMid=colMid,
                title = paste0("Agregation of differential matrices"),...
            ) + ggplot2::labs(subtitle=
                paste0("filtred by p.values, scale (auto), center(",colMid,")"))
        }else{
            plot.gp <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::annotate("text", x = 1, y = 1,
                label = "Not enough pval computed to plot a pval
                matrix (<3) or nothing significant")
        }
        plot(plot.gp)
        # FiltPval + Trim Scale + Center
        if(!is.null(attributes(aggregatedMtx)$matrices$aggDiffPvalFilt) &&
        sum(!is.na(attributes(aggregatedMtx)$matrices$pVal))>=3){
            plot.gp <- ggAPA(
                aggregatedMtx=
                attributes(aggregatedMtx)$matrices$aggDiffPvalFilt,
                colors=colors,
                trim=trim,
                colMid=colMid,
                title = paste0("Agregation of differential matrices"),...
            ) + ggplot2::labs(subtitle=paste0(
                "filtred by p.vals, scale (rm ",trim,"%), center(",colMid,")"))
        }else{
            plot.gp <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::annotate("text", x = 1, y = 1,
                label = "Not enough pval computed to plot a pval
                matrix (<3) or nothing significant")
        }
        plot(plot.gp)
        # Delta + Auto Scale + Center
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$aggDelta, 
            colors=colors,
            colMid=colMid,
            title="Differential of agregated matrices",...
        ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",colMid,")"))
        plot(plot.gp)
        # Delta + Trim Scale + Center
        if(!is.null(trim) && 0<trim){
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggDelta, 
                colors=colors,
                trim=trim,
                colMid=colMid,
                title="Differential of agregated matrices",...
            ) + ggplot2::labs(subtitle=
                paste0("scale (rm ",trim,"%), center(",colMid,")"))
            plot(plot.gp)
        }
        # Delta + Auto Scale + Center
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrectedDelta, 
            colors=colors,
            colMid=colMid,
            title="Differential of corrected agregated matrices",...
        ) + ggplot2::labs(subtitle=paste0("scale (auto), center(",colMid,")"))
        plot(plot.gp)
        # Delta + Trim Scale + Center
        if(!is.null(trim) && 0<trim){
            plot.gp <- ggAPA(
                aggregatedMtx=
                attributes(aggregatedMtx)$matrices$aggCorrectedDelta, 
                colors=colors,
                trim=trim,
                colMid=colMid,
                title="Differential of corrected agregated matrices",...
            ) + ggplot2::labs(subtitle=
                paste0("scale (rm ",trim,"%), center(",colMid,")"))
            plot(plot.gp)
        }
        # Control + Auto Scale
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
            colors=viridis(51),
            title="Agregation control",...
        ) + ggplot2::labs(subtitle="scale (auto)")
        plot(plot.gp)
        # Control + Trim Scale
        if(!is.null(trim) && 0<trim){
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
                colors=viridis(51),
                trim=trim,
                title="Agregation control",...
            ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%)"))
            plot(plot.gp)
        }
        # Control + MinMax Scale
        if(!is.null(colMinCond) || !is.null(colMaxCond)){
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
                colors=viridis(51),
                colMin=colMinCond,
                colMax=colMaxCond,
                title="Agregation control",...
            ) + ggplot2::labs(subtitle=
                paste0("scale (",colMinCond,";",colMaxCond,")"))
            plot(plot.gp)
        }
        # Condition + Auto Scale
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$agg, 
            colors=viridis(51),
            title="Agregation",...
        ) + ggplot2::labs(subtitle="scale (auto)")
        plot(plot.gp)
        # Condition + Trim Scale
        if(!is.null(trim) && 0<trim){
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$agg, 
                colors=viridis(51),
                trim=trim,
                title="Agregation",...
            ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%)"))
            plot(plot.gp)
        }
        # Condition + MinMax Scale
        if(!is.null(colMinCond) || !is.null(colMaxCond)){
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$agg, 
                colors=viridis(51),
                colMin=colMinCond,
                colMax=colMaxCond,
                title="Agregation",...
            ) + ggplot2::labs(subtitle=
                paste0("scale (",colMinCond,";",colMaxCond,")"))
            plot(plot.gp)
        }
        # Corrected condition + Auto Scale
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrected, 
            colors=viridis(51),
            title="Agregation corrected",...
        ) + ggplot2::labs(subtitle="scale (auto)")
        plot(plot.gp)
        # Corrected condition + Trim Scale
        if(!is.null(trim) && 0<trim){
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrected, 
                colors=viridis(51),
                trim=trim,
                title="Agregation corrected",...
            ) + ggplot2::labs(subtitle=paste0("scale (rm ",trim,"%)"))
            plot(plot.gp)
        }
        # Corrected condition + MinMax Scale
        if(!is.null(colMinCond) || !is.null(colMaxCond)){
            plot.gp <- ggAPA(
                aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrected, 
                colors=viridis(51),
                colMin=colMinCond,
                colMax=colMaxCond,
                title="Agregation corrected",...
            ) + ggplot2::labs(subtitle=
                paste0("scale (",colMinCond,";",colMaxCond,")"))
            plot(plot.gp)
        }
        # Grouped Scale(Condition & Control)
        colBreaks.num <- BreakVector(
            x=c(c(attributes(aggregatedMtx)$matrices$agg),
                c(attributes(aggregatedMtx)$matrices$aggCtrl)),
            n=51)
        # Control + Grouped Scale(with Condition)
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
            colors=viridis(51),
            colBreaks=colBreaks.num,
            title="Agregation control",...
        ) + ggplot2::labs(subtitle="scale (grouped with condition)")
        plot(plot.gp)
        # Condition + Grouped Scale(with Control)
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$agg, 
            colors=viridis(51),
            colBreaks=colBreaks.num,
            title="Agregation",...
        ) + ggplot2::labs(subtitle="scale (grouped with control)")
        plot(plot.gp)
        # Grouped Scale(Corrected condition & Control)
        colBreaks.num <- BreakVector(
            x=c(c(attributes(aggregatedMtx)$matrices$aggCorrected),
                c(attributes(aggregatedMtx)$matrices$aggCtrl)),
            n=51)
        # Control + Grouped Scale(with corrected condition)
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCtrl, 
            colors=viridis(51),
            colBreaks=colBreaks.num,
            title="Agregation control",...
        ) + ggplot2::labs(subtitle="scale (grouped with condition corrected)")
        plot(plot.gp)
        # Corrected condition + Grouped Scale(with Control)
        plot.gp <- ggAPA(
            aggregatedMtx=attributes(aggregatedMtx)$matrices$aggCorrected, 
            colors=viridis(51),
            colBreaks=colBreaks.num,
            title="Agregation corrected",...
        ) + ggplot2::labs(subtitle="scale (grouped with control)")
        plot(plot.gp)
        # Density

        plot.gp <- .ggDensity(
            data.lst=list(differential=stats::na.omit(c(aggregatedMtx))),
            colour.col=Hue(3)[[1]],
            title="Agregation of differential matrices density")
        plot(plot.gp)

        data.lst <- list(
            control=stats::na.omit(
                c(attributes(aggregatedMtx)$matrices$aggCtrl)),
            condition=stats::na.omit(c(attributes(aggregatedMtx)$matrices$agg))
        )
        plot.gp <- .ggDensity(
            data.lst=data.lst,
            colour.col=Hue(3)[2:3],
            title="Condition and control densities")
        plot(plot.gp)
        
        plot.gp <- .ggDensity(
            data.lst=list(deltaCorrected=
                stats::na.omit(
                    c(attributes(aggregatedMtx)$matrices$aggCorrectedDelta))),
            colour.col=Hue(3)[[1]],
            title="Differential of corrected agregated matrices density")
        plot(plot.gp)

        data.lst <- list(
            control=stats::na.omit(
                c(attributes(aggregatedMtx)$matrices$aggCtrl)),
            correctedCondition=stats::na.omit(
                c(attributes(aggregatedMtx)$matrices$aggCorrected))
        )
        plot.gp <- .ggDensity(
            data.lst=data.lst,
            colour.col=Hue(3)[2:3],
            title="Condition corrected and control densities")
        plot(plot.gp)
    }
    if(extra_info){
        # Attributes
        grid::grid.newpage()
        attr.ndx <- aggregatedMtx |>
            attributes() |>
            names() |>
            # NotIn to replace by %ni%
            NotIn(c("dim","matrices","interactions", "dimnames"))
        attr.lst <- attributes(aggregatedMtx)[attr.ndx]
        attr.lst$aggregationMethod <- function(pxl){
            pxl[is.na(pxl)]<-0;mean(pxl,na.rm=TRUE)}
        attr1.ndx <- attr.lst |>
            lapply(class) |>
            unlist() |>
            # NotIn to replace by %ni% 
            NotIn(c("matrix", "list","GInteractions","function"))
        attr1.lst <- attr.lst[attr1.ndx] |>
                    lapply(as.character) |>
                    unlist()
        attr2.ndx <- unlist(lapply(attr.lst, class)) %in% "function"
        attr2.lst <- attr.lst[attr2.ndx] |>
                    lapply(function(function.fun){
                        function.chr <- deparse(function.fun)
                        function.chr[3:(length(function.chr)-1)] |>
                        paste0(collapse=";\n")
                        return(gsub(" ","",function.chr))
                        }) |>
                    unlist()
        attr.lst <- c(attr1.lst, attr2.lst)
        attr.tbl  <- tibble::as_tibble_col(attr.lst) |>
            tibble::add_column(name=names(attr.lst)) |>
            tibble::column_to_rownames("name")
        gridExtra::grid.table(attr.tbl)
    }
}


#' Draw aggregation plots for interactions with different distances.
#' 
#' plotMultiAPA
#' @description Separates matrices based on interaction distance, performs
#'  aggregation and plots Aggregated signal for each chunk of interaction
#'  distances.
#' @param submatrices <listmatrix>: The matrices list to
#' separate using interaction distances and aggregate.
#' Chunks of distances are created with:
#' `c(0,50000*2 ^ seq(0,5,by=1))`. 
#' Other matrices with distances over 1.6 Mb are aggregated in the same
#'  final chunk.
#' @param ctrlSubmatrices <listmatrix>: The matrices list to use
#' as control condition for differential aggregation.
#' @param ... : Additional arguments to pass to
#' [Aggregation()].
#' For differential aggregation plot, `submatrices` will take the matrices
#'  of the treated condition. 
#' eg:
#' ```
#' plotMultiAPA(
#' submatrices = interactions_HS.mtx_lst,
#' ctrlSubmatrices = interactions_Ctrl.mtx_lst)```
#' @param plot.opts list of arguments to pass to [ggAPA()].
#' @return A plot with separate APAs per distance and a list of aggregated
#'  matrices as invisible output.
#' @export
#' @importFrom gridExtra grid.arrange
#' @examples
#' #' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#' data(HiC_HS.cmx_lst)
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
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for eg
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_Ctrl.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_Ctrl.cmx_lst,
#'     referencePoint = "pf"
#' )
#' interactions_HS.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_HS.cmx_lst,
#'     referencePoint = "pf"
#' )
#' interactions_Ctrl.mtx_lst <- PrepareMtxList(
#'     matrices = interactions_Ctrl.mtx_lst
#' )
#'
#' # Aggregate matrices in one matrix
#' plotMultiAPA(submatrices = interactions_Ctrl.mtx_lst)
#'
#'
#' interactions_HS.mtx_lst <- PrepareMtxList(
#'     matrices = interactions_HS.mtx_lst
#' )
#'
#' # Differential Aggregation
#' plotMultiAPA(
#'     submatrices = interactions_HS.mtx_lst,
#'     ctrlSubmatrices = interactions_Ctrl.mtx_lst,
#'     diffFun = "ratio",
#'     plot.opts = list(colors = list("blue","white","red"))
#' )
#' 
plotMultiAPA <- function(
    submatrices = NULL,
    ctrlSubmatrices=NULL,
    ...,
    plot.opts=NULL){
    # exchange variables if submatrices is null
    two_conditions <- TRUE
    if(is.null(submatrices) && !is.null(ctrlSubmatrices)){
        submatrices <- ctrlSubmatrices
        two_conditions <- FALSE
    }
    if(is.null(submatrices) && is.null(ctrlSubmatrices)){
        stop("Argument submatrices is null")
    }
    .validSubmatrices(submatrices)
    maxCol <- log2(max(attributes(submatrices)$interactions$distance)/50000)
    if(maxCol>5){
        vector_dist <- c(c(0,50000*2 ^ seq(0,5,by=1)),
            max(attributes(submatrices)$interactions$distance))
    }else{
        vector_dist <- c(0,50000*2 ^ seq(0,maxCol,by=1))
    }
    by_dist_vec_list <- list()
    noValues_vector <- c()
    n_cples <- c()
    all_cples <- attr(submatrices,"interactions")
    for(d in seq(2,length(vector_dist))){
        samples <- all_cples$name[which(all_cples$distance < vector_dist[d] & 
            all_cples$distance >= vector_dist[d-1])]
        filtered <- FilterInteractions(submatrices,targets=list(name=samples))

        if(!is.null(ctrlSubmatrices) && two_conditions){
            .validSubmatrices(ctrlSubmatrices)
            filtered_ctrl <- FilterInteractions(ctrlSubmatrices,targets=
                list(name=samples))
            if(length(filtered)>0 && length(filtered_ctrl)){
                by_dist_vec_list <- append(by_dist_vec_list,
                list(Aggregation(matrices = filtered,
                    ctrlMatrices = filtered_ctrl,...)))
                n_cples <- c(n_cples,length(attr(filtered,"interactions")))
            }else{
                noValues_vector <- c(noValues_vector,d)
            }
        }else{
            if(length(filtered)>0){
                by_dist_vec_list <- append(by_dist_vec_list,
                list(Aggregation(matrices = filtered,...)))
                n_cples <- c(n_cples,length(attr(filtered,"interactions")))
            }else{
                noValues_vector <- c(noValues_vector,d)
            }
        }
    }
    if(length(noValues_vector)>0){
        vector_dist <- vector_dist[-noValues_vector]
    }    
    plotList <- list()
    for(p in seq(1,length(by_dist_vec_list))){
        if(!is.null(plot.opts) && length(plot.opts)>0){
            chunk.plot <- list(do.call(ggAPA,c(list(
                aggregatedMtx = by_dist_vec_list[[p]]),
                as.list(plot.opts)))+
                ggplot2::labs(title=paste0(
                    "[",GenomicSystem(as.numeric(vector_dist[p])),"-",
                    GenomicSystem(as.numeric(vector_dist[p+1])),")",
                "\nn = ",n_cples[p])))
            plotList <- append(plotList,chunk.plot)
        }else{
            plotList <- append(plotList,list(ggAPA(by_dist_vec_list[[p]])+
                ggplot2::labs(title=paste0(
                    "[",GenomicSystem(as.numeric(vector_dist[p])),"-",
                    GenomicSystem(as.numeric(vector_dist[p+1])),")",
                "\nn = ",n_cples[p]))))
        }
    }
    gridExtra::grid.arrange(grobs = plotList,ncol=length(plotList),nrow=1)
    return(invisible(list("agg.matrices" = by_dist_vec_list,
    "distance.chunks"=vector_dist[-1],
    "nb.couples"=n_cples)))
}

#' preparePlotgardener
#' 
#' This function allows to obtain a dataframe that can be
#' used with plotgardener's plotHicTriangle, plotHicRectangle,
#' plotHicSquare. It is equivalent to plotgardener's readHic function.
#'
#' @param hicList hicList from which to extract data
#' @param ctrlHicList hicList to use as control
#' for a differential analysis.
#' @param submatrices list of submatrices.
#' @param submatrix.name name of submatrix
#' to focus on. use `names(submatrices)`. If submatrices is not `NULL`,
#' this argument must not be `NULL`.
#' @param diffFun The function used
#' to compute differential between counts column of hicList and ctrlHicList.
#' If the argument is a character, possible choices are:
#' \itemize{
#' \item "-", "substract" or "substraction" apply a substraction (Default)
#' \item "/" or "ratio" apply a ratio
#' \item "log2","log2-","log2/" or "log2ratio" apply a log2 on ratio
#' \item other apply a log2 on 1+ratio
#' }
#' @param which_chrom a combination
#' of chromsome names use `names(hicList)` to see the possible names.
#' @param which_range a GRanges object
#' or a string of type "chr1:1-100". see `StrToGranges()`.
#'
#' @return funcions returns a data.frame with 3 columns:
#' chrom, altchrom and counts. chrom being bin names of rows, altchrom being
#' bin names of columns and counts being the counts.
#' see also `strawr::straw()` & plotgardener::readHic()
#' @importFrom dplyr mutate select left_join join_by case_when filter
#' @importFrom GenomicRanges seqnames start
#' @importFrom checkmate assertList assert checkCharacter
#' assertChoice checkChoice
#' @export
#'
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#' data(HiC_HS.cmx_lst)
#' preparePlotgardener(
#'   hicList = HiC_HS.cmx_lst,
#'   ctrlHicList = HiC_Ctrl.cmx_lst,
#'   which_chrom = "2L_2L",
#'   diffFun = "substract")|> head()
#'   
preparePlotgardener <- function(
    hicList = NULL,
    ctrlHicList = NULL,
    submatrices = NULL,
    submatrix.name = NULL,
    diffFun = "log2ratio",
    which_chrom = NULL,
    which_range = NULL
) {
    if(is.null(submatrices) &&
        is.null(which_chrom) &&
        is.null(which_range) &&
        is.null(which_chrom)
    ){
        stop("Please provide a region to extract! You can provide:
            * 1 submatrix,
            * a combination of chromsomes eg: 'chr1_chr1',
            * a granges object
            * or a character as 'chr1:1-100'")
    }
    # Put list on correct variable
    if(!is.null(submatrices)){
        .validSubmatrices(submatrices)
        if(is.null(submatrix.name)){
        stop("provide a submatrix name to extract!
            use names(submatrices) to see the names")
        }else{
        submatrices <- FilterInteractions(
            submatrices,
            targets = list(name=submatrix.name))
        }
        checkmate::assertList(
        x = submatrices,
        max.len = 1
        )
        chromCombination <- paste(
        as.character(GenomicRanges::seqnames(InteractionSet::anchors(
        attributes(submatrices)$interactions)$first)),"_",
        as.character(GenomicRanges::seqnames(InteractionSet::anchors(
            attributes(submatrices)$interactions)$second))
        )
    }
    if(!is.null(hicList)){
        .validHicMatrices(matrices = hicList)
    } else if (!is.null(ctrlHicList) && is.null(hicList)) {
        .validHicMatrices(matrices = ctrlHicList)
        hicList <- ctrlHicList
        ctrlHicList <- NULL
    }
    if(!is.null(ctrlHicList)){
        .validHicMatrices(ctrlHicList)
    }
    if(!is.null(which_range)){
        if(is.character(which_range)){
        which_range <- StrToGRanges(which_range)
        }
        chromCombination <- paste0(
        as.character(GenomicRanges::seqnames(which_range)),
        "_",
        as.character(GenomicRanges::seqnames(which_range))
        )
    }
    if(!exists(x = "chromCombination") && 
        !is.null(which_chrom)){
        chromCombination <- which_chrom
    }
    checkmate::assert(
        checkmate::checkCharacter(
        x = chromCombination,
        null.ok = FALSE
        ),
        checkmate::checkChoice(
        x = chromCombination,
        choices = names(hicList),
        null.ok = FALSE
        ), combine = "and"
    )
    if(!is.null(ctrlHicList)){
        checkmate::assertChoice(
        x = chromCombination,
        choices = names(ctrlHicList),
        null.ok = FALSE
        )
    }
    df <- MeltSpm(hicList[[chromCombination]]@matrix) |>
        as.data.frame() |>
        dplyr::mutate(
            altchrom = 
                as.numeric(
                GenomicRanges::start(
                    hicList[[chromCombination]]@regions[.data$j])-1),
            chrom = 
                as.numeric(
                GenomicRanges::start(
                    hicList[[chromCombination]]@regions[.data$i])-1),
            counts=.data$x) |>
        dplyr::select(c("chrom","altchrom","counts"))
    
    if(!is.null(ctrlHicList)){
        
    }
    # Differential Function
    if (!is.function(diffFun) &&
        !is.null(hicList) &&
        !is.null(ctrlHicList)) {
        diffFun <- dplyr::case_when(
        tolower(diffFun) %in% c("-", "substract", "substraction") ~
            "function(mat.mtx,ctrl.mtx){mat.mtx - ctrl.mtx}",
        tolower(diffFun) %in% c("/", "ratio") ~
            "function(mat.mtx,ctrl.mtx){mat.mtx / ctrl.mtx}",
        tolower(diffFun) %in% c("log2", "log2-", "log2/", "log2ratio") ~
            "function(mat.mtx,ctrl.mtx){log2(mat.mtx)-log2(ctrl.mtx)}",
        TRUE ~
            "function(mat.mtx,ctrl.mtx){log2(mat.mtx+1)-log2(ctrl.mtx+1)}"
        )
        diffFun <- WrapFunction(diffFun)
        ctrl_df <- MeltSpm(ctrlHicList[[chromCombination]]@matrix) |>
        as.data.frame() |>
        dplyr::mutate(
            altchrom = 
            as.numeric(
                GenomicRanges::start(
                    ctrlHicList[[chromCombination]]@regions[.data$j])-1),
            chrom = 
            as.numeric(
                GenomicRanges::start(
                    ctrlHicList[[chromCombination]]@regions[.data$i])-1),
            counts = .data$x) |>
        dplyr::select(c("chrom","altchrom","counts"))
        df <- dplyr::left_join(
            x = df, y = ctrl_df,
            by = dplyr::join_by("chrom","altchrom")
        ) |> dplyr::mutate(
        counts = diffFun(.data$counts.x,.data$counts.y)
        ) |> dplyr::select(c("chrom","altchrom","counts")) |>
        dplyr::filter(!is.na(.data$counts))
    }
    return(df)
}