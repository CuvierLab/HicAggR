% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BalanceHiC.R
\name{BalanceHiC}
\alias{BalanceHiC}
\title{Compute HiC matrix-balancing.}
\usage{
BalanceHiC(
  hicLst,
  method = "ICE",
  interactionType = NULL,
  maxIter = 50,
  qtlTh = 0.15,
  cores = 1,
  verbose = FALSE
)
}
\arguments{
\item{hicLst}{<List\link[InteractionSet:ContactMatrix-class]{ContactMatrix}>:
The HiC maps list.}

\item{method}{\if{html}{\out{<character>}} : The kind of normalization method.
One of "ICE", "VC" or "VC_SQRT" (Default "ICE")}

\item{interactionType}{\if{html}{\out{<character>}} : "cis", "trans",
c("cis", "trans"), "all".
If NULL normalization is apply on cis contactMatrix then trans
contactMatrix (equivalent to c("cis", "trans")). If is "all", normalization
is apply on all contactMatrix at once. (Default NULL)}

\item{maxIter}{\if{html}{\out{<numerical>}}: The maximum iteration number.}

\item{qtlTh}{\if{html}{\out{<numerical>}}: The quantile threshold below which
the bins will be ignored. (Default 0.15)}

\item{cores}{\if{html}{\out{<numerical>}} : Number of cores to be used.
(Default 1)}

\item{verbose}{\if{html}{\out{<logical>}}: If TRUE show the progression
in console. (Default FALSE)}
}
\value{
A matrices list.
}
\description{
Apply a matrix-balancing normalization method to a
list of contacts matrix.
}
\details{
BalanceHiC
}
\examples{
data(HiC_Ctrl.cmx_lst)

HiC_Ctrl_ICE.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst,
    interactionType = "cis",
    method = "ICE"
)

HiC_Ctrl_VC.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst,
    interactionType = c("cis", "trans"),
    method = "VC"
)

HiC_Ctrl_VC_SQRT.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst,
    interactionType = "all",
    method = "VC_SQRT"
)

}
