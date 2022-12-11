#' Normalize Data of NanoString GeoMx Protein Assays
#' @description DSPNorm provides flexible options of normalization methods for raw data (only ERCC-normalized) of NanoString GeoMx protein assays for downstream performance evaluation.
#' @return a matrix or data.frame object, same class and data structure as the input of raw DSP data matrix.
#' @export
#' @importFrom EnvStats geoMean
#' @param data.raw a raw DSP data matrix that has only been ERCC-normalized. The raw data matrix may be either a matrix or data.frame object with each row for a protein and each column for a ROI.
#' @param area a numerical vector of ROI areas, must be arranged in the same order of the columns of data.raw. If missing, the normalization methods concerning area cannot be performed.
#' @param nuclei.count a numerical vector of ROI nuclei counts, must be arranged in the same order of the columns of data.raw. If missing, the normalization methods concerning nuclei count cannot be performed.
#' @param norm.instruction a list specifying scaling: a string of "None", "Area_geoMean" (in which case the argument area must be provided), or "Nuclei_geoMean"(in which case the argument nuclei.count must be provided); background removal: a string of "None", "STNR_geoMean" or "STNS_mean2SD" (in which case a character vector of negative controls, control.neg, must be provided); and if normalization to HK controls shall be applied: a string of "None" or "HK_geoMean" (in which case a character vector of housekeeping controls, control.hk, must be provided).
#' @examples
#' data(COVID.19.DSP)
#'
#' data.raw <- COVID.19.DSP$Expr.raw.DSP
#'
#' group.info <- COVID.19.DSP$Meta.DSP
#'
#' # normalize the raw data by area of each ROI:
#' data.norm <- DSPNorm(data.raw = data.raw,
#'                      area = group.info$area,
#'                      norm.instruction = list(scaling="Area_geoMean",
#'                                              bkg.removal="None",
#'                                              norm.to.hk="None"))
#'
#' # normalize the raw data by geometric means of the two housekeeping controls "S6" and "GAPDH":
#' data.norm <- DSPNorm(data.raw = data.raw,
#'                      norm.instruction = list(scaling="None",
#'                                              bkg.removal="None",
#'                                              norm.to.hk="HK_geoMean",
#'                                              control.hk = c("S6", "GAPDH")))
#'
#' # normalize the raw data first by background removal, "STNS_mean2SD", then by geometric means
#' # of three housekeeping controls "S6", "GAPDH", and "Histone H3":
#' data.norm <- DSPNorm(data.raw = data.raw,
#'                      norm.instruction = list(scaling="None",
#'                                              bkg.removal="STNS_mean2SD",
#'                                              control.neg=c("Ms IgG1", "Ms IgG2a", "Rb IgG"),
#'                                              norm.to.hk="HK_geoMean",
#'                                              control.hk = c("S6", "GAPDH", "Histone H3")))
DSPNorm <- function(data.raw, area=NA, nuclei.count=NA, norm.instruction=list(scaling="None",
                                                                              bkg.removal="None",
                                                                              control.neg=c("Ms IgG1", "Ms IgG2a", "Rb IgG"),
                                                                              norm.to.hk="HK_geoMean",
                                                                              control.hk=c("S6", "GAPDH", "Histone H3"))){

  # checking for inconsistencies:
  if(norm.instruction$scaling=="Area_geoMean" & all(is.na(area))){

    stop("for scaling = 'Area_geoMean', a numerical vector for the argument 'area' must be provided\n")

  }
  if(norm.instruction$scaling=="Nuclei_geoMean" & all(is.na(nuclei.count))){

    stop("for scaling = 'Nuclei_geoMean', a numerical vector for the argument 'nuclei.count' must be provided\n")

  }
  # if(norm.instruction$bkg.removal%in%c("STNR_geoMean", "STNS_mean2SD") & (all(is.na(norm.instruction$control.neg)))){
  #
  #   stop("for bkg.removal = 'STNR_geoMean' or 'STNS_mean2SD', a categorical vector of the negative controls for the argument 'control.neg' must be provided\n")
  #
  # }
  if(norm.instruction$bkg.removal%in%c("STNR_geoMean", "STNS_mean2SD") & (!all(norm.instruction$control.neg%in%rownames(data.raw)))){

    stop("for bkg.removal = 'STNR_geoMean' or 'STNS_mean2SD', some negative control(s) you provided for the argument 'control.neg' cannot be found in the raw data matrix\n")

  }
  if(norm.instruction$norm.to.hk=="HK_geoMean" & (!all(norm.instruction$control.hk%in%rownames(data.raw)))){

    stop("for norm.to.hk = 'HK_geoMean', some housekeeping control(s) you provided for the argument 'control.hk' cannot be found in the raw data matrix\n")

  }



  data.scaled <- DSPNorm.s(data.raw = data.raw,
                           area = area,
                           nuclei.count = nuclei.count,
                           control.hk = norm.instruction$control.hk,
                           control.neg = norm.instruction$control.neg,
                           normalization.method = norm.instruction$scaling)

  data.scaled.bkg_corrected <- DSPNorm.s(data.raw = data.scaled,
                                         area = area,
                                         nuclei.count = nuclei.count,
                                         control.hk = norm.instruction$control.hk,
                                         control.neg = norm.instruction$control.neg,
                                         normalization.method = norm.instruction$bkg.removal)


  data.scaled.bkg_corrected.hknormed <- DSPNorm.s(data.raw = data.scaled.bkg_corrected,
                                                  area = area,
                                                  nuclei.count = nuclei.count,
                                                  control.hk = norm.instruction$control.hk,
                                                  control.neg = norm.instruction$control.neg,
                                                  normalization.method = norm.instruction$norm.to.hk)

  return(data.scaled.bkg_corrected.hknormed)
}

