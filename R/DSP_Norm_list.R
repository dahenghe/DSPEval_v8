#' Normalize Data of NanoString GeoMx Protein Assays
#' @description DSPNorm.preListed provides convenience of easily choosing among the 17 pre-defined unique normalization methods for raw data (only ERCC-normalized) of NanoString GeoMx protein assays for downstream performance evaluation. Note that those pre-defined normalization methods that are not applicable due to absence of required input info will be automatically skipped.
#' @return a list object containing the following 3 slots: 1. "norm.data.list": a list of normalized DSP data matrices by selected normalization methods; 2. "methods.available": a numerical vector of those selected method IDs; 3. "methods.table.summary": a table of data.frame class manifesting the details of all the 17 methods.
#' @export
#' @importFrom EnvStats geoMean
#' @param data.raw a raw DSP data matrix that has only been ERCC-normalized. The raw data matrix may be either a matrix or data.frame object with each row for a protein and each column for a ROI.
#' @param area a numerical vector of ROI areas, must be arranged in the same order of the columns of data.raw. If missing, the normalization methods concerning area cannot be performed.
#' @param nuclei.count a numerical vector of ROI nuclei counts, must be arranged in the same order of the columns of data.raw. If missing, the normalization methods concerning nuclei count cannot be performed.
#' @param control.hk a character vector listing all the housekeeping proteins that appear in the row names of data.raw.
#' @param control.neg a character vector listing all the negative control proteins that appear in the row names of data.raw.
#' @examples
#' data(COVID.19.DSP)
#'
#' data.raw <- COVID.19.DSP$Expr.raw.DSP
#'
#' group.info <- COVID.19.DSP$Meta.DSP
#'
#' # automatically obtain normalized DSP data by the 17 built-in methods:
#' preListed <- DSPNorm.preListed(data.raw = data.raw,
#'                                area = group.info$area,
#'                                nuclei.count=NA,
#'                                control.hk=c("S6", "GAPDH", "Histone H3"),
#'                                control.neg=c("Ms IgG1", "Ms IgG2a", "Rb IgG"))
#'
#' # the methods 6, 7, and 17 are skipped due the absence of nuclei.count:
#' preListed$methods.available
#' # [1]  1  2  3  4  5  8  9 10 11 12 13 14 15 16
DSPNorm.preListed <- function(data.raw, area=NA, nuclei.count=NA, control.hk=c("S6", "GAPDH", "Histone H3"), control.neg=c("Ms IgG1", "Ms IgG2a", "Rb IgG")){

  # check if controls can be found in data:
  if(!all(control.neg%in%rownames(data.raw))){

    stop("some negative control(s) you provided for the argument 'control.neg' cannot be found in the raw data matrix\n")

  }
  if(!all(control.hk%in%rownames(data.raw))){

    stop("some housekeeping control(s) you provided for the argument 'control.hk' cannot be found in the raw data matrix\n")

  }

  methods.uniq = list()
  methods.uniq[[1]] = list(scaling="None", bkg.removal="None", control.neg="-", norm.to.hk="None", control.hk="-")

  methods.uniq[[2]] = list(scaling="None", bkg.removal="STNR_geoMean", control.neg=control.neg, norm.to.hk="None", control.hk="-")

  methods.uniq[[3]] = list(scaling="None", bkg.removal="STNS_mean2SD", control.neg=control.neg, norm.to.hk="None", control.hk="-")

  methods.uniq[[4]] = list(scaling="Area_geoMean", bkg.removal="None", control.neg="-", norm.to.hk="None", control.hk="-")

  methods.uniq[[5]] = list(scaling="Area_geoMean", bkg.removal="STNS_mean2SD", control.neg=control.neg, norm.to.hk="None", control.hk="-")

  methods.uniq[[6]] = list(scaling="Nuclei_geoMean", bkg.removal="None", control.neg="-", norm.to.hk="None", control.hk="-")

  methods.uniq[[7]] = list(scaling="Nuclei_geoMean", bkg.removal="STNS_mean2SD", control.neg=control.neg, norm.to.hk="None", control.hk="-")

  methods.uniq[[8]] = list(scaling="None", bkg.removal="None", control.neg="-", norm.to.hk="HK_geoMean", control.hk=control.hk)

  methods.uniq[[9]] = list(scaling="None", bkg.removal="None", control.neg="-", norm.to.hk="HK_geoMean", control.hk=control.hk[c(1,2)])

  methods.uniq[[10]] = list(scaling="None", bkg.removal="None", control.neg="-", norm.to.hk="HK_geoMean", control.hk=control.hk[c(1,3)])

  methods.uniq[[11]] = list(scaling="None", bkg.removal="None", control.neg="-", norm.to.hk="HK_geoMean", control.hk=control.hk[c(2,3)])

  methods.uniq[[12]] = list(scaling="None", bkg.removal="None", control.neg="-", norm.to.hk="HK_geoMean", control.hk=control.hk[1])

  methods.uniq[[13]] = list(scaling="None", bkg.removal="None", control.neg="-", norm.to.hk="HK_geoMean", control.hk=control.hk[2])

  methods.uniq[[14]] = list(scaling="None", bkg.removal="None", control.neg="-", norm.to.hk="HK_geoMean", control.hk=control.hk[3])

  methods.uniq[[15]] = list(scaling="None", bkg.removal="STNS_mean2SD", control.neg=control.neg, norm.to.hk="HK_geoMean", control.hk=control.hk)

  methods.uniq[[16]] = list(scaling="Area_geoMean", bkg.removal="STNS_mean2SD", control.neg=control.neg, norm.to.hk="HK_geoMean", control.hk=control.hk)

  methods.uniq[[17]] = list(scaling="Nuclei_geoMean", bkg.removal="STNS_mean2SD", control.neg=control.neg, norm.to.hk="HK_geoMean", control.hk=control.hk)


  norm.data.list <- list()
  methods.table <- data.frame(method.ID=as.character(1:17), scaling="-", bkg.removal="-", control.neg="-", norm.to.hk="-", control.hk="-")
  normalization.methods <- 0
  for (method_id in 1:17) {

    conditions.skipping.1 <- methods.uniq[[method_id]]$scaling=="Area_geoMean" & all(is.na(area))
    conditions.skipping.2 <- methods.uniq[[method_id]]$scaling=="Nuclei_geoMean" & all(is.na(nuclei.count))
    conditions.skipping <- conditions.skipping.1 | conditions.skipping.2

    if(conditions.skipping){

      if(conditions.skipping.1) {
        warning(paste("for the method ", method_id, ", the area info is missing, therefore the scaling by 'Area_geoMean' is not available\n", sep = ""))
      }
      if(conditions.skipping.2) {
        warning(paste("for the method ", method_id, ", the nuclei.count info is missing, therefore the scaling by 'Nuclei_geoMean' is not available\n", sep = ""))
      }

      methods.table$method.ID[method_id] <- paste(methods.table$method.ID[method_id], " (not available)", sep = "")

      methods.table$scaling[method_id] <- methods.uniq[[method_id]]$scaling
      methods.table$bkg.removal[method_id] <- methods.uniq[[method_id]]$bkg.removal
      methods.table$control.neg[method_id] <- paste(methods.uniq[[method_id]]$control.neg, collapse = "/")
      methods.table$norm.to.hk[method_id] <- methods.uniq[[method_id]]$norm.to.hk
      methods.table$control.hk[method_id] <- paste(methods.uniq[[method_id]]$control.hk, collapse = "/")

      next

    }

    norm.data.list[[method_id]] <- DSPNorm(data.raw=data.raw,
                                           area=area,
                                           nuclei.count=nuclei.count,
                                           norm.instruction=methods.uniq[[method_id]])

    normalization.methods <- c(normalization.methods, method_id)

    methods.table$scaling[method_id] <- methods.uniq[[method_id]]$scaling
    methods.table$bkg.removal[method_id] <- methods.uniq[[method_id]]$bkg.removal
    methods.table$control.neg[method_id] <- paste(methods.uniq[[method_id]]$control.neg, collapse = "/")
    methods.table$norm.to.hk[method_id] <- methods.uniq[[method_id]]$norm.to.hk
    methods.table$control.hk[method_id] <- paste(methods.uniq[[method_id]]$control.hk, collapse = "/")

  }

  methods.available <- normalization.methods[normalization.methods!=0]

  methods.table.summary <- methods.table

  contents.preListed <- list()

  contents.preListed$norm.data.list <- norm.data.list
  contents.preListed$methods.available <- methods.available
  contents.preListed$methods.table.summary <- methods.table.summary

  return(contents.preListed)
}

