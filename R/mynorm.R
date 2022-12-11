#' For internal usage of normalization
#' @description For internal usage of normalization
#' @import graphics stats
#' @importFrom EnvStats geoMean
#' @param data.raw raw DSP data matrix
#' @param area array of areas of ROIs
#' @param nuclei.count array of nuclei count
#' @param control.hk hk controls
#' @param control.neg neg controls
#' @param normalization.method normalization method
#' @noRd
DSPNorm.s <- function(data.raw, area, nuclei.count, control.neg, control.hk, normalization.method){



  if(normalization.method=="None")
  {
    data.norm.f = data.raw

  }


  # normalize data by geomean of area:
  if(normalization.method=="Area_geoMean")
  {
    data.norm.f = data.raw

    Area.geomean = geoMean(as.numeric(area), na.rm = T)
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id]*(Area.geomean/area[col_id])
    }


  }



  # normalize data by geomean of nuclei:
  if(normalization.method=="Nuclei_geoMean")
  {
    data.norm.f = data.raw

    control.Nuclei.rm0 = as.numeric(nuclei.count)
    if(sum(control.Nuclei.rm0 == 0) > 0) {
      cat("Warning: in doing Nuclei_geoMean normalization, found ROIs with 0 nuclei count, replacing 0 with 1!\n")
      control.Nuclei.rm0[control.Nuclei.rm0 == 0] = 1
    }

    Nuclei.geomean = geoMean(control.Nuclei.rm0, na.rm = T)
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id]*(Nuclei.geomean/control.Nuclei.rm0[col_id])
    }


  }




  # normalize data by geomean of HK controls:
  if(normalization.method=="HK_geoMean")
  {
    if(!all(c(control.hk)%in%row.names(data.raw)))
    {
      cat("Warning: in doing HK_geoMean normalization, some HK controls are NOT contained in your raw data!\n")
      return(NULL)
    }

    data.norm.f = data.raw

    control.HK.align = data.raw[control.hk, colnames(data.raw)]
    mean.geometric = rep(NA, dim(control.HK.align)[2])
    for(col_id in 1:dim(control.HK.align)[2])
    {
      mean.geometric[col_id] = geoMean(control.HK.align[ ,col_id], na.rm = T)
    }
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id]*(mean(mean.geometric, na.rm = T)/mean.geometric[col_id])
    }


  }





  # normalize data by geomean of negative controls (a.k.a, Signal to Noise Ratio -- "STNR_geoMean"):
  if(normalization.method=="STNR_geoMean")
  {
    if(!all(c(control.neg)%in%row.names(data.raw)))
    {
      cat("Warning: in doing STNR_geoMean normalization, some negative controls are NOT contained in your raw data!\n")
      return(NULL)
    }

    data.norm.f = data.raw

    control.SignalToNoise.align = data.raw[control.neg, colnames(data.raw)]
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id]/geoMean(control.SignalToNoise.align[ ,col_id], na.rm = T)
    }

  }



  # normalize data by subtraction of mean + 2*SD (a.k.a, "STNS_mean2SD") of negative controls:
  if(normalization.method=="STNS_mean2SD")
  {
    if(!all(c(control.neg)%in%row.names(data.raw)))
    {
      cat("Warning: in doing STNS_mean2SD normalization, some negative controls are NOT contained in your raw data!\n")
      return(NULL)
    }

    data.norm.f = data.raw

    control.SignalToNoise.align = data.raw[control.neg, colnames(data.raw)]
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id] - (mean(control.SignalToNoise.align[ ,col_id], na.rm = T) + 2*sd(control.SignalToNoise.align[ ,col_id], na.rm = T))
    }

    data.norm.f[data.norm.f <= 0] = 1.0


  }


  return(data.norm.f)
}
