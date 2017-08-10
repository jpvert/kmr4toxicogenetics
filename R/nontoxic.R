#' Non-toxic chemical compounds
#' 
#' 15 out of 106 chemical compounds in toxicity data that were shown to have no toxicity across the human cell population.
#' 
#' @return A sequence of identifiers for those 15 non-toxic chemical compounds in toxicity data.
#' 
#' @export
#' 
#' @references
#' Eduati, F., et al. "Prediction of human population responses to toxic compounds by a collaborative competition." Nature biotechnology 33.9 (2015): 933-940. \href{https://doi.org/10.1038/nbt.3299}{doi:10.1038/nbt.3299}.
#' 

nontoxic <- function()
{
  return(c("NCGC00091360.03",
           "NCGC00091837.02",
           "NCGC00016348.08",
           "NCGC00091681.02",
           "NCGC00090795.09",
           "NCGC00248438.01",
           "NCGC00091329.02",
           "NCGC00090976.02",
           "NCGC00091056.03",
           "NCGC00091546.02",
           "NCGC00091067.02",
           "NCGC00091311.03",
           "NCGC00090773.05",
           "NCGC00090742.03",
           "NCGC00091792.02"))
}
