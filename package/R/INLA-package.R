#' INLA
#'
#' @description Package to perform full Bayesian analysis on latent Gaussian models using
#'     Integrated Nested Laplace Approximations.
#'
#' See [https://www.r-inla.org/](https://www.r-inla.org/) for further details.
#' @aliases INLA
#' @seealso [inla()]
"_PACKAGE"
# This creates the appropriate INLA-package.Rd file in the way R and pkgdown likes it.
# Note: Cannot use rdname INLA, as that would clash with inla.Rd on case-insensitive platforms.
# Documentation not always available as development documentation, but is accessible
# after installing.

#' @rdname INLA-package
#' @export
INLA <- function() {
    message("Welcome to the R-INLA package!")
    utils::browseVignettes(package = "INLA")
}


# Need to import something from every Imports package...
#' @importFrom withr defer
NULL

# Special imports:
## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL




#' Bivariate Meta Analysis
#'
#' Data are taken from a meta-analysis to compare the utility of three types of
#' diagnostic imaging - lymphangiography (LAG), computed tomography (CT) and
#' magnetic resonance (MR) - to detect lymph node metastases in patients with
#' cervical cancer. The dataset consists of a total of 46 studies: the first 17
#' for LAG, the following 19 for CT and the last 10 for MR.
#'
#' @format A data frame with 92 observations on the following 9 variables.
#' \describe{
#' \item{N}{a numeric vector}
#' \item{Y}{a numeric vector}
#' \item{diid}{a numeric vector}
#' \item{lag.tp}{a numeric vector}
#' \item{lag.tn}{a numeric vector}
#' \item{ct.tp}{a numeric vector}
#' \item{ct.tn}{a numeric vector}
#' \item{mr.tp}{a numeric vector}
#' \item{mr.tn}{a numeric vector}
#' }
#' @references J. Scheidler and H. Hricak and K. K. Yu and L. Subak and M. R.
#' Segal,"Radiological evaluation of lymph node metastases in patients with
#' cervical cancer: a meta-analysis",JAMA 1997
#' @keywords datasets
#' @docType data
#' @examples
#'
#' data(BivMetaAnalysis)
"BivMetaAnalysis"





#' ~~ data name/kind ... ~~
#'
#' ~~ A concise (1-5 lines) description of the dataset. ~~
#'
#'
#' @name Cancer
#' @docType data
#' @format A data frame with 6690 observations on the following 4 variables.
#' \describe{
#' \item{Y}{Number of cases}
#' \item{N}{a numeric vector}
#' \item{Age}{a numeric vector}
#' \item{region}{a numeric vector}
#' }
#' @references Rue, H and Held, L. (2005) *Gaussian Markov Random Fields -
#' Theory and Applications* Chapman and Hall
#' @keywords datasets
NULL





#' Time series with seasonal effect
#'
#' Montly total of car drivers killed or several injuried in England from
#' January 1969 to December 1984
#'
#' NB: The last 12 lines of the data set have the first column set to
#' `NULL` since these data where not observed but we want to predict them.
#'
#'
#' @name Drivers
#' @docType data
#' @format A data frame with 204 observations on the following 4 variables.
#' \describe{
#' \item{y}{Number of deaths}
#' \item{belt}{Indicator of weather the belt was compulsory to use
#' (1) or not (0)}
#' \item{trend}{time (in months)}
#' \item{seasonal}{time (in months)}
#' }
#' @references Rue, H and Held, L. (2005) *Gaussian Markov Random Fields -
#' Theory and Applications* Chapman and Hall
#' @keywords datasets
#' @examples
#'
#' data(Drivers)
NULL





#' Repeated measures on Poisson counts
#'
#' Seizure counts in a randomised trial of anti-convulsant therpay in epilepsy
#' for 59 patients.
#'
#'
#' @name Epil
#' @docType data
#' @format A data frame with 236 observations on the following 7 variables.
#' \describe{
#' \item{y}{Number of seizures}
#' \item{Trt}{indicator for the presence of treatment}
#' \item{Base}{8-week baseline seizure counts}
#' \item{Age}{Age of the patient}
#' \item{V4}{indicator variable for the 4th visit.}
#' \item{rand}{a numeric vector}
#' \item{Ind}{indicator for the specific patient}
#' }
#' @source WinBUGS/OpenBUGS Manual Examples Vol I
#' @keywords datasets
#' @examples
#'
#' data(Epil)
NULL





#' Disease Mapping
#'
#' Cases of Oral cavity cancer in Germany from 1986-1990
#'
#'
#' @name Germany
#' @docType data
#' @format A data frame with 544 observations on the following 4 variables.
#' \describe{
#' \item{region}{Region of Germany}
#' \item{E}{Fixed quantity which accounts for number of people
#' in the district (offset)}
#' \item{Y}{Number of cases}
#' \item{x}{covariate measuring smoking consumption}
#' }
#' @references Rue, H and Held, L. (2005) *Gaussian Markov Random Fields -
#' Theory and Applications* Chapman and Hall
#' @keywords datasets
#' @examples
#'
#' data(Germany)
NULL










#' Kidney infection data
#'
#' Times of infection from the time to insertion of the catheter for 38 kindey
#' patients using portable dialysis equipment
#'
#'
#' @name Kidney
#' @docType data
#' @format A data frame with 76 observations on the following 9 variables.
#' \describe{
#' \item{time}{a numeric vector. Time to infection from the
#' insertion of catheter}
#' \item{event}{a numeric vector. 1: time of
#' infection 0: time of censuring }
#' \item{age}{a numeric vector. Age of
#' the patient at the time of infection}
#' \item{sex}{a numeric vector.
#' Sex of the patient 0: male 1:female}
#' \item{disease}{a numeric
#' vector. Type of disease}
#' \item{dis1}{a numeric vector. Dummy
#' variable to codify the disease type.}
#' \item{dis2}{a numeric vector.
#' Dummy variable to codify the disease type.}
#' \item{dis3}{a numeric
#' vector. Dummy variable to codify the disease type.}
#' \item{ID}{a
#' numeric vector. Patient code.}
#' }
#' @references McGilchrist and C.W. Aisbett (1991), Regression with frailty in
#' survival analysis, Biometrics,vol.47,pages 461--166.
#'
#' D.J. Spiegelhalter and A. Thomas and N.G. Best and W.R. Gilks (1995) BUGS:
#' Bayesian Inference Using Gibbs sampling, Version 0.50., MRC Biostatistics
#' Unit, Cambridre, England.
#' @keywords datasets
NULL





#' The Leukemia data
#'
#' This the Leukemia data from Henderson et al (2003); see source.
#'
#'
#' @name Leuk
#' @aliases Leukemia Leuk
#' @docType data
#' @format A data frame with 1043 observations on the following 9 variables.
#' \describe{
#' \item{time}{TODO}
#' \item{cens}{TODO}
#' \item{xcoord}{TODO}
#' \item{ycoord}{TODO}
#' \item{age}{TODO}
#' \item{sex}{TODO}
#' \item{wbc}{TODO}
#' \item{tpi}{TODO}
#' \item{district}{TODO}
#' }
#' @source This is the dataset from
#'
#' Henderson, R. and Shimakura, S. and Gorst, D., 2002, Modeling spatial
#' variation in leukemia survival data, JASA, 97, 460, 965--972.
#' @keywords datasets
#' @examples
#'
#' data(Leuk)
NULL





#' The Munich rent data
#'
#' The Munich rent data
#'
#'
#' @name Munich
#' @docType data
#' @format A data frame with 2035 observations on the following 17 variables.
#' \describe{
#'
#' \item{rent}{Net rent per square meter.}
#' \item{floor.size}{Size of the flat in square meters.}
#' \item{year}{Year of construction of the building in which the flat
#' is located.} \item{location}{Location index (in terms of
#' subquarters).} \item{Gute.Wohnlage}{Dummy variable for good
#' locations / good neighborhoods.} \item{Beste.Wohnlage}{Dummy
#' variable for very good locations / very good neighborhoods.}
#' \item{Keine.Wwv}{Dummy for absence of warm water supply.}
#' \item{Keine.Zh}{Dummy for absence of central heating system.}
#' \item{Kein.Badkach}{Dummy for absence of flagging in the bathroom.}
#' \item{Besond.Bad}{Dummy for special features of the bathroom.}
#' \item{Gehobene.Kueche}{Dummy for more refined kitchen equipment.}
#' \item{zim1}{Dummy for a flat with 1 room.} \item{zim2}{Dummy
#' for a flat with 2 rooms.} \item{zim3}{Dummy for a flat with 3
#' rooms.} \item{zim4}{Dummy for a flat with 4 rooms.}
#' \item{zim5}{Dummy for a flat with 5 rooms.}
#' \item{zim6}{Dummy for a flat with 6 rooms.}
#'
#' }
#' @references Rue, H and Held, L. (2005) *Gaussian Markov Random Fields -
#' Theory and Applications* Chapman and Hall
#' @source See Rue and Held (2005), Chapter 4.
#' @keywords datasets
NULL





#' The North West England map
#'
#' This map is used in association to the Leukemia data from Henderson et al
#' (2003); see source.
#'
#'
#' @name nwEngland
#' @aliases nwEngland
#' @docType data
#' @format A SpatialPolygons object.
#' @source This map are used to analyse the Leukaemia dataset from
#'
#' Henderson, R. and Shimakura, S. and Gorst, D., 2002, Modeling spatial
#' variation in leukemia survival data, JASA, 97, 460, 965--972.
#' @keywords datasets
#' @examples
#'
#' data(Leuk)
#' plot(nwEngland)
NULL





#' ~~ data name/kind ... ~~
#'
#' ~~ A concise (1-5 lines) description of the dataset. ~~
#'
#'
#' @name Oral
#' @docType data
#' @format A data frame with 544 observations on the following 3 variables.
#' \describe{ \item{region}{a numeric vector} \item{E}{a
#' numeric vector} \item{Y}{a numeric vector} }
#' @references Rue, H and Held, L. (2005) *Gaussian Markov Random Fields -
#' Theory and Applications* Chapman and Hall
#' @keywords datasets
NULL





#' The PRborder data
#'
#' A data matrix with Longitude and Latitude coordinates for the boundary of
#' Parana State.
#'
#'
#' @name PRborder
#' @aliases PRborder PRborder
#' @docType data
#' @format \describe{ \item{Longtiude}{The Longtiude coordinate}
#' \item{Latitude}{The Latitude coordinate} }
#' @seealso PRprec
#' @keywords datasets
NULL





#' The PRprec data
#'
#' A data frame with daily rainfall in the Parana State.
#'
#'
#' @name PRprec
#' @aliases PRprec PRprec
#' @docType data
#' @format A data frame .... TODO
#' \describe{
#' \item{Altitude}{TODO}
#' \item{Latitude}{TODO}
#' \item{Longitude}{TODO}
#' \item{d0101}{Daily rainfall at day "mmdd"}
#' \item{d0102}{Daily rainfall at day "mmdd"}
#' \item{d0103}{Daily rainfall at day "mmdd"}
#' \item{d0104}{Daily rainfall at day "mmdd"}
#' \item{d0105}{Daily rainfall at day "mmdd"}
#' \item{d0106}{Daily rainfall at day "mmdd"}
#' \item{d0107}{Daily rainfall at day "mmdd"} \item{d0108}{Daily
#' rainfall at day "mmdd"} \item{d0109}{Daily rainfall at day "mmdd"}
#' \item{d0110}{Daily rainfall at day "mmdd"} \item{d0111}{Daily rainfall at
#' day "mmdd"} \item{d0112}{Daily rainfall at day "mmdd"} \item{d0113}{Daily
#' rainfall at day "mmdd"} \item{d0114}{Daily rainfall at day "mmdd"}
#' \item{d0115}{Daily rainfall at day "mmdd"} \item{d0116}{Daily rainfall at
#' day "mmdd"} \item{d0117}{Daily rainfall at day "mmdd"} \item{d0118}{Daily
#' rainfall at day "mmdd"} \item{d0119}{Daily rainfall at day "mmdd"}
#' \item{d0120}{Daily rainfall at day "mmdd"} \item{d0121}{Daily rainfall at
#' day "mmdd"} \item{d0122}{Daily rainfall at day "mmdd"} \item{d0123}{Daily
#' rainfall at day "mmdd"} \item{d0124}{Daily rainfall at day "mmdd"}
#' \item{d0125}{Daily rainfall at day "mmdd"} \item{d0126}{Daily rainfall at
#' day "mmdd"} \item{d0127}{Daily rainfall at day "mmdd"} \item{d0128}{Daily
#' rainfall at day "mmdd"} \item{d0129}{Daily rainfall at day "mmdd"}
#' \item{d0130}{Daily rainfall at day "mmdd"} \item{d0131}{Daily rainfall at
#' day "mmdd"} \item{d0201}{Daily rainfall at day "mmdd"} \item{d0202}{Daily
#' rainfall at day "mmdd"} \item{d0203}{Daily rainfall at day "mmdd"}
#' \item{d0204}{Daily rainfall at day "mmdd"} \item{d0205}{Daily rainfall at
#' day "mmdd"} \item{d0206}{Daily rainfall at day "mmdd"} \item{d0207}{Daily
#' rainfall at day "mmdd"} \item{d0208}{Daily rainfall at day "mmdd"}
#' \item{d0209}{Daily rainfall at day "mmdd"} \item{d0210}{Daily rainfall at
#' day "mmdd"} \item{d0211}{Daily rainfall at day "mmdd"} \item{d0212}{Daily
#' rainfall at day "mmdd"} \item{d0213}{Daily rainfall at day "mmdd"}
#' \item{d0214}{Daily rainfall at day "mmdd"} \item{d0215}{Daily rainfall at
#' day "mmdd"} \item{d0216}{Daily rainfall at day "mmdd"} \item{d0217}{Daily
#' rainfall at day "mmdd"} \item{d0218}{Daily rainfall at day "mmdd"}
#' \item{d0219}{Daily rainfall at day "mmdd"} \item{d0220}{Daily rainfall at
#' day "mmdd"} \item{d0221}{Daily rainfall at day "mmdd"} \item{d0222}{Daily
#' rainfall at day "mmdd"} \item{d0223}{Daily rainfall at day "mmdd"}
#' \item{d0224}{Daily rainfall at day "mmdd"} \item{d0225}{Daily rainfall at
#' day "mmdd"} \item{d0226}{Daily rainfall at day "mmdd"} \item{d0227}{Daily
#' rainfall at day "mmdd"} \item{d0228}{Daily rainfall at day "mmdd"}
#' \item{d0301}{Daily rainfall at day "mmdd"} \item{d0302}{Daily rainfall at
#' day "mmdd"} \item{d0303}{Daily rainfall at day "mmdd"} \item{d0304}{Daily
#' rainfall at day "mmdd"} \item{d0305}{Daily rainfall at day "mmdd"}
#' \item{d0306}{Daily rainfall at day "mmdd"} \item{d0307}{Daily rainfall at
#' day "mmdd"} \item{d0308}{Daily rainfall at day "mmdd"} \item{d0309}{Daily
#' rainfall at day "mmdd"} \item{d0310}{Daily rainfall at day "mmdd"}
#' \item{d0311}{Daily rainfall at day "mmdd"} \item{d0312}{Daily rainfall at
#' day "mmdd"} \item{d0313}{Daily rainfall at day "mmdd"} \item{d0314}{Daily
#' rainfall at day "mmdd"} \item{d0315}{Daily rainfall at day "mmdd"}
#' \item{d0316}{Daily rainfall at day "mmdd"} \item{d0317}{Daily rainfall at
#' day "mmdd"} \item{d0318}{Daily rainfall at day "mmdd"} \item{d0319}{Daily
#' rainfall at day "mmdd"} \item{d0320}{Daily rainfall at day "mmdd"}
#' \item{d0321}{Daily rainfall at day "mmdd"} \item{d0322}{Daily rainfall at
#' day "mmdd"} \item{d0323}{Daily rainfall at day "mmdd"} \item{d0324}{Daily
#' rainfall at day "mmdd"} \item{d0325}{Daily rainfall at day "mmdd"}
#' \item{d0326}{Daily rainfall at day "mmdd"} \item{d0327}{Daily rainfall at
#' day "mmdd"} \item{d0328}{Daily rainfall at day "mmdd"} \item{d0329}{Daily
#' rainfall at day "mmdd"} \item{d0330}{Daily rainfall at day "mmdd"}
#' \item{d0331}{Daily rainfall at day "mmdd"} \item{d0401}{Daily rainfall at
#' day "mmdd"} \item{d0402}{Daily rainfall at day "mmdd"} \item{d0403}{Daily
#' rainfall at day "mmdd"} \item{d0404}{Daily rainfall at day "mmdd"}
#' \item{d0405}{Daily rainfall at day "mmdd"} \item{d0406}{Daily rainfall at
#' day "mmdd"} \item{d0407}{Daily rainfall at day "mmdd"} \item{d0408}{Daily
#' rainfall at day "mmdd"} \item{d0409}{Daily rainfall at day "mmdd"}
#' \item{d0410}{Daily rainfall at day "mmdd"} \item{d0411}{Daily rainfall at
#' day "mmdd"} \item{d0412}{Daily rainfall at day "mmdd"} \item{d0413}{Daily
#' rainfall at day "mmdd"} \item{d0414}{Daily rainfall at day "mmdd"}
#' \item{d0415}{Daily rainfall at day "mmdd"} \item{d0416}{Daily rainfall at
#' day "mmdd"} \item{d0417}{Daily rainfall at day "mmdd"} \item{d0418}{Daily
#' rainfall at day "mmdd"} \item{d0419}{Daily rainfall at day "mmdd"}
#' \item{d0420}{Daily rainfall at day "mmdd"} \item{d0421}{Daily rainfall at
#' day "mmdd"} \item{d0422}{Daily rainfall at day "mmdd"} \item{d0423}{Daily
#' rainfall at day "mmdd"} \item{d0424}{Daily rainfall at day "mmdd"}
#' \item{d0425}{Daily rainfall at day "mmdd"} \item{d0426}{Daily rainfall at
#' day "mmdd"} \item{d0427}{Daily rainfall at day "mmdd"} \item{d0428}{Daily
#' rainfall at day "mmdd"} \item{d0429}{Daily rainfall at day "mmdd"}
#' \item{d0430}{Daily rainfall at day "mmdd"} \item{d0501}{Daily rainfall at
#' day "mmdd"} \item{d0502}{Daily rainfall at day "mmdd"} \item{d0503}{Daily
#' rainfall at day "mmdd"} \item{d0504}{Daily rainfall at day "mmdd"}
#' \item{d0505}{Daily rainfall at day "mmdd"} \item{d0506}{Daily rainfall at
#' day "mmdd"} \item{d0507}{Daily rainfall at day "mmdd"} \item{d0508}{Daily
#' rainfall at day "mmdd"} \item{d0509}{Daily rainfall at day "mmdd"}
#' \item{d0510}{Daily rainfall at day "mmdd"} \item{d0511}{Daily rainfall at
#' day "mmdd"} \item{d0512}{Daily rainfall at day "mmdd"} \item{d0513}{Daily
#' rainfall at day "mmdd"} \item{d0514}{Daily rainfall at day "mmdd"}
#' \item{d0515}{Daily rainfall at day "mmdd"} \item{d0516}{Daily rainfall at
#' day "mmdd"} \item{d0517}{Daily rainfall at day "mmdd"} \item{d0518}{Daily
#' rainfall at day "mmdd"} \item{d0519}{Daily rainfall at day "mmdd"}
#' \item{d0520}{Daily rainfall at day "mmdd"} \item{d0521}{Daily rainfall at
#' day "mmdd"} \item{d0522}{Daily rainfall at day "mmdd"} \item{d0523}{Daily
#' rainfall at day "mmdd"} \item{d0524}{Daily rainfall at day "mmdd"}
#' \item{d0525}{Daily rainfall at day "mmdd"} \item{d0526}{Daily rainfall at
#' day "mmdd"} \item{d0527}{Daily rainfall at day "mmdd"} \item{d0528}{Daily
#' rainfall at day "mmdd"} \item{d0529}{Daily rainfall at day "mmdd"}
#' \item{d0530}{Daily rainfall at day "mmdd"} \item{d0531}{Daily rainfall at
#' day "mmdd"} \item{d0601}{Daily rainfall at day "mmdd"} \item{d0602}{Daily
#' rainfall at day "mmdd"} \item{d0603}{Daily rainfall at day "mmdd"}
#' \item{d0604}{Daily rainfall at day "mmdd"} \item{d0605}{Daily rainfall at
#' day "mmdd"} \item{d0606}{Daily rainfall at day "mmdd"} \item{d0607}{Daily
#' rainfall at day "mmdd"} \item{d0608}{Daily rainfall at day "mmdd"}
#' \item{d0609}{Daily rainfall at day "mmdd"} \item{d0610}{Daily rainfall at
#' day "mmdd"} \item{d0611}{Daily rainfall at day "mmdd"} \item{d0612}{Daily
#' rainfall at day "mmdd"} \item{d0613}{Daily rainfall at day "mmdd"}
#' \item{d0614}{Daily rainfall at day "mmdd"} \item{d0615}{Daily rainfall at
#' day "mmdd"} \item{d0616}{Daily rainfall at day "mmdd"} \item{d0617}{Daily
#' rainfall at day "mmdd"} \item{d0618}{Daily rainfall at day "mmdd"}
#' \item{d0619}{Daily rainfall at day "mmdd"} \item{d0620}{Daily rainfall at
#' day "mmdd"} \item{d0621}{Daily rainfall at day "mmdd"} \item{d0622}{Daily
#' rainfall at day "mmdd"} \item{d0623}{Daily rainfall at day "mmdd"}
#' \item{d0624}{Daily rainfall at day "mmdd"} \item{d0625}{Daily rainfall at
#' day "mmdd"} \item{d0626}{Daily rainfall at day "mmdd"} \item{d0627}{Daily
#' rainfall at day "mmdd"} \item{d0628}{Daily rainfall at day "mmdd"}
#' \item{d0629}{Daily rainfall at day "mmdd"} \item{d0630}{Daily rainfall at
#' day "mmdd"} \item{d0701}{Daily rainfall at day "mmdd"} \item{d0702}{Daily
#' rainfall at day "mmdd"} \item{d0703}{Daily rainfall at day "mmdd"}
#' \item{d0704}{Daily rainfall at day "mmdd"} \item{d0705}{Daily rainfall at
#' day "mmdd"} \item{d0706}{Daily rainfall at day "mmdd"} \item{d0707}{Daily
#' rainfall at day "mmdd"} \item{d0708}{Daily rainfall at day "mmdd"}
#' \item{d0709}{Daily rainfall at day "mmdd"} \item{d0710}{Daily rainfall at
#' day "mmdd"} \item{d0711}{Daily rainfall at day "mmdd"} \item{d0712}{Daily
#' rainfall at day "mmdd"} \item{d0713}{Daily rainfall at day "mmdd"}
#' \item{d0714}{Daily rainfall at day "mmdd"} \item{d0715}{Daily rainfall at
#' day "mmdd"} \item{d0716}{Daily rainfall at day "mmdd"} \item{d0717}{Daily
#' rainfall at day "mmdd"} \item{d0718}{Daily rainfall at day "mmdd"}
#' \item{d0719}{Daily rainfall at day "mmdd"} \item{d0720}{Daily rainfall at
#' day "mmdd"} \item{d0721}{Daily rainfall at day "mmdd"} \item{d0722}{Daily
#' rainfall at day "mmdd"} \item{d0723}{Daily rainfall at day "mmdd"}
#' \item{d0724}{Daily rainfall at day "mmdd"} \item{d0725}{Daily rainfall at
#' day "mmdd"} \item{d0726}{Daily rainfall at day "mmdd"} \item{d0727}{Daily
#' rainfall at day "mmdd"} \item{d0728}{Daily rainfall at day "mmdd"}
#' \item{d0729}{Daily rainfall at day "mmdd"} \item{d0730}{Daily rainfall at
#' day "mmdd"} \item{d0731}{Daily rainfall at day "mmdd"} \item{d0801}{Daily
#' rainfall at day "mmdd"} \item{d0802}{Daily rainfall at day "mmdd"}
#' \item{d0803}{Daily rainfall at day "mmdd"} \item{d0804}{Daily rainfall at
#' day "mmdd"} \item{d0805}{Daily rainfall at day "mmdd"} \item{d0806}{Daily
#' rainfall at day "mmdd"} \item{d0807}{Daily rainfall at day "mmdd"}
#' \item{d0808}{Daily rainfall at day "mmdd"} \item{d0809}{Daily rainfall at
#' day "mmdd"} \item{d0810}{Daily rainfall at day "mmdd"} \item{d0811}{Daily
#' rainfall at day "mmdd"} \item{d0812}{Daily rainfall at day "mmdd"}
#' \item{d0813}{Daily rainfall at day "mmdd"} \item{d0814}{Daily rainfall at
#' day "mmdd"} \item{d0815}{Daily rainfall at day "mmdd"} \item{d0816}{Daily
#' rainfall at day "mmdd"} \item{d0817}{Daily rainfall at day "mmdd"}
#' \item{d0818}{Daily rainfall at day "mmdd"} \item{d0819}{Daily rainfall at
#' day "mmdd"} \item{d0820}{Daily rainfall at day "mmdd"} \item{d0821}{Daily
#' rainfall at day "mmdd"} \item{d0822}{Daily rainfall at day "mmdd"}
#' \item{d0823}{Daily rainfall at day "mmdd"} \item{d0824}{Daily rainfall at
#' day "mmdd"} \item{d0825}{Daily rainfall at day "mmdd"} \item{d0826}{Daily
#' rainfall at day "mmdd"} \item{d0827}{Daily rainfall at day "mmdd"}
#' \item{d0828}{Daily rainfall at day "mmdd"} \item{d0829}{Daily rainfall at
#' day "mmdd"} \item{d0830}{Daily rainfall at day "mmdd"} \item{d0831}{Daily
#' rainfall at day "mmdd"} \item{d0901}{Daily rainfall at day "mmdd"}
#' \item{d0902}{Daily rainfall at day "mmdd"} \item{d0903}{Daily rainfall at
#' day "mmdd"} \item{d0904}{Daily rainfall at day "mmdd"} \item{d0905}{Daily
#' rainfall at day "mmdd"} \item{d0906}{Daily rainfall at day "mmdd"}
#' \item{d0907}{Daily rainfall at day "mmdd"} \item{d0908}{Daily rainfall at
#' day "mmdd"} \item{d0909}{Daily rainfall at day "mmdd"} \item{d0910}{Daily
#' rainfall at day "mmdd"} \item{d0911}{Daily rainfall at day "mmdd"}
#' \item{d0912}{Daily rainfall at day "mmdd"} \item{d0913}{Daily rainfall at
#' day "mmdd"} \item{d0914}{Daily rainfall at day "mmdd"} \item{d0915}{Daily
#' rainfall at day "mmdd"} \item{d0916}{Daily rainfall at day "mmdd"}
#' \item{d0917}{Daily rainfall at day "mmdd"} \item{d0918}{Daily rainfall at
#' day "mmdd"} \item{d0919}{Daily rainfall at day "mmdd"} \item{d0920}{Daily
#' rainfall at day "mmdd"} \item{d0921}{Daily rainfall at day "mmdd"}
#' \item{d0922}{Daily rainfall at day "mmdd"} \item{d0923}{Daily rainfall at
#' day "mmdd"} \item{d0924}{Daily rainfall at day "mmdd"} \item{d0925}{Daily
#' rainfall at day "mmdd"} \item{d0926}{Daily rainfall at day "mmdd"}
#' \item{d0927}{Daily rainfall at day "mmdd"} \item{d0928}{Daily rainfall at
#' day "mmdd"} \item{d0929}{Daily rainfall at day "mmdd"} \item{d0930}{Daily
#' rainfall at day "mmdd"} \item{d1001}{Daily rainfall at day "mmdd"}
#' \item{d1002}{Daily rainfall at day "mmdd"} \item{d1003}{Daily rainfall at
#' day "mmdd"} \item{d1004}{Daily rainfall at day "mmdd"} \item{d1005}{Daily
#' rainfall at day "mmdd"} \item{d1006}{Daily rainfall at day "mmdd"}
#' \item{d1007}{Daily rainfall at day "mmdd"} \item{d1008}{Daily rainfall at
#' day "mmdd"} \item{d1009}{Daily rainfall at day "mmdd"} \item{d1010}{Daily
#' rainfall at day "mmdd"} \item{d1011}{Daily rainfall at day "mmdd"}
#' \item{d1012}{Daily rainfall at day "mmdd"} \item{d1013}{Daily rainfall at
#' day "mmdd"} \item{d1014}{Daily rainfall at day "mmdd"} \item{d1015}{Daily
#' rainfall at day "mmdd"} \item{d1016}{Daily rainfall at day "mmdd"}
#' \item{d1017}{Daily rainfall at day "mmdd"} \item{d1018}{Daily rainfall at
#' day "mmdd"} \item{d1019}{Daily rainfall at day "mmdd"} \item{d1020}{Daily
#' rainfall at day "mmdd"} \item{d1021}{Daily rainfall at day "mmdd"}
#' \item{d1022}{Daily rainfall at day "mmdd"} \item{d1023}{Daily rainfall at
#' day "mmdd"} \item{d1024}{Daily rainfall at day "mmdd"} \item{d1025}{Daily
#' rainfall at day "mmdd"} \item{d1026}{Daily rainfall at day "mmdd"}
#' \item{d1027}{Daily rainfall at day "mmdd"} \item{d1028}{Daily rainfall at
#' day "mmdd"} \item{d1029}{Daily rainfall at day "mmdd"} \item{d1030}{Daily
#' rainfall at day "mmdd"} \item{d1031}{Daily rainfall at day "mmdd"}
#' \item{d1101}{Daily rainfall at day "mmdd"} \item{d1102}{Daily rainfall at
#' day "mmdd"} \item{d1103}{Daily rainfall at day "mmdd"} \item{d1104}{Daily
#' rainfall at day "mmdd"} \item{d1105}{Daily rainfall at day "mmdd"}
#' \item{d1106}{Daily rainfall at day "mmdd"} \item{d1107}{Daily rainfall at
#' day "mmdd"} \item{d1108}{Daily rainfall at day "mmdd"} \item{d1109}{Daily
#' rainfall at day "mmdd"} \item{d1110}{Daily rainfall at day "mmdd"}
#' \item{d1111}{Daily rainfall at day "mmdd"} \item{d1112}{Daily rainfall at
#' day "mmdd"} \item{d1113}{Daily rainfall at day "mmdd"} \item{d1114}{Daily
#' rainfall at day "mmdd"} \item{d1115}{Daily rainfall at day "mmdd"}
#' \item{d1116}{Daily rainfall at day "mmdd"} \item{d1117}{Daily rainfall at
#' day "mmdd"} \item{d1118}{Daily rainfall at day "mmdd"} \item{d1119}{Daily
#' rainfall at day "mmdd"} \item{d1120}{Daily rainfall at day "mmdd"}
#' \item{d1121}{Daily rainfall at day "mmdd"} \item{d1122}{Daily rainfall at
#' day "mmdd"} \item{d1123}{Daily rainfall at day "mmdd"} \item{d1124}{Daily
#' rainfall at day "mmdd"} \item{d1125}{Daily rainfall at day "mmdd"}
#' \item{d1126}{Daily rainfall at day "mmdd"} \item{d1127}{Daily rainfall at
#' day "mmdd"} \item{d1128}{Daily rainfall at day "mmdd"} \item{d1129}{Daily
#' rainfall at day "mmdd"} \item{d1130}{Daily rainfall at day "mmdd"}
#' \item{d1201}{Daily rainfall at day "mmdd"} \item{d1202}{Daily rainfall at
#' day "mmdd"} \item{d1203}{Daily rainfall at day "mmdd"} \item{d1204}{Daily
#' rainfall at day "mmdd"} \item{d1205}{Daily rainfall at day "mmdd"}
#' \item{d1206}{Daily rainfall at day "mmdd"} \item{d1207}{Daily rainfall at
#' day "mmdd"} \item{d1208}{Daily rainfall at day "mmdd"} \item{d1209}{Daily
#' rainfall at day "mmdd"} \item{d1210}{Daily rainfall at day "mmdd"}
#' \item{d1211}{Daily rainfall at day "mmdd"} \item{d1212}{Daily rainfall at
#' day "mmdd"} \item{d1213}{Daily rainfall at day "mmdd"} \item{d1214}{Daily
#' rainfall at day "mmdd"} \item{d1215}{Daily rainfall at day "mmdd"}
#' \item{d1216}{Daily rainfall at day "mmdd"} \item{d1217}{Daily rainfall at
#' day "mmdd"} \item{d1218}{Daily rainfall at day "mmdd"} \item{d1219}{Daily
#' rainfall at day "mmdd"} \item{d1220}{Daily rainfall at day "mmdd"}
#' \item{d1221}{Daily rainfall at day "mmdd"} \item{d1222}{Daily rainfall at
#' day "mmdd"} \item{d1223}{Daily rainfall at day "mmdd"} \item{d1224}{Daily
#' rainfall at day "mmdd"} \item{d1225}{Daily rainfall at day "mmdd"}
#' \item{d1226}{Daily rainfall at day "mmdd"} \item{d1227}{Daily rainfall at
#' day "mmdd"} \item{d1228}{Daily rainfall at day "mmdd"} \item{d1229}{Daily
#' rainfall at day "mmdd"} \item{d1230}{Daily rainfall at day "mmdd"}
#' \item{d1231}{Daily rainfall at day "mmdd"} }
#' @seealso PRborder
#' @keywords datasets
NULL





#' Extra-Poisson variation in dose-response study
#'
#' Breslow (1984) analyses some mutagenicity assay data (shown below) on
#' salmonella in which three plates have been processed at each dose i of
#' quinoline and the number of revertant colonies of TA98 Salmonella measured
#'
#'
#' @name Salm
#' @docType data
#' @format A data frame with 18 observations on the following 3 variables.
#' \describe{ \item{y}{number of salmonella bacteria}
#' \item{dose}{dose of quinoline (mg per plate)}
#' \item{rand}{indicator} }
#' @source WinBUGS/OpenBUGS manual Examples VOl.I
#' @keywords datasets
#' @examples
#'
#' data(Salm)
NULL





#' Conditional Autoregressive (CAR) model for disease mapping
#'
#'
#' The rate of lip cancer in 56 counties in Scotland is recorder. The data set
#' includes the observed and expected cases (based on the population and its
#' age and sex distribution in the country), a covariate measuring the
#' percentage of the population engaged in agricolture, fishing or forestry and
#' the "position" of each county expressed as a list of adjacent counties
#'
#'
#' @name Scotland
#' @docType data
#' @format A data frame with 56 observations on the following 4 variables.
#' \describe{ \item{Counts}{The number of lip cancer registered}
#' \item{E}{The expected number of lip cancer } \item{X}{The
#' percentage of the population engaged in agricolture, fishing or forestry }
#' \item{Region}{The county} }
#' @references
#'
#' Clayton and Kaldor (1987) and Breslow and Clayton (1993)
#' @source OpenBUGS Example manual, GeoBUGS
#' @keywords datasets
#' @examples
#'
#' data(Scotland)
NULL





#' Factorial design
#'
#'
#' Proportion of seeds that germinated on each of 21 plates arranged according
#' to a 2 by 2 factorial layout by seed and type of root extract
#'
#'
#' @name Seeds
#' @docType data
#' @format A data frame with 21 observations on the following 5 variables.
#' \describe{ \item{r}{number of germinated seeds per plate}
#' \item{n}{number of total seeds per plate} \item{x1}{seed
#' type} \item{x2}{root extracted} \item{plate}{indicator for
#' the plate} }
#' @source WinBUGS/OpenBUGS Manual Example, Vol. I
#' @keywords datasets
#' @examples
#'
#' data(Seeds)
NULL





#' toy simulated data set for the SPDE tutorial
#'
#' Simulated data set on 200 location points.  The simulation process is made
#' at the introduction of the SPDE tutorial.
#'
#'
#' @name SPDEtoy
#' @docType data
#' @format A data frame with 200 observations on the following 3 variables.
#' \describe{ \item{s1}{First element of the coordinates}
#' \item{s2}{Second element of the coordinates} \item{y}{data
#' simulated at the locations} }
#' @source SPDE tutorial
#' @keywords datasets
#' @examples
#'
#' data(SPDEtoy)
NULL





#' Surgical: Institutional ranking
#'
#' This example considers mortality rates in 12 hospitals performing cardiac
#' surgery in babies
#'
#'
#' @name Surg
#' @docType data
#' @format A data frame with 12 observations on the following 3 variables.
#' \describe{ \item{n}{Number of deaths} \item{r}{Total number
#' of cases} \item{hospital}{a factor with levels `A` `B`
#' `C` `D` `E` `F` `G` `H` `I` `J`
#' `K` `L`} }
#' @source WinBUGS/OpenBUGS Manual Examples Vol. I
#' @keywords datasets
#' @examples
#'
#' data(Surg)
NULL





#' Survival data
#'
#' Simulated data set for Weibull survival model
#'
#'
#' @name SurvSim
#' @docType data
#' @format A data frame with 100 observations on the following 3 variables.
#' \describe{ \item{y}{a numeric vector of survival times}
#' \item{cens}{a numeric vector of event indicator (0=censured
#' 1=failure)} \item{x}{a numeric vector of covariate } }
#' @keywords datasets
NULL





#' Binomial time series
#'
#' Recorded days of rain above 1 mm in Tokyo for 2 years, 1983:84
#'
#'
#' @name Tokyo
#' @docType data
#' @format A data frame with 366 observations on the following 3 variables.
#' \describe{
#' \item{y}{number of days with rain}
#' \item{n}{total number of days}
#' \item{time}{day of the year}
#' }
#' @references Rue, H and Held, L. (2005) *Gaussian Markov Random Fields -
#' Theory and Applications* Chapman and Hall
#' @source
#' <http://www.math.ntnu.no/~hrue/GMRF-book/tokyo.rainfall.data.dat>
#' @keywords datasets
#' @examples
#'
#' data(Tokyo)
NULL





#' Semiparametric regression
#'
#' Undernutrition of children in each region of Zambia is measured through a
#' score computed on the basis of some anthropometric measures. The data set
#' contains also other infomation about each child.
#'
#'
#' @name Zambia
#' @docType data
#' @format A data frame with 4847 observations on the following 10 variables.
#' \describe{ \item{hazstd}{standardised Z score of stunting}
#' \item{bmi}{body mass index of the mother} \item{agc}{age of
#' the child in months} \item{district}{district where the child lives}
#' \item{rcw}{mother employment status with categories "working" (1)
#' and "not working" (-1)} \item{edu1}{mother's educations status with
#' categories "complete primary but incomplete secondary " (`edu1=1)`,
#' "complete secondary or higher" (`edu2=1`) and "no education or
#' incomplete primary" (`edu1=edu2=-1`)} \item{edu2}{see above}
#' \item{tpr}{locality of the domicile with categories "urban" (1) and
#' "rural" (-1)} \item{sex}{gender of the child with categories "male"
#' (1) and "female" (-1)} \item{edu}{DO NOT KNOW; check source} }
#' @source BayesX Manual
#' <http://www.stat.uni-muenchen.de/~bayesx/bayesx.html>
#' @keywords datasets
#' @examples
#'
#' data(Zambia)
NULL
