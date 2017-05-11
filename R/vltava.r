#' @title Forest vegetation in deep river valley of Vltava (Czech Republic)
#' 
#' @description Datasets with species composition of forest vegetation, their species Ellenberg indicator values, species functional traits (compiled from databases only for herb species) and several measured environmental characteristics for each site. 
#' 
#' @details
#' 
#' Vegetation plots, located at even distances along transects following the steep valley slopes of Vltava river valley and collected during 2001-2003. Each transect started at the valley bottom and end up at the upper part of the valley slope. Plots are of the size 10x15 m. In each plot, all species of tree, shrub and herb layer were recorded and their abundances were estimated using 9-degree ordinal Braun-Blanquette scale (these values were consequently transformed into percentage values). At each plot, various topographical and soil factors were measured or estimated (see Table below). The dataset contains 27 transects with 97 samples.
#' 
#' For the purpose of the current dataset, species in shrub and tree layer have been merged, juveniles removed and nomenclature have been modified according to Kubat et al. (2002). Dataset has two parts: with all (tree, shrub and herb) species (\code{vltava$spe, $ell, $env} etc.) and with subset of only herb species (\code{vltava$herbs$spe, $ell, $traits etc.}). While Ellenberg indicator values are provided for both all and only herb species subset, plant functional traits are only for subset of herb species (it would perhaps not be meaningful to compare e.g. SLA or plant height for trees and herbs). Some species indicator values and trait values remain unassigned, because they are not defined in either Ellenberg or LEDA database. Exception is SLA value for \emph{Impatiens parviflora}, which in LEDA database has assigned value 148.01 mm2/mg, which seems as outlier and had been temporarily removed.
#' 
#' 
#' Environmental variables include:
#' \itemize{
#' \item ELEVATION  elevation [m a.s.l.]
#' \item SLOPE	inclination [degrees]
#' \item ASPSSW	aspect (deviation of plot aspect from 22.5, reaching the highest values in SSW orientation)
#' \item HEAT.LOAD	heat load, calculated from plot slope and aspect (McCune & Keon 2002), symetric along SW-NE axis
#' \item SURFSL	landform shape in the downslope direction (three-degree ordinal scale: -1 concave, 0 flat, 1 convex)
#' \item SURFIS	landform shape along the isohypse (three-degree ordinal scale: -1 concave, 0 flat, 1 convex)
#' \item LITHIC	lithic leptosols (shallow soils near rock outcrops)
#' \item SKELETIC	skeletic and hyperskeletic leptosols (stony soils on scree accumulatons)
#' \item CAMBISOL	cambisols (well-developed zonal soils)
#' \item FLUVISOL	fluvisols (water-influenced soils formed from alluvial deposits)
#' \item SOILDPT	soil depth [cm], measured by 0.7 m long iron rod (1.5 cm diameter) - average of 5 values measured in 5 places within the plot
#' \item pH	soil pH (measured in water solution)
#' \item COVERE3	estimated cover of tree layer [\%]
#' \item COVERE2	estimated cover of shrub layer [\%]
#' \item COVERE1	estimated cover of herb layer [\%]
#' \item COVERE0	estimated cover of moss layer [\%]
#' \item COVERE32	estimated cover of tree and shrub layer [\%] (merged tree and shrub estimations, using formula p.tree + p.shrub - p.tree*p.shrub
#' \item GROUP Classsification of the sample into one of four vegetation types using numerical classification (Ward's agglomerative clustering applied on Euclidean distances using log transformed compositional data about vltava$spe matrix with herb + merged tree and shrub species).
#' \item TBV.NO	Turboveg number - unique identifier under which the plot is stored in Czech National Phytosociological Database (http://www.sci.muni.cz/botany/vegsci/dbase.php?lang=en)
#' \item TRANSECT	transect number
#' \item LIGHT	mean Ellenberg indicator values for light, calculated as unweighted mean from data in Vltava spe (herbs + merged trees and shrubs)
#' \item TEMP	mean Ellenberg indicator values for temperature, calculated as unweighted mean from data in Vltava spe (herbs + merged trees and shrubs)
#' \item CONT	mean Ellenberg indicator values for continentality, calculated as unweighted mean from data in Vltava spe (herbs + merged trees and shrubs)
#' \item MOIST	mean Ellenberg indicator values for moisture, calculated as unweighted mean from data in Vltava spe (herbs + merged trees and shrubs)
#' \item REACT	mean Ellenberg indicator values for soil reaction calculated as unweighted mean from data in Vltava spe (herbs + merged trees and shrubs)
#' \item NUTR	mean Ellenberg indicator values for nutrients, calculated as unweighted mean from data in Vltava spe (herbs + merged trees and shrubs)
#' }
#' @usage data (vltava)
#' @format
#'  \code{vltava} is a structured list with these items:
#'  \itemize{
#'  \item \code{spe} Compositional matrix of all species (sample x species, percentage cover scale)
#'  \item \code{ell} Species Ellenberg indicator values (species x Ellenberg values for light, temperature, continentality, moisture, reaction and nutrients, compiled from Ellenberg et al. 1991).
#'  \item \code{env} Environmental variables (see Details).
#'  \item \code{spnames} Data frame with two columns: \code{Full.species.name} - original species names, and \code{Layer} - vegetation layer, in which the species occur (1 - herb layer, 23 - shrub or/and tree layer)
#'  \item \code{herbs} list with the following items, related only to the subset of herb species:
#'    \itemize{
#'    \item \code{spe} Compositional matrix of herb species (sample x species, percentage cover scale)
#'    \item \code{ell} Species Ellenberg indicator values for herb species (species x Ellenberg values for light, temperature, continentality, moisture, reaction and nutrients)
#'    \item \code{traits} Species functional traits for plant height (compiled from Czech flora, Kubat et al. 2002), specific leaf area (SLA) and seed weight (compiled from LEDA database, Kleyer et al. 2008).
#'    \item \code{spnames} Data frame with two columns: \code{Full.species.name} - original species names, and \code{Layer} - vegetation layer, in which the species occur (1 - herb layer)
#'  }
#'  \item \code{all} list with the following items, related to matrix with all species (trees, shrubs, herbs and juveniles)
#'    \itemize{
#'    \item \code{spe} Compositional matrix of all species (sample x species, percentage cover scale)
#'    \item \code{spnames} Data frame with two columns: \code{Full.species.name} - original species names, and \code{Layer} - vegetation layer, in which the species occur (3 - tree layer, 2 - shrub layer, 1 - herb layer, J - juveniles of woody species)
#'    }
#'  }

#' 
#' @name vltava
#' @docType data
#' @author David Zeleny (\email{zeleny.david@@gmail.com})
#' @references
#' 
#' Ellenberg H., Weber H.E., Dull R., Wirth V., Werner W. & Paulissen D. 1991. Zeigerwerte von Pflanzen in Mitteleuropa. Scripta Geobotanica 18: 1-248.
#' 
#' Kleyer M., Bekker R.M., Knevel I.C., Bakker J.P., Tompson K., Sonnenshein M. et al. (2008) The LEDA Traitbase: a database of life-history traits of Northwest European flora. Journal of Ecology, 96, 1266-1274.
#' 
#' Kubat K., Hrouda L., Chrtek J. Jr., Kaplan Z., Kirschner J. & Stepanek J. (eds.) (2002) Klic ke kvetene Ceske Republiky (Key to the flora of the Czech Republic). Academia, Praha, Czech Republic.
#' 
#' McCune B. & Keon D. (2002): Equations for potential annual direct incident radiation and heat load. Journal of Vegetation Sciences, 13: 603-606.
#' 
#' Zeleny D. & Chytry M. (2007): Environmental control of vegetation pattern in deep river valleys of the Bohemian Massif. Preslia, 79: 205-222.
#' 
#' 
#' @keywords vltava.spe vltava.ell vltava.env
NULL