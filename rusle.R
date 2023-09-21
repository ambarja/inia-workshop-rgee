library(rgee)
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(cptcity)

# Iniciar sesión con Google Earth Engine
ee_Initialize(user = "barja.geografo@gmail.com",drive = TRUE)

# 1. Conjunto de datos a utilizar -----------------------------------------
vegetation <- ee$ImageCollection("MODIS/061/MOD13Q1")
chirps <- ee$ImageCollection("UCSB-CHG/CHIRPS/PENTAD")
srtm <- ee$Image("USGS/SRTMGL1_003")

sand <- ee$Image("OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02")$select("b0")$rename("sand")
silt <- ee$Image("users/aschwantes/SLTPPT_I")$divide(100)$rename("silt")
clay <- ee$Image("OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02")$select("b0")$rename("clay")
org <- ee$Image("OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02")$select("b0")$rename("org")
soil <- ee$Image(c(clay, sand, silt, org))

cuenca <- ee$FeatureCollection("projects/inia-gee-399313/assets/hidrocuencas")
# 2. Parámetros a usar ----------------------------------------------------

# Años de análisis
startDate = '2021-01-01'
endDate = '2021-12-31'

# Nombre de cuenca a analizar
nombre = 'Cuenca Mantaro'
cuenca = cuenca$filter(ee$Filter$eq('nombre',nombre))
# Map$centerObject(cuenca,zoom = 7)
Map$addLayer(cuenca)

# 3. Cálculo del factor R -------------------------------------------------

pp <- chirps$filter(ee$Filter$date(startDate, endDate))

# Lista de meses (de 1 a 12).
meses <- ee$List$sequence(1, 12)

# Función de mapeo para calcular la suma de precipitación por mes.
calcularSumaPorMes <- function(mes) {
  imagenPorMes = pp$filter(ee$Filter$calendarRange(mes, mes, 'month'))$
    reduce(ee$Reducer$sum())
  return (imagenPorMes$set('month', mes))
}

pp_por_meses <- ee$ImageCollection(meses$map(ee_utils_pyfunc(calcularSumaPorMes)))$
  reduce(ee$Reducer$mean())
sum_pp <- pp_por_meses$reduce(ee$Reducer$sum())

# Factor R
factorR <- ee$Image(10)$pow(ee$Image(1.5)$multiply(pp_por_meses$pow(2)$divide(sum_pp)$log10()$subtract(-0.08188)))$multiply(1.735)
Map$addLayer(factorR$clip(cuenca),list(min = 430.1863776322073, max = 21061.884335682855))

# 4. Cálculo del factor K -------------------------------------------------
fcsand <- soil$expression(
  '0.2 + 0.3 * exp(-0.256 * sand * (1-(silt/100)))',
  list(
    'sand' = soil$select('sand'),
    'silt' = soil$select('silt')
  )
)

fcl_si <- soil$expression(
  '(silt/(clay + silt))**0.3',
  list(
    'silt' = soil$select('silt'),
    'clay' = soil$select('clay')
  )
)

forg <- soil$expression(
  '1 - ((0.25*org)/(org + exp(3.72-2.95*org)))',
  list(
    'org' = soil$select('org')
  )
)

fhyosand <- soil$expression(
  '1-((0.7*(1 - sand/100))/((1 - sand/100) + exp(-5.51 + 22.9*(1 - sand/100))))',
  list(
    'sand' = soil$select('sand')
  )
)

factores <- ee$Image(c(fcsand, fcl_si, forg, fhyosand))$rename(c("fcsand", "fcl_si", "forg", "fhyosand"))
factorK <- factores$expression(
  "fcsand * fcl_si * forg * fhyosand",
  list(
    "fcsand" = factores$select("fcsand"),
    "fcl_si" = factores$select("fcl_si"),
    "forg" = factores$select("forg"),
    "fhyosand" = factores$select("fhyosand")
  )
)

Map$addLayer(factorK$clip(cuenca),list(min =0 ,max = 0.25))

# 5. Factor LS ------------------------------------------------------------
elevation <- srtm$select('elevation')
slope <- ee$Terrain$slope(elevation)$clip(cuenca)

# De grados a %
slope <- slope$divide(180)$multiply(pi)$tan()$multiply(100)

# compute LS factor (LS = ((Q*M/22.13).pow(0.5) * (0.065 + 0.045S + 0.0065S*S), Wischmeier and Smith (1978))
# compute LS factor (LS = ((Q*M/100).pow(2) * (0.76 + 0.53S + 0.076S*S), Wischmeier and Smith (1978))
ls1 <- sqrt(500/22.13)
ls2 <- slope$multiply(0.53)
ls3 <- slope$pow(2)$multiply(0.076)
ls4 <- ls2$add(ls3)$add(0.76)
factorLS <- ls4$multiply(ls1)$rename("LS")

Map$addLayer(factorLS$clip(cuenca),list(min = 3, max= 106, palette = c('a52508','ff3818','fbff18','25cdff','2f35ff','0b2dab')))

# 6. Factor C -------------------------------------------------------------
# Factor C = 0.1x(-NDVI + 1)/2 Almagro et al (2019).
ndvi_median <- vegetation$filterDate(startDate,endDate)$mean()$multiply(0.0001)$select('NDVI')
factorC <- ee$Image(0.1)$multiply(ndvi_median)$multiply(-1)$add(1)$divide(2)       
Map$addLayer(factorC$clip(cuenca),list(min = 0, max = 0.63))

# 8. Erosión del suelo ----------------------------------------------------

erosion <- factorR$multiply(factorK)$multiply(factorLS)$multiply(factorC)$rename("erosion")
Map$addLayer(erosion)

stats <- erosion$reduceRegion(
  reducer= ee$Reducer$max(),
  geometry= cuenca$geometry(),
  scale= 5000
  )$getInfo()

viz <- c('#490eff','#12f4ff','#12ff50','#e5ff12','#ff4812')

Map$addLayer(erosion,list(min = 0.001, max= 250, palette = viz))
















