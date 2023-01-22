# LST_Landsat8_L1

# Landsat 8 Level 1 
## Land Surface Temperature

## Load library
```
library(raster)
library(sf)
```

## Set Work Directory
```
setwd(".../2022_07_19_Landsat8_L1")
```

## Load data
```
# MTL
MTL = read.delim("LC09_L1TP_192024_20220719_20220719_02_T1_MTL.txt", sep = "=")

# Band 4 (RED)
RED = raster("LC09_L1TP_192024_20220719_20220719_02_T1_B4.TIF")

# Band 5 (NIR)
NIR = raster("LC09_L1TP_192024_20220719_20220719_02_T1_B5.TIF")

# Band 10 (TIR)
TIR = raster("LC09_L1TP_192024_20220719_20220719_02_T1_B10.TIF")

# Load Cottbus Boundary
Cottbus = sf::read_sf("Cottbus_Boundary.shp")
```

## Mask and Crop to Cottbus region
```
RED = raster::mask(RED, Cottbus)
RED = raster::crop(RED, Cottbus)

NIR = raster::mask(NIR, Cottbus)
NIR = raster::crop(NIR, Cottbus)

TIR = raster::mask(TIR, Cottbus)
TIR = raster::crop(TIR, Cottbus)
```

###############################################################################################
# Start of the Calculation
###############################################################################################
## 1. Conversion to TOA Radiance

```
M = as.numeric(MTL$LANDSAT_METADATA_FILE[which(MTL$GROUP == "    RADIANCE_MULT_BAND_10 ")])
A = as.numeric(MTL$LANDSAT_METADATA_FILE[which(MTL$GROUP == "    RADIANCE_ADD_BAND_10 ")])

L = (M * TIR) + A
plot(L)
```
#####################################################
## 2. Conversion to Top of Atmosphere Brightness Temperature
```
K1 = as.numeric(MTL$LANDSAT_METADATA_FILE[which(MTL$GROUP == "    K1_CONSTANT_BAND_10 ")])
K2 = as.numeric(MTL$LANDSAT_METADATA_FILE[which(MTL$GROUP == "    K2_CONSTANT_BAND_10 ")])

T = (K2 / log((K1/L) + 1)) - 273.15
```
#####################################################
## 3. NDVI
```
NDVI = (NIR - RED) / (NIR + RED)
#plot(NDVI)
```
#####################################################
## 4. proportion of vegetation Pv
```
NDVI_min = raster::minValue(NDVI)
NDVI_max = raster::maxValue(NDVI)

Pv = ((NDVI - NDVI_min) / (NDVI_max - NDVI_min))**2
```
#####################################################
## 5. Calculate Emissivity ε
```
E = (0.004 * Pv) +0.986
```
#####################################################
## 6. Calculate Land Surface Temperature 
```
BT = T
w = 10.8 # μm

# h = Planck’s constant (6.626 * 10^-34 Js)
# s = Boltzmann constant (1.38 * 10^-23 J/K)
# c = velocity of light (2.998 * 10^8 m/s)

h = 6.626 * 10^-34
s = 1.38 * 10^-23
c = 2.998 * 10^8
p = h * (c / s)

p = 14380

E = E

LST = BT / (1 + ((w * (BT / p)) * log(E)))
LST

plot(LST)
#writeRaster(LST, ".../LST_2022_07_19_Landsat8L1.tif")
```
#####################################################

