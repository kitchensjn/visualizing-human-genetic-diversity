library(gsheet)
library(tidyverse)
library(sf)


predictions <- gsheet2tbl("docs.google.com/spreadsheets/d/1JnoYyWDmQ35xHzaHd8cVtQI-RFMrb5vOYRC1HimHmw8")


ggplot(data=predictions) +
  geom_point(aes(ACB.h, ACB.k)) +
  theme_minimal()


make_ellipse <- function(x0, y0, a, b, angle) {
  y0 <- -y0
  angle <- -angle#/360 * 2 * pi
  theta <- c(seq(0, 2 * pi, length.out = 36000), 0)
  
  crds <- cbind(a * cos(theta) * cos(angle) - b * sin(theta) * sin(angle) + x0,
                a * cos(theta) * sin(angle) + b * sin(theta) * cos(angle) + y0)
  return(sf::st_polygon(list(crds)))
}

PEL <- make_ellipse(245.9637,354.8639,1139.4186,1374.1214,1.0056)
MXL <- make_ellipse(188.7347,343.7469,1298.0515,1327.9438,-1.3532)
CEU <- make_ellipse(83.4634,437.6332,1518.3666,1150.6006,1.187)
CLM <- make_ellipse(163.3362,321.562,1247.0832,1433.2207,-0.2126)
PUR <- make_ellipse(48.8014,359.7592,1435.6344,1296.4881,-1.5072)
ASW <- make_ellipse(-186.6873,-117.9368,1753.6272,1346.1266,-1.7434)
ACB <- make_ellipse(-351.0253,-117.9368,1661.9991,1482.8142,-2.3891)

overlap <- ggplot() +
  geom_sf(data=PEL, fill=NA) +
  geom_sf(data=MXL, fill=NA) +
  geom_sf(data=CEU, fill=NA) +
  geom_sf(data=CLM, fill=NA) +
  geom_sf(data=PUR, fill=NA) +
  geom_sf(data=ASW, fill=NA) +
  geom_sf(data=ACB, fill=NA) +
  theme_void()

intersection <- st_intersection(st_intersection(st_intersection(st_intersection(st_intersection(st_intersection(PEL, MXL), CEU), CLM), PUR), ASW), ACB)

overlap +
  geom_sf(data=intersection, fill="red")

st_area(intersection)
st_area(PEL)
st_area(MXL)
st_area(CEU)
st_area(CLM)
st_area(PUR)
st_area(ASW)
st_area(ACB)


3877369


ggplot() +
  geom_sf(data=PEL, fill=NA) +
  geom_sf(data=MXL, fill=NA) +
  theme_void()

st_area(st_intersection(PEL, MXL))

