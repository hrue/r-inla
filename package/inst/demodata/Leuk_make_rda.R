
## data
fl.dat <- system.file("demodata/Leuk.dat", package = "INLA")
if(fl.dat=="")
    fl.dat <- paste0("https://raw.githubusercontent.com/hrue/r-inla/",
                     "refs/heads/devel/r-inla.org/examples/leukemia/leuk.dat")
Leuk <- read.table(fl.dat, header = TRUE)
Leuk[1:2, ]

### map coordinates and id files
fl.map <- system.file("demodata/Leuk.map", package="INLA")
fl.id  <- system.file("demodata/Leuk.names.regions", package="INLA")

## read all coordinates
allcoo <- read.table(fl.map)

## find the split indexes
length(ii0 <- which(is.na(allcoo[,1])))

### some regions contain multiple polygons
ii1 <- split(ii0, as.integer(readLines(fl.id)))
length(unique(ii1))

library(sp)

## SpatialPolygons 
nwEngland <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(
        lapply(1:length(ii1), function(i) {
            ii <- ii1[[i]]
            Polygons(
                lapply(ii, function(i) {
                    ii <- c(i + 1:allcoo[i, 2], i+1)
                    Polygon(cbind(allcoo[ii, 1], allcoo[ii, 2]), FALSE)
                }),
                ID = names(ii1)[i]
            )
        }
        )
    ),
    data = data.frame(ID = names(ii1))
)

## sf 
library(sf)

NEmap <- st_sf(
    ID = names(ii1), 
    geom = st_sfc(
        lapply(
            ii1, function(ii) 
                st_make_valid(                
                    st_polygon(
                        lapply(ii,
                               function(i) {
                                   ii <- c(i + 1:allcoo[i, 2], i+1)
                                   cbind(allcoo[ii, 1], allcoo[ii, 2])
                               }
                               )
                    )
                )
        )
    )
)

ddir <- here::here("data/")

if(FALSE)
    save(
        'Leuk', 'NEmap', 'nwEngland',
        file = paste0(ddir, "Leuk.rda")
    )
