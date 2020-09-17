## Load Libraries-----------------------------------------------------------------------------------------------------------------------------------------
spatial_libs <- c("rgdal","rayshader", "FedData","raster","sp", "automap","geoviz","sf", "ggspatial")

# general
general_libs <- c("data.table","dplyr","tidyr","ggplot2","viridisLite","data.table","grid","png","RColorBrewer", "ggnewscale")

# load libraries invisibly
invisible( lapply( c(spatial_libs,general_libs),
                   library, character.only = T) )
## Set Theme -------------------------------------------------------------------------------------------------------------------------------------------------------
theme_set(theme_bw())

# Load Data -------------------------------------------------------------------------------------------------------------------------------------------------------
# read phyloseq object
full_dat <- readRDS("data/processed/final_marine_phy.rds")

# pull out metadata
meta <-as(full_dat@sam_data, "data.frame")
meta <- as(meta, "data.table")

meta <- meta[ , c("lat","long"):= list(as.numeric(lat), as.numeric(long))]


## Set Boundaries -------------------------------------------------------------------------------------------------------------------------------------------------------

# Check geographic range of sampling points
limits <- c(
  min(meta$long),
  min(meta$lat), 
  max(meta$long),
  max(meta$lat) 
)

# Define bounding box with a cushion around sampling points

bbox <- list(
  p1 = list(long = limits[1] - 0.03, lat = limits[2] - 0.04) ,
  p2 = list(long = limits[3] + 0.03, lat = limits[4] + 0.04)
)

# add a reformatted version of this box for use with ggmap
bbox$ggmap <- unlist(bbox[1:2])
names(bbox$ggmap) <- c("left","bottom","right","top")

# create an emtpy spatial polygon template of the bounding box for pulling elevation data
bbox_extent <- polygon_from_extent(raster::extent(bbox$p1$long, bbox$p2$long, bbox$p1$lat, bbox$p2$lat),
                                   proj4string = "+proj=longlat +ellps=GRS80 +datum=WGS84 +no_defs")  



std_proj <- crs(bbox_extent)

# set up plot limits
plot_lims <- list(x=c(bbox$p1$long, bbox$p2$long), y = c(bbox$p1$lat, bbox$p2$lat))

##  Base Map -------------------------------------------------------------------------------------------------------------------------------------------------------

# hawaii coastline (state GIS, projection utm-4, nad83)
hi_coast         <- st_read("data/raw/hawaii_coastline/coast_n83.shp", crs = "epsg:26904")

# project to std proj
hi_coast <-st_transform(hi_coast,std_proj)


## Elevation Raster  -------------------------------------------------------------------------------------------------------------------------------------------------------

# get extent (federal 1 arc second elevation, lat long GRS80)
elev_extent <- spTransform(bbox_extent, CRSobj = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs ")

# pull elevation data in raster format
elev_ras_full <-
  get_ned(
    template = bbox_extent,
    label = "elev",
    res = "1",
    force.redo = F,
    raw.dir = "../data/raw/spatial",
    extraction.dir = "../data/processed/interim/"
  )

crs(elev_ras_full)


# project into standard format
elev_ras_full <- projectRaster(elev_ras_full, crs=std_proj)

crs(elev_ras_full)


## Bathymetry Rasters -------------------------------------------------------------------------------------------------------------------------------------------------------
## Bathymetry
bath_ras_50 <- raster("data/raw/mhi_mbsyn_bathyto")

bath_ras_50 <- crop(bath_ras_50, bbox_extent)

# set projection and resolution to match elevation
bath_ras_50 <- projectRaster(bath_ras_50, elev_ras)

bath_dt_50 <- as(raster::as.data.frame(bath_ras_50, xy=T), "data.table")

## Mask Bathymetry and Elevation ----------------------------------------------------------------------------------------------------------------------------

# Mask Bathymetry
bath_ras <- mask(bath_ras_50, elev_ras_full)

# convert to data table and matrix

bath_dt <- as(raster::as.data.frame(bath_ras, xy=T), "data.table")
setnames(bath_dt, 3, "depth")
bath_dt[ , alpha := ifelse(is.na(depth), 0, 1)]

# as matrix
bath_mat <- as.matrix(bath_ras) %>% t()

# test map 
bath_map <- ggplot() +
  geom_raster(data= bath_dt, mapping =aes(x = x, y = y, fill = depth))
bath_map

# Mask Elevation
elev_ras <- mask(elev_ras_full, bath_ras, inverse = T)

# and convert raster to data table and matrix 
elev_dt  <- raster::as.data.frame(elev_ras, xy=T) %>% as.data.table()
setnames(elev_dt, 3, "elev")
elev_dt[ , alpha := ifelse(is.na(elev), 0, 1)]
elev_mat <- as.matrix(elev_ras) %>% t()

# test map
elev_map <- ggplot() +
  geom_raster(data= elev_dt, mapping =aes(x = x, y = y, fill = elev))
elev_map



## -------------------------------------------------------------------------------------------------------------------------------------------------------

# function crop_sat will read in RGB landsat raster, crop it to map extent, write out as png

# use this fuction to define image size

define_image_size <- function(bbox, major_dim = 400) {
  
  # calculate aspect ration (width/height) from lat/long bounding box
  aspect_ratio <- abs((bbox$p1$long - bbox$p2$long) / (bbox$p1$lat - bbox$p2$lat))
  
  # define dimensions
  img_width <- ifelse(aspect_ratio > 1, major_dim, major_dim*aspect_ratio) %>% round()
  img_height <- ifelse(aspect_ratio < 1, major_dim, major_dim/aspect_ratio) %>% round()
  
  size_str <- paste(img_width, img_height, sep = ",")
  
  list(height = img_height, width = img_width, size = size_str)}


# calculate image dimensions based on size of elevation raster
img_size <- define_image_size(bbox = bbox, major_dim = max(dim(elev_ras)))

crop_sat <-function(rast_file = "data/raw/oahu_landsat/Oahu_Landsat_15m.jp2", extent = bbox_extent, new_proj = std_proj){
  rast <- stack(rast_file)
  
  rast_proj   <- crs(rast)  
  rast_extent <- spTransform(extent, CRSobj = rast_proj)
  
  small_rast <- crop(rast, rast_extent)
  
  message("writing out cropped png")
  
  png("data/processed/sat_basemap.png", width = img_size$width, height = img_size$height)
  par(mar = c(0,0,0,0))
  plotRGB(small_rast)
  dev.off()
  
  return(small_rast)
}

sat_ras <- crop_sat()

plotRGB(sat_ras)


sat_grob 

## Streams -----------------------------------------------------------------------------------------------------------------------------------
# check stream metadata to identify stream name
streams      <- st_read('data/raw/darstreams.kml')
waimea_river <- streams[streams$STREAM_NAM == "Waimea R",]

test <- waimea_river[10,]

## Plot 2D Map -------------------------------------------------------------------------------------------------------------------------------

elev_colors <-colorRampPalette(colors = c(rep("darkblue", 3), "lightblue",rep("darkgreen",3)))

# Plot on raster background
p1<- ggplot() + 
  # raster
  geom_raster(data= elev_dt, mapping =aes(x = x, y = y, fill = elev), alpha = elev_dt$alpha)+
  scale_fill_gradient(low = "darkgreen", high = "tan") +
  new_scale_fill()+
  geom_raster(data = bath_dt,  mapping =aes(x = x, y = y, fill = depth), alpha = bath_dt$alpha)+
  scale_fill_gradient(low = "darkblue", high = "lightblue")+

  # river
  geom_sf(data = waimea_river, color = alpha("blue",0.8)) +
  
  # points
  geom_point(aes(x = long, y = lat, color = site_name), data = point_data, shape = 16) +
  scale_color_manual(values = viridis(n=5)) +
  
  # annotations
  annotation_scale(location = "bl", width_hint = 0.4)+
  annotation_north_arrow(height = unit(0.5,"cm"),
                         width = unit(0.5,"cm"),
                         pad_y = unit(0.75,"cm"),
                         style = north_arrow_minimal("text_size = 8"))+
  
  # map settings
  coord_sf(xlim = plot_lims$x, ylim = plot_lims$y, clip = "on") +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white")) +
  labs(fill = "Depth (m)", color = "Site") 

p1

# With Oahu inset
oahu_limits <- c(-158.305366, 21.212964, -157.624682, 21.766562)

p2 <-ggplot()+
  # coastline
  geom_sf(data = hi_coast) +
  # points
  geom_point(aes(x = long, y = lat, color = site_name), data = point_data, size = 0.5, shape = 17) +
  scale_color_manual(values = c("orange","red","blue","yellow","black")) +
  # rectangle
  geom_rect(aes(xmin = plot_lims$x[1]- 0.003,
                xmax = plot_lims$x[2] + 0.003,
                ymin = plot_lims$y[1] - 0.003,
                ymax = plot_lims$y[2] + 0.003), color = "gray10", fill = NA)+
  # map settings
  theme_void()+
  theme(legend.position = "none")+
  coord_sf(xlim = c(oahu_limits[c(1,3)]), ylim= c(oahu_limits[c(2,4)])) 


p2

p1 + annotation_custom(
  grob = ggplotGrob(p2),
  ymin = plot_lims$y[1],
  ymax = plot_lims$y[1] +.02,
  xmin = plot_lims$x[2] - 0.02,
  xmax = plot_lims$x[2]
)

ggsave("output/map/ST_marine_map.png")


p3 <-  ggplot() + 
  # raster
  geom_raster(data= elev_dt, mapping =aes(x = x, y = y, fill = elev), alpha = elev_dt$alpha)+



## -------------------------------------------------------------------------------------------------------------------------------------------------------
# calculate image dimensions based on size of elevation raster
img_size <- define_image_size(bbox = bbox, major_dim = max(dim(elev_bath_ras)))

# now export the raster as a png with the same dimensions as the raster

png("test1.png", width = img_size$width, height = img_size$height)
par(mar = c(0,0,0,0))
raster::image(elev_bath_ras, axes = F, col = viridis(1000))
dev.off()



## -------------------------------------------------------------------------------------------------------------------------------------------------------
# read in custom overlay raster
elev_overlay <-  png::readPNG("test1.png")

# calculate rayshader layers
ambmat <- ambient_shade(elev_bath_mat, zscale = 30)
raymat <- ray_shade(elev_bath_mat, zscale = 30, lambert = F)
watermap <- detect_water(elev_bath_mat)

# define zscale
zscale <- 10



## -------------------------------------------------------------------------------------------------------------------------------------------------------
# plot 2D
elev_bath_mat%>%
  sphere_shade(texture = "imhof4") %>%
  add_water(watermap, color = "imhof4") %>%
  add_shadow(raymat, max_darken = 0.5) %>%
  add_shadow(ambmat, max_darken = 0.5) %>%
  add_overlay(elev_overlay, alphalayer = 0.5) %>%
  plot_map()


## -------------------------------------------------------------------------------------------------------------------------------------------------------
# clear plot
rgl::clear3d()

# plot that 3d map!
elev_bath_mat %>%
  sphere_shade(texture = "bw") %>%
  add_water(watermap, color = "imhof1") %>%
  add_shadow(raymat, max_darken = 0.5) %>%
  add_shadow(ambmat, 0) %>%
  plot_3d(
    elev_mat,
    zscale = zscale,
    windowsize = c(1200, 1000),
    water = T,
    wateralpha = 0.5,
    waterlinealpha = 0.5,
    waterdepth = 0.001,
    watercolor = "lightblue",
    theta = -120,
    phi = 50,
    zoom = 0.65,
    fov = 0
  )


Sys.sleep(0.2)
render_snapshot(filename = "output/map/test_render.png")



## -------------------------------------------------------------------------------------------------------------------------------------------------------
# pull out sample points and remove NA values
samps <- point_data
samps <- na.omit(samps)


# Define function to add points to rayshader plot
# this function came from library geoviz
# "fixed" add_gps by cutting out all the line related items, and having it plot points exclusively
# also just made it so params get written out and I can rgl:: spheres3d separately

custom_add_points <- function (raster_input, lat, long, alt, zscale, 
          colour = "red", alpha = 0.8, 
          raise_agl = 0, point_size = 20, rad = 2) {
  coords <- latlong_to_rayshader_coords(raster_input, lat, 
                                        long)
  distances_x <- coords$x
  distances_y <- coords$y

    sp_gps <- sp::SpatialPoints(cbind(long, lat), proj4string = sp::CRS("+init=epsg:4326"))
    sp_gps <- sp::spTransform(sp_gps, sp::CRS(as.character(raster::crs(raster_input))))
    gps_ground_line <- raster::extract(raster_input, sp_gps)
  

    track_altitude <- gps_ground_line + raise_agl
    
    return(list("x" = distances_x, "y" = track_altitude/zscale, "z" = -distances_y, "col" = colour, "alpha" = alpha, "size" = point_size, "radius" = rad))

}


sites_spheres <- custom_add_points( elev_ras,
                                        samps$lat,
                                        samps$long,
                                        raise_agl = 10,
                                        zscale = zscale,
                                        colour = viridis(nrow(samps)),
                                        point_size = 10,
                                        rad = 2)



## -------------------------------------------------------------------------------------------------------------------------------------------------------
## still version for checking image quality

# clear plot
rgl::clear3d()

# generate 3D map with overlay
elev_bath_mat %>%
  sphere_shade(texture = "imhof1") %>%
  add_water(watermap, color = "desert") %>%
  add_shadow(raymat, max_darken = 0.5) %>%
  add_shadow(ambmat, 0) %>%
  rayshader::plot_3d(
    elev_mat,
    zscale = zscale,
    windowsize = c(1200, 1000),
    water = T,
    wateralpha = 0.2,
    waterdepth = 0.1,
    theta = 237,
    phi = 25,
    zoom = 0.55,
    fov = 0,
    shadow = T
    )

# add points
rgl::points3d(sites_spheres$x, sites_spheres$y, sites_spheres$z, 
                   color = sites_spheres$col, alpha = sites_spheres$alpha, size = sites_spheres$size)
render_snapshot()

