# Script for processing and visualizing Gravel Bed Rivers data in the Yellowstone National Park, WY
#
# Author: Saad Tarik
# Created: Feb 23, 2019

library(lidR)
fname <- "points.las"

lidar_points <- readLAS(fname)
plot(lidar_points)

# Compute a canopy height model
# Khosravipour et al. pitfree algorithm
thr <- c(0,2,5,10,15)
edg <- c(0, 1.5)
chm <- grid_canopy(lidar_points, 1, pitfree(thr, edg))
plot(chm)

# Individual tree segmentation
lasTreeSeg <- lastrees(lidar_points, li2012())
col <- random.colors(200)
plot(lasTreeSeg, color = "treeID", colorPalette = col)

# Catalog
ctg <- catalog(".")
plot(ctg, map = TRUE)

# Digital Terrain Model
terrain <- grid_terrain(lidar_points, res = 1, algorithm = kriging(k = 10L)) # choose between tin(), kriging(k = 10L), and knnidw(k = 6L, p = 2)
plot(terrain)
plot_dtm3d(terrain)

# Run snag classification and assign classes to each point
# For the Wing2015 method, supply a matrix of snag BranchBolePtRatio conditional
# assessment thresholds (see Wing et al. 2015, Table 2, pg. 172)
bbpr_thresholds <- matrix(c(0.80, 0.80, 0.70,
                            0.85, 0.85, 0.60,
                            0.80, 0.80, 0.60,
                            0.90, 0.90, 0.55),
                          nrow =3, ncol = 4)
las_classes <- lassnags(lidar_points, wing2015(neigh_radii = c(1.5, 1, 2), BBPRthrsh_mat = bbpr_thresholds))
plot(las_classes, color="snagCls", colorPalette = rainbow(5))

# Tree-top detection
tree_tops <- tree_detection(chm, lmf(ws = 5))
plot(tree_tops)
