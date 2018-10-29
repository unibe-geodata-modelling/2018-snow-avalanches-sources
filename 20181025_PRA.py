# Script for the automated mapping of potential snow avalanche release areas
# Thalia Bertschinger
# Seminar Geodata Analysis and Modelling, Spring Semester 2018
# The model creates polygons of potential snow avalanche release areas

# imports
import arcpy

if arcpy.CheckExtension("Spatial") == "Available":
    arcpy.CheckOutExtension('Spatial')
else:
    print "no spatial analyst license available"

#**************************************************************************
# ENVIRONMENT variables - set workspace and names of input files
#**************************************************************************

# set environment and workspace and create gdb
preworkspace = "U:/Seminar_Modellieren/20181018_Test_Model"
tempdir = "C:/temp"
arcpy.env.overwriteOutput = True
# Create File GDB
gdb = "20181029_Model_PRA.gdb"
arcpy.CreateFileGDB_management(preworkspace, gdb, "CURRENT")
myworkspace = preworkspace+"/"+gdb
print "Workspace: " + myworkspace
arcpy.env.workspace = myworkspace

# input data
pre_dem = arcpy.Raster(preworkspace+"/"+"dem.tif")

# Resample dem to 5 x 5 cell size
dem_5m = myworkspace+"/"+"dem_5m"
cell_size = 5
resample_method = "BILINEAR"
arcpy.Resample_management(pre_dem, dem_5m, cell_size, resample_method)

dem = arcpy.Raster(dem_5m)
arcpy.env.cellSize = dem
arcpy.env.extent = dem
arcpy.env.snapRaster = dem

#**************************************************************************
# create temporary files if wished
#**************************************************************************

# slope = tempdir+"/"+"slope.tif" # temporary file for the slope analysis
# aspect = tempdir+"/"+"aspect.tif" # temporary file for the aspect analysis
# aspect_classes = tempdir+"/"+"aspect_classes.tif" # temporary file for the aspect classes (reclassified aspect)
# curvature = tempdir+"/"+"curvature.tif" # temporary file for the curvature analysis
plan_curvature = tempdir+"/"+"plan_curvature.tif"  # temporary file for the plan curvature
profile_curvature = tempdir+"/"+"profile_curvature.tif"  # temporary file for the plan curvature

#**************************************************************************
# start analyses
#**************************************************************************

# create slope
out_slope = "DEGREE"
z_slope = 1
outSlope = arcpy.sa.Slope(dem, out_slope, z_slope)
# outSlope.save(slope)

# create aspect
outAspect = arcpy.sa.Aspect(dem)
# outAspect.save(aspect)
# reclassify aspect to get 9 aspect classes
reclass_field = "VALUE"
range_aspect_classes = arcpy.sa.RemapRange([[-1, 0, 1],
                                            [0, 22.500000, 2],
                                            [22.500000, 67.500000, 3],
                                            [67.500000, 112.500000, 4],
                                            [112.500000, 157.500000, 5],
                                            [157.500000, 202.500000, 6],
                                            [202.500000, 247.500000, 7],
                                            [247.500000, 292.500000, 8],
                                            [292.500000, 337.500000, 9],
                                            [337.500000, 360, 2]])
outAspect_classes = arcpy.sa.Reclassify(outAspect, reclass_field, range_aspect_classes)
# outAspect_classes.save(aspect_classes)

# create curvature
z_curvature = 1
outCurvature = arcpy.sa.Curvature(dem, z_curvature, profile_curvature, plan_curvature)
# outCurvature.save(curvature)
# plan curvature is used later
planCurvature = arcpy.Raster(plan_curvature)

# identify slopes in a certain range and assign them the value 1, everything else is set to NoData
min_slope = 30
max_slope = 60
conSlope = arcpy.sa.Con((outSlope >= min_slope) & (outSlope <= max_slope), 1)
# conSlope.save(tempdir+"/"+"conslope.tif")

# set cells with a plan curvature bigger than a certain threshold to NoData (first assign them the value 2)
max_curv = 6
conCurv = arcpy.sa.Con(planCurvature > max_curv, 2, conSlope)
PRA1 = arcpy.sa.Con(conCurv == 1, 1)

# assign aspect values to the cells in PRA1
PRA1_aspect = arcpy.sa.Con(PRA1 == 1, outAspect_classes)
# PRA1_aspect.save(tempdir+"/"+"PRA1_apect.tif")

#**************************************************************************
# convert raster to polygon and simplify them
#**************************************************************************

# convert raster to polygon
PRA1_poly = myworkspace+"/"+"PRA1"  # is saved in gdb to have the area automatically added in the attribute table
arcpy.RasterToPolygon_conversion(PRA1_aspect, PRA1_poly, "SIMPLIFY", "Value")

# convert multipart to singlepart features
PRA1_single = myworkspace+"/"+"PRA1_single"
arcpy.MultipartToSinglepart_management(PRA1_poly, PRA1_single)

# eliminate polygon parts with an area smaller than a certain threshold
PRA2_1 = myworkspace+"/"+"PRA2_1"
PRA2_con = "AREA"
area_tresh = "5000 SquareMeters"
perc = 0  # Eliminate parts smaller than this percentage of a feature's total outer area.
option = "ANY"
arcpy.EliminatePolygonPart_management(PRA1_single, PRA2_1, PRA2_con, area_tresh, perc, option)

# make a feature layer and select polygons with an area smaller than 5000 m^2 and merge them with neighbouring polygons
# create a loop that is executed as long as polygons smaller than 5000 m^2 are present and can be merged
# create a list (small_polygons) where the number of polygons smaller than 5000 m2 will be added
# create a list (path_PRA2_elim) where the path of the created files will be added
small_polygons = []
path_PRA2_elim = []
count = 1
while count >= 1:
    feature_layer = "PRA_layer"+str(count)
    in_feature = myworkspace+"/"+"PRA2_"+str(count)
    arcpy.MakeFeatureLayer_management(in_feature, feature_layer)  # make a feature layer
    selection = "NEW_SELECTION"
    merge_tresh = "Shape_Area < 5000"
    # select polygons smaller than 5000 m^2
    arcpy.SelectLayerByAttribute_management(feature_layer, selection, merge_tresh)
    sel_polygons = int(arcpy.GetCount_management(feature_layer).getOutput(0))
    small_polygons.append(sel_polygons)
    print('{} polygons are smaller than 5000 m^2'.format(sel_polygons))
    while sel_polygons > 0:
        if len(small_polygons) == 1:
            count += 1
            PRA2_elim = myworkspace + "/" + "PRA2_" + str(count)
            merge_how = "LENGTH"  # The neighboring polygon is the one with the longest shared border.
            # Merges a selected polygon with a neighbouring unselected polygon by dropping the shared border.
            # merge the small polygons with neighbouring big ones
            arcpy.Eliminate_management(feature_layer, PRA2_elim, merge_how)
            path_PRA2_elim.append(PRA2_elim)
            break
        elif len(small_polygons) > 1:
            if small_polygons[- 1] != small_polygons[- 2]:
                count += 1
                PRA2_elim = myworkspace + "/" + "PRA2_"+str(count)
                merge_how = "LENGTH"  # The neighboring polygon is the one with the longest shared border.
                # Merges a selected polygon with a neighbouring unselected polygon by dropping the shared border.
                # merge the small polygons with neighbouring big ones
                arcpy.Eliminate_management(feature_layer, PRA2_elim, merge_how)
                path_PRA2_elim.append(PRA2_elim)
                break
            else:
                count = 0
                break
        else:
            count = 0
            break
    else:
        count = 0
        break

print "Loop to merge small polygons has finished."


# Create a feature class containing polygons which represent a specified minimum bounding geometry enclosing
# all polygons of the potential release areas created in the previous step
extent_PRA = myworkspace + "/" + "extent_PRA"
# take the last output of the previous step as input file (index -1 takes the last element of the list)
in_PRA2_elim = path_PRA2_elim[-1]

arcpy.MinimumBoundingGeometry_management(in_PRA2_elim, extent_PRA)

# Process: Minimum Bounding Geometry
arcpy.MinimumBoundingGeometry_management(PRA5_ohne_NoPRA__4_, PRA5_ohne_NoPRA_extent, "RECTANGLE_BY_AREA", "ALL", "", "NO_MBG_FIELDS")


#**************************************************************************
# delete all files in the temp folder
#**************************************************************************

# clear all variables except tempdir and myworkspace
for name in dir():
    if name != 'tempdir' and name != 'myworkspace':
        del globals()[name]

import os
for the_file in os.listdir(tempdir):
    file_path = os.path.join(tempdir, the_file)
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
    except Exception as e:
        print(e)

#**************************************************************************
# delete all files in the gdb except the real PRA file
#**************************************************************************

# fuktioniert noch nicht!!!! Bedingung einfügen, um ein file zu behalten!!
for the_file in os.listdir(myworkspace):
    file_path = os.path.join(myworkspace, the_file)
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
    except Exception as e:
        print(e)
