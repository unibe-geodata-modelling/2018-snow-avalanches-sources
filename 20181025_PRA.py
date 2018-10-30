# Script for the automated mapping of potential snow avalanche release areas
# Thalia Bertschinger
# Seminar Geodata Analysis and Modelling, Spring Semester 2018
# The model creates polygons of potential snow avalanche release areas

# imports
import arcpy
import os

if arcpy.CheckExtension("Spatial") == "Available":
    arcpy.CheckOutExtension('Spatial')
else:
    print "no spatial analyst license available"

# **************************************************************************
# ENVIRONMENT variables - set workspace and names of input files
# **************************************************************************

# set environment and workspace and create gdb
preworkspace = "U:/Seminar_Modellieren/20181018_Test_Model"
tempdir = "C:/temp"
arcpy.env.overwriteOutput = True
# Create File GDB
gdb = "20181030_Model_PRA.gdb"
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

# **************************************************************************
# create temporary files if wished
# **************************************************************************

# slope = tempdir+"/"+"slope.tif" # temporary file for the slope analysis
# aspect = tempdir+"/"+"aspect.tif" # temporary file for the aspect analysis
# aspect_classes = tempdir+"/"+"aspect_classes.tif" # temporary file for the aspect classes (reclassified aspect)
# curvature = tempdir+"/"+"curvature.tif" # temporary file for the curvature analysis
plan_curvature = tempdir+"/"+"plan_curvature.tif"  # temporary file for the plan curvature
profile_curvature = tempdir+"/"+"profile_curvature.tif"  # temporary file for the plan curvature

# **************************************************************************
# start analyses
# **************************************************************************

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

# **************************************************************************
# convert raster to polygon and simplify them
# **************************************************************************

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
geom_type = "RECTANGLE_BY_AREA"
group_extent = "ALL"
arcpy.MinimumBoundingGeometry_management(in_PRA2_elim, extent_PRA, geom_type, group_extent)

# unite the extent with the PRAs
PRA3_union = myworkspace + "/" + "PRA3_union"
arcpy.Union_analysis([in_PRA2_elim, extent_PRA], PRA3_union)

# convert multipart to singlepart features
PRA4_1 = myworkspace+"/"+"PRA4_1"
arcpy.MultipartToSinglepart_management(PRA3_union, PRA4_1)

# do the same thing as before to eliminate small polygons but with the area of no PRAs around
# make a feature layer and select polygons with an area smaller than 5000 m^2 and merge them with neighbouring polygons
# create a loop that is executed as long as polygons smaller than 5000 m^2 are present and can be merged
# create a list (small_polygons) where the number of polygons smaller than 5000 m2 will be added
# create a list (path_PRA2_elim) where the path of the created files will be added
small_polygons2 = []
path_PRA4_elim = []
count = 1
while count >= 1:
    feature_layer2 = "PRA_layer2"+str(count)
    in_feature2 = myworkspace+"/"+"PRA4_"+str(count)
    arcpy.MakeFeatureLayer_management(in_feature2, feature_layer2)  # make a feature layer
    selection = "NEW_SELECTION"
    merge_tresh = "Shape_Area < 5000"
    # select polygons smaller than 5000 m^2
    arcpy.SelectLayerByAttribute_management(feature_layer2, selection, merge_tresh)
    sel_polygons2 = int(arcpy.GetCount_management(feature_layer2).getOutput(0))
    small_polygons2.append(sel_polygons2)
    print('{} polygons are smaller than 5000 m^2'.format(sel_polygons2))
    while sel_polygons2 > 0:
        if len(small_polygons2) == 1:
            count += 1
            PRA4_elim = myworkspace + "/" + "PRA4_" + str(count)
            merge_how = "LENGTH"  # The neighboring polygon is the one with the longest shared border.
            # Merges a selected polygon with a neighbouring unselected polygon by dropping the shared border.
            # merge the small polygons with neighbouring big ones
            arcpy.Eliminate_management(feature_layer2, PRA4_elim, merge_how)
            path_PRA4_elim.append(PRA4_elim)
            break
        elif len(small_polygons2) > 1:
            if small_polygons2[- 1] != small_polygons2[- 2]:
                count += 1
                PRA4_elim = myworkspace + "/" + "PRA4_"+str(count)
                merge_how = "LENGTH"  # The neighboring polygon is the one with the longest shared border.
                # Merges a selected polygon with a neighbouring unselected polygon by dropping the shared border.
                # merge the small polygons with neighbouring big ones
                arcpy.Eliminate_management(feature_layer2, PRA4_elim, merge_how)
                path_PRA4_elim.append(PRA4_elim)
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

print "Loop 2 (with no PRAs) to merge small polygons has finished."

# rename the last output of the loop to get the final PRA file
rename_PRA4 = path_PRA4_elim[-1]
# PRA_final = myworkspace + "/" + "PRA_final"
PRA_final = os.path.join(myworkspace, "PRA_final")
arcpy.Rename_management(rename_PRA4, PRA_final)

# **************************************************************************
# clear all variables except tempdir and myworkspace and PRA_final
# **************************************************************************

for name in dir():
    if name != 'tempdir' and name != 'myworkspace' and name != 'PRA_final':
        del globals()[name]

# **************************************************************************
# delete all files in the temp folder
# **************************************************************************

import arcpy
import os

for the_file in os.listdir(tempdir):
    file_path = os.path.join(tempdir, the_file)
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
    except Exception as e:
        print(e)

# **************************************************************************
# delete all files in the gdb except the final PRA file
# **************************************************************************

file_list = []
for dirpath, dirnames, filenames in arcpy.da.Walk(myworkspace):
    for filename in filenames:
        file_list.append(os.path.join(dirpath, filename))

for every_file in file_list:
    try:
        if every_file != PRA_final:
            arcpy.Delete_management(every_file)
    except Exception as f:
        print (f)

print "done ..."

# **************************************************************************
# start of the validation
# **************************************************************************

# convert the final PRA file back to raseter
ras_field = "gridcode"
PRA_final_ras = myworkspace+"/"+"PRA_final_ras"
ras_method = "MAXIMUM_AREA"
ras_priority = "NONE"
ras_size = dem
arcpy.PolygonToRaster_conversion(PRA_final, ras_field, PRA_final_ras, ras_method, ras_priority, ras_size)


