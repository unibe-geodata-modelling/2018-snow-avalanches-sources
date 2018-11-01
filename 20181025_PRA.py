# Script for the automated mapping of potential snow avalanche release areas
# Thalia Bertschinger
# Seminar Geodata Analysis and Modelling, Spring Semester 2018
# The model creates polygons of potential snow avalanche release areas

# imports
import arcpy
import os
import numpy

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
gdb = "20181101_b_Model_PRA.gdb"
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

slope = tempdir+"/"+"slope.tif"  # temporary file for the slope analysis
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
outSlope.save(slope)

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

# **************************************************************************
# save the final PRA file
# **************************************************************************

# rename the last output of the loop to get the final PRA file
rename_PRA4 = path_PRA4_elim[-1]
# PRA_final = myworkspace + "/" + "PRA_final"
PRA_final = os.path.join(myworkspace, "PRA_final")
arcpy.Rename_management(rename_PRA4, PRA_final)

# add a field in the attribute table that will contain 1 for PRA and 2 for no PRA
field_PRA = "PRA"
type_PRA = "SHORT"
arcpy.AddField_management(PRA_final, field_PRA, type_PRA)
# no PRAs will get the value 2 and PRAs will get the value 1
expression_PRA = "getPRA(!gridcode!)"
code_PRA = """
def getPRA(gridcode):
    if gridcode == 0:
        return 2
    else:
        return 1"""
arcpy.CalculateField_management(PRA_final, field_PRA, expression_PRA, "PYTHON", code_PRA)

# **************************************************************************
# assign characteristical parameters to each PRA
# **************************************************************************

# add a field in the attribute table with consecutive numbers for each feature
field_PRA = "PRA_nr"
type_PRA = "SHORT"
arcpy.AddField_management(PRA_final, field_PRA, type_PRA)
# assign consecutive numbers for each feature
expression_PRA = "autoIncrement()"
code_PRA = """
rec=0 
def autoIncrement():
    global rec 
    pStart = 1  
    pInterval = 1 
    if (rec == 0):
        rec = pStart  
    else:  
        rec += pInterval  
    return rec"""
arcpy.CalculateField_management(PRA_final, field_PRA, expression_PRA, "PYTHON", code_PRA)

# calculate min, max, mean altitude and slope with zonal statistics and write each value in a new field to PRA final
# define the function add_zonal_field

from arcpy.sa import *


def add_zonal_field(PRA_final, zone_field, input_raster, field_name, field_type, stat, stat_field):

    """ Performs zonal statistics on a set of features,
    then adds and populates a new field in the feature class
    to store the results of the zonal calculations"""

    # Add zone field to features
    arcpy.AddField_management(PRA_final, field_name, field_type)

    # Clear path for temporary table
    if arcpy.Exists("zonal_table"):
        try:
                arcpy.Delete_management("zonal_table")
        except:
                arcpy.AddError("Unable to clear temp table")
                sys.exit(-1)

    # Zonal statistics
    arcpy.AddMessage("Performing zonal statistics: " + field_name + " in " + PRA_final)

    zonal_table = ZonalStatisticsAsTable(PRA_final, zone_field, input_raster, "zonal_table", "DATA", stat)

    # Digest statistics from zonal_table
    arcpy.AddMessage("Digesting " + stat + " from zonal table")

    stat_dict = {}

    with arcpy.da.SearchCursor(zonal_table, [zone_field, stat_field]) as cursor:
            for row in cursor:
                stat_dict[row[0]] = row[1]

    # update new field in feature class
    arcpy.AddMessage("Calculating " + field_name + "in " + PRA_final)

    with arcpy.da.UpdateCursor(PRA_final, [zone_field, field_name]) as cursor2:
        for row2 in cursor2:
            row2[1] = stat_dict[row2[0]]
            cursor2.updateRow(row2)


# minimum altitude
zone_field = "PRA_nr"
input_raster = arcpy.Raster(dem_5m)
field_name = "min_alt"
field_type = "FLOAT"
stat = "MINIMUM"
stat_field = "MIN"

add_zonal_field(PRA_final, zone_field, input_raster, field_name, field_type, stat, stat_field)

# maximum altitude
zone_field = "PRA_nr"
input_raster = arcpy.Raster(dem_5m)
field_name = "max_alt"
field_type = "FLOAT"
stat = "MAXIMUM"
stat_field = "MAX"

add_zonal_field(PRA_final, zone_field, input_raster, field_name, field_type, stat, stat_field)

# mean altitude
zone_field = "PRA_nr"
input_raster = arcpy.Raster(dem_5m)
field_name = "mean_alt"
field_type = "FLOAT"
stat = "MEAN"
stat_field = "MEAN"

add_zonal_field(PRA_final, zone_field, input_raster, field_name, field_type, stat, stat_field)

# minimum slope
zone_field = "PRA_nr"
input_raster = arcpy.Raster(tempdir+"/"+"slope.tif")
field_name = "min_slope"
field_type = "FLOAT"
stat = "MINIMUM"
stat_field = "MIN"

add_zonal_field(PRA_final, zone_field, input_raster, field_name, field_type, stat, stat_field)

# maximum slope
zone_field = "PRA_nr"
input_raster = arcpy.Raster(tempdir+"/"+"slope.tif")
field_name = "max_slope"
field_type = "FLOAT"
stat = "MAXIMUM"
stat_field = "MAX"

add_zonal_field(PRA_final, zone_field, input_raster, field_name, field_type, stat, stat_field)

# mean slope
zone_field = "PRA_nr"
input_raster = arcpy.Raster(tempdir+"/"+"slope.tif")
field_name = "mean_slope"
field_type = "FLOAT"
stat = "MEAN"
stat_field = "MEAN"

add_zonal_field(PRA_final, zone_field, input_raster, field_name, field_type, stat, stat_field)

# **************************************************************************
# start of the validation
# **************************************************************************

# convert the final PRA file back to raster
ras_field = "PRA"
PRA_final_ras = myworkspace+"/"+"PRA_final_ras"
ras_method = "MAXIMUM_AREA"
ras_priority = "NONE"
ras_size = dem
arcpy.PolygonToRaster_conversion(PRA_final, ras_field, PRA_final_ras, ras_method, ras_priority, ras_size)

# convert the reference data set to raster
reference = "U:/Seminar_Modellieren/20181018_Test_Model/reference_release_areas.shp"
rast_field = "RA"
reference_ras = myworkspace+"/"+"reference_ras"
arcpy.PolygonToRaster_conversion(reference, rast_field, reference_ras, ras_method, ras_priority, ras_size)

# generate variables with the PRA raster and the reference raster
PRA = arcpy.Raster(PRA_final_ras)
ref = arcpy.Raster(reference_ras)

# calculate the error matrix raster
pre_error_matrix_ras = arcpy.sa.Con(
    (PRA == 1) & (ref == 1), 1, arcpy.sa.Con(
        (PRA == 1) & (ref == 2), 2, arcpy.sa.Con(
            (PRA == 2) & (ref == 1), 3, arcpy.sa.Con(
                (PRA == 2) & (ref == 2), 4))))
pre_error_matrix_ras.save(myworkspace + "/" + "pre_error_matrix_ras")

# set noData values to 999
error_matrix_lyr = "error_matrix_lyr"
arcpy.MakeRasterLayer_management(os.path.join(myworkspace, "pre_error_matrix_ras"), error_matrix_lyr)
error_matrix_ras_final = arcpy.sa.Con(arcpy.sa.IsNull(error_matrix_lyr), 999, error_matrix_lyr)
error_matrix_ras_final.save(myworkspace + "/" + "error_matrix_ras_final")

# generate a numpy array out of the error matrix raster
error_matrix_arr = arcpy.RasterToNumPyArray(error_matrix_ras_final)

# count the number of each of the four error matrix options by looping through the whole array
rows = numpy.shape(error_matrix_arr)[0]
cols = numpy.shape(error_matrix_arr)[1]
i = 0
a_error_matrix = 0
b_error_matrix = 0
c_error_matrix = 0
d_error_matrix = 0
noData = 0
# first loop through all rows in array
while i < rows:
    j = 0
    # second loop through all columns in row
    while j < cols:
        if error_matrix_arr[i, j] == 1:
            a_error_matrix += 1
        elif error_matrix_arr[i, j] == 2:
            b_error_matrix += 1
        elif error_matrix_arr[i, j] == 3:
            c_error_matrix += 1
        elif error_matrix_arr[i, j] == 4:
            d_error_matrix += 1
        elif error_matrix_arr[i, j] == 999:
            noData += 1
        j += 1
    i += 1

if (rows * cols) == (a_error_matrix + b_error_matrix + c_error_matrix + d_error_matrix + noData):
    print "a_error_matrix: " + str(a_error_matrix)
    print "b_error_matrix: " + str(b_error_matrix)
    print "c_error_matrix: " + str(c_error_matrix)
    print "d_error_matrix: " + str(d_error_matrix)
    print "noData: " + str(noData)
else:
    print "something went wrong"

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
