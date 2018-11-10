# Automated mapping of potential snow avalanche release areas (PRAs)

Thalia Bertschinger, 12-105-607

thalia.bertschinger@students.unibe.ch

Seminar Geodata Analysis and Modelling, Spring Semester 2018

-----------------

One important advice: __DO NOT RUN THE WHOLE SCRIPT AT ONCE!__
The whole script until line 606 (including the validation) can be run as a whole.
After execution, run the last part to delete unnecessary files in the geodatabase.

-----------------


This is the explanation how to use the Python script.

After the intro part, the script has three main parts: 
1. The creation of the PRA_final file with the final PRA polygons (until line 488)
2. The validation (from line 489 until line 579)
3. The deletion of unnecessary files (from line 580 to the end)

These things have to be adjusted before the execution:
* Line 38: set a preworkspace folder
* Line 42: define the name for the file geodatabase (gdb)
* Line 49: adapt the name of the DEM (the DEM should be saved in the preworkspace folder)
* Line 67 to 74: define the parameters
* Line 80: define the reference data set (has to be polygons) (if needed, change rast_field in line 503)

The reference data set needs to be a vector file with polygons that have a field in the attribute table with 1 for PRAs and 2 for non-PRAs.
If no reference data set is available, run the script without the validation part.

If the last part of the script is used, the temporary folder and all the files in the geodatabase except the final PRA file will be deleted.
The deletion of unwanted files has two parts, one for the temporary folder and one for the geodatabase.
Sometimes an error message appears when the whole script together with this very last part is executed and PyCharm is stuck.
As it is nice to have no unnecessary files left after the analysis, but as the problem with the error message could not be solved, it is recommended to run the last part of the script separately.

The model creates polygons of potential snow avalanche release areas with characteristical parameters of each PRA in the attribute table.