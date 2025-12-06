# dataPreprocess.py
#
#
# Author: Laura Tateosian 11.20.2019
# Modified: Laura Tateosian 04.21.2023
#           Migrate from Desktop to Pro
# Modified: Laura Tateosian 04.24.2023
#           Add Color class.
# Purpose:  Create an animated gif showing the coloring the lowest income counties colored progressively from a given starting point. 
#
#
# Procedure Summary:  
#           --Take a state census tract shapefile, a county name, 
#           a low income by census tract table, indicating a binary
#           LIC value (low income or not), and a starting point.  
#           --Select the tracts within the specified county, remove the rest,
#           color the low income tracts light green, and take a screenshot.
#           --Find the tract containing the starting point, color it dark green,
#           and take a screenshot.
#           --Find the next closest low income tract, color it dark green,
#           and take a screenshot.  Repeat until all low income tracks are dark green.
#           --Compile into animated gif.
#
#
# Main Steps:
# Step 1: Tract data has the entire state.
# Get a dataset with just the current county's records.
#
# Step 2: Join the LIC data to the polygon table.
#
# Step 3: Create and initialize a flag field, named "projAdopt", for symbology
#
# Symbology is set to clear with black border, light green, and dark green
# Values: -1 for not eligible for the program;  (clear with black border)
#          0 for low income but not in the program (light green)
#          1 for low income and in the program (dark green)
# Step 4: Create the first image, with all low income tracts as light green.
#
# Step 5: Capture the image with the starting track dark green (image number 2)
#
# Step 6: Grow the program. Find the nearest eligible tract
# still not in the program, update its projAdopt value to 1,
# set it to dark green, and capture screenshot.
# Iteratively repeat until all low income tracts are in the program.
#
# Step 7: Create animation.
#
# Step 8: Create an HTML report.
#
# Software Requirements: ArcGIS Pro arcpy package and imageio package must be installed.
#
#
# Usage: base_directory state_tract_shapefile_path county_low_income_data_csv_path county_name DD_longitude_of_start_point DD_latitude_of_start_point county_FIPS_code
#
#
# Example input: C:/Users\Laura/Documents/forBart/Guilford/ C:/Users\Laura/Documents/forBart/Guilford/tl_2018_37_tract.shp C:/Users\Laura/Documents/forBart/Guilford/GuilfordLIC.csv Guilford -79.993108 35.974236 081

import arcpy
import imageio
import os
import shutil
import sys

class Color:
    """
    A class to represent an RGBA color.
    ...

    Attributes
    ----------
    R : int
        red setting (any number between 0 and 255)
    G : int
        green setting (any number between 0 and 255)
    B : int
        blue setting (any number between 0 and 255)
    A : int
        alpha setting (any number between 0 and 100)
    rgb: dict
        {'RGBA': R, G, B, A}

    Methods
    -------
    set_transparency(self, description):
        Set the transparency based on a descriptive term.
    """
    def __init__(self, red=255, green=255, blue=255, alpha=100):
        """Instantiate a color object and set its red, green, blue, and alpha attributes."""
        self.R = red
        self.G = green
        self.B = blue
        self.A = alpha
        self.rgb = {'RGB': [self.R, self.G, self.B, self.A]}

    def set_transparency(self, description):
        """Set the transparency based on a descriptive term."""
        transp_dict = {"opaque":100, "slight":95, "half":50, "mostly":5, "full":0}
        try:
            self.A = transp_dict[description.lower()]
            self.rgb = {'RGB': [self.R, self.G, self.B, self.A]}
        except KeyError:
            print(f"{description.lower()} is not an available transparency setting")
            print(f"Valid options: {transp_dict.keys()}")

def printArc(message):
    """Print message for script tool and standard output."""
    print(message)
    arcpy.AddMessage(message)


def printArgs():
    """Print user arguments."""
    printArc('Number of arguments = {0}'.format(len(sys.argv)))
    for index, arg in enumerate(sys.argv):
        printArc('Argument {0}: {1}'.format(index, arg))


def cleanADir(thePath):
    """Delete all contents of a directory."""
    contents = os.listdir(thePath)
    for f in contents:
        theFile = os.path.join(thePath, f)
        try:
            if os.path.isfile(theFile):
                os.unlink(theFile)
            elif os.path.isdir(theFile):
                shutil.rmtree(theFile)
        except Exception as e:
            print(e)


def getNumberedFileName(fileNum, theDir, word, numDigits):
    """Create a file name within a series as a word followed
    by a number this is zfilled to have numDigits digits."""
    num = str(fileNum).zfill(numDigits)
    fileName = theDir + "{0}{1}.png".format(word, num)
    return fileName


def makeADir(thePath):
    """Create the directory, if it doesn't already exist."""
    if not os.path.exists(thePath):
        os.mkdir(thePath)


def set_map_title(layoutObj, name):
    """Replace xxx in the title for the given map to the
       given name & center the title horizontally."""
    # elems = arcpy.mapping.ListLayoutElements(mapDocObj)
    mf = layoutObj.listElements("MAPFRAME_ELEMENT", "county*")[0]
    text_elems = layoutObj.listElements("TEXT_ELEMENT", "*title*")
    for e in text_elems:
        # if e.name == 'title':
        e.text = e.text.replace("xxx", name)
        e.elementPositionX = mf.elementPositionX + \
                             (mf.elementWidth * 0.5) - (e.elementWidth * 0.5)


def set_three_colors(symb, the_field, field_values=("Not eligible","Eligible","Enrolled")):
    """Set layer symbology to 3 unique colors based on the given field in the data layer"""

    # Set up Color objects for each category.
    excluded = Color()
    excluded.set_transparency("full")
    excluded_outline = Color(210, 210, 210)

    level1 = Color(161, 217, 155)
    level1.set_transparency("slight")
    level2 = Color(49, 163, 84, level1.A)

    symb.updateRenderer("UniqueValueRenderer")
    symb.renderer.fields = [the_field]

    for grp in symb.renderer.groups:
        for item in grp.items:
            value = item.values[0][0]
            if value == field_values[0]:
                item.symbol.color = excluded.rgb
                item.symbol.outlineColor = excluded_outline.rgb
            if value == field_values[1]:
                item.symbol.color = level1.rgb
            if value == field_values[2]:
                item.symbol.color = level2.rgb

    #item.label = str(value)
    return symb

# Report user arguments when running in the GUI or the IDE.
printArgs()

# Get the base directory based on the script's relative location.
scriptPath = sys.argv[0]
codePath = os.path.dirname(scriptPath)
baseDir = os.path.dirname(codePath) + "/"
# Get user input or hard-code test values.
try:
    tractPolys = sys.argv[1]
    LICtable = sys.argv[2]
    county = sys.argv[3]
    X = sys.argv[4]
    Y = sys.argv[5]
    countyFIPScode = sys.argv[6]
except IndexError:
    scriptPath = sys.argv[0]
    codePath = os.path.dirname(scriptPath)
    baseDir = os.path.dirname(codePath) + "/"
    # baseDir = r"Q:\My Drive\GIS540\sharing\environ_edu/"
    # r"C:\Users\Laura\Documents\forBart\Guilford/"

    # 2018 polygonal U.S. Census tract data from U.S. Census web interface for selecting data
    # https://www.census.gov/cgi-bin/geo/shapefiles/index.php  
    tractPolys = baseDir + "data/tl_2018_37_tract.shp"
    county = "Guilford"
    LICtable = f"{baseDir}data/{county}LIC.csv"  # "data/{0}LIC.csv".format(county)
    X = -79.993108
    Y = 35.974236
    # Data provides "Federal Information Processing Standard" (FIPS) code, not county name
    # Lookup table for FIPS codes https://www.lib.ncsu.edu/gis/countfips
    # Guilford County code: 37081 #stateFIPScode = 37
    countyFIPScode = "081"

# Specify the template project.
projectPath = baseDir + "environ_edu.aprx"

# Prepare the output workspace.
outputDir = baseDir + "output/"
imageDir = outputDir + "images/"
makeADir(outputDir)
makeADir(imageDir)
cleanADir(imageDir)

# Set geoprocessing environment.
arcpy.env.overwriteOutput = True
arcpy.env.workspace = outputDir

# ------------------
# Step 1: Tract data has the entire state.
# Get a dataset with just the current county's records.
tmp = 'tmpLayer'
sqlClause = "CountyFP = '{0}'".format(countyFIPScode)
thisCounty = outputDir + "{0}Tracts.shp".format(county)

arcpy.MakeFeatureLayer_management(in_features=tractPolys,
                                  out_layer=tmp)
arcpy.SelectLayerByAttribute_management(in_layer_or_view=tmp,
                                        where_clause=sqlClause)
arcpy.CopyFeatures_management(tmp, thisCounty)

# ------------------
# Step 2: Join the LIC data to the polygon table.
fieldObjects = arcpy.ListFields(thisCounty)
fieldNames = [f.name for f in fieldObjects]

licField = "LIC"
if licField not in fieldNames:
    LICtableView = "LIC_TV"
    # Add the LIC (low income rating) column onto the shapefile.
    # The tract values in the csv correspond to
    # the GEOID field in the census data,
    # so we use these to match the rows.
    arcpy.JoinField_management(in_data=thisCounty,
                               in_field="GEOID",
                               join_table=LICtable,
                               join_field="Tract",
                               fields=[licField])
    printArc(arcpy.GetMessages())

# ------------------
# Step 3: Create and initialize a flag field, named "projAdopt", for symbology
# Symbology is set to clear with black border, light green, and dark green
# Values: -1 for not eligible;  (clear with black border)
#          0 for eligible--low income but not in the program (light green)
#          1 for enrolled--low income and in the program (dark green)

# Add the field
flagField = "progAdopt"
categField = "Status"
if flagField not in fieldNames:
    arcpy.AddField_management(in_table=thisCounty,
                              field_name=flagField,
                              field_type="SHORT")

if categField not in fieldNames:
    arcpy.AddField_management(in_table=thisCounty,
                              field_name=categField,
                              field_type="TEXT")

# Default value for a SHORT field is zero, so if the field was just created,
# This step is not needed.
arcpy.CalculateField_management(in_table=thisCounty,
                                field=flagField,
                                expression=0)

arcpy.CalculateField_management(in_table=thisCounty,
                                field=categField,
                                expression='"Eligible"')

# If LIC is zero, set the display value to -1
# (it's not considered for the low income program).
with arcpy.da.UpdateCursor(in_table=thisCounty,
                           field_names=[flagField, categField],
                           where_clause="{0}=0".format(licField)) as uc:
    for row in uc:
        row[0] = -1
        row[1] = "Not eligible"
        uc.updateRow(row)
    del uc


# ------------------
# Step 4: Create the first image, with all low income tracts as light green.

# Get the approximate number of images needed so you can zfill the names.
tractCount = arcpy.GetCount_management(thisCounty)
digits = len(str(tractCount))

# Create the map document object.
# mxd = arcpy.mapping.MapDocument(mapPath)

# Create an ArcGISProject object
try:
    # When working inside ArcGIS Pro,
    # use the project name "CURRENT"
    aprx = arcpy.mp.ArcGISProject("CURRENT")

except OSError:
    # When working outside ArcGIS Pro,
    # use the project full path file name.
    aprx = arcpy.mp.ArcGISProject(projectPath)

### Customize the title on the map
# dfs = arcpy.mapping.ListDataFrames(mxd)
# df = dfs[0]
theMaps = aprx.listMaps()
myMap = theMaps[0]

myLayout = aprx.listLayouts()[0]

set_map_title(myLayout, county)  # customer decided against map title

# Symbology lyr created by hand by the following procedure:
# Create shapefile with projAdopt column set to mostly -1, 0, and some 1
# By hand, set the symbology to clear with outline, light green, and dark green
#   for -1, 0, and 1 respectively.
# http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
# Used colorbrewer to select colors (although customer later requested change)
# Converted hex values given in colorbrewer to rgb values.
# https://www.rapidtables.com/convert/color/hex-to-rgb.html
# By hand, set the transparency to 30%
# to set transparency, by hand, set layer Properties -> Display
# Specify the symbology template
symbLay = baseDir + "/data/symbology/polygonsForSymbology.lyr"

# Add the county's data to the map.
myLayer = myMap.addDataFromPath(thisCounty)
sym = myLayer.symbology
myLayer.symbology = set_three_colors(sym, categField)

#arcpy.ApplySymbologyFromLayer_management(in_layer=layObj,
#                                         in_symbology_layer=symbLay)

# Capture the no program image (image number 1)
# Initialize counter to track image number.
count = 1
imageName = getNumberedFileName(count, imageDir, county, digits)
lyrs = myMap.listLayers()

myLayout.exportToPNG(imageName, resolution=50)
#lyrs[0].visible = False
myMap.removeLayer(lyrs[0])


# ------------------
# Step 5: Capture the image with the starting track dark green (image number 2)

# Thanks for code from bixb0012 on
# https://community.esri.com/thread/215640-find-a-polygon-for-a-single-point-in-python-script
# to get the select layer by location tool working by first making a feature layer
# and then using the results with an update cursor.
startPoint = arcpy.Point(X, Y)  # (X,Y) is the user input starting point.
spatialRef = arcpy.Describe(thisCounty).spatialReference
polySelection = arcpy.MakeFeatureLayer_management(thisCounty, "polygonsToSelect")

inputPointGeometry = arcpy.PointGeometry(startPoint, spatialRef)
arcpy.SelectLayerByLocation_management(
    in_layer=polySelection,
    overlap_type="INTERSECT",
    select_features=inputPointGeometry,
    search_distance="",
    selection_type="NEW_SELECTION"
)

# polySelection is now the tract with the starting point.
# Set the first tract's status as enrolled and get its geoid and name.
with arcpy.da.UpdateCursor(polySelection, [flagField, categField, "geoid", "NAMELSAD"]) as uc:
    # Get the row (there should be only one, since the start point can only lie within a single tract).
    row = next(uc)

    # Update the tract's status to enrolled.
    row[0] = 1
    row[1] = "Enrolled"
    uc.updateRow(row)

    # Get the start tract GEOID to use in the near analysis.
    startTract = row[2]
    # Get the start tract's name to use in the report.
    first_tract = row[3]

    del uc

# Update the counter to track image number
# and get the image file name.
count = count + 1
imageName = getNumberedFileName(count, imageDir, county, digits)

# Add data to the map, capture screenshot, and remove data again.
myLayer = myMap.addDataFromPath(thisCounty)
sym = myLayer.symbology
myLayer.symbology = set_three_colors(sym, categField)
myLayout.exportToPNG(imageName, resolution=50)
#myLayer.visible = False
myMap.removeLayer(myLayer)

# ------------------
# Step 6: Grow the program. Find the nearest eligible tract
# still not in the program, update its projAdopt value to 1,
# set it to dark green, and capture screenshot.
# Iteratively repeat until all low income tracts are in the program.

# Count the number of remaining low income tracts.

with arcpy.da.SearchCursor(thisCounty, ["FID"], "{0}=1 AND {1}<>1".format(licField, flagField)) as sc:
    numStillNotProgram = len([row[0] for row in sc])
    del sc

arcpy.SetProgressor("step", "Creating maps...", 0, numStillNotProgram, 1)

notInProj = arcpy.MakeFeatureLayer_management(thisCounty, "notInProj")
while numStillNotProgram > 0:
    # Update the progressor message and position
    arcpy.SetProgressorLabel(f"Creating map image {count}...")
    arcpy.SetProgressorPosition(count)

    # Select the tract(s) that are low income, but not yet in the program.
    arcpy.SelectLayerByAttribute_management(
        in_layer_or_view=notInProj,
        selection_type="NEW_SELECTION",
        where_clause="{0}=1 AND {1}=0".format(licField, flagField)
    )
    # Get the closest to the starting point low income tract that's not in the program yet.
    # polySelection is the tract with the starting point.

    arcpy.Near_analysis(in_features=polySelection,
                        near_features=notInProj)

    # Get the ID of the nearest point.
    with arcpy.da.SearchCursor(polySelection, ["geoid", "NEAR_FID"]) as sc:
        for row in sc:
            theNearestOne = row[1]
        del sc

    with arcpy.da.UpdateCursor(thisCounty, [flagField, categField], "FID={0}".format(theNearestOne)) as uc:
        for row in uc:
            row[0] = 1
            row[1] = "Enrolled"
            uc.updateRow(row)
        del uc

    # Update the counter to track image number
    # and get the image file name.
    count = count + 1
    imageName = getNumberedFileName(count, imageDir, county, digits)

    # Add data to the map, capture screenshot, and remove data again.
    # By hand, set the symbology to clear with outline, light green, and dark green
    #   for -1, 0, and 1 respectively.
    myLayer = myMap.addDataFromPath(thisCounty)
    sym = myLayer.symbology
    myLayer.symbology = set_three_colors(sym, categField)

    myLayout.exportToPNG(imageName, resolution=50)
    #myLayer.visible = False
    myMap.removeLayer(myLayer)
    # Decrement the count, since we've added one more to the program.
    numStillNotProgram = numStillNotProgram - 1
    
# Delete the ArcGISProject object.
aprx.saveACopy(outputDir + "/mapCopy.aprx")
del aprx

arcpy.ResetProgressor()
arcpy.SetProgressor("default", "Creating animation...")


# ------------------
# Step 7: Create animation.
# Note, imageio library must be installed.  At command line, pip install imageio
# Thanks for code from Almar for gif creation with imageio at
# https://stackoverflow.com/questions/753190/programmatically-generate-video-or-animated-gif-in-python
# Thanks for code from Gwen for frame duration adjustment at
# https://stackoverflow.com/questions/38433425/custom-frame-duration-for-animated-gif-in-python-imageio

images = []
filenames = os.listdir(imageDir)
filenames = [imageDir + f for f in filenames if f.endswith("png")]
for filename in filenames:
    images.append(imageio.imread(filename))
gifName = outputDir + "/" + county + 'CountyMovie.gif'
imageio.mimsave(gifName, images, format='GIF', duration=0.6)
print("{} created!".format(gifName))

arcpy.SetProgressorLabel(f"Creating HTML report...")

# ------------------
# Step 8: Create an HTML report.
# Create an HTML report of the visual output and an explanation of the results.

total_tracts = int(arcpy.GetCount_management(thisCounty).getOutput(0))
summary = f"{count}/{total_tracts} ({100*count/total_tracts:.2f}%) of the tracts in {county} County are considered eligible for the program."
about = f"Mapping the program implemented first in {first_tract} and then added to its nearest neighbor, and then added to the neighbor that is nearest either of these tracts, and so forth."
html_content = f"""<!DOCTYPE html>
<html>
    <body>

    <h1 style="font-family:verdana;color:gray;">Visualizing Environmental Education</h1>
 
    <figure>
          <img src="{gifName}" alt="Animated gif of program being incrementally implemented in {county} county." style="width:425px;height:550px">
          <figcaption style=width:425px;font-family:verdana;color:gray;">{about} {summary}</figcaption>
    </figure>

    </body>
</html>"""
html_file = outputDir + f"/{county}_County_report.html"
with open(html_file, 'w') as hout:
    hout.write(html_content)

# Automatically open the prepared html report.
os.startfile(html_file)