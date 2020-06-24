# TCAGeneration

## 1. Background
This project is used for generating transit catchment areas (TCA) along road networks. Transit catchment areas are also known as transit service areas.  
Given a facility and with a cut-off distance, its catchment/service area is a buffer area along road network, where the distance from any point within the catchment area is less than or equal to the cut-off distance.  
See below for an example. 

![TCA_example](/img/illustration_of_TCA.png "A transit catchment area")

## 2. Functions and IO

### Input:

**Road network (shapefile of linestrings)**

The road network should already include topology information (for those who have no idea on how to build a road network topology, see appendix).  
The road network can be an undirected or directed road network.  
* *Undirected road network*: users need to specify the fields of "ID", "source", and "target"  
* *Directed road network*: users need to specify the fields of "ID", "source", "target " and "direction"  

|Field|Type of data|Instruction|
|----|-----|-----|
|ID|int/string|id of an edge|
|source|int/string|source node id of an edge|
|target|int/string|target node id of an edge|
|direction|int|direction of an edge: 0-double direction, 1-source to target, 2-target to source; (only required by directed road network)|

**Facility (shapefile of points)**

A facility can be represented as a point-based facility (e.g., a station represented as its station center) or a multiple-point based (MP-based) facility (e.g., a station represented as its entrances).
* *Point-based facility*: users need to specify the fields of "ID" and "cut-off"
* *MP-based facility*: users need to specify the fields of "ID", "cut-off" and "label"

|Field|Type of data|Instruction|
|----|-----|-----|
|ID|int/string|id of facility|
|cut-off|double|the cut-off distance of a facility|
|label|int/string|the label of a facility(only required by MP-based facilty, for point-based facilty, the ID equal to the label)|


**Parameters**

* *Nearest point searching radius*: the radius to search the nearest network edges of facilities
* *Unified cut-off*: if all the facilities have the same cut-off distance, then this unified distance needs to be specified

### Output:

* accessibile edges of facilities (shapefile of linestrings) 
* catchment areas of facilities (shapefile of polygons) 

## 3. Usage: standalone program

For those who have no interest in programming, we provide a simple user interface for you to run the program.
Please follow steps below:
1) Download folder "/Standalone Exe Program" and copy it to a 64-bit Windows PC.
2) Open the folder and double-click the logo "TCAGeneration.exe", you will see a GUI as below. 
Please follow the description to generate the catchment areas. 

![TCA_GUI_instruction](/img/instruction_of_the_GUI.png "the GUI instructions")


## 4. Usage: rebuild the program

For those who want to rebuild the program on their PCs. 

**External libraries:**
1) GDAL (2.3.2) prebulit library for windows: http://www.gisinternals.com/query.html?content=filelist&file=release
2) Boost (1.6.9) prebuilt library for windows: https://sourceforge.net/projects/boost/files/boost-binaries/1.69.0/
3) CGAL(4.13) https://www.cgal.org/download/windows.html
4) Qt (5.13.2) https://www.qt.io (Qt is optional, depending on if you want a GUI)

**For Windows users:**  
if you rebuild the program in visual studio, please remember to configure the 
corresponding libraries in your VS program.  
(I build the program with Visual studio 2017 on a 64-bit Windows PC.)

**For linux users:**
you need to install all the external libraries and rebuild the program.  
Unfortunately, I have not find time to test the code in a Linux PC.  

## 5. Sample data

We provide two datasets for those who want to test the program, both dataset are download
from Openstreetmap(OSM). The test area is part of Munich, Germany. 

A shp file of undirected road network: 
A shp file of point-based facilities: 

## 6. Contributors

Diao Lin, a Ph.D. candidate at Chair of Cartography, Technical University of Munich  
diao.lin@tum.de
