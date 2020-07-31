# TCAGeneration

## 1. Background
This project is used for generating transit catchment areas (TCA) along road networks. Transit catchment areas are also known as transit service areas.  
Given a facility with a cut-off distance, its catchment/service area is a buffer area along road network, where the distance from any point within the catchment area is less than or equal to the cut-off distance.  
See below for an example. 

![TCA_example](/img/illustration_of_TCA.png "A transit catchment area")

## 2. Functions and IO

`All the input shapefiles should have a local coordinate system (i.e. meter as the unit). The program does not accept shapefile with WGS84 coordinates (i.e. latitude and longitude). `

### Input:

**Road network (shapefile of linestrings)**

The road network should already include topology information (for more information on how to build a road network topology, see [Appendix](#appendix)).  
The road network can be an undirected or directed road network.  
* *Undirected road network*: users need to specify the fields of "ID", "source", and "target"  
* *Directed road network*: users need to specify the fields of "ID", "source", "target " and "direction"  

|Field|Accepted data type|Instruction|
|----|-----|-----|
|ID|int/string|id of an edge|
|source|int/string|source node id of an edge|
|target|int/string|target node id of an edge|
|direction|int|direction of an edge: 0-double direction, 1-source to target, 2-target to source; (only required by directed road network)|

**Facility (shapefile of points)**

A facility can be represented as a point-based facility (e.g., a station represented as its station center) or a multiple-point based (MP-based) facility (e.g., a station represented as its entrances).
* *Point-based facility*: users need to specify the fields of "ID" and "cut-off"
* *MP-based facility*: users need to specify the fields of "ID", "cut-off" and "label"

|Field|Accepted data type|Instruction|
|----|-----|-----|
|ID|int/string|id of facility|
|cut-off|double|the cut-off distance of a facility|
|label|int/string|the label of a facility(only required by MP-based facilty, for point-based facilty, the ID equal to the label)|


**Parameters**

* *Nearest point searching radius*: the radius to search the nearest network edges of facilities
* *Unified cut-off*: if all the facilities have the same cut-off distance, then this unified distance needs to be specified

### Output:

* accessible edges of facilities (shapefile of linestrings) 
* catchment areas of facilities (shapefile of polygons) 

## 3. Usage: standalone program

For those who have no interest in programming, we provide a simple user interface for you to run the program.
Please follow the steps below:
1) Download folder [Standalone Exe Program](./Standalone Exe Program) and copy it to a 64-bit Windows PC.
2) Open the folder and double-click the logo "TCAGeneration.exe", you will see a GUI as below. 
Please follow the description below to generate the catchment areas. 

![TCA_GUI_instruction](/img/instruction_of_the_GUI.png "the GUI instructions")

note: We have successfully tested the TCAGeneration.exe in several 64-bit window10 PCs. Nevertheless, if you get a hint of missing DLL files during the opening of TCAGeneration.exe, you can find them under this link https://drive.google.com/drive/folders/1wq4Pkzn4j4egdiXRxGFxsGFKbhbkwdlM?usp=sharing 

## 4. Usage: rebuild the program

For those who want to rebuild the program on their PCs. 

**External libraries:**
1) GDAL (2.3.2) prebuilt library for windows: http://www.gisinternals.com/query.html?content=filelist&file=release
2) Boost (1.6.9) prebuilt library for windows: https://sourceforge.net/projects/boost/files/boost-binaries/1.69.0/
3) CGAL(4.13) https://www.cgal.org/download/windows.html
4) Qt (5.13.2) https://www.qt.io (Qt is optional, depending on if you want a GUI)

**For Windows users:**  
if you rebuild the program in visual studio, please remember to configure the 
corresponding libraries in your VS program.  
(I build the program with Visual Studio 2017 on a 64-bit Windows PC.)

**For Linux users:**  
you need to install all the external libraries and rebuild the program.  
Unfortunately, I have not find time to test the code in a Linux PC.  

## 5. Sample data

We provide two datasets for those who want to test the program, both datasets are download
from Openstreetmap(OSM). The test area is a part of Munich, Germany. Please see the folder [Test data](./Test data)

*road_network.shp*: a shapefile of undirected road network  
*subway_stations.shp*: a shapefile of point-based facilities

## 6. Contributors

Diao Lin, a Ph.D. candidate at Chair of Cartography, Technical University of Munich  
diao.lin@tum.de

Lin, Diao; Zhu, Ruoxin; Yang, Jian; Meng, Liqiu: An Open-Source Framework of Generating Network-Based Transit Catchment Areas by Walking. ISPRS International Journal of Geo-Information 9 (8), 2020, 467

## Appendix

A road network with correct topology information is very important for generating correct transit catchment areas.
It is common for researchers to use the road network provided by OSM. Below are some practical suggestions.
1. download OSM road network with topology information by using [OSMnx](https://github.com/gboeing/osmnx)
2. The data quality of OSM data can be very different depending on the study area. Therefore, it is necessary for users
to check the correctness of the obtained road topology. You can use [OpenJUMP](http://www.openjump.org/) to check some common issues:  
      a. duplicated roads (some roads occur more than once, you can delete the duplicated roads using OpenJUMP)  
      b. missing intersections (two roads intersect with each other, but there are not connected)  
3. In case you find a lot of missing intersections in your data, we provide you a tool to rebuild the road network topology.
please download the tool: [TopoBuilder](https://gitlab.com/Drsulmp/topobuilder). 
Note, this tool can be used to build the undirected road network. For those who want to build a direct road network, you may need to add the missing
intersection nodes manually.




