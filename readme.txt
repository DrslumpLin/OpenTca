# TCAGeneration
## 1. Background
This project is used for generating transit catchment areas (TCA) along road networks. Transit catchment areas are also known as transit service areas. A catchment/service area is a buffer area around a facility, see below for an example.

## 2. Functions and IO

The program requires a road network and a f


### Input:

**Road network** (shp file of linestrings)
```
The road network should already include the topology information.
The road network can be a undirected or directed road network.
For undirected road network: users need to specify the fileds of "ID", source", "target"
For directed road network: users need to specify the fileds of "ID", "source", "target " and "direction"
```
**Facilities** (shp file of points)
```
The facility can be represented as point-based facility or multiple-point-based(MP-based) facility.
For point-based facility: users need to specify the fileds of "ID", "cut-off"
For Mp-based facility: users need to specify the fileds of "ID", "cut-off", "label".
```
**Parameters**
```
Nearest point searching radius: the radius to search the nearest network edges 
Unified cut-off, if you donot specify the field of "cut-off", then you need to specify a unfield cut-off for all the facilities
```

### Output:
```
1) accessibile edges of facilities (shp file of linestrings) 
2) cacthment areas of facilities (shp file of polygons) 
```

## 3. Usage:standalone program

For those who have no interest in programing, we provide a simple user interface for you to run the program.
Please follow steps below:
1. Download folder "Standalone Exe Program" and copy it to a 64-bit Windows PC.
2. Open the folder and double-click the logo "TCAGeneration.exe", you will see a GUI as below. 
Please follow the description to generate the cacthment areas. 

## 4. Usage:rebuild the program

For those who want to rebuild the program on their PCs.  

**External libraries:**
1) GDAL (2.3.2) prebulit library for windows: http://www.gisinternals.com/query.html?content=filelist&file=release
2) Boost (1.6.9) prebuilt library for windows: https://sourceforge.net/projects/boost/files/boost-binaries/1.69.0/
3) CGAL(4.13) https://www.cgal.org/download/windows.html
4) Qt (5.13.2)  https://www.qt.io

**Platform:** 
   


**For windows users:**
if you rebuild the program in visual studio, please remember to configure the 
corresponding libraries in your VS program.  
(I build the program with Visual studio s2017 on a 64-bit Windows PC.)

**For linux users:**
you need to install all the external libraries and rebuild the program.  
Unfortunately, I havenot find time to test the code in a Linux PC.

## 5. Sample data

We provide two datasets for those who want to test the program, both dataset are download
from Openstreetmap(OSM). The test area is part of Munich, Germany. 

A shp file of undirected road network: 
A shp file of point-based facilities: 

## 5. Contributors

Diao Lin, a PhD candidate at Chair of Cartography, Technical university of Munich  
diao.lin@tum.de
