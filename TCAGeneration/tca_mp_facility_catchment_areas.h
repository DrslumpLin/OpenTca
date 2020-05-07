// The abbreviation "MP" means "multiple-points", indicating the facility is represented by multiple-points facility
// An intergated interface for generating the catchment areas for multiple-points based facilities
// @author: Diao Lin
// @version: 2020.03

#pragma once
#ifndef TCA_MP_FACILITY_CATCHMENT_AREAS
#define TCA_MP_FACILITY_CATCHMENT_AREAS


#include "tca_network.h"
#include "tca_undirected_graph.h"
#include "tca_directed_graph.h"
#include "tca_build_triangulation.h"
#include "tca_tri_contour.h"
#include <chrono> 
using namespace std::chrono;

namespace tca 
{
// using boost for polygon dissolving
namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::polygon<point_t> polygon_t;

class MPFacilityCatchmentAreas
{
public:
	MPFacilityCatchmentAreas();

	~MPFacilityCatchmentAreas();
	
	// read facility and road network and generating the catchment areas
	// to variables faci_contours_ and faci_contour_labels_
	// input:
	//	@froads: path of the roads 
	//	@ffacilities: path of the facilities
	void CalculateCatchmentAreas(std::string froads, std::string ffacilities);

	// read facility and (undirected) road network and generating the accessible edges and catchment areas
	// for undirected graph and multiplepoint-based facility
	// input:
	//	@froads: path of the roads 
	//	@ffacilities: path of the facilities
	//  @netid, @source_name, @target_name: fileds for undirected graph-based network
	//  @fid,@flabel, @isunifiedcutoff,@unifiedcutoff,@cutoff_name: parameters for multiple-point based facility
	//  @searchradius: the radius for searching the nearest points for each facility
	//  @writeaccedges, @writecatchment: output accessible edges and catchment polygons
	void CalculateCatchmentAreas(const std::string froads, const std::string &netid, const std::string &source_name, const std::string &target_name,
		const std::string ffacilities, const std::string &fid, const std::string &flabel, const bool &isunifiedcutoff, const double &unifiedcutoff, const std::string &cutoff_name,
		const double &searchradius, bool writeaccedges, bool writecatchment);


	// read facility and road network and generating the accessible edges and catchment areas
	// directed road network and point-based facility
	// input:
	//  @from_facility: if the cacthment is measured from the direction of from-facility, otherwise, from to-facility direction
	//	@froads: path of the roads 
	//	@ffacilities: path of the facilities
	//  @netid, @source_name, @target_name: fileds for undirected graph-based network
	//  @fid, @isunifiedcutoff,@unifiedcutoff,@cutoff_name: parameters for point-based facility
	//  @searchradius: the radius for searching the nearest points for each facility
	//  @writeaccedges, @writecatchment: output accessible edges and catchment polygons
	void CalculateCatchmentAreasDirected(const bool &from_facility, const std::string froads,
		const std::string &netid, const std::string &source_name, const std::string &target_name, const std::string &direction_name,
		const std::string ffacilities, const std::string &fid, const std::string &flabel, const bool &isunifiedcutoff, const double &unifiedcutoff, const std::string &cutoff_name,
		const double &searchradius, bool writeaccedges, bool writecatchment);


	// write the contour as shp file of polygons
	void write_contour_polygons(const std::string &filename);

	// write the accessible edges as shp file of linestring
	void write_accessible_edges(const std::string &filename);


	//used for sorting the linestrings of a contour by the number of points of each linestring
	static bool compareBylength(const OGRLinearRing &line_a, const OGRLinearRing &line_b)
	{
		return (line_a.get_Length() >= line_b.get_Length());
	};

	//-------------------for generating the benchmark----------------------------------------------
	void OutputExtendedSPTrees(std::string froads, std::string ffacilities, std::string basic_name);

public:
	Network  network_;
	vector<Contour>  faci_contours_;
	vector<vector<int>> faci_contour_labels_;
	vector<vector<AccessibleEdge>> faci_accedges_; // the accessible edges of all the facilities
};

}// end namespace
#endif