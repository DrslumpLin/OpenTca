// 1. read the road network and facilities
//    including the undirected and directed road network and multiple-points facilities
// 2. the functions for generating the benchmark grid points are included 
// @author: Diao Lin
// @version: 2020.03

#ifndef TCA_NETWORK_H
#define TCA_NETWORK_H

// Data structures for Rtree


#include "ogrsf_frmts.h" // C++ API for GDAL
#include <iostream>
#include <math.h>       // Calulating probability
#include <iomanip>
#include <algorithm>    // Partial sort copy

#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/function_output_iterator.hpp>
#include <unordered_map>

#include "tca_types.hpp"
#include "tca_algorithm.hpp"

namespace bg = boost::geometry;
namespace bgi = bg::index;
using namespace tca::types;

namespace tca
{
// Define a type alias for the rtree for searching the nearest edges and sub-edges
typedef boost::geometry::model::point<float, 2, boost::geometry::cs::cartesian> boost_point; // Point for rtree box
typedef boost::geometry::model::box<boost_point> boost_box;                                  // box for rtree box
typedef std::pair<boost_box, Edge*> Item;											         // Item stored in rtree
typedef boost::geometry::index::rtree<Item, boost::geometry::index::quadratic<16>> Rtree;    // Rtree definition

//----for directed----
typedef std::pair<boost_box, EdgeDirected*> ItemDirected;													 // Item stored in rtree
typedef boost::geometry::index::rtree<ItemDirected, boost::geometry::index::quadratic<16>> RtreeDirected;    // Rtree definition


class Network
{

private:
	int srid_;															// Spatial reference id
	Rtree rtree_;														// Network rtree structure
	std::unordered_map<int, OGRLineString*> edge_linestring_map_;		// map between the inner id and edge linestring
	std::vector<Facility> facilities_;									// all the facilities
	
	std::vector<Edge> network_edges_;									// all edges in the network
	std::vector<Edge*> network_edge_pointers_;						    // all edge pointers in the network

	std::vector<double> deltas_;										// the cut-off distances of TCA generation	
	std::vector<std::string> labels_;									// for multiple-point facility, a label speficy the facilty a point belongs to

	//-----for directed----- 
	RtreeDirected rtree_directed_;									    // for directed road network
	std::vector<EdgeDirected> network_directed_edges_;				    // for directed road network

public:
	Network();

	~Network();	

	// read undirected road network
	void ReadRoadNetworks(const std::string filename, const std::string &id_name = "id", const std::string &source_name = "source",
		const std::string &target_name = "target");
	
	// read facilities
	void ReadFacilities(const std::string filename, std::string id_name = "id", std::string delta_name = "delta");
	
	// two combinations: <isunifiedcutoff = true, unifiedcutoff>; <isunifiedcutoff = false, cutoff_name>
	void ReadFacilities(const std::string &filename, const std::string &id_name, const bool &isunifiedcutoff,
		const double &unifiedcutoff, const std::string &cutoff_name);

	// Get the the external egde ID of an edge according to its index
	std::string GetEdgeExternalId(int internal_id);
		
	// return all the values of delta
	std::vector<double> get_deltas(void);

	// return all the labels for input points, only for multiple-point facility 
	std::vector<std::string> get_labels(void);

	// Construct a Rtree using the vector of edges
	void BuildRtreeIndex(void);
	
	//searching all the nearestpoints for input facilities by using Rtree
	NearestPoints SearchNearestPoints(double radius);
	
	// searching the nearest points under the condition of no r-tree is built
	// this is used for comparison purposes
	NearestPoints SearchNearestPointsNoRtree(double radius);

	// Extracting all the sub networks of each individual facilty 
	// Each sub-network is represented by vector<Edge*>,
	// input: 
	//	 @delta: a unified delta for all the facilities, 
	//   @facility_nearest_points: the nearest points for the entire facilities
	std::vector<std::vector<Edge*>> GetEdgesAroundDelta(double delta, NearestPoints facility_nearest_points);
	
	// Extracting all the sub networks of each individual facilty
	// in this version, each faclity has its own delta instead of a unified delta for all the facilities
	// input: @facility_nearest_points: the nearest points for the entire facilities
	std::vector<std::vector<Edge*>> GetEdgesAroundDelta(NearestPoints facility_nearest_points);

	// Get all the edge pointers of the road network
	// used for multiple-source voronoi polygons
	std::vector<Edge*> GetAllEdges();  

	// write the subegdes of all the facilities
	// for testing uses
	void write_sub_edges(std::string filename, std::vector<Edge*> sub_edges);		
	
	// comparing the distance from the projected point to the facility point of two nearest points 
	// used for sorting all the near points within the searching buffer of a facility
	static bool compareBydis(const NearestPoint &a, const NearestPoint &b)
	{
		if (a.p_2_pp != b.p_2_pp)
		{
			return a.p_2_pp < b.p_2_pp;
		}
		else
		{
			return a.edge->id < b.edge->id;
		}
	};
	
	//-------------------the functions below are used for multiple-points facility---------------------------------------------

	// read facilities that represented by multiple points, 
	// a label filed is needed to specify which facility a point belongs to
	void ReadMpFacilities(const std::string filename, std::string id_name = "id", std::string delta_name = "delta", std::string label_name = "label");

	// read facilities that represented by multiple points, 
	// a label filed is needed to specify which facility a point belongs to
	void ReadMpFacilities(const std::string filename, const std::string &id_name, std::string label_name, const bool &isunifiedcutoff,
		const double &unifiedcutoff, const std::string &cutoff_name );


	// Extracting the sub edges of a multiple-points(MP) facility
	// Input: 
	// 	@facility_nearest_points: all the points that belongs to one single facilty
	//  @delta: all the sub points of a facility have the same delta
	std::vector<Edge*> GetEdgesAroundDeltaSingleMpFacility(NearestPoints facility_nearest_points, double delta);
	// Get the sub-edges of all the multiple-points(MP) facilities (for directed roads)
	std::vector<EdgeDirected*> GetEdgesAroundDeltaSingleMpFacilityDirect(NearestPointsDirected facility_nearest_points, double delta);



	// Get the sub-edges of all the multiple-points(MP) facilities
	// Input:
	//	 @ unique_facilitis: a vector of all the MP facilities
	//   @ unique_deltas: a vector of all the deltas for the MP facilities
	std::vector<std::vector<Edge*>> GetEdgesAroundDeltaMltipleMpFacilities(std::vector<NearestPoints> unique_facilitis, std::vector<double> unique_deltas);
	
	// Get the sub-edges of all the multiple-points(MP) facilities (for directed roads)
	std::vector<std::vector<EdgeDirected*>> GetEdgesAroundDeltaMltipleMpFacilitiesDirected(std::vector<NearestPointsDirected> unique_facilitis, std::vector<double> unique_deltas);

	//-----------end for multiple points-based facility related functions --


	//-------------------the functions below are used for directed edge input---------------------------------------------
	// read road network with direction information, the defaulted version is calculated for from facilty
	void ReadDirectedRoadnetworks(std::string filename, const std::string &id_name = "id", const std::string &source_name = "source",
		const std::string &target_name = "target", const std::string &direction_name = "direction");
	
	// read road network for generating to-facility catchment area
	// the funciton reverse the direction of the input road network
	void ReadDirectedRoadnetworksToFacility(std::string filename, const std::string &id_name = "id", const std::string &source_name = "source",
		const std::string &target_name = "target", const std::string &direction_name = "direction");

	// searching the nearest points for a facility (directed network version)
	NearestPointsDirected SearchNearestPointsDirected(double radius);

	// Construct a Rtree using the vector of edges
	void BuildRtreeIndexDirected(void);

	//
	std::vector<EdgeDirected> get_edges_directed(void);
	
	//
	std::vector<std::vector<EdgeDirected*>> GetEdgesAroundDeltaDirected(double delta, NearestPointsDirected facility_nearest_points);

	// read delta from the facilty version 
	std::vector<std::vector<EdgeDirected*>> GetEdgesAroundDeltaDirected(NearestPointsDirected facility_nearest_points);

	//
	static bool compareBydisdirected(const NearestPointDirected &a, const NearestPointDirected &b)
	{
		if (a.p_2_pp != b.p_2_pp)
		{
			return a.p_2_pp < b.p_2_pp;
		}
		else
		{
			return a.edge->id < b.edge->id;
		}
	};
	//---end for directed network input related functions --


	//----------------------below are used for generating the bechmark for 10*10m grid for our article-----------------
	// a integrated function to read grids, calculate distances from each grid to the root node
	//	write out the grids with an additional filed of distances to the root node
	void BenchmarkGeneration(std::string fsub_graph, std::string fgrid_points, std::string fup_grid_points, double delta, OGRPolygon* arc_gis_CA, OGRPolygon* our_CA);

	//read the subgraph of each facilities
	std::vector<SubGraphEdge> ReadSubGraph(const std::string filename, const std::string &id_name = "id", const std::string &source_name = "source",const std::string &target_name = "target",
					  const std::string &s_distance = "s_distance", const std::string &e_distance = "e_distance");
	
	// read the grid points
	std::vector<OGRPoint>  ReadGridPoints(const std::string filename, const std::string &id_name = "id");

	// calculate the distance from each grid point to the root node
	std::vector<double> DistanceToRootNode(std::vector<SubGraphEdge> sub_edges, std::vector<OGRPoint> grid_points);

	// write out the grid points (a field of distance is added)
	void write_accessible_grid_points(std::string filename, std::vector<OGRPoint> grid_points, 
		std::vector<double> distances, const double delta, OGRPolygon* arc_gis_CA, OGRPolygon* our_CA);

	// read the CA polygons generated by ArcGIS and Ours
	std::vector<OGRPolygon*>  ReadCatchmentAreas(const std::string filename);
	//--end for benchmark generationg-related functions-----
};// end for class TCA
}// namespace tca
#endif /* TCA_NETWORK*/