// 1. build undirect graph based on the boost graph, shortest path and extended shortest path tree
// 2. generating the input for triangulation and contrianed triangulation
// 3. network voronoi construction
// @author: Diao Lin
// @version: 2020.03

#ifndef TCA_UNDIRRECTED_GRAPH_H
#define TCA_UNDIRRECTED_GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <deque>
#include <algorithm>
#include <unordered_map>

//#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include "ogrsf_frmts.h"
#include "gdal_priv.h"

#include "tca_types.hpp"
#include "tca_network.h"
#include "float.h"

// for voronoi diagram generation
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>

using namespace std;

namespace tca 
{
	// for the graph edge
	struct EdgeProperty  
	{
		std::string id;  // the external id of an edge
		double length;   // the length of an edge
	};

	//definition of undirected graph
	typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeProperty> Graph_T;
	
	// descriptors
	typedef Graph_T::vertex_descriptor vertex_descriptor;
	typedef Graph_T::edge_descriptor edge_descriptor;
	
	// iterators
	typedef boost::graph_traits<Graph_T>::vertex_iterator vertex_iterator;
	typedef boost::graph_traits<Graph_T>::out_edge_iterator out_edge_iterator;
	typedef boost::graph_traits<Graph_T>::edge_iterator edge_iterator;

	//used for driving distances
	struct found_goals {}; 

	// for storing the detected extended shortest path nodes
	struct ExtendedNode							
	{
		OGRLineString*  pLine;					//the linestring that not in the sp tree
		double break_distance_to_snode;			//the break distance to the start node of this line
		vertex_descriptor source;				// staring node of the edge
		vertex_descriptor target;				// target node of the edge
		int node_id;						    // a local id for the extended node
	};


	// definition of the dijkstra visitor for cacthment area generation 
	class CA_vistor : public boost::default_dijkstra_visitor
	{
		public:
			// Create a visitor 
			explicit CA_vistor
			(
				std::vector< vertex_descriptor > &examimed_nodes,
				std::vector<double> &distances
			) : m_nodes(examimed_nodes), m_dist(distances) {};

			template <class Graph>void examine_vertex(vertex_descriptor u, Graph &g)
			{
				m_nodes.push_back(u);
			};

		private:
			std::vector<vertex_descriptor> &m_nodes; //Precedessors
			std::vector<double> &m_dist; // Distances
	}; // CA_vistor


	// definition of the Mutiple source dijkstra visitor
	// updating the closest_facilities by using edge_relaxed
	// by using closest_facilities, each node is labeled to have a facility label
	// which representing the distance from the node to the facility is shortest among all the facilities
	// based on the discussion: 
	// http://boost.2283326.n4.nabble.com/Graph-How-to-do-multi-source-shortest-path-td4450458.html
	class MS_distance_visitor : public boost::default_dijkstra_visitor 
	{
		public:
			// Create a visitor 
			explicit MS_distance_visitor
			(
				vertex_descriptor sourcenode,
				std::vector<vertex_descriptor> &nodesInDistance,
				std::vector<vertex_descriptor> &closest_facilities
			) : m_sourcenode(sourcenode), m_nodes(nodesInDistance), m_closest_facilities(closest_facilities) {};

			template <class Graph>void examine_vertex(vertex_descriptor u, Graph &g)
			{
				m_nodes.push_back(u);
			};

			//except for the facilities nodes, the other nodes 
			template <class Graph>void edge_relaxed(edge_descriptor e, Graph &g)
			{
				vertex_descriptor  vt_s_node = boost::source(e, g);
				vertex_descriptor  vt_e_node = boost::target(e, g);
				if (m_sourcenode != vt_s_node) 
				{
					m_closest_facilities[vt_e_node] = m_closest_facilities[vt_s_node];
				}
			};
		private:
			vertex_descriptor m_sourcenode; //Delta
			std::vector<vertex_descriptor> &m_nodes; //Precedessors
			std::vector<vertex_descriptor> &m_closest_facilities; //the mark of nearest 
	}; // MS_distance_visitor


// for voronoi diagram this is for std::map index to get the inner node of virinoi diagram 


class UndirectedGraph
{

// data structures for voronoi diagram
typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel                   K;
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP>                                  VD;

typedef AT::Point_2                     Point_2;
typedef VD::Face_handle                 Face_handle;
typedef VD::Ccb_halfedge_circulator     Ccb_halfedge_circulator;

typedef CGAL::Polygon_2<K> Polygon_2;
typedef Polygon_2::Vertex_iterator polygon_vi;
typedef Polygon_2::Edge_const_iterator polygon_ei;

typedef K::Segment_2 Segment_2;
typedef K::Intersect_2 Intersect_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;


//A class to recover Voronoi diagram from stream.
//Rays, lines and segments are cropped to a rectangle
//so that only segments are stored
struct Cropped_voronoi_from_delaunay 
{
	std::list<Segment_2> m_cropped_vd;
	Iso_rectangle_2 m_bbox;

	Cropped_voronoi_from_delaunay(const Iso_rectangle_2& bbox) :m_bbox(bbox) {}

	template <class RSL>
	void crop_and_extract_segment(const RSL& rsl) 
	{
		CGAL::Object obj = CGAL::intersection(rsl, m_bbox);
		const Segment_2* s = CGAL::object_cast<Segment_2>(&obj);
		if (s) m_cropped_vd.push_back(*s);
	}

	void operator<<(const Ray_2& ray) { crop_and_extract_segment(ray); }
	void operator<<(const Line_2& line) { crop_and_extract_segment(line); }
	void operator<<(const Segment_2& seg) { crop_and_extract_segment(seg); }
};



// a point index for the convience of point operation
struct PointIndex
{
	Point_2 point;
	//PointIndex();
	PointIndex(const Point_2 & _inpoint) : point(_inpoint) {};

	bool operator<(const PointIndex& other) const
	{
		if (point.x() != other.point.x())
			return (point.x() < other.point.x());
		return (point.y() < other.point.y());
	}
};



public:
	UndirectedGraph(std::vector<Edge*> egdes);
	~UndirectedGraph();

public:
	Graph_T							udirected_g_;
	std::vector<double>				distances_map_;
	std::vector<vertex_descriptor>  predecessors_map_;
	std::vector<vertex_descriptor>  closest_facilities_map_;
	//std::vector<vertex_descriptor>  examined_vertices_;				

	std::unordered_map<std::string, OGRLineString*> edges_externalid_linestring;
	std::unordered_map<std::string, edge_descriptor> edges_externalid_discriptor;       // []	

	std::unordered_map<std::string, int> vertexs_externalid_discriptor_;
	std::vector<std::string> vertex_externalids_;									    //{inner_id,string_id}
	std::vector<OGRPoint>    vertex_coords_;										    //{inner_id,node_coords}	
	// vertex string id and vertex inner id (this inner id also equal to vertex descriptor)

	//for test
	//std::vector<OGRLineString*>   Added_linestrings_;
	
public:
	
	//----------------------------for generating the input of catchment area--------------------------------
	// insearting a nearest point to the graph
	void InsertingFacilityNode(const NearestPoint ne_point);

	// get the shortest path nodes sequence of a facility
	// return the examined nodes
	std::vector<vertex_descriptor> GetShortestPathNodes(const std::string facility_node_sid);
	
	// build the extened shortest path tree of a facility
	// input:
	//   @examined_nodes: shortest path nodes
	//   @facility_node_sid: the facility name
	std::vector<vertex_descriptor> BuildExtendedShortestPathTree(std::vector<vertex_descriptor> examined_nodes, const std::string facility_node_sid);
	
	// Get the points with a distance equal to delta to the root node	
	std::vector<OGRPoint> FindingDeltaNodes(double delta, std::vector<vertex_descriptor> ex_nodes);	
	
	// Get all the accessibile edges within a cut-off distance for a facility
	// input:
	//	  @delta: the cut-off distance
	//    @facility_id:
	// 	  @ex_nodes: the examined nodes of the extended shortest path tree for the facility
	std::vector<AccessibleEdge> GetAccessibleEdges(double delta, int facility_id, std::vector<vertex_descriptor> ex_nodes);	


	// Segmenting the extended shortest path tree edges to get the edges for constructing the constrained triangulation
	// input:
	//  @extended_sp_exaimed_nodes: the exaimed nodes of extended path tree
	// output:
	//  @edge_node_ids: the node id (local ids corresponding to the following three vectors) of all the contrianed edges
	//  @node_ids: the string-format node ids
	//  @node_distances:the distances from root node to the constrianed nodes
	//  @node_org_points: the coordinates of the contrianed nodes
	void GetEdgesForConstrianedTriangulationLineSegmentation(std::vector<vertex_descriptor> extended_sp_exaimed_nodes,
		std::vector<std::array<int, 2>>	& edge_node_ids, std::vector<string>& node_ids,
		std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points);
		
	// Inserting the dalta nodes and recalcuating the shortest path tree
	// no segmentation process is conduceted in this function
	// the parameter is the same as function: GetEdgesForConstrianedTriangulationLineSegmentation
	// except the delta is useful in this function
	void GetEdgesForConstrianedTriangulation(const std::string facility_node_sid, const double delta, 
		std::vector<std::array<int, 2>>& edge_node_ids,std::vector<string>& node_ids, 
		std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points);

	// Get the nodes for constructing triangulation, including the extended shortest path tree nodes
	// and nodes at the distance of delta
	// input:
	//  @facility_node_sid: the node id
	//  @delta: the cut-off distance
	// output:
	//  @node_ids: the string-format node ids
	//  @node_distances:the distances from root node to the constrianed nodes
	//  @node_org_points: the coordinates of the contrianed nodes
	void GetNodesForTrigulationInput(const std::string facility_node_sid, const double delta,
		std::vector<string>& node_ids, std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points);

	//--------------------------------IO functions----------------------------------------------
	void write_accessible_edges(std::string filename, double delta, int facility_id, std::vector<vertex_descriptor> ex_nodes);
	
	void wirte_sp_nodes(const std::string &filename, const std::vector<vertex_descriptor> examined_nodes);
	
	void wirte_sp_nodes_plus_delatnodes(const std::string &filename, const std::vector<vertex_descriptor> examined_nodes,
		const std::vector<OGRPoint> deltanodes, double delta);
	
	void write_sp_edges(const std::string &filename, const std::vector<vertex_descriptor> examined_nodes);
	
	//--------------------------------------basic functions---------------------------------------------------
	//
	OGRLineString GetEdgeLinestring(vertex_descriptor source, vertex_descriptor target);
	
	//
	double CalculateDistance(OGRPoint point1, OGRPoint point2);
	
	//input: point1(x1, y1) and point2(x2, y2)
	//output: (vertical_x,vertical_y)
	void GetCounterClockwiseUnitNormalVector(double x1, double y1, double x2, double y2, double & vertical_x, double & vertical_y);

	// Initialize the distance map and predecessor map
	void InitializeDistancesPredecessors();
	
	// Initialize the closest_facilities_map_ as the node itself
	void InitializeClosestFacilities(int vertexnumber);
	
	// get the shortest path node sequences of a source node
	// input:  @source
	// output: @examined_nodes
	void ShortestPathAlgorithm(const vertex_descriptor& source, std::vector<vertex_descriptor>& examined_nodes);
	
	// the multiple sources shortest path algorithm
	// input: 
	//	@source: the root node, is a virtual node
	// output: 
	//	@p_examined_nodes: the sequence of examined nodes
	//	@closest_facilities: the closest facility labels (nodes have a distance to virtual node euqal to 0) to each examined node
	void MultipleSourceShortestPathAlgorithm(const vertex_descriptor& source, 
		std::vector<vertex_descriptor>& p_examined_nodes,
		std::vector<vertex_descriptor>& closest_facilities);

	// Identifying the non-shortest path edges and the corresponding break points
	// input: @examined_nodes: the examined points of the SP tree
	// return: a vector of all the extended nodes
	std::vector<ExtendedNode> FindNonShortestPathEdges(std::vector<vertex_descriptor> examined_nodes);
	
	// Inserting the extended nodes into the graph
	void InsertingNonShortestPathEdges(const std::vector<ExtendedNode> extended_records);
	
	// inserting a extend node to a edge of the graph
	// input:
	//    @source @target: the two nodes of the edge
	//    @extend_node_id: the local id of extended node
	//	  @breakpoint
	//    @s_cutedge: the linestring from snode to break point
	//	  @e_cutedge: the linestring from break point to end node
	void InsertingExtendedNodes(vertex_descriptor source, vertex_descriptor target, int extend_node_id,
		OGRPoint breakpoint, OGRLineString* s_cutedge, OGRLineString* e_cutedge);	
	
	// Find edges with delta distances and the delta nodes
	std::vector<ExtendedNode> FindDeltaEdges(double delta, std::vector<vertex_descriptor> examined_nodes);			
	
	// Inserting the delta nodes and corresponding edges
	void InsertDeltaEdgesToGraph(const std::vector<ExtendedNode> extended_records);		
	
	// Inserting a sigle node into an edge
	void InsertingDeltaNodes(vertex_descriptor source, vertex_descriptor target, int extend_node_id,
		OGRPoint breakpoint, OGRLineString* s_cutedge, OGRLineString* e_cutedge);	

	// --------------the follwoing function is used in the article for comparison purpose------------------
	// Get network nodes as the input for triangulation
	// for both soultion 2(SP nodes) and solution 3 (Extended SP nodes)
	// this function is quite similar to "GetNodesForTrigulationInput" but, 
	// we donot explictly include delta node in this version
	// Input: @v_exmained_nodes
	// Output: @node_ids; @node_distances; @node_org_points
	void GetInputForDelaunary(std::vector<vertex_descriptor> v_exmained_nodes, vector<string>& node_ids, 
		vector<double>& node_distances, vector<OGRPoint>& node_org_points);



	//----------------------the following function is used for Mupltiple source Network Voronoi---------------------------
	
	// inserting facilities nodes, vritual nodes, and links 
	void InsertingMulltipleFacilityNodes(const NearestPoints ne_points);

	// Calculating shortest path from a vitual node
	void MultipleSourceShortestPaths(const std::string insert_node_sid);

	// the network voronoi here means each edge is assigned to its neartest facilty
	// some edges may be split into two different parts
	std::vector<NetworkVoronoiEdge> ConstructNetworkVoronoiDiagram(void);

	// the counterpart for function "ConstructNetworkVoronoiDiagram"
	std::vector<NetworkVoronoiNode> ConstructNetworkVoronoiNode(void);

	//---------------below are writing functions-----------------------
	void write_network_voronois_edges(const std::string &filename, std::vector<NetworkVoronoiEdge> vd_edges);

	void write_multiple_sources_node_distances(std::string filename);

	void write_nv_nodes(std::string filename, std::vector<NetworkVoronoiNode> nv_nodes);


	//-----------------------below are functions for cropped voronoi polygons------------------------------
	
	// build voronio diagram for the voronoi nodes and use the bounding box of these nodes to crop the voronoi diagram
	// output the voronoi diagram as segments
	// input:  @nv_nodes:the orignal voronoi nodes
	// output: @linstrings: linestrings
	// note, there is no crop at the breaking points and no labels are ouput
	void BuildCroppedVoronoiPolySegs(const std::vector<NetworkVoronoiNode> nv_nodes, vector<OGRLineString>& linstrings);

	// write the voronoi segments cropped by bounding box
	void write_multiple_sources_polygons_Segements(std::string filename, vector<OGRLineString> linstrings);

	// build voronoi polygons of the voronoi nodes including the following steps
	// 1. gengerating addtional points around the outside of the bbox of the original network voronoi points and add them as the input node for 
	//    generating voronoi diagram, the aims is try to make the original open boundary voronoi polygons into closed polygons 
	// 2. crop all the voronoi polygons by the bounding box and only keep the voronoi polygons for the orignal network voronoi points 
	// 3. if a voronoi polygon is belongs to a break point, then it needs to be split to two polygons,each sub polygon correspond
	//    to a facility 
	// The results of this program can be further processed by using dissolving in GIS software to dissoving polygons with the same label
	// Thus to obtain network voronoi based polygons 
	// input:
	//    @nv_nodes:the orignal voronoi nodes
	// output:
	//    @ploygons
	//    @labels
	void BuildCroppedVoronoiPolygon(const std::vector<NetworkVoronoiNode> nv_nodes, vector<OGRPolygon>& ploygons, vector<string>& labels);

	//write the cropped polygons with labels
	void write_multiple_sources_polygons(std::string filename, vector<OGRPolygon> ploygons, vector<string> labels);

	//  tranfer a CGAL polygon to OGRpolygon 
	void CGALToOGRPolygon(Polygon_2 bpoly, OGRPolygon & ogr_polygon);
	// transfer OGR polygon to a CGAL polygon
	void OgrPolygonToCgalPolygon(OGRPolygon ogr_polygon, Polygon_2&  cgal_poly);

	// gengerating addtional points around the outside of the bbox of the original network voronoi points
	// the aim of generating these points is to make the original open boundary voronoi polygons into closed polygons
	// input:
	//   @bbox: the bounding box of the orignal voronoi nodes
	//   @levels: defining how many addtional points are needed (|level| * 8)
	// output
	//	 @add_points: the addtional points generated
	void GenerateAddtionalNvPoints(Iso_rectangle_2 bbox, vector<double> levels, vector<Point_2>& add_points);

	// given a line (represented by a point and its slope), and a voronoi polygon (countourclockwise directions)
	// cut the voronoi into two different polygons and assign each polygon a facility label
	// the basic priciple for polygon split include:
	//	1.find the normal vector of the edge with break points
	//  2.build a line with direction as the normal vector and pass the break point
	//  3.the orginal polygon is then split into two sub polygons by the line by using the following steps
	//     1) the voronoi polygon are arranged in counterclockwise direction by default
	//     2) CCW_point:the intersected point between the counterclockwise normal vector and the voronoi polygon
	//     3) CW_point:the intersected point between the clockwise normal vector and voronoi polygon
	//     4) constructing the polygon belongs to S node(facility): by starting from CCW_point via counterclockwise direction to CW_point
	//     5) constructing the polygon belongs to T node(facility): by starting from CW_point via counterclockwise direction to CCW_point
	void SplitVoronoiPoygon(Polygon_2 poly, NetworkVoronoiNode node, OGRPolygon & ogr_polygon_1, OGRPolygon & ogr_polygon_2);	
	
};
}// end of namespace tca

#endif