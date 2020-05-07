// 1. build direct graph, shortest path and extended shortest path tree
// 2. generating the input for triangulation and contrianed triangulation
// 3. network voronoi construction 
// @author: Diao Lin
// @version: 2020.03

#pragma once
#ifndef TCA_DIRRECTED_GRAPH_H
#define TCA_DIRRECTED_GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <deque>
#include <algorithm>
#include <unordered_map>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/property_map/property_map.hpp>

#include "ogrsf_frmts.h"
#include"gdal_priv.h"

#include "tca_types.hpp"
#include "tca_network.h"
#include "float.h"




using namespace std;
namespace tca 
{
	struct EdgePropertyDrct  // for the graph edge
	{
		std::string id;  // the external id of an edge
		int direction;   // 1(for bi-edge, S-T); -1(for bi-edge, T-S); 2(for single edge, S-T); 3 (for single edge T-S);
		double length;   // the length of an edge
	};

//	struct found_goals {}; //used for driving distances

	typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgePropertyDrct> Graph_T_Directed;
	
	typedef Graph_T_Directed::vertex_descriptor vertex_descriptor_drct;
	typedef Graph_T_Directed::edge_descriptor edge_descriptor_drct;
	
	typedef boost::graph_traits<Graph_T_Directed>::vertex_iterator vertex_iterator_drct;
	typedef boost::graph_traits<Graph_T_Directed>::out_edge_iterator out_edge_iterator_drct;
	typedef boost::graph_traits<Graph_T_Directed>::in_edge_iterator in_edge_iterator;    // specialized for directed graph
	typedef boost::graph_traits<Graph_T_Directed>::edge_iterator edge_iterator_drct;

	// for storing the detected extend shortest path nodes
	struct ExtendedNodeDrct							
	{
		OGRLineString*  pLine;					//the linestring that not in the sp tree
		double break_distance_to_snode;			//the break distance to the start node of this line
		vertex_descriptor_drct source;
		vertex_descriptor_drct target;
		int direction;							// 
		int node_id;						    // start from 0
	};

	class CA_vistorDrct : public boost::default_dijkstra_visitor
	{
	public:
		// Create a visitor 
		explicit CA_vistorDrct
		(
			std::vector< vertex_descriptor_drct > &examimed_nodes,
			std::vector<double> &distances
		) : m_nodes(examimed_nodes), m_dist(distances) {};

		template <class Graph>void examine_vertex(vertex_descriptor_drct u, Graph &g)
		{
			m_nodes.push_back(u);
		};

	private:
		std::vector<vertex_descriptor_drct> &m_nodes; //Precedessors
		std::vector<double> &m_dist; // Distances
	}; // CA_vistorDrct

	// definition of the Mutiple source dijkstra visitor
	// updating the closest_facilities by using edge_relaxed
	// by using closest_facilities, each node is labeled to have a facility label
	// which representing the distance from the node to the facility is shortest among all the facilities
	// based on the discussion: 
	// http://boost.2283326.n4.nabble.com/Graph-How-to-do-multi-source-shortest-path-td4450458.html
	class MS_distance_visitor_drct : public boost::default_dijkstra_visitor
	{
	public:
		// Create a visitor 
		explicit MS_distance_visitor_drct(
			vertex_descriptor_drct sourcenode,
			std::vector<vertex_descriptor_drct> &nodesInDistance,
			std::vector<vertex_descriptor_drct> &closest_facilities
		) : m_sourcenode(sourcenode), m_nodes(nodesInDistance), m_closest_facilities(closest_facilities) {};

		template <class Graph>void examine_vertex(vertex_descriptor_drct u, Graph &g)
		{
			m_nodes.push_back(u);
		};

		//except for the facilities nodes, the other nodes 
		template <class Graph>void edge_relaxed(edge_descriptor_drct e, Graph &g)
		{
			vertex_descriptor_drct  vt_s_node = boost::source(e, g);
			vertex_descriptor_drct  vt_e_node = boost::target(e, g);
			if (m_sourcenode != vt_s_node) {
				m_closest_facilities[vt_e_node] = m_closest_facilities[vt_s_node];
			}
		};
	private:
		vertex_descriptor_drct m_sourcenode; //Delta
		std::vector<vertex_descriptor_drct> &m_nodes; //Precedessors
		std::vector<vertex_descriptor_drct> & m_closest_facilities; //the mark of nearest 
	}; // MS_distance_visitor

class DirectedGraph
{

public:

	DirectedGraph(std::vector<EdgeDirected*> egdes);
	~DirectedGraph();

	Graph_T_Directed				directed_g_;
	std::vector<double>				distances_map_;
	std::vector<vertex_descriptor_drct>  predecessors_map_;
	std::vector<vertex_descriptor_drct>  closest_facilities_map_;
	//std::vector<vertex_descriptor_drct>  examined_vertices_;				

	std::unordered_map<std::string, OGRLineString*> edges_externalid_linestring;
	std::unordered_map<std::string, edge_descriptor_drct> edges_externalid_discriptor;       // []	

	std::unordered_map<std::string, int> vertexs_externalid_discriptor_;
	std::vector<std::string> vertex_externalids_;									    //{inner_id,string_id}
	std::vector<OGRPoint>    vertex_coords_;										    //{inner_id,node_coords}	
	// vertex string id and vertex inner id (this inner id also equal to vertex descriptor)
	

public:

	//----------for generating the input of catchment area-------------
	
	// insearting a nearest point to the graph
	void InsertingFacilityNode(const NearestPointDirected ne_point);

	// get the shortest path nodes sequence of a facility
	// return the examined nodes
	std::vector<vertex_descriptor_drct> GetShortestPathNodes(const std::string facility_node_sid);
	
	// build the extened shortest path tree of a facility
	// input:
	//   @examined_nodes: shortest path nodes
	//   @facility_node_sid: the facility name
	std::vector<vertex_descriptor_drct> BuildExtendedShortestPathTree(std::vector<vertex_descriptor_drct> examined_nodes, const std::string facility_node_sid);
	
	// Get the points with a distance equal to delta to the root node	
	std::vector<OGRPoint> FindingDeltaNodes(double delta);	
	
	// Get all the accessibile edges within a cut-off distance for a facility
	// input:
	//	  @delta: the cut-off distance
	//    @facility_id:
	// 	  @ex_nodes: the examined nodes of the facility
	std::vector<AccessibleEdge> GetAccessibleEdges(double delta, int facility_id, std::vector<vertex_descriptor_drct> ex_nodes);
	

	//// A intergation of shortest path, extend sp, and find delta nodes
	//void GetNodesForTrigulationInput(const std::string facility_node_sid, const double delta,std::vector<string>& node_ids,
	//								std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points);
	

	// Segmenting the extended shortest path tree edges to get the edges for constructing the constrained triangulation
	// input:
	//  @facility_node_sid: the node id
	//  @delta: should be deleted later (no delta is inserted in this version)
	// output:
	//  @edge_node_ids: the node id (local ids corresponding to the following three vectors) of all the contrianed edges
	//  @node_ids: the string-format node ids
	//  @node_distances:the distances from root node to the constrianed nodes
	//  @node_org_points: the coordinates of the contrianed nodes
	void GetEdgesForConstrianedTriangulationLineSegmentation(const double &cutoff, std::vector<vertex_descriptor_drct> &extended_sp_exaimed_nodes,
		std::vector<std::array<int, 2>>	& edge_node_ids, std::vector<string>& node_ids,
		std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points);

	//
	// Split a line to segments and add the segments into the constrianed edges
	//input:
	// @s_dis: the distance to from root node to the starting node on edge
	// @source@target: node ids of the edge
	// @temp_line: the linstring of the edge
	//output:
	//  @edge_node_ids: constrianed edges
	//  @node_ids: the vector of node ids
	//  @node_distances: the distance from root to the node
	//  @node_org_points: the coordinate of a node
	void  SplitLineToSegments(double s_dis, int source, int target, OGRLineString* temp_line, std::vector<std::array<int, 2>>	& edge_node_ids,
		std::vector<string>& node_ids, std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points);

	//--------------------------IO functions--------------------------------------
	void write_accessible_edges(std::string filename, double delta, int facility_id, std::vector<vertex_descriptor_drct> ex_nodes);
	
	void wirte_sp_nodes(const std::string &filename, const std::vector<vertex_descriptor_drct> examined_nodes);
	
	void write_sp_edges(const std::string &filename, const std::vector<vertex_descriptor_drct> examined_nodes);
	
	void wirte_sp_nodes_plus_delatnodes(const std::string &filename, const std::vector<vertex_descriptor_drct> examined_nodes,
										const std::vector<OGRPoint> deltanodes, double delta);

	//write the edges of directed graph as shp file (for testing purposes)
	void write_for_test(std::string filename);


	//-------------------------for Mupltiple source Network Voronoi------------------------------------
	void InsertingMulltipleFacilityNodes(const NearestPointsDirected ne_points);

	void MultipleSourceShortestPaths(const std::string insert_node_sid);

	void MultipleSourceShortestPathAlgorithm(const vertex_descriptor_drct & source, std::vector<vertex_descriptor_drct>& p_examined_nodes,
		std::vector<vertex_descriptor_drct>& closest_facilities);

	std::vector<NetworkVoronoiEdge> ConstructNetworkVoronoiDiagram(void);

	void write_network_voronois_edges(const std::string &filename, std::vector<NetworkVoronoiEdge> vd_edges);

	void write_multiple_sources_node_distances(std::string filename);

private:

	OGRLineString GetEdgeLinestring(vertex_descriptor_drct source, vertex_descriptor_drct target);
	
	// Initialize the distance map and predecessor map
	void InitializeDistancesPredecessors();
	
	// Initialize the closest_facilities_map_ as the node itself
	void InitializeClosestFacilities(int vertexnumber);
	
	// get the shortest path node sequences of a source node
	// input:  @source
	// output: @examined_nodes
	void ShortestPathAlgorithm(const vertex_descriptor_drct& source, std::vector<vertex_descriptor_drct>& examined_nodes);
	
	// Identifying the non-shortest path edges and the corresponding break points
	// input: @examined_nodes: the examined points of the SP tree
	// return: a vector of all the extended nodes
	std::vector<ExtendedNodeDrct> FindNonShortestPathEdges(std::vector<vertex_descriptor_drct> examined_nodes);
	
	// Inserting the extended nodes into the graph
	void InsertingNonShortestPathEdges(const std::vector<ExtendedNodeDrct> extended_records);
	
	// inserting a extend node to a edge of the graph
	// input:
	//    @source @target: the two nodes of the edge
	//    @extend_node_id: the local id of extended node
	//    @direction: the direction of the edge
	//	  @breakpoint: 1 (for bi-edge, S-T); -1(for bi-edge, T-S); 2(for single edge, S-T); 3 (for single edge T-S);
	//    @s_cutedge: the linestring from snode to break point
	//	  @e_cutedge: the linestring from break point to end node
	void InsertingExtendedNodes(vertex_descriptor_drct source, vertex_descriptor_drct target, int extend_node_id, int direction,
		OGRPoint breakpoint, OGRLineString* s_cutedge, OGRLineString* e_cutedge);
};
}

#endif
