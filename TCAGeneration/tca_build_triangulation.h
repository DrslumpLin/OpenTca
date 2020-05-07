// build contrianed triangulation and normal Delaunay triangulation
// as preparation for generating the contour
// @author: Diao Lin
// @version: 2020.03

#pragma once
#ifndef TCA_BUILD_TRIANGULATION_H
#define TCA_BUILD_TRIANGULATION_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <unordered_map>

#include "tca_types.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>


// for voronoi diagram
//#include <CGAL/Voronoi_diagram_2.h>
//#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
//#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include "ogrsf_frmts.h" // C++ API for GDAL



using namespace std;
namespace tca
{
typedef CGAL::Exact_predicates_inexact_constructions_kernel					K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>	    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>						    Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>							    Delaunay;
typedef Delaunay::Point													    Point;

typedef CGAL::Constrained_triangulation_face_base_2<K>						      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>						      Tds_C;
typedef CGAL::Exact_predicates_tag											      Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds_C, Itag>			      CDT;
typedef CDT::Vertex_handle														  Vertex_handle;

class BuildTriangulation
{

	struct PointM
	{
		Point point;
		//PointIndex();
		PointM(const Point & _inpoint) : point(_inpoint) {};

		bool operator<(const PointM& other) const
		{
			if (point.x() != other.point.x())
				return (point.x() < other.point.x());
			return (point.y() < other.point.y());
		}
	};


private:
	// read from shapefile or pass from "tca_undirected_graph"
	// the x,y coordinates of the input nodes
	std::vector<double>							x_coords_;
	std::vector<double>							y_coords_;
	
	// z is the value for interpolation
	std::vector<double>							z_values_;

	// a pair of (node coordinates, node inner id)
	std::vector<std::pair<Point, unsigned> >	node_points_inner_ids_pair_;
	
	//the external ids of the nodes
	std::vector<std::string>					node_external_ids_;
	
	// nodes of all the edges
	std::vector<std::array<int, 2>>				edge_node_ids_;
	
	//the nodes of triangulation obtained by "BuildDelaunayTrangulation"
	//as input for generating the catchment areas
	std::vector<std::array<int, 3>>				triangle_node_ids_;
	
	// the polygons of triangulaiton 
	std::vector<OGRPolygon>						tri_polygons_;
	
	Delaunay delaunary_;
	CDT con_delau_tris_;
	
	// an unoderred map between the node inner id and its vertex handle
	std::unordered_map<int, Vertex_handle>      inner_id_vh;

	// for testing purposes
	vector<OGRLineString*> constrained_segs;
	vector<int> constrained_seg_egde_ids;
	vector<int> constrained_seg_snode_ids;
	vector<int> constrained_seg_enode_ids;

public:
	
	//defalt constructor
	BuildTriangulation();

	//initialization of node_external_ids_, z_values_, node_org_points
	BuildTriangulation(const vector<string> node_ids, const vector<double> node_distances, const vector<OGRPoint> node_org_points);
	
	~BuildTriangulation();

	// initialization from the shp files instead of constructor
	void ReadNodeCosts(const std::string & file, const std::string &id_name = "nodeid", const std::string & zvalue = "distance");

	// build normal delaunary triangulations based on the input nodes
	void BuildDelaunayTrangulation(void);
	
	// build constrainted triangulations based on the input nodes and constrained edges
	// input:
	//   @edge_node_ids:constrianed edges
	//   @node_ids: external node ids
	//   @node_distances: distances to the root node
	//   @node_org_points: coordinates of nodes
	void BuildConstrainedTrangulation(std::vector<std::array<int, 2>> edge_node_ids, const vector<string> node_ids,
		const vector<double> node_distances, const vector<OGRPoint> node_org_points);
	
	// build constrainted trangulations for multiple-points facility
	// Input is the same as that of "BuildConstrainedTrangulation"
	// compared with the point-based facility, some addtional process are needed,e.g., 
	// the excluded the virtual node during the triangulation 
	void BuildConstrainedTrangulationMultiplePointsVersion(std::vector<std::array<int, 2>> edge_node_ids, const vector<string> node_ids,
		const vector<double> node_distances, const vector<OGRPoint> node_org_points);
	
	//-------get the info for generating the contours (i.e. catchment areas)---------
	// including the x_coords; y_coords; z_values; triangle_node_ids; masks_of_triangles;
	std::vector<double> get_x_coords();
	std::vector<double> get_y_coords();
	std::vector<double> get_z_values();
	std::vector<std::array<int, 3>> get_triangle_node_ids();
	
	int get_node_num();
	int get_triangle_num();
	std::vector<bool> get_masks(); // default to be false
	
	//--------------IO fucntions------------------------------
	
	// for normal triangulations
	void write_triangles(const std::string filename);
	void write_triangles_nodes(const std::string filename);

	// for constrianed triangulations 
	void write_con_triangles_nodes(const std::string filename);
	void write_con_edges(const std::string filename);
};

}// namespace tca

#endif