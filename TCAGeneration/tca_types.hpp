// Definition of basic data types used in the TCA
// @author: Diao Lin
// @version: 2020.03

#ifndef TCA_TYPES_HPP
#define TCA_TYPES_HPP
#include <vector>
#include <string>
#include "ogrsf_frmts.h" // C++ API for GDAL

namespace tca
{
namespace types 
{
	// the data type of for the edge of undirected road netwrok
	struct Edge
	{
		int         id;
		std::string external_id; // This is the external ID attribute, which does not have to be continuous
		std::string source;      // source node ID
		std::string target;      // target node ID
		double length;	         // length of the edge polyline
		OGRLineString * line_string;     // a pointer to the edge geometry
	};
	
	// the data type for a point facility 
	struct Facility
	{
		int id;				       // This is the id, which is continuous distributed 
		std::string external_id;   // This is the external ID attribute, which does not have to be continuous
		OGRPoint * point;          // point coordinates
		double     delta;          // the cutoff cost used to generate the catchment areas
	};

	// the nearest point of a facility, i.e., the projected point of facility 
	struct NearestPoint
	{
		int		  facility_id;         // the corresponding facility ID
		double    p_2_pp;              // distance from p to pp (i.e. the projected point of p)
		double    pp_2_s;              // pp distance to the start point
		int       cut_snode_index;     // for linestring with multiple points [node_1,node_2, node_n], and cut_point located at [node_i,node_i+1], then cut_snode_index = i
		OGRPoint* point;               // projected point
		Edge *    edge;                // nearest edge
	};

	// a set of nearest points for all facilities
	typedef std::vector<NearestPoint> NearestPoints; 
	

	// an accessible node of a facility
	struct AccessibleNode
	{
		int facility_id;		// the corresponding facility id
		int node_id;			// the node id (i.e., graph node)
		int succ_nid;			// the successor node id (of the shortest path)
		int pred_nid;			// the preceding node id (of the shortest path)
		double cost;			// the cost from facility to this node
		OGRPoint* point;		// the coordinates 
	};

	// an accessibile edge of a facility 
	struct AccessibleEdge
	{
		int facility_id;
		double s_cost;			// the cost to the starting node of the edge
		double t_cost;			// the cost to the terminal node of the edge
		OGRLineString* linestring;
	};

	//the results of all accessable nodes of facilities
	typedef std::vector<AccessibleNode>  AccessibleNodes;

	// the egde for voronoi egde
	struct NetworkVoronoiEdge
	{
		std::string			facility_id;					// the facility this edge belongs to. For example, "fac1" or "fac2",...
		std::string			origanl_edge_external_id;		// the original external id of the edge
		OGRLineString*		linestring;					    // the linestrings
	};

	// the node for network voronoi 
	struct NetworkVoronoiNode
	{
		std::string			facility_id_1;					// the facility this node belongs to
		std::string         facility_id_2;                  // only under the condition that a node is a breaknode,  
															// this filed is effective, because the break node belongs to two 
															// different facilties at the same time the facilty_1 and  facilty_2 
															// correspond to the s and t node direction, respectively

		int			        inner_node_id;					// the inner node id this 
		bool                is_break_node;                  // a label mark if a point is a break point or not
		double              direction_x;                     //the x value of the unit normal vector, only used for breek node 
		double              direction_y;                    // the y value of the unit normal vector, only used for breek node
		OGRPoint		    point;					        // the coordinates of this node
		NetworkVoronoiNode() :facility_id_2("nulll"), inner_node_id(-1),is_break_node(false), direction_x(-999), direction_y(-999) {}
	};

	//-------------------special for directed road network input-------------------------------

	// for directed edge,a direction need to be added
	struct EdgeDirected
	{
		int         id;
		std::string external_id; // This is the external ID attribute, which does not have to be continuous
		std::string source;      // source node ID
		std::string target;      // target node ID
		double length;	         // length of the edge polyline
		int direction;			 // 1(for bi-edge, start-terminal(S-T) direction ); -1(for bi-edge, T-S); 2(for single edge, S-T); 3 (for single edge T-S);
		OGRLineString * line_string;     // a pointer to the edge geometry
	};
	
	// the neareastpoint for directed network input
	struct NearestPointDirected
	{
		int					facility_id;         // the corresponding facility ID
		double				p_2_pp;              // distance from p to pp (i.e. the projected point of p)
		double				pp_2_s;              // pp distance to the start point
		int					cut_snode_index;     /* for linestring with multiple points [node_1,node_2, node_n], 
												   and cut_point located at [node_i,node_i+1], then cut_snode_index = i 
												 */
		OGRPoint*			point;               // projected point
		EdgeDirected *		edge;                // nearest edge with directed info 
	};

	// a set of nearest points of all facilities for directed road network
	typedef std::vector<NearestPointDirected> NearestPointsDirected; // nearest candidates for all facilities (directed version)
	

	//-------------------for benchmark generation(grid_points and extended sp tree (subgraph))-------------------------------
	struct SubGraphEdge
	{
		int         id;
		//std::string external_id; 			// This is the external ID attribute, which does not have to be continuous
		//std::string source;      			// source node ID
		//std::string target;      			// target node ID
		double s_dis;			 			// s_distance
		//double e_dis;			 			// e_distance
		//double length;	         	  	// length of the edge polyline
		OGRLineString * line_string;     	// a pointer to the edge geometry
	};

}// namespace types
}// namespace tca

#endif