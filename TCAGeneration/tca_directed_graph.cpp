#include "tca_directed_graph.h"

namespace tca 
{


DirectedGraph::DirectedGraph(std::vector<EdgeDirected*> egdes)
{
	int current_idx = -1;
	edge_descriptor_drct e;
	bool inserted;
	int N = egdes.size();
	int source_idx = 0;
	int target_idx = 0;

	printf("Network edges :%d \n", N);

	OGRPoint spoint;
	OGRPoint epoint;

	for (int i = 0; i < N; ++i)
	{
		EdgeDirected *network_edge = egdes[i];
		network_edge->line_string->StartPoint(&spoint);
		network_edge->line_string->EndPoint(&epoint);

		// Search for source node idx
		auto search = vertexs_externalid_discriptor_.find(network_edge->source);
		if (search != vertexs_externalid_discriptor_.end())
		{
			// A node exists already
			source_idx = search->second;
		}
		else
		{
			// A new node is found
			++current_idx;
			vertexs_externalid_discriptor_.insert({ network_edge->source,current_idx });
			source_idx = current_idx;
			vertex_externalids_.push_back(network_edge->source);
			vertex_coords_.push_back(spoint);
		};

		// Search for target node idx
		search = vertexs_externalid_discriptor_.find(network_edge->target);
		if (search != vertexs_externalid_discriptor_.end())
		{
			// A node exists already
			target_idx = search->second;
		}
		else
		{
			// A new node is found
			++current_idx;
			vertexs_externalid_discriptor_.insert({ network_edge->target,current_idx });
			target_idx = current_idx;
			vertex_externalids_.push_back(network_edge->target);
			vertex_coords_.push_back(epoint);
		};

		boost::tie(e, inserted) = add_edge(source_idx, target_idx, directed_g_);
		directed_g_[e].id = network_edge->external_id;
		directed_g_[e].direction = network_edge->direction;
		directed_g_[e].length = network_edge->length;

		edges_externalid_linestring.insert({ network_edge->external_id, network_edge->line_string });
	}

	std::cout <<"number of edges" << boost::num_edges(directed_g_) << '\n';
	std::cout << "Construct graph from network edges end" << '\n';

	//edge_iterator ei, ei_end;
	//for (boost::tie(ei, ei_end) = boost::edges(directed_g_); ei != ei_end; ++ei)
	//{
	//	vertex_descriptor_drct  vt_s_node = boost::source(*ei, directed_g_);
	//	vertex_descriptor_drct  vt_e_node = boost::target(*ei, directed_g_);

	//	std::cout << "edge start node" << vt_s_node << '\n';
	//	std::cout << "edge end node" << vt_e_node << '\n';
	//}
}



DirectedGraph::~DirectedGraph()
{
	std::cout << "Cleaning the undirected graph " << '\n';
	for (auto &item : edges_externalid_linestring)
	{
		//OGRGeometryFactory::destroyGeometry(item.second);
	}
	std::cout << "Cleaning undirected graph finished" << '\n';
}


void DirectedGraph::InsertingFacilityNode(const NearestPointDirected ne_point)
{	
	OGRLineString *s_cutedge = new OGRLineString();
	OGRLineString *e_cutedge = new OGRLineString();
	algorithm::CutLinestringBasedOnBreakPoint(ne_point.edge->line_string, ne_point.cut_snode_index, ne_point.point, s_cutedge, e_cutedge);

	edge_descriptor_drct e;
	bool inserted;
	int temp_s_id = 0;
	int temp_e_id = 0;

	std::string temp_source = ne_point.edge->source;
	std::string temp_target = ne_point.edge->target;

	auto search_s = vertexs_externalid_discriptor_.find(ne_point.edge->source);
	if (search_s != vertexs_externalid_discriptor_.end())
	{
		temp_s_id = search_s->second;
	}
	auto search_t = vertexs_externalid_discriptor_.find(ne_point.edge->target);
	if (search_t != vertexs_externalid_discriptor_.end())
	{
		temp_e_id = search_t->second;
	}

	int edge_direction = abs(ne_point.edge->direction);
	
	// for single-direction edge
	if (edge_direction != 1) 
	{
		remove_edge(temp_s_id, temp_e_id, directed_g_);

		// insert one node and two edges
		int insert_node_id = boost::num_vertices(directed_g_);

		boost::tie(e, inserted) = add_edge(temp_s_id, insert_node_id, directed_g_);
		directed_g_[e].id = "fa_egde_s_" + std::to_string(ne_point.facility_id);
		directed_g_[e].length = ne_point.pp_2_s;
		directed_g_[e].direction = edge_direction; //
		edges_externalid_linestring.insert({ directed_g_[e].id, s_cutedge });

		boost::tie(e, inserted) = add_edge(insert_node_id, temp_e_id, directed_g_);
		directed_g_[e].id = "fa_egde_t_" + std::to_string(ne_point.facility_id);
		directed_g_[e].length = ne_point.edge->length - ne_point.pp_2_s;
		//added by diao, precision issue
		if (directed_g_[e].length < 0)
		{
			directed_g_[e].length = 0;
		}

		directed_g_[e].direction = edge_direction; //
		edges_externalid_linestring.insert({ directed_g_[e].id, e_cutedge });

		std::string insert_node_sid = "fac" + to_string(ne_point.facility_id);
		vertexs_externalid_discriptor_.insert({ insert_node_sid,insert_node_id });
		vertex_externalids_.push_back(insert_node_sid);
		vertex_coords_.push_back(*(ne_point.point));

	}

	else 	// edge_direction = 1 or = -1
	{
		remove_edge(temp_s_id, temp_e_id, directed_g_);
		remove_edge(temp_e_id, temp_s_id, directed_g_);

		// insert one node and two edges
		int insert_node_id = boost::num_vertices(directed_g_);

		// add the first direction edge
		boost::tie(e, inserted) = add_edge(temp_s_id, insert_node_id, directed_g_);
		directed_g_[e].id = "fa_egde_s_" + std::to_string(ne_point.facility_id);
		directed_g_[e].length = ne_point.pp_2_s;
		directed_g_[e].direction = edge_direction; //
		edges_externalid_linestring.insert({ directed_g_[e].id, s_cutedge });
		
		boost::tie(e, inserted) = add_edge(insert_node_id, temp_e_id, directed_g_);
		directed_g_[e].id = "fa_egde_t_" + std::to_string(ne_point.facility_id);
		directed_g_[e].length = ne_point.edge->length - ne_point.pp_2_s;

		//added by diao, precision issue
		if (directed_g_[e].length < 0)
		{
			directed_g_[e].length = 0;
		}

		directed_g_[e].direction = edge_direction; //
		edges_externalid_linestring.insert({ directed_g_[e].id, e_cutedge });


		// add the second (oppsite) direction edge
		boost::tie(e, inserted) = add_edge(temp_e_id, insert_node_id, directed_g_);
		directed_g_[e].id = "fa_egde_s_op_" + std::to_string(ne_point.facility_id);
		directed_g_[e].length = ne_point.edge->length - ne_point.pp_2_s;  
		directed_g_[e].direction = 0-edge_direction; //

		//added by diao, precision issue
		if (directed_g_[e].length < 0)
		{
			directed_g_[e].length = 0;
		}
		
		OGRLineString* e_cutedge_re = (OGRLineString*)e_cutedge->clone(); // here must clone once
		e_cutedge_re->reversePoints();
		edges_externalid_linestring.insert({ directed_g_[e].id,e_cutedge_re });
		
		boost::tie(e, inserted) = add_edge(insert_node_id, temp_s_id, directed_g_);
		directed_g_[e].id = "fa_egde_t_op_" + std::to_string(ne_point.facility_id);
		directed_g_[e].length = ne_point.pp_2_s;
		directed_g_[e].direction = 0-edge_direction; //

		OGRLineString* s_cutedge_re = (OGRLineString*)s_cutedge->clone(); // here must clone once
		s_cutedge_re->reversePoints();
		edges_externalid_linestring.insert({ directed_g_[e].id,s_cutedge_re });	
		
		std::string insert_node_sid = "fac" + to_string(ne_point.facility_id);
		vertexs_externalid_discriptor_.insert({ insert_node_sid,insert_node_id });
		vertex_externalids_.push_back(insert_node_sid);
		vertex_coords_.push_back(*(ne_point.point));
	}
};


/*
* Added for multiple-source dijkstra algorithm
* 1) insert all the facility nodes to graph
* 2) Add virtual node and links between virtual node and facility nodes
*/
void DirectedGraph::InsertingMulltipleFacilityNodes(const NearestPointsDirected ne_points)
{
	// inserting facility nodes
	for (int i = 0; i < ne_points.size(); i++) {
		NearestPointDirected temp_ne_point = ne_points[i];
		InsertingFacilityNode(temp_ne_point);
	}
	// inserting virtual nodes and links to virtual nodes	
	OGRPoint *vir_point = new OGRPoint(99999, 99999);
	OGRLineString *vir_edge = new OGRLineString();
	vir_edge->addPoint(vir_point);
	vir_edge->addPoint(vir_point);
	int virtual_node_id = boost::num_vertices(directed_g_);
	vertexs_externalid_discriptor_.insert({ "vir_node", virtual_node_id });
	vertex_externalids_.push_back("vir_node");
	vertex_coords_.push_back(*(vir_point));

	// Add links between virtual node and facility nodes
	for (int i = 0; i < ne_points.size(); i++)
	{
		edge_descriptor_drct e;
		bool inserted;
		std::string temp_f_node_sid = "fac" + to_string(i);
		int temp_facility_id = 0;
		auto search_f = vertexs_externalid_discriptor_.find(temp_f_node_sid);
		if (search_f != vertexs_externalid_discriptor_.end())
		{
			// A node exists already
			temp_facility_id = search_f->second;
		}
		//updating edges by adding a virtual edge
		boost::tie(e, inserted) = add_edge(virtual_node_id, temp_facility_id, directed_g_);
		directed_g_[e].id = "vir_fa_" + to_string(i);
		directed_g_[e].direction = 2;
		directed_g_[e].length = 0;
		edges_externalid_linestring.insert({ directed_g_[e].id, vir_edge });
	}
}


// calculating the shortest path from the virtual nodes
void DirectedGraph::MultipleSourceShortestPaths(const std::string insert_node_sid)
{
	std::vector<vertex_descriptor_drct> examined_nodes;
	int num_vertices = boost::num_vertices(directed_g_);
	InitializeDistancesPredecessors();
	InitializeClosestFacilities(num_vertices);

	vertex_iterator_drct vi, vend;
	for (boost::tie(vi, vend) = vertices(directed_g_); vi != vend; ++vi)
	{
		if (vertex_externalids_[*vi] == insert_node_sid)
		{
			MultipleSourceShortestPathAlgorithm(*vi, examined_nodes, closest_facilities_map_);
		}
	}
}


// Added by Diao; multiple source dijkstra shortest path algorithm
void DirectedGraph::MultipleSourceShortestPathAlgorithm(const vertex_descriptor_drct& source, std::vector<vertex_descriptor_drct>& p_examined_nodes,
	std::vector<vertex_descriptor_drct>& closest_facilities)
{
	distances_map_[source] = 0; // here is necessay, or an error will be reported, please see the info 
	double inf = std::numeric_limits<double>::max();

	dijkstra_shortest_paths_no_color_map_no_init
	(
		directed_g_,
		source,
		make_iterator_property_map(predecessors_map_.begin(), get(boost::vertex_index, directed_g_), predecessors_map_[0]),
		make_iterator_property_map(distances_map_.begin(), get(boost::vertex_index, directed_g_), distances_map_[0]),
		get(&EdgePropertyDrct::length, directed_g_),
		get(boost::vertex_index, directed_g_),
		std::less<double>(), //DistanceCompare distance_compare,
		boost::closed_plus<double>(inf),
		inf,
		0,
		MS_distance_visitor_drct(source, p_examined_nodes, closest_facilities)
	);
};// end for MS_shortest_path_algorithm 



// here is different from the undirected version
std::vector<NetworkVoronoiEdge> DirectedGraph::ConstructNetworkVoronoiDiagram(void)
{
	std::vector<NetworkVoronoiEdge> nv_edges;

	edge_iterator_drct ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(directed_g_); ei != ei_end; ++ei)
	{
		vertex_descriptor_drct  vt_s_node = boost::source(*ei, directed_g_);
		vertex_descriptor_drct  vt_e_node = boost::target(*ei, directed_g_);

		vertex_descriptor_drct s_closest_faci = closest_facilities_map_[vt_s_node];
		vertex_descriptor_drct e_closest_faci = closest_facilities_map_[vt_e_node];

		double s_dis = distances_map_[vt_s_node];
		double e_dis = distances_map_[vt_e_node];

		auto search = edges_externalid_linestring.find(directed_g_[*ei].id);
		OGRLineString* pLine = search->second;

		OGRPoint p1, p2;
		p1 = vertex_coords_[vt_s_node];

		int direction = directed_g_[*ei].direction;

		if ((direction == 2) || (direction == 3)) // for signle direction s_closest_faci and e_closest_faci must equal to each other
		{
			NetworkVoronoiEdge nv_edge;
			nv_edge.facility_id = vertex_externalids_[s_closest_faci];
			nv_edge.origanl_edge_external_id = directed_g_[*ei].id;
			nv_edge.linestring = pLine;
			nv_edges.push_back(nv_edge);
		}// end if 
		if (direction == 1) // for bi-direction edge, we only need to consider one condition
		{
			if (s_closest_faci == e_closest_faci)
			{
				NetworkVoronoiEdge nv_edge;
				nv_edge.facility_id = vertex_externalids_[s_closest_faci];
				nv_edge.origanl_edge_external_id = directed_g_[*ei].id;
				nv_edge.linestring = pLine;
				nv_edges.push_back(nv_edge);
			}
			else
			{
				// unconnected edges and virtual edges should be excluded 
				// the != 0 and != max() is not a very good style, this needs change in the future
				if ((s_dis != std::numeric_limits<double>::max())
					&& (e_dis != std::numeric_limits<double>::max())
					&& (directed_g_[*ei].length != 0.0))
				{
					double break_to_snode = (distances_map_[vt_e_node] - distances_map_[vt_s_node] + directed_g_[*ei].length) / 2;

					OGRLineString *s_cutedge = new OGRLineString();
					OGRLineString *e_cutedge = new OGRLineString();;
					OGRPoint *breakpoint = new OGRPoint();

					algorithm::CutLinestringBasedOnBreakDistance(pLine, break_to_snode, breakpoint, s_cutedge, e_cutedge);

					NetworkVoronoiEdge nv_s_edge, nv_e_edge;
					nv_s_edge.facility_id = vertex_externalids_[s_closest_faci];
					nv_s_edge.origanl_edge_external_id = directed_g_[*ei].id;
					nv_s_edge.linestring = s_cutedge;
					nv_edges.push_back(nv_s_edge);

					nv_e_edge.facility_id = vertex_externalids_[e_closest_faci];
					nv_e_edge.origanl_edge_external_id = directed_g_[*ei].id;
					nv_e_edge.linestring = e_cutedge;
					nv_edges.push_back(nv_e_edge);
				}
			}// end if 
		}// end if
	}// end for 
	return nv_edges;
}

void DirectedGraph::write_network_voronois_edges(const std::string &filename, std::vector<NetworkVoronoiEdge> vd_edges)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("network_vornoi", NULL, wkbLineString, NULL);

	// define and create filed
	OGRFieldDefn oField_edge_id("id", OFTInteger);
	oField_edge_id.SetWidth(10);
	poLayer->CreateField(&oField_edge_id);

	OGRFieldDefn oField_edge_label("label", OFTString);
	oField_edge_label.SetWidth(30);
	poLayer->CreateField(&oField_edge_label);

	OGRFieldDefn oField_edge_original_id("ori_edgeid", OFTString);
	oField_edge_original_id.SetWidth(30);
	poLayer->CreateField(&oField_edge_original_id);

	int line_nums = vd_edges.size();
	for (int i = 0; i < line_nums; i++)
	{
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());

		poFeature->SetField("id", i);
		poFeature->SetField("label", vd_edges[i].facility_id.c_str());
		poFeature->SetField("ori_edgeid", vd_edges[i].origanl_edge_external_id.c_str());
		poFeature->SetGeometry(vd_edges[i].linestring);

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(shpDataSet);
}

void DirectedGraph::write_multiple_sources_node_distances(std::string filename)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("MSnodes_labels", NULL, wkbPoint, NULL);

	// define and create filed
	OGRFieldDefn oField_distance("distance", OFTReal);
	oField_distance.SetWidth(20);
	poLayer->CreateField(&oField_distance);

	OGRFieldDefn oField_edge_label("label", OFTString);
	oField_edge_label.SetWidth(30);
	poLayer->CreateField(&oField_edge_label);

	std::pair<vertex_iterator_drct, vertex_iterator_drct > vp;
	for (vp = boost::vertices(directed_g_); vp.first != vp.second; ++vp.first)
	{
		double dis = distances_map_[*vp.first];		
		if (dis != std::numeric_limits<double>::max())
		{		
			vertex_descriptor_drct closest_faci = closest_facilities_map_[*vp.first];
			std::string label = vertex_externalids_[closest_faci];
			OGRPoint p1 = vertex_coords_[*vp.first];

			OGRFeature *poFeature;
			poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
			poFeature->SetField("distance", dis);
			poFeature->SetField("label", label.c_str());
			poFeature->SetGeometry(&p1);

			if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
			{
				printf("Failed to create feature in shapefile.\n");
				exit(1);
			}
			OGRFeature::DestroyFeature(poFeature);
		}
	}// end for
	GDALClose(shpDataSet);
}


std::vector<vertex_descriptor_drct> DirectedGraph::GetShortestPathNodes(const std::string facility_node_sid)
{
	std::vector<vertex_descriptor_drct> examined_nodes;
	int num_vertices = boost::num_vertices(directed_g_);
	InitializeDistancesPredecessors();
	vertex_iterator_drct vi, vend;
	for (boost::tie(vi, vend) = vertices(directed_g_); vi != vend; ++vi)
	{
		if (vertex_externalids_[*vi] == facility_node_sid)
		{
			ShortestPathAlgorithm(*vi, examined_nodes);
		}
	}
	return examined_nodes;
	//wirte_sp_nodes(filename, examined_nodes);
}


std::vector<vertex_descriptor_drct> DirectedGraph::BuildExtendedShortestPathTree(std::vector<vertex_descriptor_drct> examined_nodes, const std::string facility_node_sid)
{
	std::vector<ExtendedNodeDrct> extended_nodes;
	extended_nodes = FindNonShortestPathEdges(examined_nodes);
	InsertingNonShortestPathEdges(extended_nodes);

	std::vector<vertex_descriptor_drct> examined_nodes_extended_sp;
	int num_vertices = boost::num_vertices(directed_g_);
	InitializeDistancesPredecessors();

	vertex_iterator_drct vi, vend;
	for (boost::tie(vi, vend) = vertices(directed_g_); vi != vend; ++vi)
	{
		if (vertex_externalids_[*vi] == facility_node_sid)
		{
			ShortestPathAlgorithm(*vi, examined_nodes_extended_sp);
		}
	}

	return examined_nodes_extended_sp;
	//wirte_sp_nodes(filename, examined_nodes_extended_sp);
}



void DirectedGraph::ShortestPathAlgorithm(const vertex_descriptor_drct& source, std::vector<vertex_descriptor_drct>& examined_nodes)
{
	//std::vector<vertex_descriptor_drct> m_examined_nodes;
	std::vector<edge_descriptor_drct> tree_edges;
	std::vector<edge_descriptor_drct> non_tree_edges;
	//examined_nodes.push_back(source);
	distances_map_[source] = 0; // here is necessay, or an error will be reported, please see the info 
	//double delta = 1000;
	double inf = std::numeric_limits<double>::max();
	dijkstra_shortest_paths_no_color_map_no_init
	(
		directed_g_,
		source,
		make_iterator_property_map(predecessors_map_.begin(), get(boost::vertex_index, directed_g_), predecessors_map_[0]),
		make_iterator_property_map(distances_map_.begin(), get(boost::vertex_index, directed_g_), distances_map_[0]),
		get(&EdgePropertyDrct::length, directed_g_),
		get(boost::vertex_index, directed_g_),
		std::less<double>(), //DistanceCompare distance_compare,
		boost::closed_plus<double>(inf),
		inf,
		0,
		CA_vistorDrct(examined_nodes, distances_map_)
	);

	//p_examined_nodes = m_examined_nodes;

};// end for shortest_path_tree 


void DirectedGraph::wirte_sp_nodes(const std::string &filename, const std::vector<vertex_descriptor_drct> examined_nodes)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("En_sp_tree", NULL, wkbPoint, NULL);

	// define and create filed
	OGRFieldDefn oFie_faci_Id("nodeid", OFTString);
	oFie_faci_Id.SetWidth(32);
	poLayer->CreateField(&oFie_faci_Id);

	OGRFieldDefn oField_node_id("pre_node", OFTString);
	oField_node_id.SetWidth(30);
	poLayer->CreateField(&oField_node_id);

	OGRFieldDefn oField_distance("distance", OFTReal);
	oField_distance.SetWidth(30);
	poLayer->CreateField(&oField_distance);

	int number_record = examined_nodes.size();

	vertex_descriptor_drct u;
	for (int i = 0; i < number_record; i++)
	{
		u = examined_nodes[i];
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		poFeature->SetField("nodeid", vertex_externalids_[u].c_str());
		poFeature->SetField("pre_node", vertex_externalids_[predecessors_map_[u]].c_str());
		poFeature->SetField("distance", distances_map_[u]);
		poFeature->SetGeometry(&(vertex_coords_[u]));

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(shpDataSet);
};


void DirectedGraph::write_sp_edges(const std::string &filename, const std::vector<vertex_descriptor_drct> examined_nodes)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("sp_edges", NULL, wkbLineString, NULL);

	// define and create filed
	OGRFieldDefn oFie_edge_Id("edge_id", OFTInteger);
	oFie_edge_Id.SetWidth(10);
	poLayer->CreateField(&oFie_edge_Id);

	OGRFieldDefn oField_s_node("s_node", OFTString);
	oField_s_node.SetWidth(30);
	poLayer->CreateField(&oField_s_node);

	OGRFieldDefn oField_t_node("e_node", OFTString);
	oField_t_node.SetWidth(30);
	poLayer->CreateField(&oField_t_node);

	vertex_descriptor_drct u, pre_u;
	for (int i = 1; i < examined_nodes.size(); i++)
	{
		u = examined_nodes[i];
		pre_u = predecessors_map_[u];

		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		OGRLineString *temp_line = new OGRLineString();
		*temp_line = GetEdgeLinestring(pre_u, u);

		poFeature->SetField("edge_id", i);
		poFeature->SetField("s_node", vertex_externalids_[pre_u].c_str());
		poFeature->SetField("e_node", vertex_externalids_[u].c_str());
		poFeature->SetGeometry(temp_line);

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(shpDataSet);
}


// based on vertex of descriptor of source and target node to find the edge linestring
OGRLineString DirectedGraph::GetEdgeLinestring(vertex_descriptor_drct source, vertex_descriptor_drct target)
{
	OGRLineString templinestring;
	edge_descriptor_drct e;
	out_edge_iterator_drct out_i, out_end;
	OGRPoint s_point, e_point;

	for (boost::tie(out_i, out_end) = boost::out_edges(source, directed_g_); out_i != out_end; ++out_i)
	{
		s_point = vertex_coords_[source];
		e = *out_i;

		if (target == boost::target(e, directed_g_))
		{
			auto search = edges_externalid_linestring.find(directed_g_[e].id);
			templinestring = *(search->second);
		}
	}
	return templinestring;
}


void DirectedGraph::InsertingNonShortestPathEdges(const std::vector<ExtendedNodeDrct> extended_records)
{
	ExtendedNodeDrct temprecord;
	int N = extended_records.size();

	std::cout << "the number of edges not in the sp tree is : " << N << std::endl;
	std::cout << "Original graph - number of vertexs " << boost::num_vertices(directed_g_) << std::endl;
	std::cout << "Original graph - number of edges " << boost::num_edges(directed_g_) << std::endl;

	for (int i = 0; i < N; i++)
	{
		OGRLineString *s_cutedge = new OGRLineString();
		OGRLineString *e_cutedge = new OGRLineString();;
		OGRPoint *breakpoint = new OGRPoint();

		temprecord = extended_records[i];
		algorithm::CutLinestringBasedOnBreakDistance(temprecord.pLine, temprecord.break_distance_to_snode, breakpoint, s_cutedge, e_cutedge);
		InsertingExtendedNodes(temprecord.source, temprecord.target, temprecord.node_id, temprecord.direction, *breakpoint, s_cutedge, e_cutedge);
	}
	std::cout << "Updating graph - number of vertexs " << boost::num_vertices(directed_g_) << std::endl;
	std::cout << "Updating graph - number of edges " << boost::num_edges(directed_g_) << std::endl;
}


//For directed edge
std::vector<ExtendedNodeDrct> DirectedGraph::FindNonShortestPathEdges(std::vector<vertex_descriptor_drct> examined_nodes)
{
	std::vector<ExtendedNodeDrct> extend_nodes;
	int extend_node_id = 0;

	// check how many node has been exaimed
	int ex_node_num = examined_nodes.size();
	std::vector<std::string> colorMap(distances_map_.size());

	// label every node as "white"
	for (int i = 0; i < distances_map_.size(); i++)
	{
		colorMap[i] = "white";
	}
	vertex_descriptor_drct u, pre_u;

	// from the last exaimed node to the first node
	for (int i = ex_node_num - 1; i >= 0; i--)
	{
		u = examined_nodes[i];
		pre_u = predecessors_map_[u];
		colorMap[u] = "Gray";

		edge_descriptor_drct e;
		in_edge_iterator in_i, in_end;
		vertex_descriptor_drct vi_source;

		// iterating every in-edge of the current node "u"
		// if the target node of this outer edge is not the precedeer of node "u" and is not iterated
		for (boost::tie(in_i, in_end) = boost::in_edges(u, directed_g_); in_i != in_end; ++in_i)
		{
			e = *in_i;
			std::string edge_id = directed_g_[e].id;			
			int temp_direction = abs(directed_g_[e].direction);
			vi_source = boost::source(e, directed_g_);

			if ((vi_source != pre_u) && (colorMap[vi_source] == "white"))
			{
				if (temp_direction == 1) //directed edge are not considered, here only the bi-direction edges are considered
				{
					auto search = edges_externalid_linestring.find(edge_id);
					OGRLineString* pLine = search->second;
					double break_to_snode = (distances_map_[u] - distances_map_[vi_source] + directed_g_[e].length) / 2;
					
					ExtendedNodeDrct temprecord;
					temprecord.pLine = pLine;
					temprecord.break_distance_to_snode = break_to_snode;
					temprecord.source = vi_source;
					temprecord.target = u;
					temprecord.direction = temp_direction;
					temprecord.node_id = extend_node_id;

					extend_nodes.push_back(temprecord);
					extend_node_id++;
				}
			}
		}
	}
	return extend_nodes;
};


void DirectedGraph::InsertingExtendedNodes(vertex_descriptor_drct source, vertex_descriptor_drct target, int extend_node_id,
						int direction, OGRPoint breakpoint, OGRLineString* s_cutedge, OGRLineString* e_cutedge)
{
	edge_descriptor_drct e;
	bool inserted;
	
	remove_edge(source, target, directed_g_);
	remove_edge(target, source, directed_g_);

	// insert one node and two edges for each direction
	int insert_node_id = boost::num_vertices(directed_g_);
	//if (insert_node_id == 822) 
	//{
	//	std::cout << "stop here for debug" << endl;
	//}
	// add the first two direction edges
	boost::tie(e, inserted) = add_edge(source, insert_node_id, directed_g_);
	directed_g_[e].id = "ex_egde_s_" + std::to_string(extend_node_id);
	directed_g_[e].length = s_cutedge->get_Length();
	directed_g_[e].direction = direction; //
	edges_externalid_linestring.insert({ directed_g_[e].id, s_cutedge });

	boost::tie(e, inserted) = add_edge(insert_node_id, target, directed_g_);
	directed_g_[e].id = "ex_egde_t_" + std::to_string(extend_node_id);
	directed_g_[e].length = e_cutedge->get_Length();
	directed_g_[e].direction = direction; //
	edges_externalid_linestring.insert({ directed_g_[e].id, e_cutedge });
	
	// add the second two (oppsite) direction edges
	boost::tie(e, inserted) = add_edge(target, insert_node_id, directed_g_);
	directed_g_[e].id = "ex_egde_s_op_" + std::to_string(extend_node_id);
	directed_g_[e].length = e_cutedge->get_Length();
	directed_g_[e].direction = 0 - direction; //

	OGRLineString* e_cutedge_re = (OGRLineString*)e_cutedge->clone();
	e_cutedge_re->reversePoints();
	edges_externalid_linestring.insert({ directed_g_[e].id,e_cutedge_re });

	boost::tie(e, inserted) = add_edge(insert_node_id, source, directed_g_);
	directed_g_[e].id = "ex_egde_t_op_" + std::to_string(extend_node_id);
	directed_g_[e].length = s_cutedge->get_Length();
	directed_g_[e].direction = 0 - direction; //

	OGRLineString* s_cutedge_re = (OGRLineString*)s_cutedge->clone();
	s_cutedge_re->reversePoints();
	edges_externalid_linestring.insert({ directed_g_[e].id, s_cutedge_re});

	std::string insert_node_sid = "ex_node_" + to_string(extend_node_id);
	vertexs_externalid_discriptor_.insert({insert_node_sid,insert_node_id });
	vertex_externalids_.push_back(insert_node_sid);
	vertex_coords_.push_back(breakpoint);
}

// here is different from the undirected graph, beacause some single direction edge
// are not in the extend shortest path tree
std::vector<OGRPoint> DirectedGraph::FindingDeltaNodes(double delta)
{
	std::vector<OGRPoint> deltanodes;
	OGRPoint breakpoint;
	edge_iterator_drct ei, ei_end;

	for (boost::tie(ei, ei_end) = boost::edges(directed_g_); ei != ei_end; ++ei)
	{
		vertex_descriptor_drct  vt_s_node = boost::source(*ei, directed_g_);
		vertex_descriptor_drct  vt_e_node = boost::target(*ei, directed_g_);
		double s_dis = distances_map_[vt_s_node];
		double e_dis = distances_map_[vt_e_node];
		auto search = edges_externalid_linestring.find(directed_g_[*ei].id);
		OGRLineString* pLine = search->second;
		double edge_length = pLine->get_Length();
		int direction = abs(directed_g_[*ei].direction);
		
		// for unconnected edges
		if ((s_dis != std::numeric_limits<double>::max()) && (e_dis != std::numeric_limits<double>::max()))
		{
			if (direction == 1)
			{
				if ((s_dis < delta) && ( delta < e_dis))
				{
					double break_to_snode = delta - s_dis;
					breakpoint = algorithm::GetBreakPointOnLinestring(break_to_snode, pLine);
					deltanodes.push_back(breakpoint);
				}
			}
			else
			{
				if ((s_dis < delta) && (delta < (s_dis + edge_length)))
				{
					double break_to_snode = delta - s_dis;
					breakpoint = algorithm::GetBreakPointOnLinestring(break_to_snode, pLine);
					deltanodes.push_back(breakpoint);
				}
			}
		}
	}// end for

	return deltanodes;
};


// here is different from the undirected graph, beacause some single direction edge
// are not in the extend shortest path tree
std::vector<AccessibleEdge> DirectedGraph::GetAccessibleEdges(double delta, int facility_id, std::vector<vertex_descriptor_drct> ex_nodes)
{
	std::vector<AccessibleEdge> accessible_edges;
	edge_iterator_drct ei, ei_end;
	int dir_1 = 0;// for test
	int dir_s_b = 0;// for test
	for (boost::tie(ei, ei_end) = boost::edges(directed_g_); ei != ei_end; ++ei)
	{
		AccessibleEdge temp_acc_edge;
		vertex_descriptor_drct  vt_s_node = boost::source(*ei, directed_g_);
		vertex_descriptor_drct  vt_e_node = boost::target(*ei, directed_g_);
		double s_dis = distances_map_[vt_s_node];
		double e_dis = distances_map_[vt_e_node];
		auto search = edges_externalid_linestring.find(directed_g_[*ei].id);
		OGRLineString* pLine = search->second;
		double edge_length = pLine->get_Length();
		int direction = abs(directed_g_[*ei].direction);
		int ori_direction = directed_g_[*ei].direction; // for test


		OGRPoint *breakpoint = new OGRPoint();
		OGRLineString *s_cut_line = new OGRLineString();
		OGRLineString *e_cut_line = new OGRLineString();

		// for unconnected edges
		if ((s_dis != std::numeric_limits<double>::max()) && (e_dis != std::numeric_limits<double>::max()))
		{
			if (direction == 1)
			{				
				if ((s_dis < delta) && (delta < e_dis))
				{					
					double break_to_snode = delta - s_dis;
					algorithm::CutLinestringBasedOnBreakDistance(pLine, break_to_snode, breakpoint, s_cut_line, e_cut_line);
					temp_acc_edge.facility_id = facility_id;
					temp_acc_edge.s_cost = s_dis;
					temp_acc_edge.t_cost = delta;
					temp_acc_edge.linestring = s_cut_line;
					accessible_edges.push_back(temp_acc_edge);
				}
				else if ((s_dis <= delta) && (e_dis <= delta))
				{					
					if (s_dis <= e_dis) 
					{
						dir_s_b = dir_s_b + 1; // for test
						temp_acc_edge.facility_id = facility_id;
						temp_acc_edge.s_cost = s_dis;
						temp_acc_edge.t_cost = e_dis;
						temp_acc_edge.linestring = pLine;
						accessible_edges.push_back(temp_acc_edge);

						if (ori_direction == 1) {
							dir_1 = dir_1 + 1; // for test
						}
					}
				}
				else {
					continue;
				}// end if
			}
			else
			{
				if ((s_dis < delta) && (delta < (s_dis + edge_length)))
				{
					double break_to_snode = delta - s_dis;
					algorithm::CutLinestringBasedOnBreakDistance(pLine, break_to_snode, breakpoint, s_cut_line, e_cut_line);
					temp_acc_edge.facility_id = facility_id;
					temp_acc_edge.s_cost = s_dis;
					temp_acc_edge.t_cost = delta;
					temp_acc_edge.linestring = s_cut_line;
					accessible_edges.push_back(temp_acc_edge);
				}// end if 
			}// end if
		}// end if
	}// end for

	//for test
	std::cout << "start small then end: " << dir_s_b << " and the direction is s to e: "<< dir_1 << '\n';
	// start small then end: 108 and the direction is s to e: 56
	return accessible_edges;
};


void DirectedGraph::GetEdgesForConstrianedTriangulationLineSegmentation(const double &cutoff, std::vector<vertex_descriptor_drct> &extended_sp_exaimed_nodes,
	std::vector<std::array<int, 2>>	& edge_node_ids, std::vector<string>& node_ids,
	std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points)
{
	//std::vector<vertex_descriptor_drct> sp_exaimed_nodes = GetShortestPathNodes(facility_node_sid);
	//std::vector<vertex_descriptor_drct> extended_sp_exaimed_nodes = BuildExtendedShortestPathTree(sp_exaimed_nodes, facility_node_sid);
	node_ids = vertex_externalids_;
	node_distances = distances_map_;
	node_org_points = vertex_coords_;
	//int seg_node_id = vertex_coords_.size(); // started from the last node of orginal node

	int delta_id = 0;

	edge_iterator_drct ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(directed_g_); ei != ei_end; ++ei)
	{
		vertex_descriptor_drct  vt_s_node = boost::source(*ei, directed_g_);
		vertex_descriptor_drct  vt_e_node = boost::target(*ei, directed_g_);
		
		double s_dis = distances_map_[vt_s_node];
		double e_dis = distances_map_[vt_e_node];

		auto search = edges_externalid_linestring.find(directed_g_[*ei].id);
		OGRLineString* pLine = search->second;
		
		double edge_length = pLine->get_Length();
		int direction = abs(directed_g_[*ei].direction);
		
		//ex_edges.push_back(*pLine);
		//ex_snodes.push_back(int(vt_s_node));
		//ex_enodes.push_back(int(vt_e_node));

		// exclude the unconnected edge
		if ((s_dis != std::numeric_limits<double>::max()) && (e_dis != std::numeric_limits<double>::max()))
		{
			if (direction == 1) // for undirected edge
			{
				if (s_dis <= e_dis) // only need one direction,that traverse from source to target
				{
					// same procedure as undirected edge
					SplitLineToSegments(s_dis, vt_s_node, vt_e_node, pLine, edge_node_ids, node_ids, node_distances, node_org_points);
				}
			}
			else // single direction edge
			{
				if (predecessors_map_[vt_e_node] == vt_s_node) // if the edge is in the extended shortest path tree
				{
					// same procedure as undirected edge
					SplitLineToSegments(s_dis, vt_s_node, vt_e_node, pLine, edge_node_ids, node_ids, node_distances, node_org_points);
				}
				else // if the edge is not in the extended shortest path tree
				{
					if ((s_dis < cutoff) && (cutoff < (s_dis + edge_length)))
					{						
						OGRLineString *s_cutedge = new OGRLineString();
						OGRLineString *e_cutedge = new OGRLineString();
						OGRPoint *breakpoint = new OGRPoint();
						double break_to_snode = cutoff - s_dis;
						algorithm::CutLinestringBasedOnBreakDistance(pLine, break_to_snode, breakpoint, s_cutedge, e_cutedge);

						// add delta node to the node sets
						int node_idx = node_ids.size();
						std::string nodename = "delta" + std::to_string(delta_id);
						node_ids.push_back(nodename);
						node_distances.push_back(cutoff);
						node_org_points.push_back(*breakpoint);

						//split in s_cutedge
						SplitLineToSegments(s_dis, vt_s_node, node_idx, s_cutedge, edge_node_ids, node_ids, node_distances, node_org_points);
						// split in e_cut_edge
						SplitLineToSegments(cutoff, node_idx, vt_e_node, e_cutedge, edge_node_ids, node_ids, node_distances, node_org_points);

						delta_id = delta_id + 1;
					}
					else 
					{
						// same procedure as undirected edge, [delta]------s_dis--------------------e_dis-------[delta]
						// here there are some problems
						SplitLineToSegments(s_dis, vt_s_node, vt_e_node, pLine, edge_node_ids, node_ids, node_distances, node_org_points);
					}				
				}
			}
		}
	}
}



void  DirectedGraph::SplitLineToSegments(double s_dis, int source,int target, OGRLineString* temp_line,  
	std::vector<std::array<int, 2>>	& edge_node_ids, std::vector<string>& node_ids, 
	std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points)
{
	int numpoints = temp_line->getNumPoints();
	int seg_node_id = node_ids.size();
	
	if (numpoints == 2) // two nodes, added to the edge set directly
	{
		std::array<int, 2> temp_edge = { source, target };
		edge_node_ids.push_back(temp_edge);
	}
	else if (numpoints == 3) // three nodes, add one node
	{
		std::array<int, 2> temp_edge_1 = { source, seg_node_id };
		std::array<int, 2> temp_edge_2 = { seg_node_id, target };
		edge_node_ids.push_back(temp_edge_1);
		edge_node_ids.push_back(temp_edge_2);

		OGRPoint point; // middle point
		temp_line->getPoint(1, &point);

		double x1 = temp_line->getX(0);
		double y1 = temp_line->getY(0);
		double x2 = temp_line->getX(1);
		double y2 = temp_line->getY(1);
		double seg_len = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
		double node_dis = s_dis + seg_len;

		// updating the node external id, node distance, and node coords vectors
		std::string nodename = "seg_n_" + std::to_string(seg_node_id);
		node_ids.push_back(nodename);
		node_distances.push_back(node_dis);
		node_org_points.push_back(point);
	}
	else
	{
		int added_segnode_num = numpoints - 2;
		// add first edge with start node
		std::array<int, 2> temp_edge_s = { source, seg_node_id };
		edge_node_ids.push_back(temp_edge_s);

		double node_dis = s_dis;
		for (int i = 1; i < added_segnode_num; i++)
		{
			std::array<int, 2> temp_edge_s = { seg_node_id + i - 1, seg_node_id + i };
			edge_node_ids.push_back(temp_edge_s);

			OGRPoint point_i_before, point_i;
			temp_line->getPoint(i - 1, &point_i_before);
			temp_line->getPoint(i, &point_i);

			double x1 = temp_line->getX(i - 1);
			double y1 = temp_line->getY(i - 1);
			double x2 = temp_line->getX(i);
			double y2 = temp_line->getY(i);

			double seg_len = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
			node_dis = node_dis + seg_len;

			std::string nodename = "segnode" + std::to_string(seg_node_id + i - 1);
			node_ids.push_back(nodename);
			node_distances.push_back(node_dis);
			node_org_points.push_back(point_i);
		}
			
		double x1 = temp_line->getX(added_segnode_num);
		double y1 = temp_line->getY(added_segnode_num);
		double x2 = temp_line->getX(added_segnode_num - 1);
		double y2 = temp_line->getY(added_segnode_num - 1);

		double seg_len = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
		node_dis = node_dis + seg_len;

		OGRPoint point_i;
		temp_line->getPoint(added_segnode_num, &point_i);

		std::string nodename = "segnode" + std::to_string(seg_node_id + added_segnode_num - 1);
		node_ids.push_back(nodename);
		node_distances.push_back(node_dis);
		node_org_points.push_back(point_i);

		// add last edge with end node
		std::array<int, 2> temp_edge_e = { seg_node_id + added_segnode_num - 1, target};
		edge_node_ids.push_back(temp_edge_e);
	}
}


void DirectedGraph::write_for_test(std::string filename)
{

	std::vector<OGRLineString>   ex_edges;
	std::vector<int>			 ex_snodes;
	std::vector<int>			 ex_enodes;
	
	edge_iterator_drct ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(directed_g_); ei != ei_end; ++ei)
	{
		vertex_descriptor_drct  vt_s_node = boost::source(*ei, directed_g_);
		vertex_descriptor_drct  vt_e_node = boost::target(*ei, directed_g_);
		auto search = edges_externalid_linestring.find(directed_g_[*ei].id);
		OGRLineString* pLine = search->second;
		ex_edges.push_back(*pLine);
		ex_snodes.push_back(int(vt_s_node));
		ex_enodes.push_back(int(vt_e_node));
	}
	
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("ex_edges", NULL, wkbLineString, NULL);

	// define and create filed
	OGRFieldDefn oFie_faci_Id("id", OFTInteger);
	oFie_faci_Id.SetWidth(10);
	poLayer->CreateField(&oFie_faci_Id);

	OGRFieldDefn oField_s_cost("s_nodeid", OFTInteger);
	oField_s_cost.SetWidth(30);
	poLayer->CreateField(&oField_s_cost);

	OGRFieldDefn oField_t_cost("e_nodeid", OFTInteger);
	oField_t_cost.SetWidth(30);
	poLayer->CreateField(&oField_t_cost);

	int edge_num = ex_edges.size();
	for (int i = 0; i < edge_num; i++) 
	{
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		poFeature->SetField("id", i);
		poFeature->SetField("s_nodeid", ex_snodes[i]);
		poFeature->SetField("e_nodeid", ex_enodes[i]);
		poFeature->SetGeometry(&ex_edges[i]);

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}

	GDALClose(shpDataSet);
}

void DirectedGraph::write_accessible_edges(std::string filename, double delta, int facility_id, std::vector<vertex_descriptor_drct> ex_nodes)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("acc_edges", NULL, wkbLineString, NULL);

	// define and create filed
	OGRFieldDefn oFie_faci_Id("fac_id", OFTInteger);
	oFie_faci_Id.SetWidth(10);
	poLayer->CreateField(&oFie_faci_Id);

	OGRFieldDefn oField_s_cost("s_cost", OFTReal);
	oField_s_cost.SetWidth(30);
	poLayer->CreateField(&oField_s_cost);

	OGRFieldDefn oField_t_cost("t_cost", OFTReal);
	oField_t_cost.SetWidth(30);
	poLayer->CreateField(&oField_t_cost);

	edge_iterator_drct ei, ei_end;

	for (boost::tie(ei, ei_end) = boost::edges(directed_g_); ei != ei_end; ++ei)
	{
		vertex_descriptor_drct  vt_s_node = boost::source(*ei, directed_g_);
		vertex_descriptor_drct  vt_e_node = boost::target(*ei, directed_g_);
		double s_dis = distances_map_[vt_s_node];
		double e_dis = distances_map_[vt_e_node];
		auto search = edges_externalid_linestring.find(directed_g_[*ei].id);
		OGRLineString* pLine = search->second;
		double edge_length = pLine->get_Length();
		int direction = abs(directed_g_[*ei].direction);

		OGRPoint *breakpoint = new OGRPoint();
		OGRLineString *s_cut_line = new OGRLineString();
		OGRLineString *e_cut_line = new OGRLineString();
		
		if ((s_dis != std::numeric_limits<double>::max()) && (e_dis != std::numeric_limits<double>::max()))
		{
			if (direction == 1)
			{
				if ((s_dis < delta) && (delta < e_dis))
				{
					OGRFeature *poFeature;
					double break_to_snode = delta - s_dis;
					algorithm::CutLinestringBasedOnBreakDistance(pLine, break_to_snode, breakpoint, s_cut_line, e_cut_line);
					poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
					poFeature->SetField("fac_id", facility_id);
					poFeature->SetField("s_cost", s_dis);
					poFeature->SetField("t_cost", delta);
					poFeature->SetGeometry(s_cut_line);
					if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
					{
						printf("Failed to create feature in shapefile.\n");
						exit(1);
					}
					OGRFeature::DestroyFeature(poFeature);

				}
				else if ((s_dis <= delta) && (e_dis <= delta))
				{
					OGRFeature *poFeature;
					poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
					poFeature->SetField("fac_id", facility_id);
					poFeature->SetField("s_cost", s_dis);
					poFeature->SetField("t_cost", e_dis);
					poFeature->SetGeometry(pLine);
					if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
					{
						printf("Failed to create feature in shapefile.\n");
						exit(1);
					}
					OGRFeature::DestroyFeature(poFeature);
				}
				else {
					continue;
				}// end if
			}
			else
			{
				if ((s_dis < delta) && (delta < (s_dis + edge_length)))
				{
					OGRFeature *poFeature;
					double break_to_snode = delta - s_dis;

					algorithm::CutLinestringBasedOnBreakDistance(pLine, break_to_snode, breakpoint, s_cut_line, e_cut_line);
					poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
					poFeature->SetField("fac_id", facility_id);
					poFeature->SetField("s_cost", s_dis);
					poFeature->SetField("t_cost", delta);
					poFeature->SetGeometry(s_cut_line);
					if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
					{
						printf("Failed to create feature in shapefile.\n");
						exit(1);
					}
					OGRFeature::DestroyFeature(poFeature);
				}
				else if ((s_dis < delta) && ((s_dis + edge_length) <= delta))
				{
					OGRFeature *poFeature;
					poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
					poFeature->SetField("fac_id", facility_id);
					poFeature->SetField("s_cost", s_dis);
					poFeature->SetField("t_cost", e_dis);
					poFeature->SetGeometry(pLine);
					if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
					{
						printf("Failed to create feature in shapefile.\n");
						exit(1);
					}
					OGRFeature::DestroyFeature(poFeature);
				}
				else 
				{

					continue;
				}// end if
			}// end if
		}// end if
	}// end for 
	GDALClose(shpDataSet);
};



// wirte out all the nodes with distance >= delta;
void DirectedGraph::wirte_sp_nodes_plus_delatnodes(const std::string &filename, const std::vector<vertex_descriptor_drct> examined_nodes,
	const std::vector<OGRPoint> deltanodes, double delta)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("En_sp_tree_with_delta", NULL, wkbPoint, NULL);

	// define and create filed
	OGRFieldDefn oFie_faci_Id("nodeid", OFTString);
	oFie_faci_Id.SetWidth(32);
	poLayer->CreateField(&oFie_faci_Id);

	OGRFieldDefn oField_distance("distance", OFTReal);
	oField_distance.SetWidth(30);
	poLayer->CreateField(&oField_distance);

	int number_record = examined_nodes.size();

	vertex_descriptor_drct u;
	for (int i = 0; i < number_record; i++)
	{
		u = examined_nodes[i];
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		poFeature->SetField("nodeid", vertex_externalids_[u].c_str());
		poFeature->SetField("distance", distances_map_[u]);
		poFeature->SetGeometry(&(vertex_coords_[u]));

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}

	// add the delta nodes at the end of paper
	for (int j = 0; j < deltanodes.size(); j++)
	{
		std::string s_delta = "delta_";
		std::string nodeid = s_delta + to_string(j);
		OGRFeature *poFeature;

		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		poFeature->SetField("nodeid", nodeid.c_str());
		poFeature->SetField("distance", delta);

		poFeature->SetGeometry(&(deltanodes[j]));

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}

	GDALClose(shpDataSet);
}//wirte_sp_nodes_plus_delatnode


 //  Clean the distance map and predecessor map
void DirectedGraph::InitializeDistancesPredecessors()
{
	// Need initialization
	int num_vertices = boost::num_vertices(directed_g_);

	predecessors_map_ = std::vector<vertex_descriptor_drct>(num_vertices);
	distances_map_ = std::vector<double>(num_vertices);

	for (int i = 0; i < num_vertices; ++i) 
	{
		distances_map_[i] = std::numeric_limits<double>::max();
		predecessors_map_[i] = i;
	}
}


void DirectedGraph::InitializeClosestFacilities(int vertexnumber)
{
	// Need initialization
	closest_facilities_map_ = std::vector<vertex_descriptor_drct>(vertexnumber);
	for (int i = 0; i < vertexnumber; ++i) 
	{
		closest_facilities_map_[i] = i;
	}
}


}// end of namspace