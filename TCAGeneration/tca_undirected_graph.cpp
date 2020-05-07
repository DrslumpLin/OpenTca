#include "tca_undirected_graph.h"


namespace tca
{

UndirectedGraph::UndirectedGraph(std::vector<Edge*> egdes)
{
	int current_idx = -1;
	edge_descriptor e;
	bool inserted;
	int N = egdes.size();
	int source_idx = 0;
	int target_idx = 0;
	
	printf("Network edges :%d \n", N);

	OGRPoint spoint;
	OGRPoint epoint;
	
	for (int i = 0; i < N; ++i)
	{
		Edge *network_edge = egdes[i];
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
			vertexs_externalid_discriptor_.insert({ network_edge->source,current_idx});
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

		boost::tie(e, inserted) = add_edge(source_idx, target_idx, udirected_g_);
	
		udirected_g_[e].id = network_edge->external_id;
		udirected_g_[e].length = network_edge->length;
		edges_externalid_linestring.insert({ network_edge->external_id, network_edge->line_string });

		//below are clone version for sub-edges 
		//OGRLineString *line = new OGRLineString();
		//line = (OGRLineString*)(network_edge->line_string->clone());
		//edges_externalid_linestring.insert({ network_edge->external_id, line});
	}
	std::cout << "Construct graph from network edges end" << '\n';
}

UndirectedGraph::~UndirectedGraph()
{
	std::cout << "Cleaning the undirected graph " << '\n';

	// the operation of delete should be based on reference of parameter, otherwise the 
	// a copy of the parameter will be generated inside the function and delete cannot success
	// the following function still not work correctly, need to be improved
	for (auto &item : edges_externalid_linestring) 
	{
		//delete item.second; // wrong this version, when iterate for edge "fa_egde_s_0", since the OGRGeometry inside the OGRLineString is NULL		
		//if (item.second != NULL)// wrong 
		//{
		//  OGRGeometryFactory::destroyGeometry(item.second);
		//}
	}
	std::cout << "Cleaning undirected graph finished" << '\n';
}



// inserting a node to a specific edge, this is further classified into steps:
// 1) remove a edge; 2) add two edges
void UndirectedGraph::InsertingFacilityNode(const NearestPoint ne_point)
{
	OGRLineString *s_cutedge = new OGRLineString;
	OGRLineString *e_cutedge = new OGRLineString;
	algorithm::CutLinestringBasedOnBreakPoint(ne_point.edge->line_string, ne_point.cut_snode_index, ne_point.point, s_cutedge, e_cutedge);

	edge_descriptor e;
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

	// since the edge id is unique, we thus choose to not update the edge_linestring_map
	// after the following step of "remove_edge"
	remove_edge(temp_s_id, temp_e_id, udirected_g_);

	// insert one node and two edges
	int insert_node_id = boost::num_vertices(udirected_g_);
	
	boost::tie(e, inserted) = add_edge(temp_s_id, insert_node_id, udirected_g_);
	udirected_g_[e].id = "fa_egde_s_" + std::to_string(ne_point.facility_id);
	udirected_g_[e].length = ne_point.pp_2_s;
	edges_externalid_linestring.insert({udirected_g_[e].id, s_cutedge});
	
	boost::tie(e, inserted) = add_edge(insert_node_id, temp_e_id, udirected_g_);
	udirected_g_[e].id = "fa_egde_t_" + std::to_string(ne_point.facility_id);
	udirected_g_[e].length = ne_point.edge->length - ne_point.pp_2_s;

	// Ad doc solution for dealing with the precision issue
	if (udirected_g_[e].length < 0)
	{
		udirected_g_[e].length = 0;
	}

	edges_externalid_linestring.insert({udirected_g_[e].id, e_cutedge}); 
	
	std::string insert_node_sid = "fac" + to_string(ne_point.facility_id);
	vertexs_externalid_discriptor_.insert({ insert_node_sid,insert_node_id });
	vertex_externalids_.push_back(insert_node_sid);
	vertex_coords_.push_back(*(ne_point.point));
};



// Added for multiple-source dijkstra algorithm
// 1) insert all the facility nodes to graph
// 2) Add virtual node and links between virtual node and facility nodes
void UndirectedGraph::InsertingMulltipleFacilityNodes(const NearestPoints ne_points)
{
	// inserting facility nodes
	for (int i = 0; i < ne_points.size(); i++) 
	{
		NearestPoint temp_ne_point = ne_points[i];

		InsertingFacilityNode(temp_ne_point);
	}
	// inserting virtual nodes and add links to virtual nodes	
	OGRPoint *vir_point = new OGRPoint(99999, 99999);	
	OGRLineString *vir_edge = new OGRLineString();
	vir_edge->addPoint(vir_point);
	vir_edge->addPoint(vir_point);

	int virtual_node_id = boost::num_vertices(udirected_g_);
	vertexs_externalid_discriptor_.insert({"vir_node", virtual_node_id});
	vertex_externalids_.push_back("vir_node");
	vertex_coords_.push_back(*(vir_point));

	// Add links between virtual node and facility nodes
	for (int i = 0; i < ne_points.size(); i++)
	{
		edge_descriptor e;
		bool inserted;
		
		//std::string temp_f_node_sid = "fac" + to_string(i);
		std::string temp_f_node_sid = "fac" + to_string(ne_points[i].facility_id);

		int temp_facility_id = 0;
		auto search_f = vertexs_externalid_discriptor_.find(temp_f_node_sid);
		if (search_f != vertexs_externalid_discriptor_.end())
		{
			// A node exists already
			temp_facility_id = search_f->second;
		}
		//updating edges by adding a virtual edge
		boost::tie(e, inserted) = add_edge(virtual_node_id, temp_facility_id, udirected_g_);
		//udirected_g_[e].id = "vir_fa_" + to_string(i);
		udirected_g_[e].id = "vir_fa_" + to_string(ne_points[i].facility_id);
		udirected_g_[e].length = 0;
		edges_externalid_linestring.insert({udirected_g_[e].id, vir_edge});
	}
}


// calculating the shortest path from the virtual nodes
void UndirectedGraph::MultipleSourceShortestPaths(const std::string insert_node_sid)
{
	std::vector<vertex_descriptor> examined_nodes;
	int num_vertices = boost::num_vertices(udirected_g_);
	InitializeDistancesPredecessors();
	InitializeClosestFacilities(num_vertices);

	vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(udirected_g_); vi != vend; ++vi)
	{
		if (vertex_externalids_[*vi] == insert_node_sid)
		{
			MultipleSourceShortestPathAlgorithm(*vi, examined_nodes, closest_facilities_map_);
		}
	}
}

// multiple source dijkstra shortest path algorithm
void UndirectedGraph::MultipleSourceShortestPathAlgorithm(const vertex_descriptor& source, 
	std::vector<vertex_descriptor>& p_examined_nodes,
	std::vector<vertex_descriptor>& closest_facilities)
{
	distances_map_[source] = 0; // here is necessay, or an error will be reported, please see the info 
	double inf = std::numeric_limits<double>::max();

	dijkstra_shortest_paths_no_color_map_no_init
	(
		udirected_g_,
		source,
		make_iterator_property_map(predecessors_map_.begin(), get(boost::vertex_index, udirected_g_), predecessors_map_[0]),
		make_iterator_property_map(distances_map_.begin(), get(boost::vertex_index, udirected_g_), distances_map_[0]),
		get(&EdgeProperty::length, udirected_g_),
		get(boost::vertex_index, udirected_g_),
		std::less<double>(), //DistanceCompare distance_compare,
		boost::closed_plus<double>(inf),
		inf,
		0,
		MS_distance_visitor(source, p_examined_nodes, closest_facilities)
	);
};// end for MS_shortest_path_algorithm 


void UndirectedGraph::write_multiple_sources_node_distances(std::string filename)
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

	std::pair<vertex_iterator, vertex_iterator > vp;
	for (vp = boost::vertices(udirected_g_); vp.first != vp.second; ++vp.first)
	{
		double dis = distances_map_[*vp.first];
		vertex_descriptor closest_faci = closest_facilities_map_[*vp.first];
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
	}// end for
	GDALClose(shpDataSet);
}

void UndirectedGraph::write_multiple_sources_polygons_Segements(std::string filename, vector<OGRLineString> linstrings)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("MS_poly_segs", NULL, wkbLineString, NULL);

	// define and create filed
	OGRFieldDefn oField_edge_label("id", OFTInteger);
	oField_edge_label.SetWidth(10);
	poLayer->CreateField(&oField_edge_label);
	
	for (int i = 0; i < linstrings.size(); i++)
	{
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		poFeature->SetField("id", i);
		poFeature->SetGeometry(&linstrings[i]);
		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}// end for
	GDALClose(shpDataSet);

}


void UndirectedGraph::write_multiple_sources_polygons(std::string filename, vector<OGRPolygon> ploygons, vector<string> labels) 
{

	//
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("MSnodes_labels", NULL, wkbPolygon, NULL);

	// define and create filed

	OGRFieldDefn oField_edge_label("label", OFTString);
	oField_edge_label.SetWidth(30);
	poLayer->CreateField(&oField_edge_label);

	for (int i = 0; i < ploygons.size(); i++) 
	{
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		poFeature->SetField("label", labels[i].c_str());
		poFeature->SetGeometry(&ploygons[i]);
		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}// end for
	GDALClose(shpDataSet);
}


std::vector<NetworkVoronoiEdge> UndirectedGraph::ConstructNetworkVoronoiDiagram(void)
{
	std::vector<NetworkVoronoiEdge> nv_edges;

	edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(udirected_g_); ei != ei_end; ++ei)
	{
		vertex_descriptor  vt_s_node = boost::source(*ei, udirected_g_);
		vertex_descriptor  vt_e_node = boost::target(*ei, udirected_g_);

		vertex_descriptor s_closest_faci = closest_facilities_map_[vt_s_node];
		vertex_descriptor e_closest_faci = closest_facilities_map_[vt_e_node];

		double s_dis = distances_map_[vt_s_node];
		double e_dis = distances_map_[vt_e_node];
		
		auto search = edges_externalid_linestring.find(udirected_g_[*ei].id);
		// here cannot use the = for two pointer directly, this will change the original linestring of an edge
		// under the condition of reverse coordinates
		OGRLineString* pLine = search->second;

		OGRPoint p1, p2;
		p1 = vertex_coords_[vt_s_node];
				
		if (s_closest_faci == e_closest_faci)
		{
			NetworkVoronoiEdge nv_edge;
			nv_edge.facility_id = vertex_externalids_[s_closest_faci];
			nv_edge.origanl_edge_external_id = udirected_g_[*ei].id;
			nv_edge.linestring = pLine;			
			nv_edges.push_back(nv_edge);
		}
		else
		{
			// unconnected edges and virtual edges should be excluded 
			// the != 0 and != max() is not a very good style, this needs change in the future
			if ((s_dis != std::numeric_limits<double>::max())
				&& (e_dis != std::numeric_limits<double>::max())
				&& (udirected_g_[*ei].length != 0.0))
			{
				double break_to_snode = (distances_map_[vt_e_node] - distances_map_[vt_s_node] + udirected_g_[*ei].length) / 2;
				pLine->StartPoint(&p2);
				
				//if (p1 != p2)
				// precision issue, since the topology generated by qgis, have two coords with same node string id
				// but different coords, when add this short edge into the graph, the node coords might different from the 
				// start node of another linestring
				double distance = CalculateDistance(p1, p2);

				if (distance > 0.03) 
				{
					pLine->reversePoints(); // may be its better to new a pline£¬because reverse may change the sequence of the orginal *pointer
				}

				OGRLineString *s_cutedge = new OGRLineString();
				OGRLineString *e_cutedge = new OGRLineString();
				OGRPoint *breakpoint = new OGRPoint();

				algorithm::CutLinestringBasedOnBreakDistance(pLine, break_to_snode, breakpoint, s_cutedge, e_cutedge);

				NetworkVoronoiEdge nv_s_edge, nv_e_edge;
				nv_s_edge.facility_id = vertex_externalids_[s_closest_faci];
				nv_s_edge.origanl_edge_external_id = udirected_g_[*ei].id;
				nv_s_edge.linestring = s_cutedge;
				nv_edges.push_back(nv_s_edge);
				
				nv_e_edge.facility_id = vertex_externalids_[e_closest_faci];
				nv_e_edge.origanl_edge_external_id = udirected_g_[*ei].id;
				nv_e_edge.linestring = e_cutedge;		
				nv_edges.push_back(nv_e_edge);
				
				// Pointer Change mark
				//OGRGeometryFactory::destroyGeometry(breakpoint);
			}
			//else
			//{
			//	std::cout << "distance to s node: " << s_dis << "distance to e node: " << e_dis << std::endl;
			//	std::cout << "edge_length: " << udirected_g_[*ei].length << std::endl;
			//}								
		}// end if
	}// end for 
	return nv_edges;
}


std::vector<NetworkVoronoiNode> UndirectedGraph::ConstructNetworkVoronoiNode(void)
{
	std::vector<NetworkVoronoiNode> nv_nodes;

	edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(udirected_g_); ei != ei_end; ++ei)
	{
		vertex_descriptor  vt_s_node = boost::source(*ei, udirected_g_);
		vertex_descriptor  vt_e_node = boost::target(*ei, udirected_g_);

		vertex_descriptor s_closest_faci = closest_facilities_map_[vt_s_node];
		vertex_descriptor e_closest_faci = closest_facilities_map_[vt_e_node];

		double s_dis = distances_map_[vt_s_node];
		double e_dis = distances_map_[vt_e_node];

		auto search = edges_externalid_linestring.find(udirected_g_[*ei].id);
		OGRLineString* pLine = search->second;

		OGRPoint p_s, p_e;
		OGRPoint line_s_point;
		p_s = vertex_coords_[vt_s_node];
		p_e = vertex_coords_[vt_e_node];
		// unconnected edges and virtual edges should be excluded 
        // the != 0 and != max() is not a very good style, this needs change in the future
		if ((s_dis != std::numeric_limits<double>::max())
			&& (e_dis != std::numeric_limits<double>::max())
			&& (udirected_g_[*ei].length != 0.0))
		{
		// add the start and end nodes respectively
			NetworkVoronoiNode nv_node_s, nv_node_e;
			nv_node_s.facility_id_1 = vertex_externalids_[s_closest_faci];
			nv_node_s.inner_node_id = int(vt_s_node);
			nv_node_s.point = p_s;

			nv_node_e.facility_id_1 = vertex_externalids_[e_closest_faci];
			nv_node_e.inner_node_id = int(vt_e_node);
			nv_node_e.point = p_e;

			nv_nodes.push_back(nv_node_s);
			nv_nodes.push_back(nv_node_e);

			if (s_closest_faci != e_closest_faci) // the two ends nodes with different nodes
			{
				// make sure that the retrived linestring is the same direction from (Snode to Nnode)
				double break_to_snode = (distances_map_[vt_e_node] - distances_map_[vt_s_node] + udirected_g_[*ei].length) / 2;
				pLine->StartPoint(&line_s_point);
				double distance = CalculateDistance(p_s, line_s_point); 
				
				OGRLineString *s_cutedge = new OGRLineString();
				OGRLineString *e_cutedge = new OGRLineString();;
				OGRPoint *breakpoint = new OGRPoint();

				if (distance > 0.03) // some points with 
				{
					OGRLineString* pLine_re = new OGRLineString();// Pointer Change mark
					*pLine_re = *(pLine);// Pointer Change mark
					pLine_re->reversePoints(); // Pointer Change mark
					// change the direction of the pline
					algorithm::CutLinestringBasedOnBreakDistance(pLine_re, break_to_snode, breakpoint, s_cutedge, e_cutedge);
				}
				else
				{
					algorithm::CutLinestringBasedOnBreakDistance(pLine, break_to_snode, breakpoint, s_cutedge, e_cutedge);
				}
				NetworkVoronoiNode nv_break_node;
				nv_break_node.facility_id_1 = vertex_externalids_[s_closest_faci];
				nv_break_node.facility_id_2 = vertex_externalids_[e_closest_faci];
				nv_break_node.is_break_node = true;
				nv_break_node.point = *breakpoint;
				GetCounterClockwiseUnitNormalVector(p_s.getX(), p_s.getY(), p_e.getX(), p_e.getY(), nv_break_node.direction_x, nv_break_node.direction_y);
				nv_nodes.push_back(nv_break_node);			
			}
		}
	}
	return nv_nodes;
}


void UndirectedGraph::write_nv_nodes(std::string filename, std::vector<NetworkVoronoiNode> nv_nodes)
{
	std::map<PointIndex, int> points_ids; // used for speed index of nv_nodes, thus to judge if a nv_node is a break node
	std::vector<NetworkVoronoiNode> clean_nodes;

	for (int i = 0; i < nv_nodes.size(); i++)
	{
		NetworkVoronoiNode node = nv_nodes[i];		
		auto result = points_ids.insert({ PointIndex(Point_2(node.point.getX(), node.point.getY())),i });
		if (result.second) clean_nodes.push_back(node);
	}

	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("network_vornoi", NULL, wkbPoint, NULL);

	// define and create filed
	OGRFieldDefn oField_node_id("innerid", OFTInteger);
	oField_node_id.SetWidth(10);
	poLayer->CreateField(&oField_node_id);
	
	OGRFieldDefn oField_fac_s("fac_s", OFTString);
	oField_fac_s.SetWidth(30);
	poLayer->CreateField(&oField_fac_s);
	
	OGRFieldDefn oField_fac_e("fac_e", OFTString);
	oField_fac_e.SetWidth(30);
	poLayer->CreateField(&oField_fac_e);

	OGRFieldDefn oField_isbreak("isBreak", OFTInteger);
	oField_isbreak.SetWidth(30);
	poLayer->CreateField(&oField_isbreak);

	OGRFieldDefn oField_dirX("dirX", OFTReal);
	oField_dirX.SetWidth(30);
	poLayer->CreateField(&oField_dirX);

	OGRFieldDefn oField_dirY("dirY", OFTReal);
	oField_dirY.SetWidth(30);
	poLayer->CreateField(&oField_dirY);

	int node_nums = clean_nodes.size();
	for (int i = 0; i < node_nums; i++)
	{
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());

		poFeature->SetField("innerid", clean_nodes[i].inner_node_id);

		poFeature->SetField("fac_s", clean_nodes[i].facility_id_1.c_str());
		poFeature->SetField("fac_e", clean_nodes[i].facility_id_2.c_str());
		poFeature->SetField("isBreak", int(clean_nodes[i].is_break_node));
		poFeature->SetField("dirX", clean_nodes[i].direction_x);
		poFeature->SetField("dirY", clean_nodes[i].direction_y);
		poFeature->SetGeometry(&clean_nodes[i].point);

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(shpDataSet);
}

void UndirectedGraph::GenerateAddtionalNvPoints(Iso_rectangle_2 bbox, vector<double> levels, vector<Point_2>& add_points)
{
	/*
			7------6--------5
			|    ____       |
			|	|	 |      |
			8   |____|      4
			|               |
			|               |
			|1----2--------3
	*/
	double xmin = bbox.xmin();
	double xmax = bbox.xmax();
	double ymin = bbox.ymin();
	double ymax = bbox.ymax();
	double dx = xmax - xmin;
	double dy = ymax - xmin;


	int level = levels[0];
	Point_2 p1((xmin - level * dx), (ymin - level * dy));
	Point_2 p2((xmin + xmax) / 2, (ymin - level * dy));
	Point_2 p3((xmax + level * dx), (ymin - level * dy));
	Point_2 p4((xmax + level * dx), (ymax + ymin) / 2);
	Point_2 p5((xmax + level * dx), (ymax + level * dy));
	Point_2 p6((xmin + xmax) / 2, (ymax + level * dy));
	Point_2 p7((xmin - level * dx), (ymax + level * dy));
	Point_2 p8((xmin - level * dx), (ymax + ymin) / 2);
	add_points.push_back(p1);
	add_points.push_back(p2);
	add_points.push_back(p3);
	add_points.push_back(p4);
	add_points.push_back(p5);
	add_points.push_back(p6);
	add_points.push_back(p7);
	add_points.push_back(p8);

	level = levels[1];
	p1 = Point_2 ((xmin - level * dx), (ymin - level * dy));
	p2 = Point_2 ((xmin + xmax)/2.1, (ymin - level * dy));
	p3 = Point_2 ((xmax + level * dx), (ymin - level * dy));
	p4 = Point_2 ((xmax + level * dx), (ymax + ymin) / 2.1);
	p5 = Point_2 ((xmax + level * dx), (ymax + level * dy));
	p6 = Point_2 ((xmin + xmax) / 2.1, (ymax + level * dy));
	p7 = Point_2 ((xmin - level * dx), (ymax + level * dy));
	p8 = Point_2 ((xmin - level * dx), (ymax + ymin) / 2.1);
	add_points.push_back(p1);
	add_points.push_back(p2);
	add_points.push_back(p3);
	add_points.push_back(p4);
	add_points.push_back(p5);
	add_points.push_back(p6);
	add_points.push_back(p7);
	add_points.push_back(p8);

}

// the code is adapted  from https://stackoverflow.com/questions/26364212/cgal-voronoi-diagram-link-input-sites-to-faces
// this version is supper slow (the part of reconstruction trigulations), need to be improved here
void UndirectedGraph::BuildCroppedVoronoiPolygon(const std::vector<NetworkVoronoiNode> nv_nodes, vector<OGRPolygon>& ploygons, vector<string>& labels)
{
	VD vd;
	vector<Point_2> points;
	std::map<PointIndex, int> points_ids; // used for speed index of nv_nodes, thus to judge if a nv_node is a break node	
	for (int i = 0; i < nv_nodes.size(); i++)
	{
		NetworkVoronoiNode node = nv_nodes[i];
		auto result = points_ids.insert({ PointIndex(Point_2(node.point.getX(), node.point.getY())),i });
		// some duplicated nodes, beacause we get the nv_nodes by iterate over the graph edges	
		if (result.second) points.push_back(Point_2(node.point.getX(), node.point.getY()));
	}

	// a bounding box of points to crop the open (boundary) voronoi polygons,
	const K::Iso_rectangle_2 bbox = CGAL::bounding_box(points.begin(), points.end());
	CGAL::Polygon_2<K> bpoly;
	bpoly.push_back(K::Point_2(bbox.xmin(), bbox.ymin()));
	bpoly.push_back(K::Point_2(bbox.xmax(), bbox.ymin()));
	bpoly.push_back(K::Point_2(bbox.xmax(), bbox.ymax()));
	bpoly.push_back(K::Point_2(bbox.xmin(), bbox.ymax()));
	
	vector<double> levels = {1,2};
	vector<Point_2> add_points;
	GenerateAddtionalNvPoints(bbox, levels, add_points);	
	for (int i = 0; i < add_points.size(); i++) points.push_back(add_points[i]);
	
	// add points to generate voronoi diagram
	vd.insert(points.begin(), points.end());
	

	int fac_num = 0;
	for (VD::Face_iterator fit = vd.faces_begin(); fit != vd.faces_end(); ++fit)
	{
		if (!fit->is_unbounded()) 
		{
			if (fac_num % 1000 == 0) std::cout << "finishing processing the fac" << fac_num << endl;
			Polygon_2 pgon;
			// the dual point here is the center node of 
			Point_2 dual_point = fit->dual()->point();
			//auto search = points_ids.find(PointIndex(dual_point.x(), dual_point.y()));
			auto search = points_ids.find(PointIndex(dual_point));
			int inner_id = -1;
		
			// only output the polygons corresponds to the original points
			if (search != points_ids.end()) 
			{
				inner_id = search->second;
				
				//----------------below is used to find the Cropped polygon------------------
				//Edge circulators traverse endlessly around a face. Make a note of the
				//starting point so we know when to quit.
				Ccb_halfedge_circulator ec_start = fit->ccb();
				Ccb_halfedge_circulator ec = ec_start;				
				do
				{	
					pgon.push_back(ec->source()->point());
				} while (++ec != ec_start); //Loop until we get back to the beginning
				
				OGRPolygon ogr_bpoly;
				OGRPolygon ogr_polygon;
				
				CGALToOGRPolygon(bpoly, ogr_bpoly);   //the OGR representation of the bbox (the outer ring of ogr polygons are clockwise)
				
				CGALToOGRPolygon(pgon, ogr_polygon);  // the OGR representation of the voronoi polygon
				
				// using the intersection methods to keep the parts inside the bbox of each voronoi polygon
				OGRPolygon* pTemp = static_cast<OGRPolygon*>(ogr_bpoly.Intersection(&ogr_polygon));
				ogr_polygon = *pTemp;
				
				// if a voronoi polygon is belongs to a break point 
				// then it needs to be split to two polygons,each sub polygon correspond
				// to a facility
				if (!nv_nodes[inner_id].is_break_node) // not a break node
				{
					ploygons.push_back(ogr_polygon);
					labels.push_back(nv_nodes[inner_id].facility_id_1);
				}
				else// a break node, then need to split the polygons 
				{	
					OGRPolygon ogr_polygon_1, ogr_polygon_2;
					Polygon_2 cgal_poly;
					OgrPolygonToCgalPolygon(ogr_polygon, cgal_poly);
					SplitVoronoiPoygon(cgal_poly, nv_nodes[inner_id], ogr_polygon_1, ogr_polygon_2);
					ploygons.push_back(ogr_polygon_1);
					labels.push_back(nv_nodes[inner_id].facility_id_1);
					ploygons.push_back(ogr_polygon_2);
					labels.push_back(nv_nodes[inner_id].facility_id_2);
				}
			}
		fac_num = fac_num + 1;
		}
	}
}


// this one can successfully generate the segments
// one possbile post-processing for this is to use virtual OGRGeometry * Polygonize () const CPL_WARN_UNUSED_RESULT
// https://www.gdal.org/classOGRGeometry.html#aadb34b556c52aef3e93f03cf65d3f4bc
void UndirectedGraph::BuildCroppedVoronoiPolySegs(const std::vector<NetworkVoronoiNode> nv_nodes, vector<OGRLineString>& linstrings)
{
	DT dt2;
	vector<Point_2> points;
	std::map<PointIndex, int> points_ids; // used for speed index of nv_nodes, thus to judge if a nv_node is a break node	
	for (int i = 0; i < nv_nodes.size(); i++)
	{
		NetworkVoronoiNode node = nv_nodes[i];
		auto result = points_ids.insert({ PointIndex(Point_2(node.point.getX(), node.point.getY())),i });
		// some duplicated nodes, beacause we get the nv_nodes by iterate over the graph edges	
		if (result.second) points.push_back(Point_2(node.point.getX(), node.point.getY()));
	}
	//
	dt2.insert(points.begin(), points.end());

	// a bounding box of points to crop the voronoi polygons
	const K::Iso_rectangle_2 bbox = CGAL::bounding_box(points.begin(), points.end());

	Cropped_voronoi_from_delaunay vor(bbox);
	//extract the cropped Voronoi diagram
	dt2.draw_dual(vor);
	
	for (auto it = vor.m_cropped_vd.begin(); it != vor.m_cropped_vd.end(); ++it) 
	{
		Segment_2 temp = *it;
		Point_2 source = temp.source();
		Point_2 target = temp.target();
		OGRLineString temp_line;
		temp_line.addPoint(source.x(), source.y());
		temp_line.addPoint(target.x(), target.y());
		linstrings.push_back(temp_line);
	}
}

//  tranfer a OGRpolygon to CGAL polygon 
// the interctions  between a voronoi polygon and bbox must be a polygon without inner ring 
void  UndirectedGraph::OgrPolygonToCgalPolygon( OGRPolygon ogr_polygon, Polygon_2&  cgal_poly)
{
	//OGRLineRing;	
	OGRLinearRing *linering = ogr_polygon.getExteriorRing();
	int pnum = linering->getNumPoints();

	// here we start from second last coord, to generate a non-closed clockwise CGAL polygon  
	for (int i = pnum - 1; i > 0; i--)
	{
		cgal_poly.push_back(Point_2(linering->getX(i), linering->getY(i)));
	}
}

//  tranfer a CGAL polygon to OGRpolygon 
void  UndirectedGraph::CGALToOGRPolygon(Polygon_2 bpoly, OGRPolygon & ogr_polygon)
{
	//OGRLineRing;
	OGRLinearRing linering;
	// based on the info here: https://doc.cgal.org/latest/Voronoi_diagram_2/index.html
	// the orinigal vertex is arranged in a contourclockwise direction
	for (auto v = bpoly.vertices_begin(); v != bpoly.vertices_end(); v++) 
		//linering.addPoint(CGAL::to_double(v->x()), CGAL::to_double(v->y()));
		linering.addPoint(v->x(), v->y());
	linering.closeRings();	
	if (!linering.isClockwise()) linering.reverseWindingOrder(); // the or
	ogr_polygon.addRing(&linering);
}

// this one already pass test
// given a line (represented by a point and its slope), and a voronoi polygon
// cut the voronoi  into two different polygons

void  UndirectedGraph::SplitVoronoiPoygon(Polygon_2 poly, NetworkVoronoiNode node, OGRPolygon & ogr_polygon_1, OGRPolygon & ogr_polygon_2)
{
	Polygon_2 pgon_s, pgon_t;
	// get a line first
	double max_distance = 0;
	vector<Point_2> poly_points;
	
	for (polygon_vi vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi) 
	{
		poly_points.push_back(Point_2(vi->x(), vi->y()));
		for (polygon_vi vi_1 = poly.vertices_begin(); vi_1 != poly.vertices_end(); ++vi_1) 
		{
			double x = vi->x();
			double y = vi->y();
			double x1 = vi_1->x();
			double y1 =vi_1->y();
			double distance = std::sqrt((x-x1)*(x - x1) + (y-y1)*(y-y1));			
			if (distance > max_distance) max_distance = distance;
		}
	}
	poly_points.push_back(poly_points.front());//  closed the polygon in counterclockwise direction
	
	// construct a line, given 
	double point_x_counterclockwise = node.point.getX() + node.direction_x * max_distance;
	double point_y_counterclockwise = node.point.getY() + node.direction_y * max_distance;

	double point_x_clockwise = node.point.getX() + (-node.direction_x) * max_distance;
	double point_y_clockwise = node.point.getY() + (-node.direction_y) * max_distance;

	Segment_2 seg_counterclockwise(Point_2(node.point.getX(), node.point.getY()), Point_2(point_x_counterclockwise, point_y_counterclockwise));
	Segment_2 seg_clockwise(Point_2(node.point.getX(), node.point.getY()), Point_2(point_x_clockwise, point_y_clockwise));

	Point_2 point_counterclockwise, point_clockwise;
	int counterclockwise_vi, clockwise_vi;

	// find the labels and points of two diretions intersect points between polygon segments and normal vector
	for (int i= 0; i< poly_points.size()-1; i++)
	{
		Segment_2 temp_seg(poly_points[i], poly_points[i+1]);	//Segment_2 tmp_segs[] = { seg_opp };				
		CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type result_1 = intersection(seg_counterclockwise, temp_seg);
		if (result_1)
		{
			point_counterclockwise = (*boost::get<Point_2 >(&*result_1));
			counterclockwise_vi = i + 1;
		}
		// 
		CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type result_2 = intersection(seg_clockwise, temp_seg);
		if (result_2)
		{
			point_clockwise = (*boost::get<Point_2 >(&*result_2));
			clockwise_vi = i + 1;
		}
	}
	
	// note all the nodes of the polygon are in contourclockwise direction
	if (counterclockwise_vi < clockwise_vi) 
	{
		// constructed ploygon for Snode side (i.e. facility 1)
		pgon_s.push_back(point_counterclockwise);
		for ( int i= counterclockwise_vi; i < clockwise_vi; i++)		
			pgon_s.push_back(poly_points[i]);
		pgon_s.push_back(point_clockwise);

		// constricted polygon for Tnode side (i.e. facility 2)
		pgon_t.push_back(point_clockwise);			
		for (int i = clockwise_vi; i < (counterclockwise_vi + poly_points.size()); i++)
			pgon_t.push_back(poly_points[(i%poly_points.size())]);
		pgon_t.push_back(point_counterclockwise);
	}
	else
	{
		// construct polygon for Snode side (i.e. facility 2)
		pgon_s.push_back(point_counterclockwise);
		for (int i = counterclockwise_vi; i < (clockwise_vi + poly_points.size()); i++)
			pgon_s.push_back(poly_points[(i%poly_points.size())]);
		pgon_s.push_back(point_clockwise);
		
		// constructed ploygon for Tnode side (i.e. facility 1)
		pgon_t.push_back(point_clockwise);
		for (int i = clockwise_vi; i < counterclockwise_vi; i++)
			pgon_t.push_back(poly_points[i]);
		pgon_t.push_back(point_counterclockwise);
	}


	CGALToOGRPolygon(pgon_s, ogr_polygon_1);
	CGALToOGRPolygon(pgon_t, ogr_polygon_2);
}


// calculate the unit normal vector of a segment
void UndirectedGraph::GetCounterClockwiseUnitNormalVector(double x1, double y1, double x2, double y2, double & vertical_x, double & vertical_y)
{
	vertical_x = -(y2-y1)/ std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
	vertical_y = (x2-x1)/ std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}


void UndirectedGraph::write_network_voronois_edges(const std::string &filename, std::vector<NetworkVoronoiEdge> vd_edges)
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


std::vector<vertex_descriptor> UndirectedGraph::GetShortestPathNodes(const std::string facility_node_sid)
{
	std::vector<vertex_descriptor> examined_nodes;
	int num_vertices = boost::num_vertices(udirected_g_);
	// for test
	//std:string test_node = "39069";
	InitializeDistancesPredecessors();
	vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(udirected_g_); vi != vend; ++vi)
	{
		if (vertex_externalids_[*vi] == facility_node_sid)
		{
			ShortestPathAlgorithm(*vi, examined_nodes);
		}
	}
	return examined_nodes;
	//wirte_sp_nodes(filename, examined_nodes);
}


std::vector<vertex_descriptor> UndirectedGraph::BuildExtendedShortestPathTree(std::vector<vertex_descriptor> examined_nodes, const std::string facility_node_sid)
{	
	std::vector<ExtendedNode> extended_nodes;
	extended_nodes = FindNonShortestPathEdges(examined_nodes);
	InsertingNonShortestPathEdges(extended_nodes);

	std::vector<vertex_descriptor> examined_nodes_extended_sp;
	int num_vertices = boost::num_vertices(udirected_g_);
	InitializeDistancesPredecessors();
	
	vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(udirected_g_); vi != vend; ++vi)
	{
		if (vertex_externalids_[*vi] == facility_node_sid)
		{
			ShortestPathAlgorithm(*vi, examined_nodes_extended_sp);
		}
	}

	return examined_nodes_extended_sp;
	//wirte_sp_nodes(filename, examined_nodes_extended_sp);
}


std::vector<ExtendedNode> UndirectedGraph::FindNonShortestPathEdges(std::vector<vertex_descriptor> examined_nodes)
{
	std::vector<ExtendedNode> extend_nodes;
	int extend_node_id = 0;

	// check how many node has been exaimed
	int N_size = examined_nodes.size();
	std::vector<std::string> colorMap(distances_map_.size());

	// label every node as "white"
	for (int i = 0; i < distances_map_.size(); i++)
	{
		colorMap[i] = "white";
	}
	vertex_descriptor u, pu;
	OGRPoint p1, p2;
	// from the last exaimed node to the first node
	for (int i = N_size - 1; i >= 0; i--)
	{
		u = examined_nodes[i];
		p1 = vertex_coords_[u]; // the coords of node "u"

		pu = predecessors_map_[u];
		colorMap[u] = "Gray";

		edge_descriptor e;
		out_edge_iterator out_i, out_end;
		vertex_descriptor vi_target;

		// iterating every outer-edge of the current node "u"
		// if the target node of this outer edge is not the precedeer of node "u" and is not iterated
		for (boost::tie(out_i, out_end) = boost::out_edges(u, udirected_g_); out_i != out_end; ++out_i)
		{
			e = *out_i;
			std::string edge_id = udirected_g_[e].id;
			vi_target = boost::target(e, udirected_g_);
			//std::cout << "stop here for debuging" << vertex_id_vec[vi_target] << std::endl;

			if ((vi_target != pu) && (colorMap[vi_target] == "white"))
			{	
				//std::cout << "--node id: " << i << std::endl;

				auto search = edges_externalid_linestring.find(edge_id); // question here, not matter a node [1,2] or [2, 1] 
				OGRLineString* pLine = search->second;
				double break_to_snode = (distances_map_[vi_target] - distances_map_[u] + udirected_g_[e].length) / 2;
				// this was added for undirected graph
				pLine->StartPoint(&p2);
				double distance = CalculateDistance(p1, p2);

				ExtendedNode temprecord;
				if (distance > 0.03) // only under the condition of reverse, a new linestring need to be created
				//if (p1 != p2) // due to the precision, we can not use this criteria
				{
					// Pointer Change mark
					OGRLineString* pLine_re = new OGRLineString();// Pointer Change mark
					*pLine_re = *(search->second);                // Pointer Change mark
					pLine_re->reversePoints();                    // Pointer Change mark
					temprecord.pLine = pLine_re;                 // Pointer Change mark
				}
				else // Pointer Change mark
				{
					temprecord.pLine = pLine;// Pointer Change mark
				}		
				temprecord.break_distance_to_snode = break_to_snode;
				temprecord.source = u;
				temprecord.target = vi_target;
				temprecord.node_id = extend_node_id;
				extend_nodes.push_back(temprecord);
				extend_node_id++;
			}
		}
	}
	return extend_nodes;
};


void UndirectedGraph::InsertingNonShortestPathEdges(const std::vector<ExtendedNode> extended_records)
{
	ExtendedNode temprecord;
	int N = extended_records.size();

	std::cout << "the number of edges not in the sp tree is : " << N << std::endl;
	std::cout << "Original graph - number of vertexs " << boost::num_vertices(udirected_g_) << std::endl;
	std::cout << "Original graph - number of edges " << boost::num_edges(udirected_g_) << std::endl;

	for (int i = 0; i < N; i++)
	{
		OGRLineString *s_cutedge = new OGRLineString();
		OGRLineString *e_cutedge = new OGRLineString();
		OGRPoint *breakpoint = new OGRPoint();
		temprecord = extended_records[i];
		algorithm::CutLinestringBasedOnBreakDistance(temprecord.pLine, temprecord.break_distance_to_snode, breakpoint, s_cutedge, e_cutedge);
		InsertingExtendedNodes(temprecord.source, temprecord.target, temprecord.node_id, *breakpoint, s_cutedge, e_cutedge);
	}
	std::cout << "Updating graph - number of vertexs " << boost::num_vertices(udirected_g_) << std::endl;
	std::cout << "Updating graph - number of edges " << boost::num_edges(udirected_g_) << std::endl;
}

/*
* Inserting nodes into the extended shortest path tree
* for each edge condcut following steps:
*
* 1) remove a edge
* 2) inserting edge [snode, breakpoint]
* 3) inserting a node
* 4) inserting edge [snode, breakpoint]
* 5) inserting a node
* 6)
* 7)
*/

void UndirectedGraph::InsertingExtendedNodes(vertex_descriptor source, vertex_descriptor target, int extend_node_id,
	OGRPoint breakpoint, OGRLineString* s_cutedge, OGRLineString* e_cutedge)
{
	//if (extend_node_id == 56) 
	//{
	//	std::cout << "stop here for debuging" << std::endl;
	//}

	edge_descriptor e;
	bool inserted;

	// remove the edge
	remove_edge(source, target, udirected_g_);

	// insert two nodes with different id but with same point 
	// to make sure each node are included in the extended SP tree
	int num_vertices = boost::num_vertices(udirected_g_);
	int insert_node_sid = num_vertices;

	boost::tie(e, inserted) = add_edge(source, insert_node_sid, udirected_g_);

	//edge
	std::string tempSedge = "ex_E_s_" + to_string(extend_node_id);
	udirected_g_[e].id = tempSedge;
	udirected_g_[e].length = s_cutedge->get_Length();
	edges_externalid_linestring.insert({tempSedge, s_cutedge });
	//node
	std::string snode_string = "ex_N_s" + to_string(extend_node_id);
	vertexs_externalid_discriptor_.insert({ snode_string,insert_node_sid });
	vertex_externalids_.push_back(snode_string);
	vertex_coords_.push_back(breakpoint);	
	//edge
	int insert_node_eid = num_vertices + 1;
	std::string tempEedge = "ex_E_e_" + to_string(extend_node_id);
	boost::tie(e, inserted) = add_edge(insert_node_eid, target, udirected_g_);
	udirected_g_[e].id = tempEedge;
	udirected_g_[e].length = e_cutedge->get_Length();
	edges_externalid_linestring.insert({ tempEedge, e_cutedge});	
	//node
	std::string enode_string = "ex_N_e" + to_string(extend_node_id);
	vertexs_externalid_discriptor_.insert({ enode_string,insert_node_eid });
	vertex_externalids_.push_back(enode_string);
	vertex_coords_.push_back(breakpoint);
};


// newly added for generating the input for constrianted delaunary trigulations
void UndirectedGraph::GetEdgesForConstrianedTriangulation(const std::string facility_node_sid, const double delta,
	std::vector<std::array<int, 2>>	& edge_node_ids, std::vector<string>& node_ids,
	std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points)
{
	std::vector<vertex_descriptor> sp_exaimed_nodes = GetShortestPathNodes(facility_node_sid);
	std::vector<vertex_descriptor> extended_sp_exaimed_nodes = BuildExtendedShortestPathTree(sp_exaimed_nodes, facility_node_sid);
	std::vector<ExtendedNode>  delta_extended_nodes = FindDeltaEdges(delta, extended_sp_exaimed_nodes);
	InsertDeltaEdgesToGraph(delta_extended_nodes);
	// re-calculating the shortest path tree
	GetShortestPathNodes(facility_node_sid);

	edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(udirected_g_); ei != ei_end; ++ei)
	{
		auto search = edges_externalid_linestring.find(udirected_g_[*ei].id);
		OGRLineString* pLine = search->second;
		vertex_descriptor  vt_s_node = boost::source(*ei, udirected_g_);
		vertex_descriptor  vt_e_node = boost::target(*ei, udirected_g_);
		double s_dis = distances_map_[vt_s_node];
		double e_dis = distances_map_[vt_e_node];

		// un-connnected edges should be excluded
		if ((s_dis != std::numeric_limits<double>::max()) && (e_dis != std::numeric_limits<double>::max()))
		{
			std::array<int, 2> temp_edge = { int(vt_s_node), int(vt_e_node) };
			edge_node_ids.push_back(temp_edge);
		}
	}

	node_ids = vertex_externalids_;
	node_distances = distances_map_;
	node_org_points = vertex_coords_;
}


// here all the nodes in sub-gragh + intermidiate points are included
// note: some nodes that not in extended sp tree are included,i.e.,disatance[u] = DOUBLE_MAX are included
// but will be filtered at the constrained trianglation construction period
void UndirectedGraph::GetEdgesForConstrianedTriangulationLineSegmentation(std::vector<vertex_descriptor> extended_sp_exaimed_nodes,
	std::vector<std::array<int, 2>>	& edge_node_ids, std::vector<string>& node_ids,
	std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points)
{
	//std::string fextend_sp_edges = "C://Users//dlint//source//repos//TCA_nets//data_input//undirected//select//extend_shortest_sp_trees//"+ facility_node_sid + ".shp";
	//write_sp_edges(fextend_sp_edges, sp_exaimed_nodes);
	//std::string sp_nodes = "C://Users//dlint//source//repos//TCA_nets//data_input//PHD_Extra_Experiments//test//sp_node_1.shp";
	//std::string fextend_sp_edges = "C://Users//dlint//source//repos//TCA_nets//data_input//PHD_Extra_Experiments//test//sp_edges_label_1_" + to_string(delta) +".shp";
	//write_sp_edges(fextend_sp_edges, extended_sp_exaimed_nodes);
	
	node_ids = vertex_externalids_;
	node_distances = distances_map_;
	node_org_points = vertex_coords_;
	int seg_node_id = vertex_coords_.size(); // started from the last node of orginal node

	vertex_descriptor u, pre_u;

	for (int i = 1; i < extended_sp_exaimed_nodes.size(); i++)
	{
		OGRLineString *temp_line = new OGRLineString();
		u = extended_sp_exaimed_nodes[i];
		pre_u = predecessors_map_[u];
		*temp_line = GetEdgeLinestring(pre_u, u);
		int numpoints = temp_line->getNumPoints();

		if (numpoints == 2) // two nodes, added to the edge set directly
		{
			std::array<int, 2> temp_edge = { int(pre_u), int(u) };
			edge_node_ids.push_back(temp_edge);
		}
		else if (numpoints == 3) // three nodes, add one node
		{
			std::array<int, 2> temp_edge_1 = { int(pre_u), seg_node_id };
			std::array<int, 2> temp_edge_2 = { seg_node_id, int(u) };
			edge_node_ids.push_back(temp_edge_1);
			edge_node_ids.push_back(temp_edge_2);

			OGRPoint point; // middle point
			temp_line->getPoint(1, &point);

			double x1 = temp_line->getX(0);
			double y1 = temp_line->getY(0);
			double x2 = temp_line->getX(1);
			double y2 = temp_line->getY(1);
			double seg_len = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
			double node_dis = distances_map_[pre_u] + seg_len;

			// updating the node external id, node distance, and node coords vectors
			std::string nodename = "seg_n_" + std::to_string(seg_node_id);
			node_ids.push_back(nodename);
			node_distances.push_back(node_dis);
			node_org_points.push_back(point);

			// 
			seg_node_id = seg_node_id + 1;
		}
		else // more than three nodes
		{
			int added_segnode_num = numpoints - 2;

			// add first edge with start node
			std::array<int, 2> temp_edge_s = { int(pre_u), seg_node_id };
			edge_node_ids.push_back(temp_edge_s);

			double node_dis = distances_map_[pre_u];

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
			double x2 = temp_line->getX(added_segnode_num-1);
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
			std::array<int, 2> temp_edge_e = { seg_node_id + added_segnode_num - 1, int(u) };
			edge_node_ids.push_back(temp_edge_e);
			seg_node_id = seg_node_id + added_segnode_num;
		}
	}
}

std::vector<ExtendedNode> UndirectedGraph::FindDeltaEdges(double delta, std::vector<vertex_descriptor> examined_nodes)
{
	std::vector<ExtendedNode> extend_deltanodes;

	vertex_descriptor u, pu;
	
	OGRPoint breakpoint;
	int delta_node_id = 0;

	for (int i = 1; i < examined_nodes.size(); i++)
	{
		u = examined_nodes[i];
		pu = predecessors_map_[u];
		if (distances_map_[u] > delta)
		{
			if (distances_map_[pu] < delta)
			{
				OGRLineString *temp_line = new OGRLineString();
				double break_to_snode = delta - distances_map_[pu];
				*temp_line = GetEdgeLinestring(pu, u); // need to be changed of pointer, there is no need to new a linestring here
				
				ExtendedNode temprecord;
				temprecord.pLine = temp_line;
				temprecord.break_distance_to_snode = break_to_snode;
				temprecord.source = pu;
				temprecord.target = u;
				temprecord.node_id = delta_node_id;
				extend_deltanodes.push_back(temprecord);
				delta_node_id++;
			}
		}
	}
	return extend_deltanodes;
}

void UndirectedGraph::InsertDeltaEdgesToGraph(const std::vector<ExtendedNode> extended_records)
{
	ExtendedNode temprecord;
	int N = extended_records.size();

	std::cout << "the number of edges not in the sp tree is : " << N << std::endl;
	std::cout << "Original extended graph - number of vertexs " << boost::num_vertices(udirected_g_) << std::endl;
	std::cout << "Original extended graph - number of edges " << boost::num_edges(udirected_g_) << std::endl;

	for (int i = 0; i < N; i++)
	{
		OGRLineString *s_cutedge = new OGRLineString();
		OGRLineString *e_cutedge = new OGRLineString();
		OGRPoint *breakpoint = new OGRPoint();
		temprecord = extended_records[i];
		algorithm::CutLinestringBasedOnBreakDistance(temprecord.pLine, temprecord.break_distance_to_snode, breakpoint, s_cutedge, e_cutedge);
		InsertingDeltaNodes(temprecord.source, temprecord.target, temprecord.node_id, *breakpoint, s_cutedge, e_cutedge);
	}
	std::cout << "after adding delta egde - graph - number of vertexs " << boost::num_vertices(udirected_g_) << std::endl;
	std::cout << "after adding delta egde - graph - number of edges " << boost::num_edges(udirected_g_) << std::endl;

}

void UndirectedGraph::InsertingDeltaNodes(vertex_descriptor source, vertex_descriptor target, int extend_node_id,
	OGRPoint breakpoint, OGRLineString* s_cutedge, OGRLineString* e_cutedge)
{
	edge_descriptor e;
	bool inserted;

	// remove the edge
	remove_edge(source, target, udirected_g_);

	// different from InsertingExtendedNodes, here we only need to insert one delta node
	// insert two edges
	int num_vertices = boost::num_vertices(udirected_g_);
	int insert_node_sid = num_vertices;

	// inserting two edges
	boost::tie(e, inserted) = add_edge(source, insert_node_sid, udirected_g_);
	std::string tempSedge = "delta_E_s_" + to_string(extend_node_id);
	udirected_g_[e].id = tempSedge;
	udirected_g_[e].length = s_cutedge->get_Length();
	edges_externalid_linestring.insert({ tempSedge, s_cutedge });

	boost::tie(e, inserted) = add_edge(insert_node_sid, target, udirected_g_);
	std::string tempEedge = "delta_E_e_" + to_string(extend_node_id);	
	udirected_g_[e].id = tempEedge;
	udirected_g_[e].length = e_cutedge->get_Length();
	edges_externalid_linestring.insert({ tempEedge, e_cutedge });
	
	// inserting one node (delta)
	std::string node_string = "delta_N_" + to_string(extend_node_id);
	vertexs_externalid_discriptor_.insert({ node_string,insert_node_sid });
	vertex_externalids_.push_back(node_string);
	vertex_coords_.push_back(breakpoint);
};



// based on vertex of descriptor of source and target node to find the edge linestring
// since the extended_nodes are already inserted, thereby no need to compare the edge length
OGRLineString UndirectedGraph::GetEdgeLinestring(vertex_descriptor source, vertex_descriptor target)
{
	OGRLineString templinestring;
	edge_descriptor e;
	out_edge_iterator out_i, out_end;
	OGRPoint s_point, e_point;
	
	for (boost::tie(out_i, out_end) = boost::out_edges(source, udirected_g_); out_i != out_end; ++out_i)
	{
		s_point = vertex_coords_[source];
		e = *out_i;

		if (target == boost::target(e, udirected_g_))
		{
			auto search = edges_externalid_linestring.find(udirected_g_[e].id); 
			templinestring = *(search->second);
			templinestring.StartPoint(&e_point);
			double distance = CalculateDistance(s_point, e_point);
			if (distance > 0.03)
			{
				templinestring.reversePoints();                       // for undirected graph, here is very important
			}
		}
	}
	return templinestring;
}


double UndirectedGraph::CalculateDistance(OGRPoint point1, OGRPoint point2) 
{
	double distance = std::sqrt((point1.getX() - point2.getX())*(point1.getX() - point2.getX()) + 
		(point1.getY() - point2.getY())*(point1.getY() - point2.getY()));

	return distance;
}


std::vector<OGRPoint> UndirectedGraph::FindingDeltaNodes(double delta, std::vector<vertex_descriptor> ex_nodes)
{
	std::vector<OGRPoint> deltanodes;
	vertex_descriptor u, pu;
	OGRLineString temp_line;
	OGRPoint breakpoint;
	OGRLineString ptemplinestring;

	for (int i = 1; i < ex_nodes.size(); i++)
	{
		u = ex_nodes[i];
		pu = predecessors_map_[u];
		if (distances_map_[u] > delta)
		{
			if (distances_map_[pu] < delta)
			{
				double break_to_snode = delta - distances_map_[pu];
				ptemplinestring = GetEdgeLinestring(pu, u);
				breakpoint = algorithm::GetBreakPointOnLinestring(break_to_snode, &ptemplinestring);
				deltanodes.push_back(breakpoint);
			}
		}
	}
	return deltanodes;
};


std::vector<AccessibleEdge> UndirectedGraph::GetAccessibleEdges(double delta,int facility_id, std::vector<vertex_descriptor> ex_nodes)
{
	std::vector<AccessibleEdge> accessible_edges;
	vertex_descriptor u, pre_u;
	OGRLineString temp_line;	

	for (int i = 1; i < ex_nodes.size(); i++)
	{
		AccessibleEdge temp_acc_edge;

		OGRPoint *breakpoint = new OGRPoint();
		OGRLineString *ptemplinestring = new OGRLineString();
		OGRLineString *s_cut_line = new OGRLineString();
		OGRLineString *e_cut_line = new OGRLineString();
		u = ex_nodes[i];
		pre_u = predecessors_map_[u];

		if ((distances_map_[pre_u] < delta) && ( delta < distances_map_[u]))
		{
			double break_to_snode = delta - distances_map_[pre_u];
			*ptemplinestring = GetEdgeLinestring(pre_u, u);
			algorithm::CutLinestringBasedOnBreakDistance(ptemplinestring, break_to_snode, breakpoint, s_cut_line, e_cut_line); 
			
			temp_acc_edge.facility_id = facility_id;
			temp_acc_edge.s_cost = distances_map_[pre_u];
			temp_acc_edge.t_cost = delta;
			temp_acc_edge.linestring = s_cut_line;
			accessible_edges.push_back(temp_acc_edge);
		}
		else if ((distances_map_[pre_u] <= delta) && (distances_map_[u] <= delta))
		{
			*ptemplinestring = GetEdgeLinestring(pre_u, u);
			temp_acc_edge.facility_id = facility_id;
			temp_acc_edge.s_cost = distances_map_[pre_u];
			temp_acc_edge.t_cost = distances_map_[u];
			temp_acc_edge.linestring = ptemplinestring;
			accessible_edges.push_back(temp_acc_edge);
		}
		else 
		{
			continue;
		}		
	}
	return accessible_edges;
};



void InsertDeltaNodesToGraph(double delta, std::vector<vertex_descriptor> ex_nodes) 
{


}



void UndirectedGraph::ShortestPathAlgorithm(const vertex_descriptor& source, std::vector<vertex_descriptor>& examined_nodes)
{
	//std::vector<vertex_descriptor> m_examined_nodes;
	std::vector<edge_descriptor> tree_edges;
	std::vector<edge_descriptor> non_tree_edges;
	//examined_nodes.push_back(source);
	distances_map_[source] = 0; // here is necessay, or an error will be reported, please see the info 
	//double delta = 1000;
	double inf = std::numeric_limits<double>::max();
	dijkstra_shortest_paths_no_color_map_no_init
	(
		udirected_g_,
		source,
		make_iterator_property_map(predecessors_map_.begin(), get(boost::vertex_index, udirected_g_), predecessors_map_[0]),
		make_iterator_property_map(distances_map_.begin(), get(boost::vertex_index, udirected_g_), distances_map_[0]),
		get(&EdgeProperty::length, udirected_g_),
		get(boost::vertex_index, udirected_g_),
		std::less<double>(), //DistanceCompare distance_compare,
		boost::closed_plus<double>(inf),
		inf,
		0,
		CA_vistor(examined_nodes, distances_map_)
	);

	//p_examined_nodes = m_examined_nodes;

};// end for shortest_path_tree 


/*
   Clean the distance map and predecessor map
 */
void UndirectedGraph::InitializeDistancesPredecessors()
{
	// Need initialization
	int num_vertices = boost::num_vertices(udirected_g_);

	predecessors_map_ = std::vector<vertex_descriptor>(num_vertices);
	distances_map_ = std::vector<double>(num_vertices);

	for (int i = 0; i < num_vertices; ++i) {
		distances_map_[i] = std::numeric_limits<double>::max();
		predecessors_map_[i] = i;
	}
}

void UndirectedGraph::InitializeClosestFacilities(int vertexnumber)
{
	// Need initialization
	closest_facilities_map_ = std::vector<vertex_descriptor>(vertexnumber);
	for (int i = 0; i < vertexnumber; ++i) {
		closest_facilities_map_[i] = i;
	}
}

//
void UndirectedGraph::write_accessible_edges(std::string filename, double delta, int facility_id, std::vector<vertex_descriptor> ex_nodes)
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

	vertex_descriptor u, pre_u;
	OGRLineString temp_line;
	
	for (int i = 1; i < ex_nodes.size(); i++)
	{
		OGRPoint *breakpoint = new OGRPoint();
		OGRLineString* ptemplinestring = new OGRLineString();
		OGRLineString* s_cut_line = new OGRLineString();
		OGRLineString* e_cut_line = new OGRLineString();
		u = ex_nodes[i];
		pre_u = predecessors_map_[u];

		if ((distances_map_[pre_u] < delta) && (delta < distances_map_[u]))
		{
			OGRFeature *poFeature;
			double break_to_snode = delta - distances_map_[pre_u];
			*ptemplinestring = GetEdgeLinestring(pre_u, u);
			algorithm::CutLinestringBasedOnBreakDistance(ptemplinestring, break_to_snode, breakpoint, s_cut_line, e_cut_line);

			poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
			poFeature->SetField("fac_id", facility_id);
			poFeature->SetField("s_cost", distances_map_[pre_u]);
			poFeature->SetField("t_cost", delta);
			poFeature->SetGeometry(s_cut_line);
			
			if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
			{
				printf("Failed to create feature in shapefile.\n");
				exit(1);
			}
			OGRFeature::DestroyFeature(poFeature);
			//delete breakpoint;
			//delete e_cut_line;
			//e_cut_line = NULL;
			//breakpoint = NULL;
		}
		else if ((distances_map_[pre_u] <= delta) && (distances_map_[u] <= delta))
		{
			OGRFeature *poFeature;
			*ptemplinestring = GetEdgeLinestring(pre_u, u);
			
			poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
			poFeature->SetField("fac_id", facility_id);
			poFeature->SetField("s_cost", distances_map_[pre_u]);
			poFeature->SetField("t_cost", distances_map_[u]);
			poFeature->SetGeometry(ptemplinestring);

			if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
			{
				printf("Failed to create feature in shapefile.\n");
				exit(1);
			}
			OGRFeature::DestroyFeature(poFeature);
			//delete breakpoint;
			//delete s_cut_line; 
			//delete e_cut_line;
			//s_cut_line = NULL;
			//e_cut_line = NULL;
			//breakpoint = NULL;
		}
		else 
		{
			//delete ptemplinestring;
			//delete breakpoint;
			//delete s_cut_line;
			//delete e_cut_line;
			//ptemplinestring = NULL;
			//s_cut_line = NULL;
			//e_cut_line = NULL;
			//breakpoint = NULL;
			continue;
		}
	}
	GDALClose(shpDataSet);
};


void UndirectedGraph::wirte_sp_nodes(const std::string &filename, const std::vector<vertex_descriptor> examined_nodes)
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

	vertex_descriptor u;
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


void UndirectedGraph::write_sp_edges(const std::string &filename, const std::vector<vertex_descriptor> examined_nodes) 
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

	OGRFieldDefn oField_s_dis("s_distance", OFTReal);
	oField_s_dis.SetWidth(30);
	oField_s_dis.SetPrecision(5);
	poLayer->CreateField(&oField_s_dis);

	OGRFieldDefn oField_e_dis("e_distance", OFTReal);
	oField_e_dis.SetWidth(30);
	oField_e_dis.SetPrecision(5);
	poLayer->CreateField(&oField_e_dis);


	vertex_descriptor u, pre_u;
	double s_dis, e_dis;
	// note here we start from 1 to exclude the root node
	for (int i = 1; i < examined_nodes.size(); i++)
	{
		u = examined_nodes[i];
		pre_u = predecessors_map_[u];

		s_dis = distances_map_[pre_u];
		e_dis = distances_map_[u];
		
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		OGRLineString *temp_line = new OGRLineString();
		*temp_line = GetEdgeLinestring(pre_u, u);

		poFeature->SetField("edge_id", i);
		poFeature->SetField("s_node", vertex_externalids_[pre_u].c_str());
		poFeature->SetField("e_node", vertex_externalids_[u].c_str());
		
		poFeature->SetField("s_distance", s_dis);
		poFeature->SetField("e_distance", e_dis);

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


// wirte out all the nodes with distance >= delta;
void UndirectedGraph::wirte_sp_nodes_plus_delatnodes(const std::string &filename, const std::vector<vertex_descriptor> examined_nodes, 
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

	vertex_descriptor u;
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


// this is an integrated version of generating nodes information for generating the delaunary trigulations 
void UndirectedGraph::GetNodesForTrigulationInput(const std::string facility_node_sid, const double delta, 
	std::vector<string>& node_ids, std::vector<double>& node_distances, std::vector<OGRPoint>& node_org_points)
{
	std::vector<vertex_descriptor> sp_exaimed_nodes = GetShortestPathNodes(facility_node_sid);
	std::vector<vertex_descriptor> extended_sp_exaimed_nodes = BuildExtendedShortestPathTree(sp_exaimed_nodes, facility_node_sid);
	std::vector<OGRPoint> delta_nodes = FindingDeltaNodes(delta, extended_sp_exaimed_nodes);
	
	int num_ex_nodes = extended_sp_exaimed_nodes.size();
	int num_records = num_ex_nodes + delta_nodes.size();

	node_ids = std::vector<string>(num_records);
	node_distances = std::vector<double>(num_records);
	node_org_points = std::vector<OGRPoint>(num_records);
			
	vertex_descriptor vd_node;
	for (int i = 0; i < extended_sp_exaimed_nodes.size(); i++)
	{
		vd_node = extended_sp_exaimed_nodes[i];
		node_ids[i] = vertex_externalids_[vd_node];
		node_distances[i] = distances_map_[vd_node];
		node_org_points[i] = vertex_coords_[vd_node];
	}

	for (int j = 0; j < delta_nodes.size(); j++)
	{	
		string basic_name = "delta";
		string delta_node_id = basic_name + to_string(j);
		node_ids[(j+ num_ex_nodes)] = delta_node_id;
		node_distances[(j + num_ex_nodes)] = delta;
		node_org_points[(j + num_ex_nodes)] = delta_nodes[j];
	}
} //GetNodesForTrigulationInput


// this function can be used for both solution 2 and solution 3
void UndirectedGraph::GetInputForDelaunary(std::vector<vertex_descriptor> v_exmained_nodes, vector<string>& node_ids, vector<double>& node_distances, vector<OGRPoint>& node_org_points)
{
	int num_records = v_exmained_nodes.size();
	node_ids = std::vector<string>(num_records);
	node_distances = std::vector<double>(num_records);
	node_org_points = std::vector<OGRPoint>(num_records);

	vertex_descriptor vd_node;
	for (int i = 0; i < v_exmained_nodes.size(); i++)
	{
		vd_node = v_exmained_nodes[i];
		node_ids[i] = vertex_externalids_[vd_node];
		node_distances[i] = distances_map_[vd_node];
		node_org_points[i] = vertex_coords_[vd_node];
	}
}
}// end of namespace tca