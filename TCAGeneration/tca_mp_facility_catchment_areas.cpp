#include "tca_mp_facility_catchment_areas.h"

namespace tca 
{

	MPFacilityCatchmentAreas::MPFacilityCatchmentAreas()
{


}

	MPFacilityCatchmentAreas::~MPFacilityCatchmentAreas()
{

}


void MPFacilityCatchmentAreas::CalculateCatchmentAreas(std::string froads, std::string ffacilities)
{
	
	network_.ReadRoadNetworks(froads);

	// the three fields need to be set
	network_.ReadMpFacilities(ffacilities, "ObjectID", "Tq7550", "station_id");

	network_.BuildRtreeIndex();
	std::vector<double>  deltas = network_.get_deltas();
	std::vector<std::string> labels = network_.get_labels();

	std::vector<int> unique_label_indices;               // the local indices of unique facility label
	std::vector<int> duptimes;                          // how many sub-points a facility have 
	std::unordered_map<std::string, int> unique_labels; // (facility label, facility local id)

	// insert the first element
	unique_labels.insert(std::make_pair(labels[0], 0));
	unique_label_indices.push_back(0);

	int dup_time = 1; // indicate how many times a unique value is duplicated
	// iterate from the second label
	for (int i = 1; i < labels.size(); i++)
	{
		auto iter_pair = unique_labels.insert(std::make_pair(labels[i], i));
		if (iter_pair.second) // successfully inserted means a new MP facility is found
		{
			duptimes.push_back(dup_time);
			unique_label_indices.push_back(i); // local index of the unique values
			dup_time = 1;
		}
		else //failed to insert means more sub-points of a facilty is found
		{
			dup_time = dup_time + 1;
		}
	}
	duptimes.push_back(dup_time); // the duplicated times of the last label 
	
	// search the nearest points for all the facilities
	// the paramter should be adapted based on the need
	NearestPoints fn_ps = network_.SearchNearestPoints(500);

	// the number of facility euqal to the number of unique labels
	int facility_num = unique_label_indices.size();

	// a vector of facilities
	std::vector<NearestPoints> unique_fn_ps(facility_num);
	std::vector<double> unique_deltas(facility_num);
	for (int i = 0; i < facility_num; i++)
	{
		int start_id = unique_label_indices[i];
		int end_id = unique_label_indices[i] + duptimes[i];
		unique_fn_ps[i].assign(fn_ps.begin() + start_id, fn_ps.begin() + end_id); // [first, last)
		unique_deltas[i] = deltas[start_id];		
	}

	// Initialization 
	faci_contours_ = vector<Contour>(facility_num);
	faci_contour_labels_ = vector<vector<int>>(facility_num);
	std::vector<std::vector<Edge*>> facility_sub_edges = network_.GetEdgesAroundDeltaMltipleMpFacilities(unique_fn_ps, unique_deltas);

	//int start_minis_id = 0;
	for (int i = 0; i < facility_num; i++)
	{
		double delta = unique_deltas[i];
		std::vector<Edge*> sub_edges = facility_sub_edges[i];

		UndirectedGraph ungraph(sub_edges);
		
		//start_minis_id = start_minis_id + duptimes[i];
		ungraph.InsertingMulltipleFacilityNodes(unique_fn_ps[i]);
				
		std::string  fac_name = "vir_node";
		std::vector<string> node_ids;
		std::vector<double> node_distances;
		std::vector<OGRPoint> node_org_points;
		std::vector<std::array<int, 2>>	edge_node_ids;

		std::vector<vertex_descriptor> sp_exaimed_nodes = ungraph.GetShortestPathNodes(fac_name);
		std::vector<vertex_descriptor> extended_sp_exaimed_nodes = ungraph.BuildExtendedShortestPathTree(sp_exaimed_nodes, fac_name);
		ungraph.GetEdgesForConstrianedTriangulationLineSegmentation(extended_sp_exaimed_nodes, edge_node_ids, node_ids, node_distances, node_org_points);
				
		// build the constrianed trigulations
		BuildTriangulation tribuilder;
		tribuilder.BuildConstrainedTrangulationMultiplePointsVersion(edge_node_ids, node_ids, node_distances, node_org_points);

		std::vector<double> xcoords = tribuilder.get_x_coords();
		std::vector<double> ycoords = tribuilder.get_y_coords();
		std::vector<double> distances = tribuilder.get_z_values();
		std::vector<std::array<int, 3>> delaunay_traiangles = tribuilder.get_triangle_node_ids();
		std::vector<bool> masks = tribuilder.get_masks();
		int triangle_num = tribuilder.get_triangle_num();
		std::vector<std::array<int, 3>> neighbor_array(triangle_num);
		std::vector<std::array<int, 2>> edge_array;
		int correct = 0;

		Triangulation triangles(xcoords, ycoords, delaunay_traiangles, masks, edge_array, neighbor_array, correct);
		TriContourGenerator tri_countour_gtor(triangles, distances);
		
		// this label is used to indicate if the countour line is generated by following boudary lines
		// or following the inner edges
		vector<int> contour_labels;
		Contour nonfilled_countor = tri_countour_gtor.create_contour(delta, contour_labels);
		faci_contours_[i] = nonfilled_countor;
		faci_contour_labels_[i] = contour_labels;
	}
}


void MPFacilityCatchmentAreas::CalculateCatchmentAreas(const std::string froads, const std::string &netid, const std::string &source_name, const std::string &target_name,
	const std::string ffacilities, const std::string &fid, const std::string &flabel, const bool &isunifiedcutoff, const double &unifiedcutoff, const std::string &cutoff_name,
	const double &searchradius, bool writeaccedges, bool writecatchment) 
{
	network_.ReadRoadNetworks(froads, netid, source_name, target_name);
	network_.ReadMpFacilities(ffacilities, fid, flabel, isunifiedcutoff, unifiedcutoff, cutoff_name);
	network_.BuildRtreeIndex();
	std::vector<double>  deltas = network_.get_deltas();
	std::vector<std::string> labels = network_.get_labels();

	std::vector<int> unique_label_indices;               // the local indices of unique facility label
	std::vector<int> duptimes;                          // how many sub-points a facility have 
	std::unordered_map<std::string, int> unique_labels; // (facility label, facility local id)

	// insert the first element
	unique_labels.insert(std::make_pair(labels[0], 0));
	unique_label_indices.push_back(0);

	int dup_time = 1; // indicate how many times a unique value is duplicated
	// iterate from the second label
	for (int i = 1; i < labels.size(); i++)
	{
		auto iter_pair = unique_labels.insert(std::make_pair(labels[i], i));
		if (iter_pair.second) // successfully inserted means a new MP facility is found
		{
			duptimes.push_back(dup_time);
			unique_label_indices.push_back(i); // local index of the unique values
			dup_time = 1;
		}
		else //failed to insert means more sub-points of a facilty is found
		{
			dup_time = dup_time + 1;
		}
	}
	duptimes.push_back(dup_time); // the duplicated times of the last label 

	// search the nearest points for all the facilities
	// the paramter should be adapted based on the need
	NearestPoints fn_ps = network_.SearchNearestPoints(searchradius);

	// the number of facility euqal to the number of unique labels
	int facility_num = unique_label_indices.size();
	// a vector of facilities
	std::vector<NearestPoints> unique_fn_ps(facility_num);
	std::vector<double> unique_deltas(facility_num);
	for (int i = 0; i < facility_num; i++)
	{
		int start_id = unique_label_indices[i];
		int end_id = unique_label_indices[i] + duptimes[i];
		unique_fn_ps[i].assign(fn_ps.begin() + start_id, fn_ps.begin() + end_id); // [first, last)
		unique_deltas[i] = deltas[start_id];
	}
	
	// Initialization 
	std::vector<std::vector<Edge*>> facility_sub_edges = network_.GetEdgesAroundDeltaMltipleMpFacilities(unique_fn_ps, unique_deltas);
	
	if (writeaccedges)
	{
		faci_accedges_ = vector<vector<AccessibleEdge>>(facility_num);
	}
	if (writecatchment)
	{
		faci_contours_ = vector<Contour>(facility_num);
		faci_contour_labels_ = vector<vector<int>>(facility_num);
	}

	//int start_minis_id = 0;
	for (int i = 0; i < facility_num; i++)
	{
		double delta = unique_deltas[i];
		std::vector<Edge*> sub_edges = facility_sub_edges[i];
		UndirectedGraph ungraph(sub_edges);

		//start_minis_id = start_minis_id + duptimes[i];
		ungraph.InsertingMulltipleFacilityNodes(unique_fn_ps[i]);

		std::string  fac_name = "vir_node";
		std::vector<string> node_ids;
		std::vector<double> node_distances;
		std::vector<OGRPoint> node_org_points;
		std::vector<std::array<int, 2>>	edge_node_ids;

		std::vector<vertex_descriptor> sp_exaimed_nodes = ungraph.GetShortestPathNodes(fac_name);
		std::vector<vertex_descriptor> extended_sp_exaimed_nodes = ungraph.BuildExtendedShortestPathTree(sp_exaimed_nodes, fac_name);
		
		if (writeaccedges) // if output the accessible edges
		{
			vector<AccessibleEdge> siglefa_edges = ungraph.GetAccessibleEdges(delta, i, extended_sp_exaimed_nodes);
			faci_accedges_[i] = siglefa_edges;
		}
		if (writecatchment) 
		{
			ungraph.GetEdgesForConstrianedTriangulationLineSegmentation(extended_sp_exaimed_nodes, edge_node_ids, node_ids, node_distances, node_org_points);

			// build the constrianed trigulations
			BuildTriangulation tribuilder;
			tribuilder.BuildConstrainedTrangulationMultiplePointsVersion(edge_node_ids, node_ids, node_distances, node_org_points);

			std::vector<double> xcoords = tribuilder.get_x_coords();
			std::vector<double> ycoords = tribuilder.get_y_coords();
			std::vector<double> distances = tribuilder.get_z_values();
			std::vector<std::array<int, 3>> delaunay_traiangles = tribuilder.get_triangle_node_ids();
			std::vector<bool> masks = tribuilder.get_masks();
			int triangle_num = tribuilder.get_triangle_num();
			std::vector<std::array<int, 3>> neighbor_array(triangle_num);
			std::vector<std::array<int, 2>> edge_array;
			int correct = 0;

			Triangulation triangles(xcoords, ycoords, delaunay_traiangles, masks, edge_array, neighbor_array, correct);
			TriContourGenerator tri_countour_gtor(triangles, distances);

			// this label is used to indicate if the countour line is generated by following boudary lines
			// or following the inner edges
			vector<int> contour_labels;
			Contour nonfilled_countor = tri_countour_gtor.create_contour(delta, contour_labels);
			faci_contours_[i] = nonfilled_countor;
			faci_contour_labels_[i] = contour_labels;

		}
	}
}




void MPFacilityCatchmentAreas::CalculateCatchmentAreasDirected(const bool &from_facility, const std::string froads,
	const std::string &netid, const std::string &source_name, const std::string &target_name, const std::string &direction_name,
	const std::string ffacilities, const std::string &fid, const std::string &flabel,const bool &isunifiedcutoff, const double &unifiedcutoff, const std::string &cutoff_name,
	const double &searchradius, bool writeaccedges, bool writecatchment) 
{
	if (from_facility)
	{
		network_.ReadDirectedRoadnetworks(froads, netid, source_name, target_name, direction_name);
	}
	else
	{
		network_.ReadDirectedRoadnetworksToFacility(froads, netid, source_name, target_name, direction_name);
	}
	
	network_.ReadMpFacilities(ffacilities, fid, flabel, isunifiedcutoff, unifiedcutoff, cutoff_name);
	network_.BuildRtreeIndexDirected();					//directed
	std::vector<double>  deltas = network_.get_deltas();
	std::vector<std::string> labels = network_.get_labels();

	std::vector<int> unique_label_indices;               // the local indices of unique facility label
	std::vector<int> duptimes;                          // how many sub-points a facility have 
	std::unordered_map<std::string, int> unique_labels; // (facility label, facility local id)

	// insert the first element
	unique_labels.insert(std::make_pair(labels[0], 0));
	unique_label_indices.push_back(0);

	int dup_time = 1; // indicate how many times a unique value is duplicated
	// iterate from the second label
	for (int i = 1; i < labels.size(); i++)
	{
		auto iter_pair = unique_labels.insert(std::make_pair(labels[i], i));
		if (iter_pair.second) // successfully inserted means a new MP facility is found
		{
			duptimes.push_back(dup_time);
			unique_label_indices.push_back(i); // local index of the unique values
			dup_time = 1;
		}
		else //failed to insert means more sub-points of a facilty is found
		{
			dup_time = dup_time + 1;
		}
	}
	duptimes.push_back(dup_time); // the duplicated times of the last label 

	NearestPointsDirected fn_ps = network_.SearchNearestPointsDirected(searchradius);

	// the number of facility euqal to the number of unique labels
	int facility_num = unique_label_indices.size();
	// a vector of facilities
	std::vector<NearestPointsDirected> unique_fn_ps(facility_num);
	std::vector<double> unique_deltas(facility_num);
	for (int i = 0; i < facility_num; i++)
	{
		int start_id = unique_label_indices[i];
		int end_id = unique_label_indices[i] + duptimes[i];
		unique_fn_ps[i].assign(fn_ps.begin() + start_id, fn_ps.begin() + end_id); // [first, last)
		unique_deltas[i] = deltas[start_id];
	}

	// Initialization 
	std::vector<std::vector<EdgeDirected*>> facility_sub_edges = network_.GetEdgesAroundDeltaMltipleMpFacilitiesDirected(unique_fn_ps, unique_deltas);
	
	if (writeaccedges)
	{
		faci_accedges_ = vector<vector<AccessibleEdge>>(facility_num);
	}
	if (writecatchment)
	{
		faci_contours_ = vector<Contour>(facility_num);
		faci_contour_labels_ = vector<vector<int>>(facility_num);
	}
	
	//int start_minis_id = 0;
	for (int i = 0; i < facility_num; i++)
	{
		double delta = unique_deltas[i];
		std::vector<EdgeDirected*> sub_edges = facility_sub_edges[i];
		DirectedGraph directedgraph(sub_edges);

		//start_minis_id = start_minis_id + duptimes[i];
		directedgraph.InsertingMulltipleFacilityNodes(unique_fn_ps[i]);

		std::string  fac_name = "vir_node";
		std::vector<string> node_ids;
		std::vector<double> node_distances;
		std::vector<OGRPoint> node_org_points;
		std::vector<std::array<int, 2>>	edge_node_ids;

		std::vector<vertex_descriptor_drct> sp_exaimed_nodes = directedgraph.GetShortestPathNodes(fac_name);
		std::vector<vertex_descriptor_drct> extended_sp_exaimed_nodes = directedgraph.BuildExtendedShortestPathTree(sp_exaimed_nodes, fac_name);

		if (writeaccedges) // if output the accessible edges
		{
			vector<AccessibleEdge> siglefa_edges = directedgraph.GetAccessibleEdges(delta, i, extended_sp_exaimed_nodes);
			faci_accedges_[i] = siglefa_edges;
		}
		if (writecatchment)
		{
			directedgraph.GetEdgesForConstrianedTriangulationLineSegmentation(delta, extended_sp_exaimed_nodes, 
				edge_node_ids, node_ids, node_distances, node_org_points);

			// build the constrianed trigulations
			BuildTriangulation tribuilder;
			tribuilder.BuildConstrainedTrangulationMultiplePointsVersion(edge_node_ids, node_ids, node_distances, node_org_points);

			std::vector<double> xcoords = tribuilder.get_x_coords();
			std::vector<double> ycoords = tribuilder.get_y_coords();
			std::vector<double> distances = tribuilder.get_z_values();
			std::vector<std::array<int, 3>> delaunay_traiangles = tribuilder.get_triangle_node_ids();
			std::vector<bool> masks = tribuilder.get_masks();
			int triangle_num = tribuilder.get_triangle_num();
			std::vector<std::array<int, 3>> neighbor_array(triangle_num);
			std::vector<std::array<int, 2>> edge_array;
			int correct = 0;

			Triangulation triangles(xcoords, ycoords, delaunay_traiangles, masks, edge_array, neighbor_array, correct);
			TriContourGenerator tri_countour_gtor(triangles, distances);

			// this label is used to indicate if the countour line is generated by following boudary lines
			// or following the inner edges
			vector<int> contour_labels;
			Contour nonfilled_countor = tri_countour_gtor.create_contour(delta, contour_labels);
			faci_contours_[i] = nonfilled_countor;
			faci_contour_labels_[i] = contour_labels;

		}
	}

}



void MPFacilityCatchmentAreas::OutputExtendedSPTrees(std::string froads, std::string ffacilities, std::string basic_name)
{
	// read the entire road network and the facilities	
	network_.ReadRoadNetworks(froads);
	network_.ReadFacilities(ffacilities, "id", "d1200");
	// network_.ReadFacilities(ffacilities, "id", "delta")
	//	network_.ReadFacilities(ffacilities, "ObjectID", "Breaks_len");
	network_.BuildRtreeIndex();

	std::vector<double>  delats = network_.get_deltas();
	// search the nearest points for all the facilities
	NearestPoints fn_ps = network_.SearchNearestPoints(500);
	// NearestPoints fn_ps = network_.SearchNearestPoints(100);
	// get all the sub-edges around delta buffer for each facility
	std::vector<std::vector<Edge*>> facility_sub_edges = network_.GetEdgesAroundDelta(fn_ps);

	int facility_num = facility_sub_edges.size();
	faci_contours_ = vector<Contour>(facility_num);
	faci_contour_labels_ = vector<vector<int>>(facility_num);

	//std::string basic_name = "C://Users//dlint//source//repos//TCA_Nets//data_input//acticle_experiment//output//Extende_SP_graphs//extended_sp_graph_";
	std::string fextend_sp_edges;
	for (int i = 0; i < facility_num; i++)
	{
		std::string  fac_name = "fac" + to_string(i);
		double delta = delats[i];

		UndirectedGraph ungraph(facility_sub_edges[i]);
		ungraph.InsertingFacilityNode(fn_ps[i]);
		std::vector<string> node_ids;
		std::vector<double> node_distances;
		std::vector<OGRPoint> node_org_points;
		std::vector<std::array<int, 2>>	edge_node_ids;

		std::vector<vertex_descriptor> sp_exaimed_nodes = ungraph.GetShortestPathNodes(fac_name);
		std::vector<vertex_descriptor> extended_sp_exaimed_nodes = ungraph.BuildExtendedShortestPathTree(sp_exaimed_nodes, fac_name);
		fextend_sp_edges = basic_name + std::to_string(i) + ".shp";
		ungraph.write_sp_edges(fextend_sp_edges, extended_sp_exaimed_nodes);
	}
	//code for generating the node distances + generating the catchment areas
}



void MPFacilityCatchmentAreas::write_contour_polygons(const std::string &filename)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("con_ploys", NULL, wkbPolygon, NULL);

	// define and create filed
	OGRFieldDefn oField_tri_id("id", OFTInteger);
	oField_tri_id.SetWidth(10);
	poLayer->CreateField(&oField_tri_id);
	int contour_num = faci_contours_.size();

	for (int k = 0; k < contour_num; k++)
	{
		Contour contour_lines = faci_contours_[k];
		int line_num = contour_lines.size();

		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		OGRPolygon poly;
		for (int i = 0; i < line_num; i++)
		{
			OGRLinearRing  linering;
			if (contour_lines[i].size() != 1) // some contour line is not closed with only one point, thus need to be excluded
			{
				for (int j = 0; j < contour_lines[i].size(); j++)
				{
					linering.addPoint(contour_lines[i][j].x, contour_lines[i][j].y);
				}
				linering.closeRings();
			}
			poly.addRing(&linering);
		}

		poFeature->SetField("id", k);		
		poFeature->SetGeometry(&poly);

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	
	GDALClose(shpDataSet);
}


void MPFacilityCatchmentAreas::write_accessible_edges(const std::string &filename)
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

	// iterate every facility 
	for (int i = 0; i < faci_accedges_.size(); i++)
	{
		vector<AccessibleEdge> single_ficility = faci_accedges_[i];
		// iterate every accessible edge of each facility
		for (int j = 0; j < single_ficility.size(); j++)
		{
			AccessibleEdge accedge = single_ficility[j];
			OGRFeature *poFeature;
			poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
			poFeature->SetField("fac_id", accedge.facility_id);
			poFeature->SetField("s_cost", accedge.s_cost);
			poFeature->SetField("t_cost", accedge.t_cost);
			poFeature->SetGeometry(accedge.linestring);
			if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
			{
				printf("Failed to create feature in shapefile.\n");
				exit(1);
			}
			OGRFeature::DestroyFeature(poFeature);
		}

	}
	GDALClose(shpDataSet);
}



}// end of namespace tca