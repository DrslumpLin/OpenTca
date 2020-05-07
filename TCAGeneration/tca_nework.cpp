
#include "tca_network.h"

namespace tca 
{

Network::Network()
{

}

Network::~Network()
{
	std::cout << "Cleaning network" << '\n';
	for (auto &item : network_edges_)
	{
		OGRGeometryFactory::destroyGeometry(item.line_string);
	}
	std::cout << "Cleaning network finished" << '\n';
}

void Network::ReadRoadNetworks(const std::string filename, const std::string &id_name, const std::string &source_name,
	const std::string &target_name)
{
	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(),GDAL_OF_VECTOR, NULL, NULL, NULL));
	if (poDS == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}
	OGRLayer *ogrlayer = poDS->GetLayer(0);

	// Initialize network edges
	OGRFeatureDefn *ogrFDefn = ogrlayer->GetLayerDefn();
	OGRFeature *ogrFeature;
	int num_features = ogrlayer->GetFeatureCount();
	network_edges_ = std::vector<Edge>(num_features);
	
	// Fetch the field index given field name.
	int id_idx = ogrFDefn->GetFieldIndex(id_name.c_str());
	int source_idx = ogrFDefn->GetFieldIndex(source_name.c_str());
	int target_idx = ogrFDefn->GetFieldIndex(target_name.c_str());


	while ((ogrFeature = ogrlayer->GetNextFeature()) != NULL)
	{
		int id = ogrFeature->GetFID();
		Edge e;

		e.id = id;
		e.external_id = std::string(ogrFeature->GetFieldAsString(id_idx));
		e.source = std::string(ogrFeature->GetFieldAsString(source_idx));
		e.target = std::string(ogrFeature->GetFieldAsString(target_idx));
		OGRGeometry *rawlinegeometry = ogrFeature->GetGeometryRef();
		e.line_string = (OGRLineString*)rawlinegeometry->clone();
		e.length = e.line_string->get_Length();
		network_edges_[id] = e;
		edge_linestring_map_.insert({ e.id, e.line_string }); //edge linestring map, id is a inner id and continous
		OGRFeature::DestroyFeature(ogrFeature);
	}
	GDALClose(poDS);

	std::cout << "Read network finish." << '\n';
	std::cout << "\tTotal number of edges read " << network_edges_.size() << '\n';
}



std::string Network::GetEdgeExternalId(int internal_id)
{
	return network_edges_[internal_id].external_id;
}


void Network::BuildRtreeIndex()
{
	std::cout << "Start to construct boost rtree" << '\n';
	for (std::size_t i = 0; i < network_edges_.size(); ++i)
	{
		Edge *edge = &network_edges_[i];
		double x1, y1, x2, y2;
		
		algorithm::GetLinestringBoundingbox(edge->line_string, &x1, &y1, &x2, &y2);
		boost_box b(boost_point(x1, y1), boost_point(x2, y2));
		rtree_.insert(std::make_pair(b, edge));
	}
	std::cout << "Finish construct boost rtree" << '\n';
}

NearestPoints Network::SearchNearestPoints(double radius)
{
	int num_facilities = facilities_.size();
	NearestPoints ficility_nearest_points(num_facilities);
	
	for (int i = 0; i < num_facilities; ++i)
	{
		//if (i == 11) {std::cout << "test";}
		std::vector<NearestPoint> local_nearest_points;

		OGRPoint* point = facilities_[i].point;
		boost_point pt = boost_point(point->getX(), point->getY());

		std::vector<Item> temp;
		boost_box b(boost_point(point->getX() - radius, point->getY() - radius), boost_point(point->getX() + radius, point->getY() + radius));
		rtree_.query(bgi::intersects(b), std::back_inserter(temp));
		
		for (Item &j : temp)  // Rtree can only detect intersect with a the bounding box of the geometry stored.
		{
			Edge *edge = j.second;
			float break_distance;
			double virtical_distance;
			int cut_seg_index;
			
			OGRPoint* p_break_point = new OGRPoint;
			// here the offset distance is measured to the starting node
			algorithm::ProjectPointOnLinestring(point, edge->line_string, &cut_seg_index, &virtical_distance, &break_distance);
			*p_break_point = algorithm::GetBreakPointOnLinestring(break_distance, edge->line_string);
			NearestPoint np = NearestPoint{i,virtical_distance,break_distance,cut_seg_index, p_break_point, edge};			
			local_nearest_points.push_back(np);
		}

		// below are used to fetch the point with smallest distance to the facility
		if (local_nearest_points.empty())
		{
			return NearestPoints();
			std::cout << "no nearest edges are achieved, consider enlarge the search radius\n";
			
		};
		if (local_nearest_points.size() == 1)
		{
			ficility_nearest_points[i] = local_nearest_points[0];
		}
		else
		{
			std::sort(local_nearest_points.begin(), local_nearest_points.end(), compareBydis);
			ficility_nearest_points[i] = local_nearest_points[0];
		}
	}
	return ficility_nearest_points;
}



// This is used for comparision purposes,when no r-tree is built
NearestPoints Network::SearchNearestPointsNoRtree(double radius)
{
	int num_facilities = facilities_.size();
	NearestPoints ficility_nearest_points(num_facilities);
	
	for (int i = 0; i < num_facilities; ++i)
	{
		OGRPoint* point = facilities_[i].point;
		double nearest_vertical_distance =  9999999999999999;
		
		float break_distance;
		double virtical_distance;
		int cut_seg_index;
		NearestPoint np;
		
		for (int j = 0; j < network_edges_.size(); ++j)
		{
			Edge edge = network_edges_[j];
			OGRPoint* p_break_point = new OGRPoint;
			// here the offset distance is measured to the starting node
			algorithm::ProjectPointOnLinestring(point, edge.line_string, &cut_seg_index, &virtical_distance, &break_distance);

			if (virtical_distance < nearest_vertical_distance)
			{
				*p_break_point = algorithm::GetBreakPointOnLinestring(break_distance, edge.line_string);
				np = NearestPoint{j,virtical_distance,break_distance,cut_seg_index, p_break_point, &edge};
			}
		}
		ficility_nearest_points[i] = np;
	}
	return ficility_nearest_points;
}


std::vector<std::vector<Edge*>> Network::GetEdgesAroundDelta(double delta, NearestPoints  facility_nearest_points)
{
	int num_facilities = facility_nearest_points.size();
	std::vector<std::vector<Edge*>> facility_sub_networks(num_facilities);

	for (int i = 0; i < num_facilities; ++i)
	{
		std::vector<Edge*> temp_sub_network;
		// cut-off point
		OGRPoint* cut_point = facility_nearest_points[i].point;

		std::vector<Item> temp;
		// build a boudning box based on delata and 
		boost_box b(boost_point(cut_point->getX() - delta, cut_point->getY() - delta), boost_point(cut_point->getX() + delta, cut_point->getY() + delta));
		rtree_.query(bgi::intersects(b), std::back_inserter(temp));

		// Rtree can only detect intersect with a the bounding box of the geometry stored.
		for (Item &j : temp)
		{
			Edge *edge = j.second;
			temp_sub_network.push_back(edge);
		}
		facility_sub_networks[i] = temp_sub_network;
	}
	return facility_sub_networks;
}


std::vector<std::vector<Edge*>> Network::GetEdgesAroundDelta(NearestPoints  facility_nearest_points)
{
	int num_facilities = facility_nearest_points.size();
	std::vector<std::vector<Edge*>> facility_sub_networks(num_facilities);

	for (int i = 0; i < num_facilities; ++i)
	{
		double delta = facilities_[i].delta;
		std::vector<Edge*> temp_sub_network;
		// cut-off point
		OGRPoint* cut_point = facility_nearest_points[i].point;

		std::vector<Item> temp;
		// build a boudning box based on delata and 
		boost_box b(boost_point(cut_point->getX() - delta, cut_point->getY() - delta), boost_point(cut_point->getX() + delta, cut_point->getY() + delta));
		rtree_.query(bgi::intersects(b), std::back_inserter(temp));

		// Rtree can only detect intersection with the bounding box of the geometry stored.
		for (Item &j : temp)
		{
			Edge *edge = j.second;			
			temp_sub_network.push_back(edge);
		}
		facility_sub_networks[i] = temp_sub_network;
	}
	return facility_sub_networks;
}

// This is to get the sub edges of a single facility that represented by multiple points
std::vector<Edge*>  Network::GetEdgesAroundDeltaSingleMpFacility(NearestPoints  facility_nearest_points, double delta)
{
	std::vector<Edge*> sub_edges;
	int num_multiple_points = facility_nearest_points.size();
	double box_x_min = DBL_MAX;
	double box_y_min = DBL_MAX;
	double box_x_max = DBL_MIN;
	double box_y_max = DBL_MIN;
	//double delta;
	// get the bounding box of the multiple points
	for (int i = 0; i < num_multiple_points; ++i)
	{
		//delta = facilities_[i].delta; // here we assume that the deltas of all the points are equal
		OGRPoint* projected_point = facility_nearest_points[i].point;		
		double x = projected_point->getX();
		double y = projected_point->getY();
		if (x < box_x_min) box_x_min = x;
		if (y < box_y_min) box_y_min = y;
		if (x > box_x_max) box_x_max = x;
		if (y > box_y_max) box_y_max = y;
	}
	// create the bbox for sub-edges and extract corresponding sub edges
	std::vector<Item> temp;
	boost_box b(boost_point(box_x_min - delta, box_y_min - delta), boost_point(box_x_max + delta, box_y_max + delta));
	rtree_.query(bgi::intersects(b), std::back_inserter(temp));	
	for (Item &j : temp)
	{
		Edge *edge = j.second;
		sub_edges.push_back(edge);
	}
	return sub_edges;	
}


std::vector<EdgeDirected*> Network::GetEdgesAroundDeltaSingleMpFacilityDirect(NearestPointsDirected facility_nearest_points, double delta) 
{
	std::vector<EdgeDirected*> sub_edges;
	int num_multiple_points = facility_nearest_points.size();
	double box_x_min = DBL_MAX;
	double box_y_min = DBL_MAX;
	double box_x_max = DBL_MIN;
	double box_y_max = DBL_MIN;
	//double delta;
	// get the bounding box of the multiple points
	for (int i = 0; i < num_multiple_points; ++i)
	{
		//delta = facilities_[i].delta; // here we assume that the deltas of all the points are equal
		OGRPoint* projected_point = facility_nearest_points[i].point;
		double x = projected_point->getX();
		double y = projected_point->getY();
		if (x < box_x_min) box_x_min = x;
		if (y < box_y_min) box_y_min = y;
		if (x > box_x_max) box_x_max = x;
		if (y > box_y_max) box_y_max = y;
	}
	// create the bbox for sub-edges and extract corresponding sub edges
	std::vector<ItemDirected> temp;
	boost_box b(boost_point(box_x_min - delta, box_y_min - delta), boost_point(box_x_max + delta, box_y_max + delta));
	rtree_directed_.query(bgi::intersects(b), std::back_inserter(temp));
	for (ItemDirected &j : temp)
	{
		EdgeDirected *edge = j.second;
		sub_edges.push_back(edge);
	}
	return sub_edges;
}


std::vector<std::vector<Edge*>> Network::GetEdgesAroundDeltaMltipleMpFacilities(std::vector<NearestPoints> unique_facilitis, std::vector<double> unique_deltas)
{

	int facility_num = unique_facilitis.size();
	std::vector<std::vector<Edge*>> facility_sub_networks(facility_num);

	for (int i = 0; i < facility_num; i++)
	{
		double delta = unique_deltas[i];
		std::vector<Edge*> sub_edges = GetEdgesAroundDeltaSingleMpFacility(unique_facilitis[i], delta);		
		facility_sub_networks[i] = sub_edges;
	}
	return facility_sub_networks;
}

std::vector<std::vector<EdgeDirected*>> Network::GetEdgesAroundDeltaMltipleMpFacilitiesDirected(std::vector<NearestPointsDirected> unique_facilitis, std::vector<double> unique_deltas)
{

	int facility_num = unique_facilitis.size();
	std::vector<std::vector<EdgeDirected*>> facility_sub_networks(facility_num);

	for (int i = 0; i < facility_num; i++)
	{
		double delta = unique_deltas[i];
		std::vector<EdgeDirected*> sub_edges = GetEdgesAroundDeltaSingleMpFacilityDirect(unique_facilitis[i], delta);
		facility_sub_networks[i] = sub_edges;
	}
	return facility_sub_networks;
}


std::vector<Edge*> Network::GetAllEdges() 
{
	return network_edge_pointers_;
}


//----------------------------------for directed road network input------------------------------------------------------
void Network::ReadDirectedRoadnetworks(std::string filename, const std::string &id_name, const std::string &source_name, const std::string &target_name, const std::string &direction_name)
{
	//"C://Users//dlint//source//repos//TCA_Nets//data_input//directed//select//edges_select_dire.shp"
	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(),GDAL_OF_VECTOR, NULL, NULL, NULL));
	if (poDS == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}
	OGRLayer *ogrlayer = poDS->GetLayer(0);

	// Initialize network edges
	OGRFeatureDefn *ogrFDefn = ogrlayer->GetLayerDefn();
	OGRFeature *ogrFeature;

	// Fetch the field index given field name.
	int id_idx = ogrFDefn->GetFieldIndex(id_name.c_str());
	int source_idx = ogrFDefn->GetFieldIndex(source_name.c_str());
	int target_idx = ogrFDefn->GetFieldIndex(target_name.c_str());
	int direction_idx = ogrFDefn->GetFieldIndex(direction_name.c_str());

	int inner_id = 0;
	while ((ogrFeature = ogrlayer->GetNextFeature()) != NULL)
	{
		int direction = ogrFeature->GetFieldAsInteger(direction_idx);

		if (direction == 0)
		{
			EdgeDirected e1, e2; // e1 means source - target; e2 means traget to source
			e1.id = inner_id;
			e1.external_id = std::string(ogrFeature->GetFieldAsString(id_idx));
			e1.source = std::string(ogrFeature->GetFieldAsString(source_idx));
			e1.target = std::string(ogrFeature->GetFieldAsString(target_idx));
			e1.direction = 1;
			OGRGeometry *rawlinegeometry1 = ogrFeature->GetGeometryRef();
			e1.line_string = (OGRLineString*)rawlinegeometry1->clone();
			e1.length = e1.line_string->get_Length();

			e2.id = inner_id + 1;
			e2.external_id = std::string(ogrFeature->GetFieldAsString(id_idx)) + "_op";  // the external id here need to be different with the original one
			e2.source = std::string(ogrFeature->GetFieldAsString(target_idx));
			e2.target = std::string(ogrFeature->GetFieldAsString(source_idx));
			e2.direction = -1;
			OGRGeometry *rawlinegeometry2 = ogrFeature->GetGeometryRef();
			e2.line_string = ((OGRLineString*)rawlinegeometry2->clone());
			e2.length = e2.line_string->get_Length();
			e2.line_string->reversePoints();         // the points of linestring need to be reversed here

			network_directed_edges_.push_back(e1);
			network_directed_edges_.push_back(e2);

			edge_linestring_map_.insert({ e1.id, e1.line_string });
			edge_linestring_map_.insert({ e2.id, e2.line_string });

			inner_id = inner_id + 2;
		}
		else if (direction == 1)
		{
			EdgeDirected e1; // e1 means source - target
			e1.id = inner_id;
			e1.external_id = std::string(ogrFeature->GetFieldAsString(id_idx));
			e1.source = std::string(ogrFeature->GetFieldAsString(source_idx));
			e1.target = std::string(ogrFeature->GetFieldAsString(target_idx));
			e1.direction = 2;              // the direction here infact is no use any more 
			OGRGeometry *rawlinegeometry1 = ogrFeature->GetGeometryRef();
			e1.line_string = (OGRLineString*)rawlinegeometry1->clone();
			e1.length = e1.line_string->get_Length();

			network_directed_edges_.push_back(e1);
			edge_linestring_map_.insert({ e1.id, e1.line_string });
			inner_id = inner_id + 1;
		}
		else
		{
			EdgeDirected e2; // e1 means target - source
			e2.id = inner_id;
			e2.external_id = std::string(ogrFeature->GetFieldAsString(id_idx));
			e2.source = std::string(ogrFeature->GetFieldAsString(target_idx));
			e2.target = std::string(ogrFeature->GetFieldAsString(source_idx));
			e2.direction = 3;                   // the direction here infact is no use any more 
			OGRGeometry *rawlinegeometry2 = ogrFeature->GetGeometryRef();
			e2.line_string = ((OGRLineString*)rawlinegeometry2->clone());
			e2.length = e2.line_string->get_Length();
			e2.line_string->reversePoints();         // the points of linestring need to be reversed here

			network_directed_edges_.push_back(e2);
			edge_linestring_map_.insert({ e2.id, e2.line_string });
			inner_id = inner_id + 1;
		}
		OGRFeature::DestroyFeature(ogrFeature);
	}
	GDALClose(poDS);

	std::cout << "Read --directed-- network finish." << '\n';
	std::cout << "\tTotal number of edges read " << network_edges_.size() << '\n';
}

// this version is for calculting the to-facilty catchment area and network voronoi
// we change the direction of all directed edge into its oppsite direction
void Network::ReadDirectedRoadnetworksToFacility(std::string filename, const std::string &id_name, const std::string &source_name, const std::string &target_name, const std::string &direction_name)
{
	//"C://Users//dlint//source//repos//TCA_Nets//data_input//directed//select//edges_select_dire.shp"
	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(),GDAL_OF_VECTOR, NULL, NULL, NULL));
	if (poDS == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}
	OGRLayer *ogrlayer = poDS->GetLayer(0);

	// Initialize network edges
	OGRFeatureDefn *ogrFDefn = ogrlayer->GetLayerDefn();
	OGRFeature *ogrFeature;

	// Fetch the field index given field name.
	int id_idx = ogrFDefn->GetFieldIndex(id_name.c_str());
	int source_idx = ogrFDefn->GetFieldIndex(source_name.c_str());
	int target_idx = ogrFDefn->GetFieldIndex(target_name.c_str());
	int direction_idx = ogrFDefn->GetFieldIndex(direction_name.c_str());

	int inner_id = 0;
	while ((ogrFeature = ogrlayer->GetNextFeature()) != NULL)
	{
		int direction = ogrFeature->GetFieldAsInteger(direction_idx);

		if (direction == 0) // no change for 
		{
			EdgeDirected e1, e2; // e1 means source - target; e2 means traget to source
			e1.id = inner_id;
			e1.external_id = std::string(ogrFeature->GetFieldAsString(id_idx));
			e1.source = std::string(ogrFeature->GetFieldAsString(source_idx));
			e1.target = std::string(ogrFeature->GetFieldAsString(target_idx));
			e1.direction = 1;
			OGRGeometry *rawlinegeometry1 = ogrFeature->GetGeometryRef();
			e1.line_string = (OGRLineString*)rawlinegeometry1->clone();
			e1.length = e1.line_string->get_Length();

			e2.id = inner_id + 1;
			e2.external_id = std::string(ogrFeature->GetFieldAsString(id_idx)) + "_op";  // the external id here need to be different with the original one
			e2.source = std::string(ogrFeature->GetFieldAsString(target_idx));
			e2.target = std::string(ogrFeature->GetFieldAsString(source_idx));
			e2.direction = -1;
			OGRGeometry *rawlinegeometry2 = ogrFeature->GetGeometryRef();
			e2.line_string = ((OGRLineString*)rawlinegeometry2->clone());
			e2.length = e2.line_string->get_Length();
			e2.line_string->reversePoints();         // the points of linestring need to be reversed here

			network_directed_edges_.push_back(e1);
			network_directed_edges_.push_back(e2);

			edge_linestring_map_.insert({ e1.id, e1.line_string });
			edge_linestring_map_.insert({ e2.id, e2.line_string });

			inner_id = inner_id + 2;
		}
		else if (direction == 1)// orignal is source - target, we need to transfer to target - source and change the direction of linestring
		{
			EdgeDirected e1; 
			e1.id = inner_id;
			e1.external_id = std::string(ogrFeature->GetFieldAsString(id_idx));
			
			OGRGeometry *rawlinegeometry1 = ogrFeature->GetGeometryRef();
			e1.line_string = (OGRLineString*)rawlinegeometry1->clone();
			e1.length = e1.line_string->get_Length();				
			e1.line_string->reversePoints();  // the points of linestring need to be reversed here																 
			e1.source = std::string(ogrFeature->GetFieldAsString(target_idx));								
			e1.target = std::string(ogrFeature->GetFieldAsString(source_idx));
			
			e1.direction = 2;  // the direction here infact is no use any more 
			network_directed_edges_.push_back(e1);
			edge_linestring_map_.insert({ e1.id, e1.line_string });
			inner_id = inner_id + 1;
		}
		else // orinigal is target - source, we need to transfer to source - target, then there is no need to change the direciton of linestring 
		{
			EdgeDirected e2;
			e2.id = inner_id;
			e2.external_id = std::string(ogrFeature->GetFieldAsString(id_idx));
			
			OGRGeometry *rawlinegeometry2 = ogrFeature->GetGeometryRef();
			e2.line_string = ((OGRLineString*)rawlinegeometry2->clone());
			e2.length = e2.line_string->get_Length();

			e2.source = std::string(ogrFeature->GetFieldAsString(source_idx));				
			e2.target = std::string(ogrFeature->GetFieldAsString(target_idx));
			
			e2.direction = 3; 
			network_directed_edges_.push_back(e2);
			edge_linestring_map_.insert({ e2.id, e2.line_string });
			inner_id = inner_id + 1;
		}
		OGRFeature::DestroyFeature(ogrFeature);
	}
	GDALClose(poDS);

	std::cout << "Read --directed-- network finish." << '\n';
	std::cout << "\tTotal number of edges read " << network_edges_.size() << '\n';
}

//------------read the facilities points-------------------------------
void Network::ReadFacilities(const std::string filename, std::string id_name, std::string delta_name)
{
	std::cout << "Read the facility list" << '\n';
	GDALDataset *poDS1 = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(),GDAL_OF_VECTOR, NULL, NULL, NULL));
	OGRLayer *polayer = poDS1->GetLayer(0);
	int num_points = polayer->GetFeatureCount();
	std::cout << "\tNumber of points in file: " << num_points << '\n';

	deltas_ = std::vector<double>(num_points);
	facilities_ = std::vector<Facility>(num_points);

	OGRFeatureDefn *poiFDefn = polayer->GetLayerDefn();
	OGRFeature *poiFeature;

	int id_idx = poiFDefn->GetFieldIndex(id_name.c_str());
	int delat_idx = poiFDefn->GetFieldIndex(delta_name.c_str());

	while ((poiFeature = polayer->GetNextFeature()) != NULL)
	{
		int id = poiFeature->GetFID();
		Facility *temp_facility = &facilities_[id];		
		temp_facility->id = id;
		temp_facility->external_id = std::string(poiFeature->GetFieldAsString(id_idx));
		temp_facility->delta = poiFeature->GetFieldAsDouble(delat_idx);
		OGRGeometry *rawpointgeometry = poiFeature->GetGeometryRef();
		temp_facility->point = (OGRPoint*)rawpointgeometry->clone();

		deltas_[id] = poiFeature->GetFieldAsDouble(delat_idx);
	}
	GDALClose(poDS1);

	std::cout << "finish reading the facilities list" << '\n';
	std::cout << "\tTotal number of facilities " << num_points << '\n';
}


//------------read the facilities points-------------------------------
void Network::ReadFacilities(const std::string &filename, const std::string &id_name, const bool &isunifiedcutoff,
	const double &unifiedcutoff, const std::string &cutoff_name)
{
	std::cout << "Read the facility list" << '\n';
	GDALDataset *poDS1 = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	OGRLayer *polayer = poDS1->GetLayer(0);
	int num_points = polayer->GetFeatureCount();
	std::cout << "\tNumber of points in file: " << num_points << '\n';

	deltas_ = std::vector<double>(num_points);
	facilities_ = std::vector<Facility>(num_points);

	OGRFeatureDefn *poiFDefn = polayer->GetLayerDefn();
	OGRFeature *poiFeature;

	int id_idx = poiFDefn->GetFieldIndex(id_name.c_str());
	
	if (isunifiedcutoff) // is unified cut-off cost, then the "unifiedcutoff" is used 
	{
		while ((poiFeature = polayer->GetNextFeature()) != NULL)
		{
			int id = poiFeature->GetFID();
			Facility *temp_facility = &facilities_[id];
			temp_facility->id = id;
			temp_facility->external_id = std::string(poiFeature->GetFieldAsString(id_idx));
			temp_facility->delta = unifiedcutoff;
			OGRGeometry *rawpointgeometry = poiFeature->GetGeometryRef();
			temp_facility->point = (OGRPoint*)rawpointgeometry->clone();
			deltas_[id] = unifiedcutoff;
		}
	}
	else // is individual cut-off cost, then the cut-off is read from file based on the provided "cutoff_name"
	{
		int delat_idx = poiFDefn->GetFieldIndex(cutoff_name.c_str());
		while ((poiFeature = polayer->GetNextFeature()) != NULL)
		{
			int id = poiFeature->GetFID();
			Facility *temp_facility = &facilities_[id];
			temp_facility->id = id;
			temp_facility->external_id = std::string(poiFeature->GetFieldAsString(id_idx));
			temp_facility->delta = poiFeature->GetFieldAsDouble(delat_idx);
			OGRGeometry *rawpointgeometry = poiFeature->GetGeometryRef();
			temp_facility->point = (OGRPoint*)rawpointgeometry->clone();
			deltas_[id] = poiFeature->GetFieldAsDouble(delat_idx);
		}		
	}
	GDALClose(poDS1);	
	std::cout << "finish reading the facilities list" << '\n';
	std::cout << "\tTotal number of facilities " << num_points << '\n';
}


void Network::ReadMpFacilities(const std::string filename, std::string id_name, std::string delta_name , std::string label_name)
{
	std::cout << "Read the facility list" << '\n';
	GDALDataset *poDS1 = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	OGRLayer *polayer = poDS1->GetLayer(0);
	int num_points = polayer->GetFeatureCount();
	std::cout << "\tNumber of points in file: " << num_points << '\n';

	deltas_ = std::vector<double>(num_points);
	// for multiple-point facility, a label is needed for a point to mark its facility label
	labels_ = std::vector<std::string>(num_points);
	facilities_ = std::vector<Facility>(num_points);

	OGRFeatureDefn *poiFDefn = polayer->GetLayerDefn();
	OGRFeature *poiFeature;

	int id_idx = poiFDefn->GetFieldIndex(id_name.c_str());
	int delat_idx = poiFDefn->GetFieldIndex(delta_name.c_str());
	int label_idx = poiFDefn->GetFieldIndex(label_name.c_str());

	
	while ((poiFeature = polayer->GetNextFeature()) != NULL)
	{
		int id = poiFeature->GetFID();
		Facility *temp_facility = &facilities_[id];
		temp_facility->id = id;
		temp_facility->external_id = std::string(poiFeature->GetFieldAsString(id_idx));
		temp_facility->delta = poiFeature->GetFieldAsDouble(delat_idx);
		OGRGeometry *rawpointgeometry = poiFeature->GetGeometryRef();
		temp_facility->point = (OGRPoint*)rawpointgeometry->clone();

		deltas_[id] = poiFeature->GetFieldAsDouble(delat_idx);
		// for multiple-point facility, a label is needed for a point to mark its facility label
		labels_[id] = poiFeature->GetFieldAsString(label_idx);
	}
	GDALClose(poDS1);

	std::cout << "finish reading the facilities list" << '\n';
	std::cout << "\tTotal number of facilities " << num_points << '\n';
}



void Network::ReadMpFacilities(const std::string filename, const std::string &id_name, std::string label_name, const bool &isunifiedcutoff,
	const double &unifiedcutoff, const std::string &cutoff_name)
{
	std::cout << "Read the facility list" << '\n';
	GDALDataset *poDS1 = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	OGRLayer *polayer = poDS1->GetLayer(0);
	int num_points = polayer->GetFeatureCount();
	std::cout << "\tNumber of points in file: " << num_points << '\n';

	deltas_ = std::vector<double>(num_points);
	// for multiple-point facility, a label is needed for a point to mark its facility label
	labels_ = std::vector<std::string>(num_points);
	facilities_ = std::vector<Facility>(num_points);
	
	if (isunifiedcutoff) 
	{
		OGRFeatureDefn *poiFDefn = polayer->GetLayerDefn();
		OGRFeature *poiFeature;
		int id_idx = poiFDefn->GetFieldIndex(id_name.c_str());
		int label_idx = poiFDefn->GetFieldIndex(label_name.c_str());

		while ((poiFeature = polayer->GetNextFeature()) != NULL)
		{
			int id = poiFeature->GetFID();
			Facility *temp_facility = &facilities_[id];
			temp_facility->id = id;
			temp_facility->external_id = std::string(poiFeature->GetFieldAsString(id_idx));
			OGRGeometry *rawpointgeometry = poiFeature->GetGeometryRef();
			temp_facility->point = (OGRPoint*)rawpointgeometry->clone();
			labels_[id] = poiFeature->GetFieldAsString(label_idx);
			temp_facility->delta = unifiedcutoff;
			deltas_[id] = unifiedcutoff;
		}
	}
	else 
	{
		OGRFeatureDefn *poiFDefn = polayer->GetLayerDefn();
		OGRFeature *poiFeature;
		int id_idx = poiFDefn->GetFieldIndex(id_name.c_str());
		int label_idx = poiFDefn->GetFieldIndex(label_name.c_str());		
		int delat_idx = poiFDefn->GetFieldIndex(cutoff_name.c_str());// add an addtional cutoff idx
		while ((poiFeature = polayer->GetNextFeature()) != NULL)
		{
			int id = poiFeature->GetFID();
			Facility *temp_facility = &facilities_[id];
			temp_facility->id = id;
			temp_facility->external_id = std::string(poiFeature->GetFieldAsString(id_idx));			
			OGRGeometry *rawpointgeometry = poiFeature->GetGeometryRef();
			temp_facility->point = (OGRPoint*)rawpointgeometry->clone();
			labels_[id] = poiFeature->GetFieldAsString(label_idx);			
			temp_facility->delta = poiFeature->GetFieldAsDouble(delat_idx);
			deltas_[id] = poiFeature->GetFieldAsDouble(delat_idx);
			
		}
	}

	GDALClose(poDS1);
	std::cout << "finish reading the facilities list" << '\n';
	std::cout << "\tTotal number of facilities " << num_points << '\n';
}



std::vector<double> Network::get_deltas(void)
{
	return deltas_;
}


std::vector<std::string> Network::get_labels(void)
{
	return labels_;
}

// for directed road network input
void Network::BuildRtreeIndexDirected()
{
	std::cout << "Start to construct boost rtree for --directed-- network" << '\n';
	for (std::size_t i = 0; i < network_directed_edges_.size(); ++i)
	{
		EdgeDirected *edge = &network_directed_edges_[i];
		double x1, y1, x2, y2;

		algorithm::GetLinestringBoundingbox(edge->line_string, &x1, &y1, &x2, &y2);
		boost_box b(boost_point(x1, y1), boost_point(x2, y2));
		rtree_directed_.insert(std::make_pair(b, edge));
	}
	std::cout << "Finish construct boost rtree for --directed-- network" << '\n';
}

// for directed road network input
NearestPointsDirected Network::SearchNearestPointsDirected(double radius)
{
	int num_facilities = facilities_.size();
	NearestPointsDirected ficility_nearest_points(num_facilities);

	for (int i = 0; i < num_facilities; ++i)
	{
		std::vector<NearestPointDirected> local_nearest_points;

		OGRPoint* point = facilities_[i].point;
		boost_point pt = boost_point(point->getX(), point->getY());

		std::vector<ItemDirected> temp;
		boost_box b(boost_point(point->getX() - radius, point->getY() - radius), boost_point(point->getX() + radius, point->getY() + radius));
		rtree_directed_.query(bgi::intersects(b), std::back_inserter(temp));

		for (ItemDirected &j : temp)  // Rtree can only detect intersect with a the bounding box of the geometry stored.
		{
			EdgeDirected *edge = j.second;
			float break_distance;
			double virtical_distance;
			int cut_seg_index;

			OGRPoint* p_break_point = new OGRPoint;
			// here the offset distance is measured to the starting node
			algorithm::ProjectPointOnLinestring(point, edge->line_string, &cut_seg_index, &virtical_distance, &break_distance);
			*p_break_point = algorithm::GetBreakPointOnLinestring(break_distance, edge->line_string);
			NearestPointDirected np = NearestPointDirected{ i,virtical_distance,break_distance,cut_seg_index, p_break_point, edge };
			local_nearest_points.push_back(np);
		}

		// below are used to fetch the point with smallest distance to the facility
		if (local_nearest_points.empty())
		{
			return NearestPointsDirected();
		};
		if (local_nearest_points.size() == 1)
		{
			ficility_nearest_points[i] = local_nearest_points[0];
		}
		else
		{
			std::sort(local_nearest_points.begin(), local_nearest_points.end(), compareBydisdirected);
			ficility_nearest_points[i] = local_nearest_points[0];
		}
	}
	return ficility_nearest_points;
}


std::vector<EdgeDirected> Network::get_edges_directed(void) 
{
	return network_directed_edges_;
}

std::vector<std::vector<EdgeDirected*>> Network::GetEdgesAroundDeltaDirected(double delta, NearestPointsDirected facility_nearest_points)
{
	int num_facilities = facility_nearest_points.size();
	std::vector<std::vector<EdgeDirected*>> facility_sub_networks(num_facilities);

	for (int i = 0; i < num_facilities; ++i)
	{
		std::vector<EdgeDirected*> temp_sub_network;
		// cut-off point
		OGRPoint* cut_point = facility_nearest_points[i].point;

		std::vector<ItemDirected> temp;
		// build a boudning box based on delata and 
		boost_box b(boost_point(cut_point->getX() - delta, cut_point->getY() - delta), boost_point(cut_point->getX() + delta, cut_point->getY() + delta));

		rtree_directed_.query(bgi::intersects(b), std::back_inserter(temp));

		// Rtree can only detect intersect with a the bounding box of the geometry stored.
		for (ItemDirected &j : temp)
		{
			EdgeDirected *edge = j.second;
			temp_sub_network.push_back(edge);
		}
		facility_sub_networks[i] = temp_sub_network;
	}
	return facility_sub_networks;
	
}


std::vector<std::vector<EdgeDirected*>> Network::GetEdgesAroundDeltaDirected(NearestPointsDirected facility_nearest_points)
{
	int num_facilities = facility_nearest_points.size();
	std::vector<std::vector<EdgeDirected*>> facility_sub_networks(num_facilities);
	
	for (int i = 0; i < num_facilities; ++i)
	{
		double delta = facilities_[i].delta;
		std::vector<EdgeDirected*> temp_sub_network;
		// cut-off point
		OGRPoint* cut_point = facility_nearest_points[i].point;

		std::vector<ItemDirected> temp;
		// build a boudning box based on delata and 
		boost_box b(boost_point(cut_point->getX() - delta, cut_point->getY() - delta), boost_point(cut_point->getX() + delta, cut_point->getY() + delta));
		rtree_directed_.query(bgi::intersects(b), std::back_inserter(temp));

		// Rtree can only detect intersect with a the bounding box of the geometry stored.
		for (ItemDirected &j : temp)
		{
			EdgeDirected *edge = j.second;
			temp_sub_network.push_back(edge);
		}
		facility_sub_networks[i] = temp_sub_network;
	}
	return facility_sub_networks;
}


void Network::write_sub_edges(std::string filename, std::vector<Edge*> sub_edges)
{
	//create ESRI shp
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("sub_edges", NULL, wkbLineString, NULL);

	// define and create filed
	OGRFieldDefn oFie_edge_Id("edge_id", OFTString);
	oFie_edge_Id.SetWidth(30);
	poLayer->CreateField(&oFie_edge_Id);

	OGRFieldDefn oField_s_node("s_node", OFTString);
	oField_s_node.SetWidth(30);
	poLayer->CreateField(&oField_s_node);

	OGRFieldDefn oField_t_node("e_node", OFTString);
	oField_t_node.SetWidth(30);
	poLayer->CreateField(&oField_t_node);


	for (int i = 1; i < sub_edges.size(); i++)
	{
		std::string ex_id = sub_edges[i]->external_id;
		std::string source = sub_edges[i]->source;
		std::string target = sub_edges[i]->target;
		
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		poFeature->SetField("edge_id", ex_id.c_str());
		poFeature->SetField("s_node", source.c_str());
		poFeature->SetField("e_node", target.c_str());
		poFeature->SetGeometry(sub_edges[i]->line_string);

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(shpDataSet);

}


/*-----------below functions are used for generating the benchmark, i.e., the accessible grid points--------*/

void Network::BenchmarkGeneration(std::string fsub_graph, std::string fgrid_points, std::string fup_grid_points, double delta, OGRPolygon* arc_gis_CA, OGRPolygon* our_CA)
{
	std::vector<SubGraphEdge> sub_edges = ReadSubGraph(fsub_graph);
	std::vector<OGRPoint> grid_points = ReadGridPoints(fgrid_points);
	std::vector<double>  pdistances = DistanceToRootNode(sub_edges, grid_points);
	write_accessible_grid_points(fup_grid_points, grid_points, pdistances, delta, arc_gis_CA, our_CA);
}


std::vector<SubGraphEdge> Network::ReadSubGraph(const std::string filename, const std::string &id_name, const std::string &source_name,
						   const std::string &target_name, const std::string &s_distance, const std::string &e_distance)
{
	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	if (poDS == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}
	OGRLayer *ogrlayer = poDS->GetLayer(0);

	// Initialize network edges
	OGRFeatureDefn *ogrFDefn = ogrlayer->GetLayerDefn();
	OGRFeature *ogrFeature;
	int num_features = ogrlayer->GetFeatureCount();
	std::vector<SubGraphEdge> subgraph_edges(num_features);

	// Fetch the field index given field name.
	int id_idx = ogrFDefn->GetFieldIndex(id_name.c_str());
	int source_idx = ogrFDefn->GetFieldIndex(source_name.c_str());
	int target_idx = ogrFDefn->GetFieldIndex(target_name.c_str());
	int s_dis_idx = ogrFDefn->GetFieldIndex(s_distance.c_str());
	int e_dis_idx = ogrFDefn->GetFieldIndex(e_distance.c_str());

	while ((ogrFeature = ogrlayer->GetNextFeature()) != NULL)
	{
		int id = ogrFeature->GetFID();
		SubGraphEdge e;	
		e.id = id;
		//e.external_id = std::string(ogrFeature->GetFieldAsString(id_idx));
		//e.source = std::string(ogrFeature->GetFieldAsString(source_idx));
		//e.target = std::string(ogrFeature->GetFieldAsString(target_idx));
		e.s_dis = ogrFeature->GetFieldAsDouble(s_dis_idx);
		//e.e_dis = std::double(ogrFeature->GetFieldAsString(e_dis_idx));
		OGRGeometry *rawlinegeometry = ogrFeature->GetGeometryRef();		
		e.line_string = (OGRLineString*)rawlinegeometry->clone();
		//e.length = e.line_string->get_Length();
		subgraph_edges[id] = e;
		//edge_linestring_map_.insert({e.id, e.line_string }); //edge linestring map, id is a inner id and continous
		OGRFeature::DestroyFeature(ogrFeature);
	}
	GDALClose(poDS);

	std::cout << "Read network finish." << '\n';
	std::cout << "\tTotal number of edges read " << num_features << '\n';

	return subgraph_edges;
}

std::vector<OGRPolygon*>  Network::ReadCatchmentAreas(const std::string filename) 
{
	GDALAllRegister();
	std::cout << "Read the catchment list" << '\n';
	GDALDataset *poDS1 = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	OGRLayer *polayer = poDS1->GetLayer(0);
	int num_CAs = polayer->GetFeatureCount();

	std::vector<OGRPolygon*> CAs(num_CAs);

	OGRFeatureDefn *poiFDefn = polayer->GetLayerDefn();
	OGRFeature *poiFeature;

	while ((poiFeature = polayer->GetNextFeature()) != NULL)
	{
		int id = poiFeature->GetFID();
		OGRGeometry *rawpointgeometry = poiFeature->GetGeometryRef();
		CAs[id] = (OGRPolygon*)rawpointgeometry->clone();
	}
	GDALClose(poDS1);
	std::cout << "finish reading the catchment area" << '\n';
	std::cout << "\tTotal number of  catchment area " << num_CAs << '\n';

	return CAs;
}

std::vector<OGRPoint>  Network::ReadGridPoints(const std::string filename, const std::string &id_name)
{
	std::cout << "Read the grid point list" << '\n';
	GDALDataset *poDS1 = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	OGRLayer *polayer = poDS1->GetLayer(0);
	int num_points = polayer->GetFeatureCount();

	std::vector<OGRPoint> grid_points(num_points);

	OGRFeatureDefn *poiFDefn = polayer->GetLayerDefn();
	OGRFeature *poiFeature;

	while ((poiFeature = polayer->GetNextFeature()) != NULL)
	{
		int id = poiFeature->GetFID();
		OGRGeometry *rawpointgeometry = poiFeature->GetGeometryRef();
		//grid_points[id] = (OGRPoint*)rawpointgeometry->clone();
		grid_points[id] = *(OGRPoint*)rawpointgeometry;
	}
	GDALClose(poDS1);
	std::cout << "finish reading the grid point list" << '\n';
	std::cout << "\tTotal number of  grid point " << num_points << '\n';

	return grid_points;
}

std::vector<double>  Network::DistanceToRootNode(std::vector<SubGraphEdge> sub_edges, std::vector<OGRPoint> grid_points)
{
	int num_grid_points = grid_points.size();
	std::vector<double> distances(num_grid_points);

	for (int i = 0; i < num_grid_points; ++i)
	{
		OGRPoint point = grid_points[i];
		double distance_by_nearest_edge = 0;
		float break_distance;
		double virtical_distance;
		int cut_seg_index;
		double nearest_vertical_distance = 9999999999999999;
		
		for (int j = 0; j < sub_edges.size(); ++j)
		{
			SubGraphEdge edge = sub_edges[j];
			OGRPoint* p_break_point = new OGRPoint;
			// here the offset distance is measured to the starting node
			algorithm::ProjectPointOnLinestring(&point, edge.line_string, &cut_seg_index, &virtical_distance, &break_distance);
			if (virtical_distance < nearest_vertical_distance)
			{
				nearest_vertical_distance = virtical_distance;
				distance_by_nearest_edge = edge.s_dis + break_distance + nearest_vertical_distance;
			}
		}
		distances[i] = distance_by_nearest_edge;
	}
	return distances;
}

void Network::write_accessible_grid_points(std::string filename, std::vector<OGRPoint> grid_points, std::vector<double> distances, 
	const double delta, OGRPolygon* arc_gis_CA, OGRPolygon* our_CA)
{
	//创建ESRI shp文件
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("sub_edges", NULL, wkbPoint, NULL);

	// define and create filed
	OGRFieldDefn oFie_point_Id("id", OFTString);
	oFie_point_Id.SetWidth(10);
	poLayer->CreateField(&oFie_point_Id);

	OGRFieldDefn oField_dis ("distance", OFTReal);
	oField_dis.SetWidth(10);
	oField_dis.SetPrecision(5);
	poLayer->CreateField(&oField_dis);

	OGRFieldDefn oField_acclabel("accLabel", OFTInteger);
	oField_acclabel.SetWidth(2);
	poLayer->CreateField(&oField_acclabel);

	OGRFieldDefn oField_inacrgis("InArcgis", OFTInteger);
	oField_inacrgis.SetWidth(2);
	poLayer->CreateField(&oField_inacrgis);

	OGRFieldDefn oField_inours("InOurs", OFTInteger);
	oField_inours.SetWidth(2);
	poLayer->CreateField(&oField_inours);

	
	for (int i = 1; i < grid_points.size(); i++)
	{
		int acc_label = 1;
		int within_acrgis = 0;
		int within_ours = 0;
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		poFeature->SetField("id", std::to_string(i).c_str());
		poFeature->SetField("distance", distances[i]);
		if (distances[i] > delta) acc_label = 0;
		poFeature->SetField("accLabel", acc_label);
		
		//if (grid_points[i]->Within(arc_gis_CA)) within_acrgis = 1;
		//if (grid_points[i]->Within(our_CA)) within_ours = 1;
		if (grid_points[i].Within(arc_gis_CA)) within_acrgis = 1;
		if (grid_points[i].Within(our_CA)) within_ours = 1;

		poFeature->SetField("InArcgis", within_acrgis);
		poFeature->SetField("InOurs", within_ours);
		//poFeature->SetGeometry(grid_points[i]);
		poFeature->SetGeometry(&grid_points[i]);

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(shpDataSet);
}


}// end for namespace tca