#include "tca_build_triangulation.h"

namespace tca 
{

BuildTriangulation::BuildTriangulation()
{
	
}

BuildTriangulation::BuildTriangulation(const vector<string> node_ids, const vector<double> node_distances, const vector<OGRPoint> node_org_points)
	:node_external_ids_(node_ids), z_values_(node_distances)
{
	x_coords_ = std::vector<double>(node_ids.size());
	y_coords_ = std::vector<double>(node_ids.size());
	node_points_inner_ids_pair_ = std::vector< std::pair<Point, unsigned> >(node_ids.size());

	for (int i = 0; i < node_ids.size(); ++i) 
	{		
		x_coords_[i] = node_org_points[i].getX();
		y_coords_[i] = node_org_points[i].getY();
		Point cgal_point(node_org_points[i].getX(), node_org_points[i].getY());
		node_points_inner_ids_pair_[i] = std::make_pair(cgal_point, i);
	}
}

BuildTriangulation::~BuildTriangulation()
{

}

void BuildTriangulation::ReadNodeCosts(const std::string & file, const std::string &id_name, const std::string & zvalue)
{
	GDALAllRegister();
	GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(file.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
	if (poDS == NULL)
	{
		printf("Open failed.\n");
		exit(1);
	}
	
	OGRLayer *ogrlayer = poDS->GetLayer(0);
	int num_nodes = ogrlayer->GetFeatureCount();
	std::cout << "\tNumber of nodes for generating the triangulation: " << num_nodes << '\n';

	x_coords_ = std::vector<double>(num_nodes);
	y_coords_ = std::vector<double>(num_nodes);
	z_values_ = std::vector<double>(num_nodes);
	node_points_inner_ids_pair_ = std::vector<std::pair<Point, unsigned>>(num_nodes);
	node_external_ids_ = std::vector<std::string>(num_nodes);

	OGRFeatureDefn *ogrFDefn = ogrlayer->GetLayerDefn();
	OGRFeature *ogrFeature;

	// Fetch the field index given field name.
	int id_idx = ogrFDefn->GetFieldIndex(id_name.c_str());
	int zvalue_idx = ogrFDefn->GetFieldIndex(zvalue.c_str());

	std::string temp_id;
	double temp_z;
	int id;

	while ((ogrFeature = ogrlayer->GetNextFeature()) != NULL)
	{
		id = ogrFeature->GetFID();
		temp_id = std::string(ogrFeature->GetFieldAsString(id_idx));
		temp_z = (ogrFeature->GetFieldAsDouble(zvalue_idx));
		OGRGeometry *rawpointgeometry = ogrFeature->GetGeometryRef();
		Point cgal_point(((OGRPoint*)rawpointgeometry->clone())->getX(), ((OGRPoint*)rawpointgeometry->clone())->getY());

		x_coords_[id] = ((OGRPoint*)rawpointgeometry->clone())->getX();
		y_coords_[id] = ((OGRPoint*)rawpointgeometry->clone())->getY();
		z_values_[id] = temp_z;
		node_external_ids_[id] = temp_id;			
		node_points_inner_ids_pair_[id] = std::make_pair(cgal_point, id);

		OGRFeature::DestroyFeature(ogrFeature);
	}
	GDALClose(poDS);
}// ReadNodeCosts


//https://doc.cgal.org/latest/Triangulation_2/Triangulation_2_2info_insert_with_pair_iterator_2_8cpp-example.html#_a1
void BuildTriangulation::BuildDelaunayTrangulation(void)
{
	delaunary_.insert(node_points_inner_ids_pair_.begin(), node_points_inner_ids_pair_.end());

	for (Delaunay::Finite_faces_iterator fit = delaunary_.finite_faces_begin();fit != delaunary_.finite_faces_end(); ++fit) 
	{
		Delaunay::Face_handle face = fit;
		std::array<int, 3> temp_tr = { face->vertex(0)->info(),face->vertex(1)->info(),face->vertex(2)->info() };
		triangle_node_ids_.push_back(temp_tr);

		OGRLinearRing *ring = new OGRLinearRing();
		OGRPolygon *poly = new OGRPolygon();
		ring->addPoint(delaunary_.triangle(face)[0].x(), delaunary_.triangle(face)[0].y());
		ring->addPoint(delaunary_.triangle(face)[1].x(), delaunary_.triangle(face)[1].y());
		ring->addPoint(delaunary_.triangle(face)[2].x(), delaunary_.triangle(face)[2].y());
		ring->addPoint(delaunary_.triangle(face)[0].x(), delaunary_.triangle(face)[0].y());
		//ring->closeRings();
		poly->addRing(ring);
		tri_polygons_.push_back(*poly);
	}
}// BuildDelaunayTrangulation

void BuildTriangulation::BuildConstrainedTrangulation(std::vector<std::array<int, 2>> edge_node_ids, const vector<string> node_ids,
	const vector<double> node_distances, const vector<OGRPoint> node_org_points)
{
	std::map<PointM, Vertex_handle> points_handles_map;
	node_external_ids_ = node_ids;
	z_values_ = node_distances;
	x_coords_ = std::vector<double>(node_ids.size());
	y_coords_ = std::vector<double>(node_ids.size());

	for (int i = 0; i < node_ids.size(); ++i)
	{
		x_coords_[i] = node_org_points[i].getX();
		y_coords_[i] = node_org_points[i].getY();
		if (node_distances[i] != std::numeric_limits<double>::max()) // the nodes not in the extended sp tree are exclude
		{
			Point cgal_point(node_org_points[i].getX(), node_org_points[i].getY());
			node_points_inner_ids_pair_.push_back(std::make_pair(cgal_point, i));
		}
	}

	// build delaunary trigulation
	con_delau_tris_.insert(node_points_inner_ids_pair_.begin(), node_points_inner_ids_pair_.end());

	CDT::Finite_vertices_iterator vit;
	for (vit = con_delau_tris_.finite_vertices_begin(); vit != con_delau_tris_.finite_vertices_end(); ++vit)
	{
		Vertex_handle handle = vit;
		inner_id_vh.insert({ handle->info(), handle });
		// since there are some duplicated points without handle, we thereby add a map to 
		// index the handle of the vertex with same point
		points_handles_map.insert({ PointM(handle->point()),handle });
	}
	//  the number of handle is smaller than the number of node_points_inner_ids_pair_
	//   e.g., 710 VS 873 VS, because of two reasons 
	//  1) for undirected graph ,the break points have been stored two times in the extend SP tree
	//  2) if the distance between two nodes are too small the cgal will automatically keep one of them (e.g., regard the two points as the same)
	//  Therefore, we build a "points_handles_map" to solve the first problem,
	//  for the second problem, it is unsolved currently
	Vertex_handle s_node_handle, e_node_handle;
	int edge_size = edge_node_ids.size();
	for (int i = 0; i < edge_size; i++)
	{
		bool insert_con_seg = true;

		int s_node_id = edge_node_ids[i][0];
		int e_node_id = edge_node_ids[i][1];

		auto search_s = inner_id_vh.find(s_node_id);
		auto search_e = inner_id_vh.find(e_node_id);

		if ((search_s != inner_id_vh.end()) && (search_e != inner_id_vh.end()))
		{
			s_node_handle = search_s->second;
			e_node_handle = search_e->second;
		}
		else // this is added to confirm that all 
		{
			if (search_s == inner_id_vh.end())
			{
				Point cgal_point_s(node_org_points[s_node_id].getX(), node_org_points[s_node_id].getY());
				auto ser_s = points_handles_map.find(PointM(cgal_point_s));
				if (ser_s != points_handles_map.end())// here, we can make sure that this handle can be found
				{
					s_node_handle = ser_s->second;
				}
				else 
				{
					insert_con_seg = false; // sometimes if the distance between two points is too small (but not excatly same), 
											 //one of them will be ignored during the construction of denaulary 
				}			
			}
			else 
			{
				s_node_handle = search_s->second;
			}

			if (search_e == inner_id_vh.end())
			{
				Point cgal_point_e(node_org_points[e_node_id].getX(), node_org_points[e_node_id].getY());
				auto ser_e = points_handles_map.find(PointM(cgal_point_e));
				if (ser_e != points_handles_map.end()) 
				{ 
					e_node_handle = ser_e->second; 
				}	// here, we can make sure that this a handle can be found
				else 
				{
					insert_con_seg = false;
				}
			}
			else 
			{ 
				e_node_handle = search_e->second; 
			}
		}
		if (insert_con_seg) // only if the handles of the nodes for constrained segs found, we add this as constrians
		{
			OGRPoint * spoint = new OGRPoint();
			spoint->setX(s_node_handle->point().x());
			spoint->setY(s_node_handle->point().y());

			OGRPoint * epoint = new OGRPoint();
			epoint->setX(e_node_handle->point().x());
			epoint->setY(e_node_handle->point().y());

			if (*epoint != *spoint) // when two points are same a error CGAL error: precondition violation!
			 //Expression : vaa != vbb will be reported
			{
				con_delau_tris_.insert_constraint(s_node_handle, e_node_handle);// insert the constrained graph edges
			}

			OGRLineString * seg = new OGRLineString();

			seg->addPoint(spoint);
			seg->addPoint(epoint);
			constrained_segs.push_back(seg);

			// below are for test
			int s_node = s_node_handle->info();
			int e_node = e_node_handle->info();
			constrained_seg_egde_ids.push_back(i);
			constrained_seg_snode_ids.push_back(s_node);
			constrained_seg_enode_ids.push_back(e_node);
		}
		//else // the handle of the edge is not found
		//{
		//	cout << "stop here for debug" << endl;
		//}
	}

	// when two constrianed edges intersect with each other
	// the intersection points will be added automatically during the constrained triangulation
	// under such case, the error "point index out of bounds" will occur during the contour generation,as following
	// Assertion failed: point >= 0 && point < _triangulation.get_npoints() && "Point index out of bounds."
	// ...\tca_tri_contour.cpp, line 1067
	// Therefore, we filter triangles including such newly added intersection points
	for (CDT::Finite_faces_iterator fit = con_delau_tris_.finite_faces_begin(); fit != con_delau_tris_.finite_faces_end(); ++fit)
	{
		CDT::Face_handle face = fit;
		int verx_1 = face->vertex(0)->info();
		int verx_2 = face->vertex(1)->info();
		int verx_3 = face->vertex(2)->info();	
		int max_nodeID = node_ids.size();
		if (verx_1 >= 0 && verx_1 < max_nodeID && verx_2 >= 0 && verx_2 < max_nodeID && verx_3 >= 0 && verx_3 < max_nodeID )
		{
			std::array<int, 3> temp_tr = {verx_1,verx_2,verx_3 };
			triangle_node_ids_.push_back(temp_tr);
			OGRLinearRing *ring = new OGRLinearRing();
			OGRPolygon *poly = new OGRPolygon();
			ring->addPoint(con_delau_tris_.triangle(face)[0].x(), con_delau_tris_.triangle(face)[0].y());
			ring->addPoint(con_delau_tris_.triangle(face)[1].x(), con_delau_tris_.triangle(face)[1].y());
			ring->addPoint(con_delau_tris_.triangle(face)[2].x(), con_delau_tris_.triangle(face)[2].y());
			ring->addPoint(con_delau_tris_.triangle(face)[0].x(), con_delau_tris_.triangle(face)[0].y());
			poly->addRing(ring);
			tri_polygons_.push_back(*poly);
		}
	}
}// BuildDelaunayTrangulation


// this is for multiple point version, since we added a "virtual node" during the construction
// this point should not be included as an input point for triangulation
// correspondingly, the edge between the projected points and virtual node should not be added as the constrained segments.
void BuildTriangulation::BuildConstrainedTrangulationMultiplePointsVersion(std::vector<std::array<int, 2>> edge_node_ids, const vector<string> node_ids,
	const vector<double> node_distances, const vector<OGRPoint> node_org_points)
{
	std::map<PointM, Vertex_handle> points_handles_map;
	node_external_ids_ = node_ids;
	z_values_ = node_distances;

	x_coords_ = std::vector<double>(node_ids.size());
	y_coords_ = std::vector<double>(node_ids.size());

	for (int i = 0; i < node_ids.size(); ++i)
	{
		x_coords_[i] = node_org_points[i].getX();
		y_coords_[i] = node_org_points[i].getY();
		if (node_distances[i] != std::numeric_limits<double>::max()) // the nodes not in the extended sp tree are excluded
		{
			if (!((node_distances[i] == 0) && (node_ids[i] == "vir_node"))) // excluding the virtual node
			{
				Point cgal_point(node_org_points[i].getX(), node_org_points[i].getY());
				node_points_inner_ids_pair_.push_back(std::make_pair(cgal_point, i));
			}
		}
	}

	// build delaunary trigulation
	con_delau_tris_.insert(node_points_inner_ids_pair_.begin(), node_points_inner_ids_pair_.end());

	CDT::Finite_vertices_iterator vit;
	for (vit = con_delau_tris_.finite_vertices_begin(); vit != con_delau_tris_.finite_vertices_end(); ++vit)
	{
		Vertex_handle handle = vit;
		inner_id_vh.insert({ handle->info(), handle });
		// since there are some duplicated points without handle, we thereby add a map to 
		// index the handle of the vertex with same point
		points_handles_map.insert({ PointM(handle->point()),handle });
	}

	Vertex_handle s_node_handle, e_node_handle;
	int edge_size = edge_node_ids.size();
	for (int i = 0; i < edge_size; i++)
	{
		bool insert_con_seg = true;

		int s_node_id = edge_node_ids[i][0];
		int e_node_id = edge_node_ids[i][1];

		auto search_s = inner_id_vh.find(s_node_id);
		auto search_e = inner_id_vh.find(e_node_id);
		//bool edge_handle_found = true;

		if ((search_s != inner_id_vh.end()) && (search_e != inner_id_vh.end()))
		{
			s_node_handle = search_s->second;
			e_node_handle = search_e->second;
		}
		else
		{
			if (search_s == inner_id_vh.end())
			{
				Point cgal_point_s(node_org_points[s_node_id].getX(), node_org_points[s_node_id].getY());
				auto ser_s = points_handles_map.find(PointM(cgal_point_s));

				// different from the version of signle point facilty, we can NOT make sure that this a handle can be found here
				if (ser_s != points_handles_map.end()) s_node_handle = ser_s->second;
				else insert_con_seg = false;
			}
			else s_node_handle = search_s->second;

			if (search_e == inner_id_vh.end())
			{
				Point cgal_point_e(node_org_points[e_node_id].getX(), node_org_points[e_node_id].getY());
				auto ser_e = points_handles_map.find(PointM(cgal_point_e));
				// different from the version of signle point facilty, we can NOT make sure that this a handle can be found here
				if (ser_e != points_handles_map.end())e_node_handle = ser_e->second;
				else insert_con_seg = false;
			}
			else e_node_handle = search_e->second;
		}
		// only under the condition that the handle of an edge is succesfully found
		if (insert_con_seg)
		{
			OGRPoint * spoint = new OGRPoint();
			spoint->setX(s_node_handle->point().x());
			spoint->setY(s_node_handle->point().y());

			OGRPoint * epoint = new OGRPoint();
			epoint->setX(e_node_handle->point().x());
			epoint->setY(e_node_handle->point().y());

			if (*epoint != *spoint) // this is added during the test of directed edge
			{
				con_delau_tris_.insert_constraint(s_node_handle, e_node_handle);// insert the constrained graph edges
			}

			OGRLineString * seg = new OGRLineString();
			seg->addPoint(spoint);
			seg->addPoint(epoint);
			constrained_segs.push_back(seg);

			// below are for test
			int s_node = s_node_handle->info();
			int e_node = e_node_handle->info();

			constrained_seg_egde_ids.push_back(i);
			constrained_seg_snode_ids.push_back(s_node);
			constrained_seg_enode_ids.push_back(e_node);
		}
	}

	for (CDT::Finite_faces_iterator fit = con_delau_tris_.finite_faces_begin(); fit != con_delau_tris_.finite_faces_end(); ++fit)
	{
		CDT::Face_handle face = fit;

		int verx_1 = face->vertex(0)->info();
		int verx_2 = face->vertex(1)->info();
		int verx_3 = face->vertex(2)->info();
		//if (verx_1 >= 0 and verx_2 >= 0 and verx_3 >= 0)
		if (verx_1 >= 0 && verx_2 >= 0 && verx_3 >= 0)
		{
			std::array<int, 3> temp_tr = { verx_1,verx_2,verx_3 };
			triangle_node_ids_.push_back(temp_tr);
			OGRLinearRing *ring = new OGRLinearRing();
			OGRPolygon *poly = new OGRPolygon();
			ring->addPoint(con_delau_tris_.triangle(face)[0].x(), con_delau_tris_.triangle(face)[0].y());
			ring->addPoint(con_delau_tris_.triangle(face)[1].x(), con_delau_tris_.triangle(face)[1].y());
			ring->addPoint(con_delau_tris_.triangle(face)[2].x(), con_delau_tris_.triangle(face)[2].y());
			ring->addPoint(con_delau_tris_.triangle(face)[0].x(), con_delau_tris_.triangle(face)[0].y());
			//ring->closeRings();
			poly->addRing(ring);

			tri_polygons_.push_back(*poly);

		}

	}
}


std::vector<double> BuildTriangulation::get_x_coords() {
	return x_coords_;
}

std::vector<double> BuildTriangulation::get_y_coords() {
	return y_coords_;
}

std::vector<double> BuildTriangulation::get_z_values() {
	return z_values_;
}

std::vector<std::array<int, 3>> BuildTriangulation::get_triangle_node_ids() {
	return triangle_node_ids_;
}

int BuildTriangulation::get_node_num() {
	return y_coords_.size();
}

int BuildTriangulation::get_triangle_num() {
	return triangle_node_ids_.size();
}

std::vector<bool>  BuildTriangulation::get_masks() {
	std::vector<bool> masks(triangle_node_ids_.size());
	for (int i = 0; i < get_triangle_num(); ++i) {
		masks[i] = false;
	}
	return masks;
}


void BuildTriangulation::write_triangles(const std::string filename)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("Delaunary_triangles", NULL, wkbPolygon, NULL);

	// define and create filed
	OGRFieldDefn oField_tri_id("id", OFTInteger);
	oField_tri_id.SetWidth(10);
	poLayer->CreateField(&oField_tri_id);

	OGRFieldDefn oField_node_id_1("node1", OFTInteger);
	oField_node_id_1.SetWidth(10);
	poLayer->CreateField(&oField_node_id_1);

	OGRFieldDefn oField_node_id_2("node2", OFTInteger);
	oField_node_id_2.SetWidth(10);
	poLayer->CreateField(&oField_node_id_2);

	OGRFieldDefn oField_node_id_3("node3", OFTInteger);
	oField_node_id_3.SetWidth(10);
	poLayer->CreateField(&oField_node_id_3);


	OGRFieldDefn oField_mask("mask", OFTInteger);
	oField_mask.SetWidth(10);
	poLayer->CreateField(&oField_mask);

	int tri_number = triangle_node_ids_.size();
	int t_mask = 0;

	for (int i = 0; i < tri_number; i++)
	{
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());

		poFeature->SetField("id", i);
		poFeature->SetField("node1", triangle_node_ids_[i][0]);
		poFeature->SetField("node2", triangle_node_ids_[i][1]);
		poFeature->SetField("node3", triangle_node_ids_[i][2]);
		poFeature->SetField("mask", t_mask);

		poFeature->SetGeometry(&(tri_polygons_[i]));

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(shpDataSet);
} // write_triangles



void BuildTriangulation::write_triangles_nodes(const std::string filename)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("Delaunary_triangles_nodes", NULL, wkbPoint, NULL);

	// define and create filed
	OGRFieldDefn oField_node_id("inner_id", OFTInteger);
	oField_node_id.SetWidth(10);
	poLayer->CreateField(&oField_node_id);

	OGRFieldDefn oField_node_distance("distance", OFTReal);
	oField_node_distance.SetWidth(10);
	poLayer->CreateField(&oField_node_distance);

	int node_number = node_points_inner_ids_pair_.size();

	for (int i = 0; i < node_number; i++)
	{
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());

		int inner_id = node_points_inner_ids_pair_[i].second;

		poFeature->SetField("inner_id", inner_id);
		poFeature->SetField("distance", z_values_[inner_id]);

		OGRPoint * point = new OGRPoint();
		point->setX(node_points_inner_ids_pair_[i].first.x());
		point->setY(node_points_inner_ids_pair_[i].first.y());
		poFeature->SetGeometry(point);

		if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
		{
			printf("Failed to create feature in shapefile.\n");
			exit(1);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(shpDataSet);

}//write_triangles_nodes



void BuildTriangulation::write_con_triangles_nodes(const std::string filename)
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("Delaunary_triangles_nodes", NULL, wkbPoint, NULL);

	// define and create filed
	OGRFieldDefn oField_node_id("id", OFTInteger);
	oField_node_id.SetWidth(10);
	poLayer->CreateField(&oField_node_id);

	OGRFieldDefn oField_node_distance("distance", OFTReal);
	oField_node_distance.SetWidth(10);
	poLayer->CreateField(&oField_node_distance);

	CDT::Finite_vertices_iterator vit;
	for (vit = con_delau_tris_.finite_vertices_begin(); vit != con_delau_tris_.finite_vertices_end(); ++vit)
	{
		Vertex_handle handle = vit;
		int id = handle->info(); 
		if (id >= 0 && id < z_values_.size()) 
		{
			OGRFeature *poFeature;
			poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
			poFeature->SetField("id", id);
			poFeature->SetField("distance", z_values_[id]);
			OGRPoint * point = new OGRPoint();
			point->setX(vit->point().x());
			point->setY(vit->point().y());
			poFeature->SetGeometry(point);
			if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
			{
				printf("Failed to create feature in shapefile.\n");
				exit(1);
			}
			OGRFeature::DestroyFeature(poFeature);
		}
	}
	GDALClose(shpDataSet);
}//write_triangles_nodes


void BuildTriangulation::write_con_edges(const std::string filename) 
{
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("constrined_segs", NULL, wkbLineString, NULL);

	// define and create filed
	OGRFieldDefn oField_seg_id("id", OFTInteger);
	oField_seg_id.SetWidth(10);
	poLayer->CreateField(&oField_seg_id);
	
	OGRFieldDefn oField_edge_id("edgeid", OFTInteger);
	oField_edge_id.SetWidth(10);
	poLayer->CreateField(&oField_edge_id);

	OGRFieldDefn oField_snode_id("sn_id", OFTInteger);
	oField_snode_id.SetWidth(10);
	poLayer->CreateField(&oField_snode_id);

	OGRFieldDefn oField_enode_id("en_id", OFTInteger);
	oField_enode_id.SetWidth(10);
	poLayer->CreateField(&oField_enode_id);


	int seg_number = constrained_segs.size();
	for (int i = 0; i < seg_number; i++)
	{
		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());

		poFeature->SetField("id", i);
		poFeature->SetField("edgeid", constrained_seg_egde_ids[i]);
		poFeature->SetField("sn_id", constrained_seg_snode_ids[i]);
		poFeature->SetField("en_id", constrained_seg_enode_ids[i]);
		poFeature->SetGeometry(constrained_segs[i]);

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

