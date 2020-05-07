#include "tca_tri_contour.h"
#include <algorithm>
#include <set>

namespace tca 
{

TriEdge::TriEdge()
	: tri(-1), edge(-1)
{}

TriEdge::TriEdge(int tri_, int edge_)
	: tri(tri_), edge(edge_)
{}

bool TriEdge::operator<(const TriEdge& other) const
{
	if (tri != other.tri)
		return tri < other.tri;
	else
		return edge < other.edge;
}

bool TriEdge::operator==(const TriEdge& other) const
{
	return tri == other.tri && edge == other.edge;
}

bool TriEdge::operator!=(const TriEdge& other) const
{
	return !operator==(other);
}

std::ostream& operator<<(std::ostream& os, const TriEdge& tri_edge)
{
	return os << tri_edge.tri << ' ' << tri_edge.edge;
}


XY::XY()
{}

XY::XY(const double& x_, const double& y_)
	: x(x_), y(y_)
{}

double XY::angle() const
{
	return atan2(y, x);
}

double XY::cross_z(const XY& other) const
{
	return x * other.y - y * other.x;
}

bool XY::is_right_of(const XY& other) const
{
	if (x == other.x)
		return y > other.y;
	else
		return x > other.x;
}

bool XY::operator==(const XY& other) const
{
	return x == other.x && y == other.y;
}

bool XY::operator!=(const XY& other) const
{
	return x != other.x || y != other.y;
}

XY XY::operator*(const double& multiplier) const
{
	return XY(x*multiplier, y*multiplier);
}

const XY& XY::operator+=(const XY& other)
{
	x += other.x;
	y += other.y;
	return *this;
}

const XY& XY::operator-=(const XY& other)
{
	x -= other.x;
	y -= other.y;
	return *this;
}

XY XY::operator+(const XY& other) const
{
	return XY(x + other.x, y + other.y);
}

XY XY::operator-(const XY& other) const
{
	return XY(x - other.x, y - other.y);
}

std::ostream& operator<<(std::ostream& os, const XY& xy)
{
	return os << '(' << xy.x << ' ' << xy.y << ')';
}



XYZ::XYZ(const double& x_, const double& y_, const double& z_)
	: x(x_), y(y_), z(z_)
{}

XYZ XYZ::cross(const XYZ& other) const
{
	return XYZ(y*other.z - z * other.y,
		z*other.x - x * other.z,
		x*other.y - y * other.x);
}

double XYZ::dot(const XYZ& other) const
{
	return x * other.x + y * other.y + z * other.z;
}

double XYZ::length_squared() const
{
	return x * x + y * y + z * z;
}

XYZ XYZ::operator-(const XYZ& other) const
{
	return XYZ(x - other.x, y - other.y, z - other.z);
}

std::ostream& operator<<(std::ostream& os, const XYZ& xyz)
{
	return os << '(' << xyz.x << ' ' << xyz.y << ' ' << xyz.z << ')';
}



BoundingBox::BoundingBox()
	: empty(true)
{}

void BoundingBox::add(const XY& point)
{
	if (empty) {
		empty = false;
		lower = upper = point;
	}
	else {
		if (point.x < lower.x) lower.x = point.x;
		else if (point.x > upper.x) upper.x = point.x;

		if (point.y < lower.y) lower.y = point.y;
		else if (point.y > upper.y) upper.y = point.y;
	}
}

void BoundingBox::expand(const XY& delta)
{
	if (!empty) {
		lower -= delta;
		upper += delta;
	}
}


ContourLine::ContourLine()
	: std::vector<XY>()
{}

void ContourLine::insert_unique(iterator pos, const XY& point)
{
	if (empty() || pos == end() || point != *pos)
		std::vector<XY>::insert(pos, point);
}

void ContourLine::push_back(const XY& point)
{
	if (empty() || point != back())
		std::vector<XY>::push_back(point);
}

void ContourLine::write() const
{
	std::cout << "ContourLine of " << size() << " points:";
	for (const_iterator it = begin(); it != end(); ++it)
		std::cout << ' ' << *it;
	std::cout << std::endl;
}


void write_contour(const Contour& contour)
{
	std::cout << "Contour of " << contour.size() << " lines." << std::endl;
	for (Contour::const_iterator it = contour.begin(); it != contour.end(); ++it)
		it->write();
}

//
Triangulation::Triangulation(const CoordinateArray& x,
	const CoordinateArray& y,
	const TriangleArray& triangles,
	const MaskArray& mask,
	const EdgeArray& edges,
	const NeighborArray& neighbors,
	int correct_triangle_orientations)
	: _x(x),
	_y(y),
	_triangles(triangles),
	_mask(mask),
	_edges(edges),
	_neighbors(neighbors)
{
	if (correct_triangle_orientations)
		correct_triangles();

	//Added by Diao
	calculate_neighbors_at_first();
}


Triangulation::~Triangulation()
{


}



void Triangulation::calculate_boundaries()
{
	get_neighbors();  // Ensure _neighbors has been created.

	// Create set of all boundary TriEdges, which are those which do not
	// have a neighbor triangle.
	typedef std::set<TriEdge> BoundaryEdges;
	BoundaryEdges boundary_edges;
	for (int tri = 0; tri < get_ntri(); ++tri) 
	{
		if (!is_masked(tri)) 
		{
			for (int edge = 0; edge < 3; ++edge) 
			{
				if (get_neighbor(tri, edge) == -1) 
				{
					boundary_edges.insert(TriEdge(tri, edge));
				}
			}
		}
	}

	// Take any boundary edge and follow the boundary until return to start
	// point, removing edges from boundary_edges as they are used.  At the same
	// time, initialise the _tri_edge_to_boundary_map.
	while (!boundary_edges.empty()) 
	{
		// Start of new boundary.
		BoundaryEdges::iterator it = boundary_edges.begin();
		int tri = it->tri;
		int edge = it->edge;
		_boundaries.push_back(Boundary());
		Boundary& boundary = _boundaries.back();

		while (true) 
		{
			boundary.push_back(TriEdge(tri, edge));
			boundary_edges.erase(it);
			_tri_edge_to_boundary_map[TriEdge(tri, edge)] = BoundaryEdge(_boundaries.size() - 1, boundary.size() - 1);

			// Move to next edge of current triangle.
			edge = (edge + 1) % 3;

			// Find start point index of boundary edge.
			int point = get_triangle_point(tri, edge);

			// Find next TriEdge by traversing neighbors until find one
			// without a neighbor.
			while (get_neighbor(tri, edge) != -1) 
			{
				tri = get_neighbor(tri, edge);
				edge = get_edge_in_triangle(tri, point);
			}

			if (TriEdge(tri, edge) == boundary.front())
				break;  // Reached beginning of this boundary, so finished it.
			else
				it = boundary_edges.find(TriEdge(tri, edge));
		}
	}
}

void Triangulation::calculate_edges()
{
	assert(_edges.empty() && "Expected empty edges array");

	// Create set of all edges, storing them with start point index less than
	// end point index.
	typedef std::set<Edge> EdgeSet;
	EdgeSet edge_set;
	for (int tri = 0; tri < get_ntri(); ++tri) {
		if (!is_masked(tri)) {
			for (int edge = 0; edge < 3; edge++) {
				int start = get_triangle_point(tri, edge);
				int end = get_triangle_point(tri, (edge + 1) % 3);
				edge_set.insert(start > end ? Edge(start, end) : Edge(end, start));
			}
		}
	}

	// initialization for _edge
	int i = 0;
	for (EdgeSet::const_iterator it = edge_set.begin(); it != edge_set.end(); ++it) {
		_edges[i][0] = it->start;
		_edges[i++][1] = it->end;
	}
}

// Added by Diao 
void Triangulation::calculate_neighbors_at_first()
{
	//assert(_neighbors.empty() && "Expected empty neighbors array");

	int tri, edge;
	for (tri = 0; tri < get_ntri(); ++tri) {
		for (edge = 0; edge < 3; ++edge)
			_neighbors[tri][edge] = -1;
	}

	// For each triangle edge (start to end point), find corresponding neighbor
	// edge from end to start point.  Do this by traversing all edges and
	// storing them in a map from edge to TriEdge.  If corresponding neighbor
	// edge is already in the map, don't need to store new edge as neighbor
	// already found.
	typedef std::map<Edge, TriEdge> EdgeToTriEdgeMap;
	EdgeToTriEdgeMap edge_to_tri_edge_map;
	for (tri = 0; tri < get_ntri(); ++tri) 
	{
		if (!is_masked(tri)) 
		{
			for (edge = 0; edge < 3; ++edge) 
			{
				int start = get_triangle_point(tri, edge);
				int end = get_triangle_point(tri, (edge + 1) % 3);

				EdgeToTriEdgeMap::iterator it = edge_to_tri_edge_map.find(Edge(end, start));
				
				if (it == edge_to_tri_edge_map.end()) 
				{
					// No neighbor edge exists in the edge_to_tri_edge_map, so
					// add this edge to it.
					edge_to_tri_edge_map[Edge(start, end)] = TriEdge(tri, edge);
				}
				else 
				{
					// Neighbor edge found, set the two elements of _neighbors
					// and remove edge from edge_to_tri_edge_map.
					_neighbors[tri][edge] = it->second.tri;
					_neighbors[it->second.tri][it->second.edge] = tri;
					edge_to_tri_edge_map.erase(it);
				}
			}
		}
	}

	// Note that remaining edges in the edge_to_tri_edge_map correspond to
	// boundary edges, but the boundaries are calculated separately elsewhere.
}




void Triangulation::calculate_neighbors()
{
	assert(_neighbors.empty() && "Expected empty neighbors array");

	int tri, edge;
	for (tri = 0; tri < get_ntri(); ++tri) {
		for (edge = 0; edge < 3; ++edge)
			_neighbors[tri][edge] = -1;
	}

	// For each triangle edge (start to end point), find corresponding neighbor
	// edge from end to start point.  Do this by traversing all edges and
	// storing them in a map from edge to TriEdge.  If corresponding neighbor
	// edge is already in the map, don't need to store new edge as neighbor
	// already found.
	typedef std::map<Edge, TriEdge> EdgeToTriEdgeMap;
	EdgeToTriEdgeMap edge_to_tri_edge_map;
	for (tri = 0; tri < get_ntri(); ++tri) {
		if (!is_masked(tri)) {
			for (edge = 0; edge < 3; ++edge) {
				int start = get_triangle_point(tri, edge);
				int end = get_triangle_point(tri, (edge + 1) % 3);
				EdgeToTriEdgeMap::iterator it =
					edge_to_tri_edge_map.find(Edge(end, start));
				if (it == edge_to_tri_edge_map.end()) {
					// No neighbor edge exists in the edge_to_tri_edge_map, so
					// add this edge to it.
					edge_to_tri_edge_map[Edge(start, end)] = TriEdge(tri, edge);
				}
				else {
					// Neighbor edge found, set the two elements of _neighbors
					// and remove edge from edge_to_tri_edge_map.
					_neighbors[tri][edge] = it->second.tri;
					_neighbors[it->second.tri][it->second.edge] = tri;
					edge_to_tri_edge_map.erase(it);
				}
			}
		}
	}

	// Note that remaining edges in the edge_to_tri_edge_map correspond to
	// boundary edges, but the boundaries are calculated separately elsewhere.
}



Triangulation::TwoCoordinateArray Triangulation::calculate_plane_coefficients(
	const CoordinateArray& z)
{
	Triangulation::TwoCoordinateArray planes(get_ntri());
	
	int point;
	for (int tri = 0; tri < get_ntri(); ++tri) {
		if (is_masked(tri)) {
			planes[tri][0] = 0.0;
			planes[tri][1] = 0.0;
			planes[tri][2] = 0.0;
		}
		else {
			// Equation of plane for all points r on plane is r.normal = p
			// where normal is vector normal to the plane, and p is a
			// constant.  Rewrite as
			// r_x*normal_x + r_y*normal_y + r_z*normal_z = p
			// and rearrange to give
			// r_z = (-normal_x/normal_z)*r_x + (-normal_y/normal_z)*r_y +
			//       p/normal_z
			point = _triangles[tri][0];
			
			XYZ point0(_x[point], _y[point], z[point]);
						
			point = _triangles[tri][1];
			XYZ side01 = XYZ(_x[point], _y[point], z[point]) - point0;

			point = _triangles[tri][2];
			XYZ side02 = XYZ(_x[point], _y[point], z[point]) - point0;

			XYZ normal = side01.cross(side02);

			if (normal.z == 0.0) {
				// Normal is in x-y plane which means triangle consists of
				// colinear points. To avoid dividing by zero, we use the
				// Moore-Penrose pseudo-inverse.
				double sum2 = (side01.x*side01.x + side01.y*side01.y +
					side02.x*side02.x + side02.y*side02.y);
				double a = (side01.x*side01.z + side02.x*side02.z) / sum2;
				double b = (side01.y*side01.z + side02.y*side02.z) / sum2;
				planes[tri][0] = a;
				planes[tri][1] = b;
				planes[tri][2] = point0.z - a * point0.x - b * point0.y;
			}
			else {
				planes[tri][0] = -normal.x / normal.z;           // x
				planes[tri][1] = -normal.y / normal.z;           // y
				planes[tri][2] = normal.dot(point0) / normal.z;  // constant
			}
		}
	}

	return planes;
}



void Triangulation::correct_triangles()
{
	for (int tri = 0; tri < get_ntri(); ++tri) {
		XY point0 = get_point_coords(_triangles[tri][0]);
		XY point1 = get_point_coords(_triangles[tri][1]);
		XY point2 = get_point_coords(_triangles[tri][2]);
		if ((point1 - point0).cross_z(point2 - point0) < 0.0) {
			// Triangle points are clockwise, so change them to anticlockwise.
			std::swap(_triangles[tri][1], _triangles[tri][2]);
			if (!_neighbors.empty())
				std::swap(_neighbors[tri][1], _neighbors[tri][2]);
		}
	}
}


const Triangulation::Boundaries& Triangulation::get_boundaries() const
{
	if (_boundaries.empty())
		const_cast<Triangulation*>(this)->calculate_boundaries();
	return _boundaries;
}


void Triangulation::get_boundary_edge(const TriEdge& triEdge,
	int& boundary,
	int& edge) const
{
	get_boundaries();  // Ensure _tri_edge_to_boundary_map has been created.
	TriEdgeToBoundaryMap::const_iterator it =
		_tri_edge_to_boundary_map.find(triEdge);
	assert(it != _tri_edge_to_boundary_map.end() &&
		"TriEdge is not on a boundary");
	boundary = it->second.boundary;
	edge = it->second.edge;
}


int Triangulation::get_edge_in_triangle(int tri, int point) const
{
	assert(tri >= 0 && tri < get_ntri() && "Triangle index out of bounds");
	assert(point >= 0 && point < get_npoints() && "Point index out of bounds.");
	for (int edge = 0; edge < 3; ++edge) {
		if (_triangles[tri][edge] == point)
			return edge;
	}
	return -1;  // point is not in triangle.
}


Triangulation::EdgeArray& Triangulation::get_edges()
{
	if (_edges.empty())
		calculate_edges();
	return _edges;
}


int Triangulation::get_neighbor(int tri, int edge) const
{
	assert(tri >= 0 && tri < get_ntri() && "Triangle index out of bounds");
	assert(edge >= 0 && edge < 3 && "Edge index out of bounds");
	if (_neighbors.empty())
		const_cast<Triangulation&>(*this).calculate_neighbors();
	return _neighbors[tri][edge];
}

TriEdge Triangulation::get_neighbor_edge(int tri, int edge) const
{
	int neighbor_tri = get_neighbor(tri, edge);
	if (neighbor_tri == -1)
		return TriEdge(-1, -1);
	else
		return TriEdge(neighbor_tri,
			get_edge_in_triangle(neighbor_tri,
				get_triangle_point(tri,
				(edge + 1) % 3)));
}



Triangulation::NeighborArray& Triangulation::get_neighbors()
{
	if (_neighbors.empty())
		calculate_neighbors();
	return _neighbors;
}

int Triangulation::get_npoints() const
{
	return _x.size();
}



int Triangulation::get_ntri() const
{
	return _triangles.size();
}

XY Triangulation::get_point_coords(int point) const
{
	assert(point >= 0 && point < get_npoints() && "Point index out of bounds.");
	return XY(_x[point], _y[point]);
}

int Triangulation::get_triangle_point(int tri, int edge) const
{
	assert(tri >= 0 && tri < get_ntri() && "Triangle index out of bounds");
	assert(edge >= 0 && edge < 3 && "Edge index out of bounds");
	return _triangles[tri][edge];
}


int Triangulation::get_triangle_point(const TriEdge& tri_edge) const
{
	return get_triangle_point(tri_edge.tri, tri_edge.edge);
}

bool Triangulation::is_masked(int tri) const
{
	assert(tri >= 0 && tri < get_ntri() && "Triangle index out of bounds.");
	return !_mask.empty() && _mask[tri];
}

void Triangulation::write_boundaries() const
{
    const Boundaries& bs = get_boundaries();
    std::cout << "Number of boundaries: " << bs.size() << std::endl;
    for (Boundaries::const_iterator it = bs.begin(); it != bs.end(); ++it) {
        const Boundary& b = *it;
        std::cout << "  Boundary of " << b.size() << " points: ";
        for (Boundary::const_iterator itb = b.begin(); itb != b.end(); ++itb) {
            std::cout << *itb << ", ";
        }
        std::cout << std::endl;
    }
}


TriContourGenerator::TriContourGenerator(Triangulation& triangulation,
	const CoordinateArray& z)
	: _triangulation(triangulation),
	_z(z),
	_interior_visited(2 * _triangulation.get_ntri()),
	_boundaries_visited(0),
	_boundaries_used(0)
{}

TriContourGenerator::~TriContourGenerator() 
{
	
}


void TriContourGenerator::clear_visited_flags(bool include_boundaries)
{
	// Clear _interiorVisited.
	std::fill(_interior_visited.begin(), _interior_visited.end(), false);

	if (include_boundaries) {
		if (_boundaries_visited.empty()) {
			const Boundaries& boundaries = get_boundaries();

			// Initialise _boundaries_visited.
			_boundaries_visited.reserve(boundaries.size());
			for (Boundaries::const_iterator it = boundaries.begin();
				it != boundaries.end(); ++it)
				_boundaries_visited.push_back(BoundaryVisited(it->size()));

			// Initialise _boundaries_used.
			_boundaries_used = BoundariesUsed(boundaries.size());
		}

		// Clear _boundaries_visited.
		for (BoundariesVisited::iterator it = _boundaries_visited.begin();
			it != _boundaries_visited.end(); ++it)
			std::fill(it->begin(), it->end(), false);

		// Clear _boundaries_used.
		std::fill(_boundaries_used.begin(), _boundaries_used.end(), false);
	}
}



// Changed by Diao 
Contour TriContourGenerator::create_contour(const double& level, std::vector<int> & contour_labels)
{
	clear_visited_flags(false);
	Contour contour;

	find_boundary_lines(contour, level, contour_labels);
	find_interior_lines(contour, level, false, false, contour_labels);

	return contour;
}



Contour TriContourGenerator::create_filled_contour(const double& lower_level, const double& upper_level, std::vector<int> & contour_labels)
{
	clear_visited_flags(true);
	Contour contour;
	
	find_boundary_lines_filled(contour, lower_level, upper_level, contour_labels);
	find_interior_lines(contour, lower_level, false, true, contour_labels);
	find_interior_lines(contour, upper_level, true, true, contour_labels);
	return contour;
}

XY TriContourGenerator::edge_interp(int tri, int edge, const double& level)
{
	return interp(_triangulation.get_triangle_point(tri, edge),
		_triangulation.get_triangle_point(tri, (edge + 1) % 3),
		level);
}

void TriContourGenerator::find_boundary_lines(Contour& contour,
	const double& level, std::vector<int> & contour_labels)
{
	// Traverse boundaries to find starting points for all contour lines that
	// intersect the boundaries.  For each starting point found, follow the
	// line to its end before continuing.
	const Triangulation& triang = _triangulation;
	const Boundaries& boundaries = get_boundaries();
	for (Boundaries::const_iterator it = boundaries.begin(); it != boundaries.end(); ++it) 
	{
		const Boundary& boundary = *it;
		bool startAbove, endAbove = false;
		for (Boundary::const_iterator itb = boundary.begin(); itb != boundary.end(); ++itb) 
		{
			if (itb == boundary.begin())
				startAbove = get_z(triang.get_triangle_point(*itb)) >= level;
			else
				startAbove = endAbove;
			endAbove = get_z(triang.get_triangle_point(itb->tri,
				(itb->edge + 1) % 3)) >= level;
			if (startAbove && !endAbove) 
			{
				// This boundary edge is the start point for a contour line,
				// so follow the line.
				contour.push_back(ContourLine());
				ContourLine& contour_line = contour.back();
				contour_labels.push_back(1); // added by Diao, contour label
				TriEdge tri_edge = *itb;
				follow_interior(contour_line, tri_edge, true, level, false);
			}
		}
	}
}

// 
void TriContourGenerator::find_boundary_lines_filled(Contour& contour,
	const double& lower_level,
	const double& upper_level,
	std::vector<int> & contour_labels)
{
	// Traverse boundaries to find starting points for all contour lines that
	// intersect the boundaries.  For each starting point found, follow the
	// line to its end before continuing.
	const Triangulation& triang = _triangulation;
	const Boundaries& boundaries = get_boundaries();
	for (Boundaries::size_type i = 0; i < boundaries.size(); ++i) 
	{
		const Boundary& boundary = boundaries[i];
		for (Boundary::size_type j = 0; j < boundary.size(); ++j) 
		{
			if (!_boundaries_visited[i][j]) 
			{
				// z values of start and end of this boundary edge.
				double z_start = get_z(triang.get_triangle_point(boundary[j]));
				double z_end = get_z(triang.get_triangle_point(
					boundary[j].tri, (boundary[j].edge + 1) % 3));

				// Does this boundary edge's z increase through upper level
				// and/or decrease through lower level?
				bool incr_upper = (z_start < upper_level && z_end >= upper_level);
				bool decr_lower = (z_start >= lower_level && z_end < lower_level);

				if (decr_lower || incr_upper) 
				{
					// Start point for contour line, so follow it.
					contour.push_back(ContourLine());
					
					contour_labels.push_back(1); // added by Diao; contour label 

					ContourLine& contour_line = contour.back();
					TriEdge start_tri_edge = boundary[j];
					TriEdge tri_edge = start_tri_edge;

					// Traverse interior and boundaries until return to start.
					bool on_upper = incr_upper;
					do {
						follow_interior(contour_line, tri_edge, true,
							on_upper ? upper_level : lower_level, on_upper);
						on_upper = follow_boundary(contour_line, tri_edge,
							lower_level, upper_level, on_upper);
					} while (tri_edge != start_tri_edge);

					// Filled contour lines must not have same first and last
					// points.
					if (contour_line.size() > 1 &&
						contour_line.front() == contour_line.back())
						contour_line.pop_back();
				}
			}
		}
	}

	// Add full boundaries that lie between the lower and upper levels.  These
	// are boundaries that have not been touched by an internal contour line
	// which are stored in _boundaries_used.
	for (Boundaries::size_type i = 0; i < boundaries.size(); ++i) {
		if (!_boundaries_used[i]) {
			const Boundary& boundary = boundaries[i];
			double z = get_z(triang.get_triangle_point(boundary[0]));
			if (z >= lower_level && z < upper_level) {
				contour.push_back(ContourLine());
				ContourLine& contour_line = contour.back();
				for (Boundary::size_type j = 0; j < boundary.size(); ++j)
					contour_line.push_back(triang.get_point_coords(
						triang.get_triangle_point(boundary[j])));
			}
		}
	}
}

void TriContourGenerator::find_interior_lines(Contour& contour,
	const double& level,
	bool on_upper,
	bool filled, std::vector<int> & contour_labels)
{
	const Triangulation& triang = _triangulation;
	int ntri = triang.get_ntri();
	for (int tri = 0; tri < ntri; ++tri) 
	{
		int visited_index = (on_upper ? tri + ntri : tri);

		if (_interior_visited[visited_index] || triang.is_masked(tri))
			continue;  // Triangle has already been visited or is masked.

		_interior_visited[visited_index] = true;

		// Determine edge via which to leave this triangle.
		int edge = get_exit_edge(tri, level, on_upper);
		assert(edge >= -1 && edge < 3 && "Invalid exit edge");
		if (edge == -1)
			continue;  // Contour does not pass through this triangle.

		// Found start of new contour line loop.
		contour.push_back(ContourLine());
		ContourLine& contour_line = contour.back();
		
		contour_labels.push_back(0);// added by Diao, contour label; 1 means boudary followers, 0 means interior
		
		TriEdge tri_edge = triang.get_neighbor_edge(tri, edge);
		follow_interior(contour_line, tri_edge, false, level, on_upper);

		if (!filled)
			// Non-filled contour lines must be closed.
			contour_line.push_back(contour_line.front());
		else if (contour_line.size() > 1 &&
			contour_line.front() == contour_line.back())
			// Filled contour lines must not have same first and last points.
			contour_line.pop_back();
	}
}

bool TriContourGenerator::follow_boundary(ContourLine& contour_line,
	TriEdge& tri_edge,
	const double& lower_level,
	const double& upper_level,
	bool on_upper)
{
	const Triangulation& triang = _triangulation;
	const Boundaries& boundaries = get_boundaries();

	// Have TriEdge to start at, need equivalent boundary edge.
	int boundary, edge;
	triang.get_boundary_edge(tri_edge, boundary, edge);
	_boundaries_used[boundary] = true;

	bool stop = false;
	bool first_edge = true;
	double z_start, z_end = 0;
	while (!stop)
	{
		assert(!_boundaries_visited[boundary][edge] && "Boundary already visited");
		_boundaries_visited[boundary][edge] = true;

		// z values of start and end points of boundary edge.
		if (first_edge)
			z_start = get_z(triang.get_triangle_point(tri_edge));
		else
			z_start = z_end;
		z_end = get_z(triang.get_triangle_point(tri_edge.tri,
			(tri_edge.edge + 1) % 3));

		if (z_end > z_start) {  // z increasing.
			if (!(!on_upper && first_edge) &&
				z_end >= lower_level && z_start < lower_level) {
				stop = true;
				on_upper = false;
			}
			else if (z_end >= upper_level && z_start < upper_level) {
				stop = true;
				on_upper = true;
			}
		}
		else {  // z decreasing.
			if (!(on_upper && first_edge) &&
				z_start >= upper_level && z_end < upper_level) {
				stop = true;
				on_upper = true;
			}
			else if (z_start >= lower_level && z_end < lower_level) {
				stop = true;
				on_upper = false;
			}
		}

		first_edge = false;

		if (!stop) {
			// Move to next boundary edge, adding point to contour line.
			edge = (edge + 1) % (int)boundaries[boundary].size();
			tri_edge = boundaries[boundary][edge];
			contour_line.push_back(triang.get_point_coords(
				triang.get_triangle_point(tri_edge)));
		}
	}

	return on_upper;
}

void TriContourGenerator::follow_interior(ContourLine& contour_line,
	TriEdge& tri_edge,
	bool end_on_boundary,
	const double& level,
	bool on_upper)
{
	int& tri = tri_edge.tri;
	int& edge = tri_edge.edge;

	// Initial point.
	contour_line.push_back(edge_interp(tri, edge, level));

	while (true) 
	{
		int visited_index = tri;
		if (on_upper)
			visited_index += _triangulation.get_ntri();

		// Check for end not on boundary.
		if (!end_on_boundary && _interior_visited[visited_index])
			break;  // Reached start point, so return.

		// Determine edge by which to leave this triangle.
		edge = get_exit_edge(tri, level, on_upper);
		assert(edge >= 0 && edge < 3 && "Invalid exit edge");

		_interior_visited[visited_index] = true;

		// Append new point to point set.
		assert(edge >= 0 && edge < 3 && "Invalid triangle edge");
		contour_line.push_back(edge_interp(tri, edge, level));

		// Move to next triangle.
		TriEdge next_tri_edge = _triangulation.get_neighbor_edge(tri, edge);

		// Check if ending on a boundary.
		if (end_on_boundary && next_tri_edge.tri == -1)
			break;

		tri_edge = next_tri_edge;
		assert(tri_edge.tri != -1 && "Invalid triangle for internal loop");
	}
}

const TriContourGenerator::Boundaries& TriContourGenerator::get_boundaries() const
{
	return _triangulation.get_boundaries();
}

int TriContourGenerator::get_exit_edge(int tri,
	const double& level,
	bool on_upper) const
{
	assert(tri >= 0 && tri < _triangulation.get_ntri() &&
		"Triangle index out of bounds.");

	unsigned int config =
		(get_z(_triangulation.get_triangle_point(tri, 0)) >= level) |
		(get_z(_triangulation.get_triangle_point(tri, 1)) >= level) << 1 |
		(get_z(_triangulation.get_triangle_point(tri, 2)) >= level) << 2;

	if (on_upper) config = 7 - config;

	switch (config) {
	case 0: return -1;
	case 1: return  2;
	case 2: return  0;
	case 3: return  2;
	case 4: return  1;
	case 5: return  1;
	case 6: return  0;
	case 7: return -1;
	default: assert(0 && "Invalid config value"); return -1;
	}
}

const double& TriContourGenerator::get_z(int point) const
{
	assert(point >= 0 && point < _triangulation.get_npoints() &&
		"Point index out of bounds.");
	return _z[point];
}

XY TriContourGenerator::interp(int point1,
	int point2,
	const double& level) const
{
	assert(point1 >= 0 && point1 < _triangulation.get_npoints() &&
		"Point index 1 out of bounds.");
	assert(point2 >= 0 && point2 < _triangulation.get_npoints() &&
		"Point index 2 out of bounds.");
	assert(point1 != point2 && "Identical points");
	double fraction = (get_z(point2) - level) / (get_z(point2) - get_z(point1));
	return _triangulation.get_point_coords(point1)*fraction +
		_triangulation.get_point_coords(point2)*(1.0 - fraction);
}


// Added by Diao
void TriContourGenerator::write_contour_lines(const std::string &filename, const Contour contour_lines) {
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	OGRRegisterAll();
	GDALDriver *shpDriver;
	GDALDataset *shpDataSet;

	shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	shpDataSet = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
	OGRLayer *poLayer = shpDataSet->CreateLayer("contourlines", NULL, wkbLineString, NULL);
	
	// define and create filed
	OGRFieldDefn oField_tri_id("id", OFTInteger);
	oField_tri_id.SetWidth(10);
	poLayer->CreateField(&oField_tri_id);
	
	int line_num = contour_lines.size();
		
	for (int i = 0; i < line_num; i++)
	{
		if (contour_lines[i].size() != 1)
		{
			OGRFeature *poFeature;
			poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
			OGRLineString *linestring = new OGRLineString();
			
			for (int j = 0; j < contour_lines[i].size(); j++) 
			{
				linestring->addPoint(contour_lines[i][j].x, contour_lines[i][j].y);
			}			
			poFeature->SetField("id", i);
			poFeature->SetGeometry(linestring);

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

// write countour polygons
void TriContourGenerator::write_contour_polygons(const std::string &filename, const Contour contour_lines, int facility_id) 
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

	int line_num = contour_lines.size();

	OGRFeature *poFeature;
	poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
	OGRPolygon poly;
	for (int i = 0; i < line_num; i++)
	{
		OGRLinearRing  linering;
		if (contour_lines[i].size() != 1)
		{					
			for (int j = 0; j < contour_lines[i].size(); j++)
			{
				linering.addPoint(contour_lines[i][j].x, contour_lines[i][j].y);
			}
			linering.closeRings();						
		}
		poly.addRing(&linering);
	}

	poFeature->SetField("id", facility_id);
	poFeature->SetGeometry(&poly);

	if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
	{
		printf("Failed to create feature in shapefile.\n");
		exit(1);
	}
	OGRFeature::DestroyFeature(poFeature);
	GDALClose(shpDataSet);
}
}// namespace tca
