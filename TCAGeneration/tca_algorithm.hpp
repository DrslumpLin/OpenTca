/**
 * Definitions of some basic algorithms that frequetly used in TCA Project
 * @author: Diao Lin
 * @version: 2020.03
 */


#ifndef TCA_ALGORITHM_HPP
#define TCA_ALGORITHM_HPP
#include <cmath>
#include "ogrsf_frmts.h" // C++ API for GDAL

namespace tca 
{
namespace algorithm 
{
	
	// Compute the boundary of an OGRLineString
	// Input: an OGRLineString pointer, @linestring
	// output: @box_x_min, @box_y_min, @box_x_max, @box_y_max

	inline	void GetLinestringBoundingbox(OGRLineString *linestring, double * box_x_min, double * box_y_min,
		double *box_x_max, double *box_y_max)
	{
		int num_points = linestring->getNumPoints();
		double x, y;

		*box_x_min = DBL_MAX;
		*box_y_min = DBL_MAX;
		*box_x_max = DBL_MIN;
		*box_y_max = DBL_MIN;

		for (int i = 0; i < num_points; ++i)
		{
			x = linestring->getX(i);
			y = linestring->getY(i);
			if (x < *box_x_min) *box_x_min = x;
			if (y < *box_y_min) *box_y_min = y;
			if (x > *box_x_max) *box_x_max = x;
			if (y > *box_y_max) *box_y_max = y;
		};
	}; // GetLinestringBoundingbox


	// Project a point to a segement
	// input: a point (@input_point_x, @input_point_y) and a segment [(@seg_x1, @seg_y1), (@seg_x2, @seg_y2)]
	// output: 
	//		  @vertical_distance: the distance between the point and its projected point
	//		  @break_distance: the distance from the start node (of the segment) to the projected point
	inline	void ProjectPointOnSegment(double input_point_x, double input_point_y,
		double seg_x1, double seg_y1, double seg_x2, double seg_y2, double * vertical_distance, double * break_distance)
	{
		double seg_length = (seg_x2 - seg_x1)*(seg_x2 - seg_x1) + (seg_y2 - seg_y1)*(seg_y2 - seg_y1);

		if (seg_length == 0.0)
		{
			*vertical_distance = std::sqrt((input_point_x - seg_x1)*(input_point_x - seg_x1) + (input_point_y - seg_y1)*(input_point_y - seg_y1));
			*break_distance = 0.0;
			return;
		}

		double x1_x = input_point_x - seg_x1;
		double y1_y = input_point_y - seg_y1;

		double x1_x2 = seg_x2 - seg_x1;
		double y1_y2 = seg_y2 - seg_y1;

		double ratio = (x1_x*x1_x2 + y1_y * y1_y2) / seg_length;

		ratio = (ratio > 1) ? 1 : ratio;
		ratio = (ratio < 0) ? 0 : ratio;

		double prj_x = seg_x1 + ratio * (x1_x2);
		double prj_y = seg_y1 + ratio * (y1_y2);

		*vertical_distance = std::sqrt((prj_x - input_point_x)*(prj_x - input_point_x) + (prj_y - input_point_y)*(prj_y - input_point_y));
		*break_distance = std::sqrt((prj_x - seg_x1)*(prj_x - seg_x1) + (prj_y - seg_y1)*(prj_y - seg_y1));

	}; // ProjectPointOnSegment

	 // Projecting a point to a linestring
	 // Input: @point, @linestring
	 // output:
	 // 	@break_seg_idx: idx of seg with minimum vertical distance
	 // 	@vertical_distance: minimum distance to the linestring
	 // 	@ break_distance: break distance from projected point to the start point of the linestring
	inline	void ProjectPointOnLinestring(OGRPoint *point, OGRLineString *linestring, int * break_seg_idx, double * vertical_distance, float *break_distance)
	{
		int numpoints = linestring->getNumPoints();

		double min_vertical_distance = DBL_MAX;
		double final_break_distance = 0;

		double length_parsed = 0;
		int min_seg_index = 0;

		double x = point->getX();
		double y = point->getY();

		int i = 0;
		while (i < numpoints - 1)
		{
			double x1 = linestring->getX(i);
			double y1 = linestring->getY(i);
			double x2 = linestring->getX(i + 1);
			double y2 = linestring->getY(i + 1);

			double temp_min_vertical_distance;
			double temp_break_distance;

			ProjectPointOnSegment(x, y, x1, y1, x2, y2, &temp_min_vertical_distance, &temp_break_distance);

			if (temp_min_vertical_distance < min_vertical_distance)
			{
				min_seg_index = i;
				min_vertical_distance = temp_min_vertical_distance;
				final_break_distance = length_parsed + temp_break_distance;
			}

			length_parsed += std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));

			++i;
		};

		*break_seg_idx = min_seg_index;
		*vertical_distance = min_vertical_distance;
		*break_distance = final_break_distance;

	}; // linear_referencing


	//Given a break distance (to start node) and segment([x1, y1], [x2, y2]), return the break point on linestring
	inline	OGRPoint GetBreakPointOnSegment(double seg_x1, double seg_y1, double seg_x2, double seg_y2, double break_distance_to_snode)
	{
		OGRPoint temppoint;
	
		double seg_length = std::sqrt((seg_x2 - seg_x1)*(seg_x2 - seg_x1) + (seg_y2 - seg_y1)*(seg_y2 - seg_y1));
		double ratio = break_distance_to_snode / seg_length;

		double break_x = seg_x1 + ratio * (seg_x2 - seg_x1);
		double break_y = seg_y1 + ratio * (seg_y2 - seg_y1);
	
		temppoint.setX(break_x);
		temppoint.setY(break_y);
	
		return temppoint;
	};//GetBreakPointOnSegment


	//Given a break distance to the start of a (OGR)linestring, return the break point
	inline	OGRPoint GetBreakPointOnLinestring(double break_distance_to_snode, OGRLineString * linestring)
	{	
		int i = 0;
		double length_parsed = 0;
		int	   min_seg_index = 0; // the seg index equal to its start node idx		
		double seg_len;	
		double local_break_distance;
		double break_point_x, break_point_yd;
		OGRPoint tempbreakpoint;

		int numpoints = linestring->getNumPoints();

		while (i < numpoints - 1)
		{
			double x1 = linestring->getX(i);
			double y1 = linestring->getY(i);
			double x2 = linestring->getX(i + 1);
			double y2 = linestring->getY(i + 1);
			seg_len = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
			length_parsed += seg_len;
		
			// due to the precision issue, the length_parsed can not be accurate represented
			// we manually add a "0.0001" at the end of "length_parsed", otherwise an error will occur.
			// this is Ad hoc solution which needs to be improved later
			if (break_distance_to_snode <= length_parsed + 0.0001) 
			{
				min_seg_index = i;
				local_break_distance = break_distance_to_snode - length_parsed + seg_len;
				tempbreakpoint = GetBreakPointOnSegment(x1, y1, x2, y2, local_break_distance);
				break;
			}
			++i;
		};

		return tempbreakpoint;
	}; // GetBreakPointOnLinestring

	
	 // Given a (OGR)linestring and breakpoint at the linestring, and the seg idx correspond to the break point 
	 // cut the linestring into two (OGR)linestrings at the breakpoint, 
	 //return:
	 // @s_cutline: s_node to break point
	 // @e_cutline: break point to end node
	inline	void CutLinestringBasedOnBreakPoint(OGRLineString *linestring, int min_seg_index, OGRPoint* break_point,
										OGRLineString *s_cutline, OGRLineString * e_cutline) 
	{
		int numpoints = linestring->getNumPoints();
		
		//OGRPoint temppoint;
		OGRLineString s_cutedge;
		OGRLineString e_cutedge;

		for (int m = 0; m < (min_seg_index + 1); m++) 
		{
			double x1 = linestring->getX(m);
			double y1 = linestring->getY(m);
			s_cutedge.addPoint(x1,y1);
		}
		s_cutedge.addPoint(break_point);
		
		e_cutedge.addPoint(break_point);
		for (int m = (min_seg_index + 1); m < numpoints; m++) 
		{
			double x1 = linestring->getX(m);
			double y1 = linestring->getY(m);
			e_cutedge.addPoint(x1, y1);
		}
		*s_cutline = s_cutedge;
		*e_cutline = e_cutedge;
	}

	// Given a linestring and the break distance to the start node of linestring
	// cut the linestring into two lines at the breakpoint
	 //return:
	 // @s_cutline: s_node to break point
	 // @e_cutline: break point to end node
	inline	void CutLinestringBasedOnBreakDistance(OGRLineString *linestring, double break_distance_to_snode, OGRPoint *breakpoint,
											OGRLineString *s_cutline, OGRLineString *e_cutline)
	{
		int numpoints = linestring->getNumPoints();
		double length_parsed = 0;
		int i = 0;
		double seg_len;
		int min_seg_index = 0;
		double local_break_distanc;
		OGRPoint tempbreakpoint;
		double cutx, cuty;

		while (i < numpoints - 1)
		{
			double x1 = linestring->getX(i);
			double y1 = linestring->getY(i);

			double x2 = linestring->getX(i + 1);
			double y2 = linestring->getY(i + 1);
		
			seg_len = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
			length_parsed += seg_len;

			if (break_distance_to_snode <= length_parsed)
			{
				min_seg_index = i;
				local_break_distanc = break_distance_to_snode - length_parsed + seg_len;
				tempbreakpoint = GetBreakPointOnSegment(x1, y1, x2, y2, local_break_distanc);
				break;
			}
			++i;
		};

		*breakpoint = tempbreakpoint;

		CutLinestringBasedOnBreakPoint(linestring, min_seg_index, breakpoint, s_cutline, e_cutline);
	}

	// calculating the unit normal vector of two points
	// here no direct are specified
	inline	OGRPoint GetNormolizedVerticalVector(OGRPoint* pointA, OGRPoint* pointB)
	{
		OGRPoint normal_vec;

		double x1 = pointA->getX();
		double y1 = pointA->getY();
		double x2 = pointB->getX();
		double y2 = pointB->getY();

		double xd = x1 - x2;
		double yd = y1 - y2;

		if (xd == 0)
		{
			normal_vec.setX(1);
			normal_vec.setY(0);
		}
		else {
			double yvec = 1;
			// normal vector:  n = (a,b)  ax+by = 0; a = -(by/x), here the b is seting to 1
			double xvec = -yd / xd;
			double veclen = std::sqrt(xvec*xvec + 1);
			double n_xvec = xvec / veclen;
			double n_yvec = 1 / veclen;
			normal_vec.setX(n_xvec);
			normal_vec.setY(n_yvec);
		}
		return normal_vec;
	};

}// namespace algorithm
}// namespace tca

#endif /* TCA_ALGORITHM*/
