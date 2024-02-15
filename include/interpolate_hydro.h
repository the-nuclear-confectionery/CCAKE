#ifndef INTERPOLATE_HYDRO_H
#define INTERPOLATE_HYDRO_H

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "delaunay/table_delaunay.h"

using namespace std;

constexpr double EPS = 1e-10;
constexpr double INF = 1e10;

class Delaunay2D
{
	private:
		vector<vector<double> > vertices, faces;
		vector<vector<int> > face_indices;
		vector<vector<double> > value_fields;
		ostream & out;
		ostream & err;

		inline double d2( double x, double y, double x0, double y0 )
		{
			return ( ( x - x0 )*( x - x0 ) + ( y - y0 )*( y - y0 ) );
		}

		inline double sign( const vector<double> & p1, const vector<double> & p2, const vector<double> & p3 )
		{
		    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);
		}
		
		bool pointInTriangle ( const vector<double> & pt, const vector<double> & v1,
		                       const vector<double> & v2, const vector<double> & v3 )
		{
		    double d1 = sign(pt, v1, v2);
		    double d2 = sign(pt, v2, v3);
		    double d3 = sign(pt, v3, v1);
		
		    bool has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
		    bool has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
		
		    return !(has_neg && has_pos);
		}
		
		// this function evaluates the plane passing through v1, v2, v3 at the point p
		double basic_interp( const vector<double> & p, const vector<double> & v1,
		                     const vector<double> & v2, const vector<double> & v3 )
		{
			double det_x = (v2[1]-v1[1])*(v3[2]-v1[2]) - (v3[1]-v1[1])*(v2[2]-v1[2]);
			double det_y = (v2[2]-v1[2])*(v3[0]-v1[0]) - (v3[2]-v1[2])*(v2[0]-v1[0]);
			double det_z = (v2[0]-v1[0])*(v3[1]-v1[1]) - (v3[0]-v1[0])*(v2[1]-v1[1]);
			return ( v1[2] - ( det_x*(p[0] - v1[0]) + det_y*(p[1] - v1[1]) ) / (det_z+1e-100) );
		}

		void initialize_all( const vector<vector<double> > & vertices_in,
                             const vector<vector<double> > & value_fields_in )
		{
			vertices = vertices_in;
			value_fields = value_fields_in;

			// generate grid itself
			int success = generate_delaunay_2d( vertices, faces );

			// associate vertices to the faces that contain them
			face_indices.resize( vertices.size() );
			for ( int iFace = 0 ; iFace < faces.size(); iFace++ )
			for ( int iVertex = 0 ; iVertex < 3; iVertex++ )
				face_indices[ faces[iFace][iVertex] ].push_back( iFace );


			return;
		}

	public:
		Delaunay2D( const vector<vector<double> > & vertices_in,
					const vector<vector<double> > & value_fields_in,
					ostream & out_stream = std::cout,
					ostream & err_stream = std::cerr )
					:
					out(out_stream),
					err(err_stream)
					{ initialize_all( vertices_in, value_fields_in ); };


		void interpolate( vector<double> & outPoint )
		{
			double x0 = outPoint[0];
			double y0 = outPoint[1];

			// find closest vertex to chosen point
			int iClosestVertex = 0;
			double closestDistance2 = d2( vertices[0][0], vertices[0][1], x0, y0 );
			for ( int iVertex = 0; iVertex < vertices.size(); iVertex++ )
			{
				double distance2 = d2( vertices[iVertex][0], vertices[iVertex][1], x0, y0 );
				if ( distance2 < closestDistance2 )
				{
					iClosestVertex = iVertex;
					closestDistance2 = distance2;
				}
			}

			// if chosen point *is* the closest vertex, just return that function value
			if ( sqrt(closestDistance2) < EPS )
			{
				outPoint[2] = value_fields[iClosestVertex][0];
				outPoint[3] = value_fields[iClosestVertex][1];
				outPoint[4] = value_fields[iClosestVertex][2];
				outPoint[5] = value_fields[iClosestVertex][3];
			}
			else
			{
				// otherwise, point must be in face containing closest vertex or outside
				// convex hull entirely
				bool no_solution_found_yet = true;
				vector<double> point = { x0, y0 };
				for ( const auto & face_index : face_indices[ iClosestVertex ] )
				{
					// look at next face
					const auto & face = faces[ face_index ];
					vector<double> v0 = vertices[face[0]];
					vector<double> v1 = vertices[face[1]];
					vector<double> v2 = vertices[face[2]];
					if ( pointInTriangle( point, v0, v1, v2 ) )
					{
						v0.resize(3);
						v1.resize(3);
						v2.resize(3);

						for (int iField = 0; iField < 4; iField++)
						{
							v0[2] = value_fields[ face[0] ][iField];
							v1[2] = value_fields[ face[1] ][iField];
							v2[2] = value_fields[ face[2] ][iField];

							outPoint[iField+2] = basic_interp( point, v0, v1, v2 );
						}
		
						no_solution_found_yet = false;
						break;
					}
				}

				if ( no_solution_found_yet )
				{
					//cout << "Result is out of bounds!" << endl;
					outPoint[2] = 0.0;
					outPoint[3] = 0.0;
					outPoint[4] = 0.0;
					outPoint[5] = 0.0;
				}
			}

			return;
		}
		
};

void read_in_data(vector<vector<double> > & points, string filename, int nHeaderLines = 0)
{
	points.clear();

	ifstream infile(filename.c_str());

	if (infile.is_open())
	{
		string line;
		int count = 0;
		while ( getline (infile, line) )
		{
			if ( count++ < nHeaderLines ) continue;

			istringstream iss(line);
			vector<double> point(6);
			for (int iCol = 0; iCol < 6; iCol++)
				iss >> point[iCol];

			// to hold distances
			point.push_back( 0.0 );
			
			points.push_back( point );
		}
	}

	infile.close();
	return;
}


bool sortByDistance( const vector<double> & v1, const vector<double> & v2 ) {
  return ( v1[6] < v2[6] );
}


void interpolate(vector<vector<double> > & points, vector<double> & outPoint)
{
	double x = outPoint[0];
	double y = outPoint[1];

	//the first col is x axis, the second col is y axis
	for ( auto & point : points )
	{
		double dx = point[0]-x, dy = point[1]-y;
		//point[6] = sqrt(dx*dx+dy*dy);	//distance
		point[6] = dx*dx+dy*dy;	//distance^2
		if (point[2] < EPS || point[3] < 0.15)
		{
			point[6] = INF; //if energy = 0, set distance to be large enough to eliminate the influence
			continue;
		}
	}
	
	// partial sort only smallest 8 distances
	//sort(points.begin(), points.end(), sortByDistance);
	std::nth_element( points.begin(), points.begin() + 7, points.end(), sortByDistance );
	
	int num = 8;
	double totalrate   = 0;
	double totaltemp   = 0;
	double totalvinx   = 0;
	double totalviny   = 0;
	double totalenergy = 0;
	
	for( int i = 0; i < num; i++)
	{
		const auto & point = points[i];

		// if the distance is 0, which means they are the same point,
		// copy all the info instead of doing interpolation
		if (point[6] < EPS)
		{
			totaltemp   = point[3];
			totalvinx   = point[4];
			totalviny   = point[5];
			totalenergy = point[2];
			break;
		}
		else
		{
			// f += pow(1.0 / point[3], r);
			double rate  = 1.0 / point[6];
			totalrate   += rate;
			totalenergy += point[2] * rate;
			totaltemp   += point[3] * rate;
			totalvinx   += point[4] * rate;
			totalviny   += point[5] * rate;
		}
	}
	
	totaltemp   /= totalrate;
	totalvinx   /= totalrate;
	totalviny   /= totalrate;
	totalenergy /= totalrate;
	
	/*std::cout <<"Interpolation temp " <<totaltemp << '\n';
	std::cout <<"Interpolation velocity in x  " <<totalvinx << '\n';
	std::cout <<"Interpolation velocity in y  " <<totalviny << '\n';
	std::cout <<"Interpolation energydensity " <<totalenergy << '\n';*/

	outPoint[2] = totalenergy;
	outPoint[3] = totaltemp;
	outPoint[4] = totalvinx;
	outPoint[5] = totalviny;
	
	return;
}

void interpolate_hydro_driver(
		string filename, vector<vector<double> > & outputGrid,
		const vector<double> & xGrid, const vector<double> & yGrid,
		int interpolation_mode, bool verbose = false )
{
	// this function conducts the interpolation
	// interpolation_mode: 	0 - Steven Li's interpolation algorithm (weighted average of 8 nearby points)
	// 						1 - Christopher Plumberg's modified Delaunay (bilinear) interpolation routine

	if (verbose) cout << endl << " - Reading in data from " << filename << endl;
	vector<vector<double> > points;
	read_in_data(points, filename, 1);

	const int xGridSize = xGrid.size();
	const int yGridSize = yGrid.size();
	
	if (verbose)
		cout << endl << " - Performing interpolation (interpolation_mode == "
			 << interpolation_mode << ")" << endl;

	// set up output grid and do the interpolation
	if ( interpolation_mode == 0 )
	{
		outputGrid = vector<vector<double> >( xGridSize*yGridSize, vector<double> ( 6 ) );
		int idx = 0;
		for ( int ix = 0; ix < xGridSize; ix++ )
		for ( int iy = 0; iy < yGridSize; iy++ )
		{
			outputGrid[idx][0] = xGrid[ix];
			outputGrid[idx][1] = yGrid[iy];
			interpolate( points, outputGrid[idx] );
			idx++;
		}
	}
	else if ( interpolation_mode == 1 )
	{
		// set Delaunay grid
		if (verbose) cout << endl << " - Constructing Delaunay grid" << endl;
		vector<vector<double> > vertices ( points.size(), vector<double>(2) );
		vector<vector<double> > value_fields ( points.size(), vector<double>(4) );
		for (int iVertex = 0; iVertex < vertices.size(); iVertex++)
		{
			vertices[iVertex][0]     = points[iVertex][0];
			vertices[iVertex][1]     = points[iVertex][1];
			value_fields[iVertex][0] = points[iVertex][2];
			value_fields[iVertex][1] = points[iVertex][3];
			value_fields[iVertex][2] = points[iVertex][4];
			value_fields[iVertex][3] = points[iVertex][5];
		}
	
		// create Delaunay interpolation object
		Delaunay2D interpolationGrid( vertices, value_fields );
	
		// set up output grid and do the interpolation
		if (verbose) cout << endl << " - Performing interpolation" << endl;
		outputGrid = vector<vector<double> >( xGridSize*yGridSize, vector<double> ( 6 ) );
		int idx = 0;
		for ( int ix = 0; ix < xGridSize; ix++ )
		for ( int iy = 0; iy < yGridSize; iy++ )
		{
			outputGrid[idx][0] = xGrid[ix];
			outputGrid[idx][1] = yGrid[iy];
			interpolationGrid.interpolate( outputGrid[idx] );
			idx++;
		}
	}
	else
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Invalid choice of interpolation_mode!" << endl;
		exit;
	}

	return;
}

#endif
