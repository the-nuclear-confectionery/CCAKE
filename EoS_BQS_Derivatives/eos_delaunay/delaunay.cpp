#include "libqhullcpp/PointCoordinates.h"
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/QhullUser.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"

#include <cstdio>   /* for printf() of help message */
#include <fstream>
#include <iomanip> // setw
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::vector;

using orgQhull::PointCoordinates;
using orgQhull::Qhull;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullFacetListIterator;
using orgQhull::QhullFacetSet;
using orgQhull::QhullFacetSetIterator;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;
using orgQhull::QhullPointsIterator;
using orgQhull::QhullQh;
using orgQhull::QhullUser;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;
using orgQhull::QhullVertexSet;
using orgQhull::QhullVertexSetIterator;
using orgQhull::RboxPoints;

int compute_delaunay(double * arr, const size_t pc_dimension, const size_t n_pts,
					 vector<vector<size_t> > & simplices)
{
	Qhull qhull;
	qhull.runQhull("", pc_dimension, n_pts, arr, "d Qbb Qt");
	
	/*qhull.outputQhull();
	if(qhull.hasQhullMessage())
	{
		cerr << "\nResults of qhull\n" << qhull.qhullMessage();
		qhull.clearQhullMessage();
	}*/

	for ( const auto & facet : qhull.facetList().toStdVector() )
	{
		vector<size_t> simplex;
		for ( const auto & vertex : facet.vertices().toStdVector() )
			simplex.push_back( vertex.point().id() );
		simplices.push_back( simplex );
	}

    return 0;
}


