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

int main(int argc, char **argv);
int user_eg3(int argc, char **argv);

/*--------------------------------------------
-user_eg3-  main procedure of user_eg3 application
*/
/*int main(int argc, char **argv){

    QHULL_LIB_CHECK

    try{
        return user_eg3(argc, argv);
    }catch(QhullError &e){
        cerr << e.what() << std::endl;
        return e.errorCode();
    }
}//main*/

double gaussian(double x, double y)
{
	return exp(-x*x-y*y);
}

int compute_delaunay(double * arr, const int pc_dimension, const int n_pts)
{
Qhull qhull;
qhull.runQhull("", pc_dimension, n_pts, arr, "d Qbb Qt");

QhullFacetList facets = qhull.facetList();

std::vector<QhullFacet> facetVector = facets.toStdVector();

for ( int ifacet = 0; ifacet < facetVector.size(); ifacet++ )

//for ( int ifacet = 0; ifacet < facetVector.size(); ifacet++ )
//{
//	cout << "facet #" << ifacet+1 << " of " << facetVector.size() << endl;
//	QhullVertexSet verts = facetVector[ifacet].vertices();
//	std::vector<QhullVertex> vertVector = verts.toStdVector();
//	for ( int ivertex = 0; ivertex < vertVector.size(); ivertex++ )
//	{
//		cout << "   id() = " << vertVector[ivertex].id() << ": "
//			<< ivertex+1 << "   " << vertVector[ivertex]
//			<< "; " << (vertVector[ivertex].point()).id() << "   "
//			<< (vertVector[ivertex].point()).size() << "   "
//			<< (vertVector[ivertex].point()).count() << endl;
//		const size_t pointID = (vertVector[ivertex].point()).id();
//		cout << "   Check: " << arr[2*pointID] << "   " << arr[2*pointID+1] << endl;
//	}
//}


    return 0;
}//user_eg3

