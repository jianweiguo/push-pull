/*
 * Conflict-Coverage Optimization in periodic domain.
 * Created by Abdalla Gafar
 * 2015-06-17 (this file, the algorithm is a few months older)
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex>

#include "getopt.c"
#include "drand48.c"

#define VL(x) sqrt((x).squared_length())

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                                  Point;
typedef K::Vector_2                                                 Vector;
typedef std::vector<Point>                                          Points;
typedef std::vector<Vector>                                         Vectors;
struct VInfo {                                                                  // Information stored in vertices in DT.
    int id;                                                                     // Index of point in t-map
};
struct FInfo {                                                                  // Information stored in faces in DT.
    Point c;                                                                    // Circumcenter, to avoid repeated calculation
};
typedef CGAL::Triangulation_vertex_base_with_info_2<VInfo, K>       Tvb;
typedef CGAL::Triangulation_face_base_with_info_2<FInfo, K>         Tfb;
typedef CGAL::Triangulation_data_structure_2<Tvb,Tfb>               Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      DT;
typedef DT::Vertex_handle                                           VH;
typedef DT::Face_handle                                             FH;
typedef DT::Vertex_circulator                                       VC;
typedef DT::Face_circulator                                         FC;

unsigned *shuffle(const unsigned N) {                                           // Return a randomly ordered list.
    unsigned *list = new unsigned[N];
    for (int i = 0; i < N; i++) list[i] = i;
    for (unsigned i = 0; i < N-1; i++) {
        unsigned r = i + rand() % (N-1 - i);
        std::swap(list[i], list[r]);
    }
    return list;
}

inline double rnd(double &max) { return drand48() * max; }

const double TWO_PI = 6.28318530717959;
const double epsilon = 1e-12;
double maxNorm = 0.03;                                                      // From BNOT


struct Statistics {                                                             // This structure is adapted from psa code "https://code.google.com/p/psa/"
    double mindist;                                                             // Global minimum nearest neighbor distance
    double avgmindist;                                                          // Average minimum nearest neighbor distance
    double orientorder;                                                         // Bond orientation order
    double effnyquist;                                                          // Effective Nyquist frequency
    double oscillations;                                                        // Oscillations, please refer to Heck et al. for meaning of these measures
    double coverageRadius;                                                      // We added this to monitor the size of holes
    double N4, N5, N6, N7, N8;                                                  // We added this to track number of neighbors percentage
    double sdA;                                                                 // Normalized standard deviation of areas of Voronoi cells
    double centroidDist;                                                        // Average distance of points to centroids of their cells
    double positionNorm;
    Statistics() : mindist(0), avgmindist(0), orientorder(0),
    effnyquist(0), oscillations(0), coverageRadius(0),
    N4(0), N5(0), N6(0), N7(0), N8(0) {};
};

inline bool isInRect(const Point &p, const Point &BL, const Point &TR) {        // Check if a point is inside a given rectangle
    return BL.x() <= p.x() && BL.y() <= p.y() &&
    p.x() < TR.x() && p.y() < TR.y();
}

inline double triangleType(Point &p1, Point &p2, Point &p3) {
    double sq1 = (p3 - p2).squared_length();
    double sq2 = (p1 - p3).squared_length();
    double sq3 = (p2 - p1).squared_length();
    if (sq1 < sq2) std::swap(sq1, sq2);
    if (sq1 < sq3) std::swap(sq1, sq3);
    return sq1 - (sq2 + sq3);                                                   // 0 for right, < 0 for acute, > 0 for obtuse
}

inline double triangleType(FC &fc) {
    return triangleType(
        fc->vertex(0)->point(),
        fc->vertex(1)->point(),
        fc->vertex(2)->point()
    );
}

class CPointSet {
private:
    int n;                                                                      // number of points in one period on an AABitmap
    double dhex;                                                                // Reference spacing of hexagonal packing.
    double rel_dmin;                                                     // Relative minimum distance between points; twice the conflict radius
    double rel_rc;                                                       // Relative maximum coverage radius.
    double sdA;                                                   // Target standard deviation of cell areas; default is from Schlomer thesis p 64.
    bool allStable;                                                             // To implement termination criteria
    double ONE;                                                           // Make it possible to use another size for torroidal domain.
    double HALF;
    DT dt;                                                                      // Will maintain a Delaunay triangulation for relaxation and queries
    double maxShift;
    void updateFaceInfo() {
        DT::Finite_faces_iterator fit = dt.finite_faces_begin();
        for ( ; fit != dt.finite_faces_end(); fit++) {                          // Iterate through all (finite) faces in triangulation
            fit->info().c = dt.circumcenter(fit);                               // Circumcenter of face is one end of a Voronoi edge
        }
    };
    struct TSite {                                                              // This is to handle the Delaunay triangulation.
        Point p;
        VH vh[9];                                                               // Handles to up to 9 replicas to maintain toroidal domain
        Vector force;                                                           // Resultant attraction/repulsion force exerted on this point by neighbors
        bool isStable, becomeStable;                                            // Used during optimization to skip processing a point if neither it nor neighbors move.
        inline double x() { return p.x(); };
        inline double y() { return p.y(); };
    };
    std::vector<TSite> s;                                                       // The final coordinates of points.
    inline double toroidallinearDist(double x1, double x2) const {              // 1D Nearest distance between replicas of two points
        double dx = x1 - x2;                                                    // Find distance in primary period
        while (dx > HALF) dx -= ONE;                                            // If larger than half the period length another replica is certainly nearer
        while (dx < -HALF) dx += ONE;                                           // Same, but opposite ordering of points
        return dx;
    };
    inline double toroidalSqDist(Point &p1, Point &p2) const {                  // 2D Nearest distance; square to avoid costly square root operation
        double dx = toroidallinearDist(p1.x(), p2.x());
        double dy = toroidallinearDist(p1.y(), p2.y());
        return dx * dx + dy * dy;
    };
    inline double toroidalDist(Point &p1, Point &p2) const {                    // 2D Nearest distance; square to avoid costly square root operation
        return sqrt( toroidalSqDist(p1, p2) );
    };
    inline Point mainReplica(Point &p) {
        double x = p.x(), y = p.y();
        while (x < 0) x += ONE;
        while (x >= ONE) x -= ONE;
        while (y < 0) y += ONE;
        while (y >= ONE) y -= ONE;
        return Point(x, y);
    };

    inline Point replica(Point &p, int i) {                                     // Find one of the 9 replicas of a point
        i = (i+4) % 9;                                                          // We make the middle replica at index 0
        double x = p.x() + (i%3 - 1) * ONE;                                     // Add -ONE, 0, or ONE to x
        double y = p.y() + (i/3 - 1) * ONE;                                     // Same for y
        return Point(x, y);
    };
    Point marginBL, marginTR;                                                   // Points within these margins affect points in the main replica
    Point setSite(int index, Point p);
	void moveSite(int index, Point p);

    inline void moveSite(int index, Vector shift) {                             // A handy function to shift points relative to their current location
        if (shift.squared_length() > epsilon) {
            moveSite(index, s[index].p + shift);
        }
    };
    // Relaxation forces:
    typedef Vector (CPointSet::*forceFunction) (int);                           // A pointer to a member function for calculating forces & return shift vectors
    Vector centroid(int index);
    Vector conflict(int index);
    Vector coverage(int index);
	Vector capacitySerial(int index);

    inline Point normalize(Point p) {                                           // Scale to unit & wrap around if points were shifted outside primary period
        p = mainReplica(p);
        double x = p.x() / ONE;
        double y = p.y() / ONE;
        return Point(x, y);
    };
    void initRandom();
    void initDarts();
    void initJittered();
    void initGrid();

public:
	std::string outputPath;                                                // Path for output files; could be set to a folder in /tmp
    int stableCount;
    CPointSet(int number_of_points, int initType = 0);                          // 0: random, 1: darts, 2: jittered grid, 3: regular grid

	double PPO_serial(std::string seq, double scale = 1);

    void setOutputPath(std::string path) { outputPath = path; };                // Let the user choose output folder. May be /tmp/??
    void setdmin(double d) { rel_dmin = d; };                                   // Set target NND for spring().
    void setRc(double r) { rel_rc = r; };
    void set_sdA(double sd) { sdA = sd; }
    void plotEPS(std::string fileName);                                         // Plot EPS file of points
    void printText(std::string fileName);                                       // Generate a text printout
    bool isAllStable() { return allStable; };
    double getMaxShift() { return sqrt(maxShift); }
    Statistics GetStatistics();
    void setAllUnstable() {
        for (int i = 0; i < n; i++) s[i].isStable = false;
        allStable = false;
    }

	void statis_multi_sets(std::string filename_base);
};

CPointSet::CPointSet(int number_of_points, int initType) {

	rel_dmin = 0.87;  
	rel_rc = 0.65;    
	sdA = 0.038600518;          
	//ONE = 1.0;        
	//HALF = 0.5;
	outputPath = "";

    n = number_of_points;                                                   // Number of points in one period
    ONE = sqrt(n);
    HALF = 0.5 * ONE;
    dhex = ONE * sqrt(2/(sqrt(3) * n));                                     // Maximum packing distance
    double margin = std::min(10/sqrt(n), 1.0);                              // Margin for toroidal domain. Heuristically, use 10 layers of points.
    marginBL = Point(-margin * ONE, -margin * ONE);                              // Bottom-left of primary period + margin
    marginTR = Point((1+margin) * ONE, (1+margin) * ONE);                        // Top-right. In our convention BL is included, TR is excluded
    s.reserve(n);                                                           // The final location including shifts
    switch (initType) {
        case 0: initRandom(); break;
        case 1: initDarts(); break;
        case 2: initJittered(); break;
        case 3: initGrid(); break;
        default:
            fprintf(
                stderr, "Undefined initialization option %d\n", initType);
            exit(1);
    }
}

/*
generate intial point set randomly
*/
void CPointSet::initRandom() {
    for (int i = 0; i < n; i++) {
        Point p (rnd(ONE), rnd(ONE));
        setSite(i, p);
    }
}

/*
generate intial point set by Poisson-disk darts
*/
void CPointSet::initDarts() {
    const int MAX = 10000;                                                  // To avoid infinite loop
    double r = 0.75 * dhex;
    double rr = r * r;
    for (int i = 0; i < n; i++) {
        Point p;
        bool accepted;
        do {
            p = Point(rnd(ONE), rnd(ONE));
            accepted = true;
            int counter = 0;
            for (int j = 0; (j < i) && accepted; j++) {
                accepted = accepted && (toroidalSqDist(p, s[j].p) > rr);
            }
            if (counter++ >= MAX) {
                fprintf(stderr, "Darts failed at %d\n", i);
                exit(1);
            }
        } while (!accepted);
        setSite(i, p);
    }
}

/*
generate intial point set by jitter
*/
void CPointSet::initJittered() {
    int m = sqrt(n);
    for (int i = 0; i < n; i++) {
        double x = i % m + 2 * drand48();
        double y = i / m + 2 * drand48();
        setSite(i, Point(x, y));
    }
}

/*
generate intial point set regularly
*/
void CPointSet::initGrid() {
    int m = sqrt(n);
    for (int i = 0; i < n; i++) {
        double x = i % m + 0.5;
        double y = i / m + 0.5;
        setSite(i, Point(x, y));
    }
}

Vector CPointSet::centroid(int index) {                                         // See http://en.wikipedia.org/wiki/Centroid
    double a = 0, cx = 0, cy = 0;                                               // Cell area (actually twice the area) and centroid coordinates
    double XProduct;                                                            // Cross product of vectors to adjacent vertices
    FC fc = dt.incident_faces(s[index].vh[0]), done(fc);
    do {
        Point p1 = dt.circumcenter(fc);
        Point p2 = dt.circumcenter(++fc);
        XProduct = p1.x() * p2.y() - p1.y() * p2.x();
        a += XProduct;                                                          // Accumulate areas
        cx += (p1.x() + p2.x()) * XProduct;
        cy += (p1.y() + p2.y()) * XProduct;
    } while (fc != done);
    cx /= 3.0 * a;
    cy /= 3.0 * a;
    return Point(cx, cy) - s[index].p;                                          // Return shift from current position to centroid
};

Point CPointSet::setSite(int index, Point p) {                                         // Set location of the indexed point (in t-domain) and insert it in triangulation
	p = mainReplica(p);
	s[index].p = p;                                                         // Save a handy copy of point coordinates
	s[index].isStable = false;
	for (int i = 0; i < 9; i++) {                                           // We loop through the 9 replicas,
		if (isInRect(replica(p, i), marginBL, marginTR)) {                  // if the location of a replica is within margin
			s[index].vh[i] = dt.insert(replica(p, i));                      // insert replica in triangulation and keep handle to it
			s[index].vh[i]->info().id = index;                              // Point the DT point back to map entry
		}
		else s[index].vh[i] = NULL;
	}
	return p;
};

void CPointSet::moveSite(int index, Point p) {                                         // Adjust location of indexed point (in t-domain) and update in triangulation
	double l = (p - s[index].p).squared_length();
	maxShift = std::max(maxShift, l);
	s[index].p = p;                                                         // Save a handy copy of updated point coordinates
	for (int i = 0; i < 9; i++) {
		if (s[index].vh[i] != NULL)
			s[index].vh[i] = dt.move(s[index].vh[i], replica(p, i));
	}
	s[index].becomeStable = false;                                          // Mark the point instable
	VC vc = dt.incident_vertices(s[index].vh[0]), done(vc);
	do {                                                                    // Mark neighbors instable
		s[vc->info().id].becomeStable = false;
	} while (++vc != done);
	allStable = false;                                                      // Mark the whole point set instable
};

Vector CPointSet::capacitySerial(int i) {                                       // Immediately update a site and neighbors
    double d[20], el[20], pressure;
    double sum_w = 0;
    double a = 0;                                                               // Area of Voronoi cell
    double XProduct;
    FC fc2 = dt.incident_faces(s[i].vh[0]), fc1(fc2++);                         // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
    VC vc = dt.incident_vertices(s[i].vh[0], fc2), done(vc);                    // The vertex sharing fc1 anf fc2 with v[i].vh
    int m = 0;                                                                  // Number of neighbors
    Vector dir[20];                                                             // Direction vectors to neighbors
    int id[20];                                                                 // Id's of neighbors. We can't use the circulator for updating
    do {
        Point c1 = dt.circumcenter(fc1), c2 = dt.circumcenter(fc2);             // Circumcenters of faces are endpoints of Voronoi cell edge
        XProduct = c1.x() * c2.y() - c1.y() * c2.x();
        a += XProduct;                                                          // Accumulate areas
        el[m] = sqrt((c2 - c1).squared_length());                               // Length of Voronoi edge
        dir[m] = (vc->point() - s[i].p);
        d[m] = sqrt(dir[m].squared_length());                                   // Distance to neighbor (= 2 x distance to Voronoi edge)
        dir[m] = dir[m] / d[m];                                                 // Normalize direction vector
        id[m] = vc->info().id;
        ++fc1;
        ++fc2;
        ++m;
    } while (++vc != done);
    a /= 2;
    double dA = a - 1;                                                          // Required expansion or contraction
    if (fabs(dA) > sdA) {
        for (int j = 0; j < m; j++) {
            sum_w += el[j] * el[j];
        }
        pressure = -2 * dA / sum_w;
        for (int j = 0; j < m; j++) {                                               // Loop again through neighbors to give each an appropriately sized nudge
            Vector force = pressure * el[j] * dir[j];
            moveSite(id[j], force);
        }
    }
    return Vector(0, 0);
}

// test conflict resolution by pushing neighbors
Vector CPointSet::conflict(int index) {
    double dmin = rel_dmin * dhex;
    Point &p = s[index].p;
    bool conflict[30];
    Vector shift[30];
    int id[30];
    int m = 0;
    VC vc = dt.incident_vertices(s[index].vh[0]), done(vc);
    do {
        Vector edge = vc->point() - p;
        double l = VL(edge);
        if (l < dmin) {
            conflict[m] = true;
            shift[m] = (1.001 * dmin/l - 1) * edge;
            id[m] = vc->info().id;
        } else conflict[m] = false;
        m++;
    } while (++vc != done);
    for (int i = 0; i < m; i++) {
        if (conflict[i]) moveSite(id[i], shift[i]);
    }
    return Vector(0, 0);
}

// A coverage routine which pulls neighbors
Vector CPointSet::coverage (int index) {                                        // Apply push and pull shifts to optimze conflict and coverage, respectivly
    double rc = rel_rc * dhex;
    double dmin = rel_dmin * dhex;
    int id[30];
    double scale[30];
    Vector edge[30];
    int m = 0;
    FC fc = dt.incident_faces(s[index].vh[0]), done(fc);
    VC vc = dt.incident_vertices(s[index].vh[0], fc);
    do {
        vc++;
        edge[m] = s[index].p - vc->point();
        id[m] = vc->info().id;
        if (triangleType (fc) <= 0) {
            Point c = dt.circumcenter(fc);
            double l = VL(c - s[index].p);
            scale[m] = rc / l;
        } else scale[m] = 2;                                                    // > 1
        m++;
    } while (++fc != done);
    for (int i = 0; i < m; i++) {
        double scl = std::min(scale[i], scale[(i+1)%m]);
        if (scl < 1) {
            Vector shift = (1 - scl) * edge[i];
            moveSite(id[i], shift);
        }
    }
    return Vector(0, 0);
}

double CPointSet::PPO_serial(std::string seq, double scale) {
    maxShift = 0;
    for (int i = 0; i < n; i++) s[i].becomeStable = true;                       // Assume all points will be found stable.
    allStable = true;                                                           // Assume the point set will be found stable.
    unsigned *order = shuffle(n);
    for (int i = 0; i < n; i++) {                                               // Iterate through numbers up to maximum index
        int index = order[i];                                                   // Translate index according to the random order
        if (s[index].isStable) continue;                                        // Skip stable points
        for (int k = 0; k < seq.length(); k++) {
            switch (seq[k]) {
                case '0': coverage(index); break;
                case '1': conflict(index); break;
                case '2': capacitySerial(index); break;
            }
        }
    }
    delete[] order;
    stableCount = 0;
    for (int i = 0; i < n; i++) {
        s[i].isStable = s[i].becomeStable;
        if (s[i].becomeStable) stableCount++;
    }
    return sqrt(maxShift);                                                      // Return maximum shift which can be used to monitor convergence
}

void CPointSet::printText(std::string fileName) {                               // Generate a text printout
	std::string fname = outputPath + fileName;
	const char *fullFileName = fname.c_str();
    FILE *file = fopen(fullFileName, "w");
    if (!file) {
        fprintf(stderr, "Failed to open %s\n", fullFileName);
        exit(1);
    }
    fprintf(file, "%d\n", n);                                                   // Print number of points in first line; a common convention
    for (int i = 0; i < n; i++) {                                               // Loop through all points
        Point p = normalize(s[i].p);                                            // Normalize and wrap back to unit torus
        fprintf(file, "%0.12f %0.12f\n", p.x(), p.y());
    }
    fclose(file);
}

Statistics CPointSet::GetStatistics()  {                                        // Return statistics useful to monitor progress of optimization
    Statistics stats;                                                           // This function is adapted from PSA code.
    stats.mindist = ONE;                                                        // The minimum distance can't go above ONE :)
    stats.avgmindist = 0;                                                       // Reset this to 0 to accumulate for averaging
    stats.orientorder = 0;                                                      // ~
    std::complex<double> acc = 0;
    unsigned long nacc = 0;
    for (int i = 0; i < n; i++) {
        double localmd = ONE;
        std::complex<double> localacc = 0;
        VC vc = dt.incident_vertices(s[i].vh[0]),
            done(vc), next;
        do {
            next = vc; ++next;
            const Point &v1 = vc->point();
            const Point &v2 = next->point();
            // Local mindist
            double dist = CGAL::squared_distance(s[i].p, v1);
            localmd = std::min(localmd, (double) dist);
            // Orientational order
            std::complex<double> c1(v1.x(), v1.y());
            std::complex<double> c2(v2.x(), v2.y());
            localacc += std::polar(1.0, 6.0 * arg(c1 - c2));                    // This is how it's done in PSA (divide only at end); but I feel it needs review
            ++nacc;
        } while (++vc != done);
        stats.mindist = std::min(stats.mindist, localmd);
        stats.avgmindist += sqrtf(localmd);
        acc += abs(localacc);
    }
    // -------- Coverage Radius
    double Rc = 0;
    Point BL(0, 0), TR(ONE, ONE);                                               // We will consider only faces in one period
    DT::Finite_faces_iterator it;
    for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++) {     // Iterate through all finite faces in triangulation
        if (isInRect(it->vertex(0)->point(), BL, TR)) {                         // It suffices to have one vertex within main period to consider the face
            Point c = dt.circumcenter(it);                                      // Farthest point in current face
            double r = CGAL::squared_distance(c, it->vertex(0)->point());       // Farthest distance in current face
            if (r > Rc) { Rc = r; }                                             // If larger than current farthest distance update the latter
        }
    }
    // -------- Voronoi cell N-Gon types:                                   // Before orientation order this used to be a measure; see paper by Balzer et al.
    int histogram[10] = {0,0,0,0,0,0,0,0,0,0};                              // A simple histogram, where index 0 is interpreted as 3
    for (int i = 0; i < n; i++) {                                          // Iterate through all points
        VC vc = dt.incident_vertices(s[i].vh[0]), done(vc);                 // We did not call neighbors to save the overhead of allocating Points
        int n = 0; do { n++; } while (++vc != done);                        // Count number of neighbors
        ++histogram[n - 3];                                                 // Increment the respective histogram bin; triangle is the smallest possible cell
    }
    // -------- capacity variations
    updateFaceInfo();                                                       // Update circumcenters
    double aa = 0;
    for (int i = 0; i < n; i++) {                                           // First we need to calculate pressure of each cell
        double a = 0;                                                       // Area of Voronoi cell
        double XProduct;
        FC fc2 = dt.incident_faces(s[i].vh[0]), fc1(fc2++), done(fc1);      // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
        do {
            Point &c1 = fc1->info().c, &c2 = fc2->info().c;                 // Circumcenters of faces are endpoints of Voronoi cell edge
            XProduct = c1.x() * c2.y() - c1.y() * c2.x();
            a += XProduct;                                                  // Accumulate areas
            ++fc2;
        } while (++fc1 != done);
        a /= 2;
        aa += a * a;                                                        // We track convergence by monitoring variance in areas
    }
    double areaVariance = (aa / n) - 1;                                     // E(a^2) - (E(a)) ^ 2; where E(a) = 1;

    // position gradient; according to BNOT:
    double norm = 0;
    for (int i = 0; i < n; i++) {
        Vector dist = centroid(i);
        norm += dist.squared_length();
    }

    // -------- Aggregate statistics:
    stats.mindist = sqrtf(stats.mindist) / dhex;                            // Nearest Neighbor Distance (NND) should
    stats.avgmindist /= n * dhex;                                           // be relative to dhex, the NND of hexagonal lattice
    stats.orientorder = abs(acc) / nacc;                                    // I don't quite get why we should average this once rather than average averages
    stats.coverageRadius = sqrtf(Rc) / dhex;                                // Size of largest hole (cf dmin); there should be an average too (cf davg)!
    stats.N4 = (double)histogram[1] / n;                                    // Ratios of various Voronoi cell types
    stats.N5 = (double)histogram[2] / n;                                    // 3 and 9+ are quite unlikely after relaxation
    stats.N6 = (double)histogram[3] / n;
    stats.N7 = (double)histogram[4] / n;
    stats.N8 = (double)histogram[5] / n;
    stats.sdA = sqrt(areaVariance);
    stats.positionNorm = sqrt(norm/n);
    return stats;
}

void CPointSet::plotEPS(std::string fileName) {                                        // Plot EPS file of points
	std::string fname = outputPath + fileName;
	const char *fullFileName = fname.c_str();

    FILE *epsfile = fopen(fullFileName, "w");
    if (!epsfile) {
        fprintf(stderr, "Failed to open %s\n", fullFileName);
        exit(1);
    }
    double radius = 0.1 * dhex;
    fprintf(epsfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(epsfile, "%%%%BoundingBox: 0 0 1000 1000\n");                   // Plot a single period
    fprintf(epsfile, "/ONE %f def\n", ONE);
    fprintf(epsfile, "1000 dup scale\n");                                   // We scale down to unit torus to follow the BNOT convention
    fprintf(epsfile, "/r %f def\n", radius / ONE);
    fprintf(epsfile, "/p {r 0 360 arc fill} def\n");                        // A typical PS function for point plotting
    for (int i = 0; i < n; i++) {                                           // Loop through all points
        Point p = normalize(s[i].p);                                        // Normalize and wrap back to unit torus
        fprintf(epsfile, "%f %f p\n", p.x(), p.y());                        // Print out
    }
    fprintf(epsfile, "showpage\n");                                         // Some viewer (e.g. new version of evince) may need this to display properly
    fclose(epsfile);                                                        // We need to close the file explicitly coz we are generating so many files
}

const char *USAGE_MESSAGE = "Usage: %s [options] <n>\n"
    "Options:\n"
    "-I <initialization type: 0: random, 1: darts, 2: jittered grid, 3: regular grid>\n"
    "-i <number of iterations>\n"
    "-p <path/prefix of output files>\n"
    "-s <scale applied to forces>\n"
    "-d <target dmin>\n"
	"-r <target coverage radius>\n"
	"-v <normalized maximum deviation of cell areas>\n"
    "-e order eps plots\n"
    "-t order text printouts\n"
    "-m order printout of measurements\n"
    "-q <forces sequence>\n"
    "-R <rounds>\n"
;

int main(int argc,char **argv) {
    srand(time(NULL)); srand48(time(NULL));                                     // Random seeds to random number generators (int and real)
    int opt;                                                                    // For use by getopt, the command line options utility
    int iterations = 1;                                                         // Number of iterations to apply
    char *fileNamePrefix = (char *)"";                                          // Current working directory is default output directory
    bool printText = false;                                                     // Whether to make text printouts
    double scale = 1.0;                                                         // Scale applied to forces. From experience We prefer damping the forces (~5%).
    bool serial = false;                                                        // Whether to apply shifts serially (FPO style), we found it to improve isotropy
    bool printEPS = false;                                                      // Whether to make EPS plots
    bool printMeasurments = false;                                              // Whether to print measurements after each iteration. Do not set if profiling speed.
    double dmin = 0.87;                                                         // Target dmin for use by our spring relaxation
    double rc = 0.65;                                                           // Target rc for use by coverage relaxation
    int forceType = 0;                                                          // The force to use; default is cco.
    int initialization = 0;
    double sdA = -1;                                                            // Target sdA.
    std::string seq = "012";
    int rounds = 1;
    while ((opt = getopt(argc, argv, "I:i:p:k:s:d:r:v:q:x:R:etSm")) != -1) {        // Modify default settings with command line options
        switch (opt) {
            case 'I': initialization = atoi(optarg); break;
            case 'i': iterations = atoi(optarg); break;
            case 'p': fileNamePrefix = optarg; break;
            case 's': scale = atof(optarg); break;
            case 'd': dmin = atof(optarg); break;
            case 'r': rc = atof(optarg); break;
            case 'e': printEPS = true; break;
            case 't': printText = true; break;
            case 'm': printMeasurments = true; break;
            case 'v': sdA = atof(optarg); break;
            case 'q': seq = optarg; break;
            case 'R': rounds = atoi(optarg); break;
            default: fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
        }
    }
    if (optind > argc - 1) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
    }
    int n = atoi(argv[optind]);
    fprintf(stderr, "Initializing ..");
    CPointSet ps(n, initialization);                                            // Create the point set
    fprintf(stderr, ". done\n");

    ps.setOutputPath(fileNamePrefix);                                           // Set path to output files (typically a temporary folder)
	ps.setdmin(dmin);                                                           // Set target dmin
	ps.setRc(rc);                                                               // Set target coverage radius
	if (sdA >= 0) ps.set_sdA(sdA);

	
    if (printEPS) { ps.plotEPS("original.eps"); }                               // If EPS plots are ordered, make a plot prior to optimization
    if (printText) { ps.printText("original.txt"); }                            // If text prints are ordered, make a printout prior to optimization

	int i;
	clock_t tOpt0 = clock();
	ps.setAllUnstable();
	for (i = 0; i < iterations; i++) {                                          // The main optimization loop
		int ft = forceType;
		std::string q = seq;
		
		ps.PPO_serial(seq);

		double iterationMax = ps.getMaxShift();
		double ratio_stable = (double)ps.stableCount / n;

		fprintf(
			stderr,
			"%4d - stable = %6d, stable-ratio = %10.8f, max-step = %10.8f",
			i, ps.stableCount, ratio_stable, iterationMax
			);

		fprintf(stderr, "\n");

		if (ps.isAllStable()) break;
	}

    clock_t tOpt1 = clock();
    double totalTime = (double)(tOpt1 - tOpt0) / CLOCKS_PER_SEC;
    fprintf(stderr, "Total optimization time: %10.6fs\n", totalTime);

	if (printMeasurments) {
		Statistics stats = ps.GetStatistics();                              // Read statistics after this iteration
		fprintf(stderr,                                                     // Print stats. It could be convenient to pipeline stderr to file then "tail -f"
			"Statistics Info: sdA = %8f, "
			"N4-N8%% = {%.4f,%.4f,%.4f,%.4f,%.4f}, "
			"Rc = %.3f, Belta = %.3f, G-MD = %.3f, A-MD = %.3f, BOO = %.3f, "
			"Xnorm = %.6f\n",
			stats.sdA, 100 * stats.N4, 100 * stats.N5,
			100 * stats.N6, 100 * stats.N7, 100 * stats.N8,
			stats.coverageRadius, stats.coverageRadius / stats.mindist,
			stats.mindist, stats.avgmindist,
			stats.orientorder,
			stats.positionNorm
			);
	}

	char fName[200];                                                            // A buffer to compile output file names
	if (printEPS) { sprintf(fName, "pushpull2d_%.3f_%.3f_%02d.eps", dmin, rc, rounds); ps.plotEPS(fName); }     // Make EPS plot is ordered
	if (printText) { sprintf(fName, "pushpull2d_%.3f_%.3f_%02d.txt", dmin, rc, rounds); ps.printText(fName); }  // Make text printout if ordered

}
