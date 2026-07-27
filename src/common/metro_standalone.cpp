// metro_standalone.cpp — Standalone Hausdorff distance between triangle meshes
// Self-contained reimplementation of VCGlib's metro tool. Zero external dependencies.
// Supports PLY (ascii + binary), OFF, OBJ mesh formats.
// Usage: metro_standalone file1 file2 [options]

// ============================================================================
// Section 1: Standard Library Includes
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cfloat>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <numeric>
#include <climits>

// ============================================================================
// Section 2: Vec3d — 3D Vector (double)
// ============================================================================
struct Vec3d {
    double x, y, z;
    Vec3d() : x(0), y(0), z(0) {}
    Vec3d(double a, double b, double c) : x(a), y(b), z(c) {}
    explicit Vec3d(double v) : x(v), y(v), z(v) {}

    double& operator[](int i) { return (&x)[i]; }
    double  operator[](int i) const { return (&x)[i]; }

    Vec3d operator+(const Vec3d& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3d operator-(const Vec3d& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3d operator*(double s) const { return {x*s, y*s, z*s}; }
    Vec3d operator/(double s) const { double inv=1.0/s; return {x*inv, y*inv, z*inv}; }
    Vec3d operator-() const { return {-x, -y, -z}; }

    Vec3d& operator+=(const Vec3d& o) { x+=o.x; y+=o.y; z+=o.z; return *this; }
    Vec3d& operator-=(const Vec3d& o) { x-=o.x; y-=o.y; z-=o.z; return *this; }
    Vec3d& operator*=(double s) { x*=s; y*=s; z*=s; return *this; }

    double dot(const Vec3d& o) const { return x*o.x + y*o.y + z*o.z; }
    Vec3d cross(const Vec3d& o) const {
        return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x};
    }
    double squaredNorm() const { return x*x + y*y + z*z; }
    double norm() const { return sqrt(squaredNorm()); }
    Vec3d normalized() const { double n=norm(); return n>0 ? *this*(1.0/n) : Vec3d(); }
};

inline Vec3d operator*(double s, const Vec3d& v) { return v*s; }

inline double Distance(const Vec3d& a, const Vec3d& b) { return (a-b).norm(); }
inline double SquaredDistance(const Vec3d& a, const Vec3d& b) { return (a-b).squaredNorm(); }

// ============================================================================
// Section 3: Box3d — Axis-Aligned Bounding Box (double)
// ============================================================================
struct Box3d {
    Vec3d min, max;
    bool null;

    Box3d() : null(true) {}

    void SetNull() { null = true; }
    bool IsNull() const { return null; }

    void Add(const Vec3d& p) {
        if (null) { min = max = p; null = false; }
        else {
            if (p.x < min.x) min.x = p.x; if (p.x > max.x) max.x = p.x;
            if (p.y < min.y) min.y = p.y; if (p.y > max.y) max.y = p.y;
            if (p.z < min.z) min.z = p.z; if (p.z > max.z) max.z = p.z;
        }
    }
    void Add(const Box3d& b) {
        if (b.null) return;
        Add(b.min); Add(b.max);
    }
    void Offset(double d) {
        if (null) return;
        min.x -= d; min.y -= d; min.z -= d;
        max.x += d; max.y += d; max.z += d;
    }
    double Diag() const { return null ? 0 : (max - min).norm(); }
    Vec3d Dim() const { return null ? Vec3d() : max - min; }
    bool IsInEx(const Vec3d& p) const {
        return !null && p.x>=min.x && p.x<max.x && p.y>=min.y && p.y<max.y && p.z>=min.z && p.z<max.z;
    }
    void Intersect(const Box3d& b) {
        if (null || b.null) { SetNull(); return; }
        if (min.x < b.min.x) min.x = b.min.x;
        if (min.y < b.min.y) min.y = b.min.y;
        if (min.z < b.min.z) min.z = b.min.z;
        if (max.x > b.max.x) max.x = b.max.x;
        if (max.y > b.max.y) max.y = b.max.y;
        if (max.z > b.max.z) max.z = b.max.z;
        if (min.x > max.x || min.y > max.y || min.z > max.z) SetNull();
    }
};

// ============================================================================
// Section 4-5: Point3i / Box3i — Integer Grid Types
// ============================================================================
struct Point3i {
    int v[3];
    Point3i() { v[0]=v[1]=v[2]=0; }
    Point3i(int a, int b, int c) { v[0]=a; v[1]=b; v[2]=c; }
    int& operator[](int i) { return v[i]; }
    int  operator[](int i) const { return v[i]; }
    Point3i operator-(const Point3i& o) const { return {v[0]-o.v[0], v[1]-o.v[1], v[2]-o.v[2]}; }
    Point3i operator+(const Point3i& o) const { return {v[0]+o.v[0], v[1]+o.v[1], v[2]+o.v[2]}; }
};

struct Box3i {
    Point3i min, max;
    bool null;
    Box3i() : null(true) {}
    Box3i(const Point3i& a, const Point3i& b) : min(a), max(b), null(false) {}
    void SetNull() { null = true; }
    bool IsNull() const { return null; }
    void Add(const Point3i& p) {
        if (null) { min = max = p; null = false; }
        else {
            for (int i=0;i<3;i++) { if(p[i]<min[i]) min[i]=p[i]; if(p[i]>max[i]) max[i]=p[i]; }
        }
    }
    void Intersect(const Box3i& b) {
        if (null || b.null) { SetNull(); return; }
        for (int i=0;i<3;i++) { if(min[i]<b.min[i]) min[i]=b.min[i]; if(max[i]>b.max[i]) max[i]=b.max[i]; }
        if (min[0]>max[0]||min[1]>max[1]||min[2]>max[2]) SetNull();
    }
};

// ============================================================================
// Section 6: Mesh Data Structures
// ============================================================================
struct Triangle { int v[3]; };

struct Mesh {
    std::vector<Vec3d> verts;
    std::vector<Triangle> faces;
    Box3d bbox;
    std::vector<float> vertQuality;  // per-vertex quality for error saving

    int vn() const { return (int)verts.size(); }
    int fn() const { return (int)faces.size(); }

    void computeBBox() {
        bbox.SetNull();
        for (size_t i = 0; i < verts.size(); i++)
            bbox.Add(verts[i]);
    }

    double faceDoubleArea(int fi) const {
        const Vec3d& v0 = verts[faces[fi].v[0]];
        const Vec3d& v1 = verts[faces[fi].v[1]];
        const Vec3d& v2 = verts[faces[fi].v[2]];
        return (v1 - v0).cross(v2 - v0).norm();
    }

    double computeArea() const {
        double area = 0;
        for (int i = 0; i < fn(); i++)
            area += faceDoubleArea(i);
        return area / 2.0;
    }
};

// ============================================================================
// Section 7: PointTriangleDist — Voronoi Region Method (Ericson)
// ============================================================================
// Returns distance from point pt to triangle (a,b,c). Sets closestPt.
inline double PointTriangleDist(const Vec3d& pt, const Vec3d& a, const Vec3d& b,
                                const Vec3d& c, Vec3d& closestPt)
{
    Vec3d ab = b - a;
    Vec3d ac = c - a;
    Vec3d ap = pt - a;

    double d1 = ab.dot(ap);
    double d2 = ac.dot(ap);
    if (d1 <= 0 && d2 <= 0) { closestPt = a; return (pt - a).norm(); }

    Vec3d bp = pt - b;
    double d3 = ab.dot(bp);
    double d4 = ac.dot(bp);
    if (d3 >= 0 && d4 <= d3) { closestPt = b; return (pt - b).norm(); }

    double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        closestPt = a + ab * v;
        return (pt - closestPt).norm();
    }

    Vec3d cp = pt - c;
    double d5 = ab.dot(cp);
    double d6 = ac.dot(cp);
    if (d6 >= 0 && d5 <= d6) { closestPt = c; return (pt - c).norm(); }

    double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        closestPt = a + ac * w;
        return (pt - closestPt).norm();
    }

    double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        closestPt = b + (c - b) * w;
        return (pt - closestPt).norm();
    }

    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;
    closestPt = a + ab * v + ac * w;
    return (pt - closestPt).norm();
}

// ============================================================================
// Section 8: Histogram
// ============================================================================
class Histogram {
    std::vector<double> H;  // bin counts
    std::vector<double> R;  // bin ranges
    double minv, maxv;
    int n;                   // number of valid intervals
    double cnt, sum, rmsVal;

public:
    Histogram() : minv(0), maxv(1), n(0), cnt(0), sum(0), rmsVal(0) {}

    void SetRange(double _minv, double _maxv, int _n) {
        Clear();
        minv = _minv; maxv = _maxv; n = _n;
        H.resize(n + 2, 0);
        R.resize(n + 3);
        R[0] = -DBL_MAX;
        R[n + 2] = DBL_MAX;
        double delta = maxv - minv;
        for (int i = 0; i <= n; i++)
            R[i + 1] = minv + delta * (double)i / n;
    }

    void Clear() {
        H.clear(); R.clear();
        cnt = sum = rmsVal = 0;
        n = 0; minv = 0; maxv = 1;
    }

    int BinIndex(double val) const {
        auto it = std::lower_bound(R.begin(), R.end(), val);
        if (it == R.begin()) return 0;
        if (it == R.end()) return (int)H.size() - 1;
        int pos = (int)(it - R.begin()) - 1;
        return pos;
    }

    void Add(double v) {
        int pos = BinIndex(v);
        if (pos < 0) pos = 0;
        if (pos >= (int)H.size()) pos = (int)H.size() - 1;
        H[pos] += 1.0;
        cnt += 1.0;
        sum += v;
        rmsVal += v * v;
    }

    double Avg() const { return cnt > 0 ? sum / cnt : 0; }
    double RMS() const { return cnt > 0 ? sqrt(rmsVal / cnt) : 0; }
    double Cnt() const { return cnt; }

    double Percentile(double frac) const {
        if (H.empty() || R.empty()) return 0;
        double s = 0;
        for (size_t i = 0; i < H.size(); i++) s += H[i];
        double target = s * frac;
        double partsum = 0;
        size_t i;
        for (i = 0; i < H.size(); i++) {
            partsum += H[i];
            if (partsum >= target) break;
        }
        if (i >= H.size()) i = H.size() - 1;
        return R[i + 1];
    }

    void FileWrite(const std::string& filename) const {
        FILE* fp = fopen(filename.c_str(), "w");
        if (!fp) return;
        for (size_t i = 0; i < H.size(); i++)
            fprintf(fp, "%12.8lf , %12.8lf \n", R[i], cnt > 0 ? H[i] / cnt : 0.0);
        fclose(fp);
    }
};

// ============================================================================
// Section 9: BestDim — Grid Dimension Computation
// ============================================================================
static void BestDim(long long elems, const Vec3d& size, Point3i& dim) {
    const long long mincells = 1;
    double diag = size.norm();
    double eps = diag * 1e-4;

    long long ncell = elems;
    if (ncell < mincells) ncell = mincells;

    dim[0] = dim[1] = dim[2] = 1;

    if (size[0] > eps) {
        if (size[1] > eps) {
            if (size[2] > eps) {
                double k = pow((double)(ncell / (size[0]*size[1]*size[2])), 1.0/3.0);
                dim[0] = (int)(size[0] * k); if (dim[0]<1) dim[0]=1;
                dim[1] = (int)(size[1] * k); if (dim[1]<1) dim[1]=1;
                dim[2] = (int)(size[2] * k); if (dim[2]<1) dim[2]=1;
            } else {
                dim[0] = (int)sqrt(ncell * size[0] / size[1]); if (dim[0]<1) dim[0]=1;
                dim[1] = (int)sqrt(ncell * size[1] / size[0]); if (dim[1]<1) dim[1]=1;
            }
        } else {
            if (size[2] > eps) {
                dim[0] = (int)sqrt(ncell * size[0] / size[2]); if (dim[0]<1) dim[0]=1;
                dim[2] = (int)sqrt(ncell * size[2] / size[0]); if (dim[2]<1) dim[2]=1;
            } else {
                dim[0] = (int)ncell; if (dim[0]<1) dim[0]=1;
            }
        }
    } else {
        if (size[1] > eps) {
            if (size[2] > eps) {
                dim[1] = (int)sqrt(ncell * size[1] / size[2]); if (dim[1]<1) dim[1]=1;
                dim[2] = (int)sqrt(ncell * size[2] / size[1]); if (dim[2]<1) dim[2]=1;
            } else {
                dim[1] = (int)ncell; if (dim[1]<1) dim[1]=1;
            }
        } else if (size[2] > eps) {
            dim[2] = (int)ncell; if (dim[2]<1) dim[2]=1;
        }
    }
}

// ============================================================================
// Section 10: UniformGrid — Spatial Index + Closest Face Query
// ============================================================================
class UniformGrid {
    struct Link {
        int faceIdx;
        int cellIdx;
        Link() : faceIdx(-1), cellIdx(0) {}
        Link(int f, int c) : faceIdx(f), cellIdx(c) {}
        bool operator<(const Link& o) const { return cellIdx < o.cellIdx; }
    };

    std::vector<Link> links;
    std::vector<int> gridStart; // gridStart[i] = index into links where cell i begins
    Box3d bbox;
    Point3i siz;
    Vec3d voxel;
    Vec3d dim;
    Mesh* mesh;

    // marker system: generation counter to avoid re-testing faces
    std::vector<int> faceMarker;
    int currentMark;

    void PToIP(const Vec3d& p, Point3i& ip) const {
        Vec3d t = p - bbox.min;
        ip[0] = (int)(t.x / voxel.x);
        ip[1] = (int)(t.y / voxel.y);
        ip[2] = (int)(t.z / voxel.z);
    }

    void BoxToIBox(const Box3d& b, Box3i& ib) const {
        PToIP(b.min, ib.min);
        PToIP(b.max, ib.max);
        ib.null = false;
    }

    int CellIndex(int x, int y, int z) const {
        return x + siz[0] * (y + siz[1] * z);
    }

public:
    UniformGrid() : mesh(nullptr), currentMark(0) {}

    void Build(Mesh& m) {
        mesh = &m;

        // Compute bbox from faces
        Box3d fbox;
        for (int fi = 0; fi < m.fn(); fi++) {
            for (int j = 0; j < 3; j++)
                fbox.Add(m.verts[m.faces[fi].v[j]]);
        }

        // Inflate slightly
        double infl = fbox.Diag() / m.fn();
        fbox.min -= Vec3d(infl);
        fbox.max += Vec3d(infl);

        bbox = m.bbox; // use the already-set combined bbox
        // But ensure it contains all faces
        bbox.Add(fbox);

        // Compute grid dimensions
        Vec3d sz = bbox.Dim();
        BestDim((long long)m.fn(), sz, siz);

        dim = sz;
        voxel = Vec3d(dim.x/siz[0], dim.y/siz[1], dim.z/siz[2]);

        int totalCells = siz[0] * siz[1] * siz[2];

        // Build links: for each face, compute its AABB in grid coords, insert into all overlapping cells
        links.clear();
        for (int fi = 0; fi < m.fn(); fi++) {
            Box3d fb;
            for (int j = 0; j < 3; j++)
                fb.Add(m.verts[m.faces[fi].v[j]]);
            fb.Intersect(bbox);
            if (fb.IsNull()) continue;

            Box3i ib;
            BoxToIBox(fb, ib);
            // Clamp
            for (int a=0;a<3;a++) {
                if (ib.min[a] < 0) ib.min[a] = 0;
                if (ib.max[a] >= siz[a]) ib.max[a] = siz[a] - 1;
            }

            for (int z = ib.min[2]; z <= ib.max[2]; z++)
                for (int y = ib.min[1]; y <= ib.max[1]; y++)
                    for (int x = ib.min[0]; x <= ib.max[0]; x++)
                        links.push_back(Link(fi, CellIndex(x, y, z)));
        }

        // Sentinel
        links.push_back(Link(-1, totalCells));

        // Sort by cell index
        std::sort(links.begin(), links.end());

        // Build grid start array
        gridStart.resize(totalCells + 1);
        int li = 0;
        for (int ci = 0; ci <= totalCells; ci++) {
            gridStart[ci] = li;
            while (li < (int)links.size() && links[li].cellIdx == ci)
                li++;
        }

        // Init markers
        faceMarker.assign(m.fn(), 0);
        currentMark = 1;
    }

    // Find closest face to point p. Returns distance, updates closestPt.
    // maxDist is the search upper bound; minDist is set to the found distance.
    int GetClosest(const Vec3d& p, double maxDist, double& minDist, Vec3d& closestPt) {
        minDist = maxDist;
        int winner = -1;
        currentMark++;

        double voxelNorm = voxel.norm();
        double newradius = voxelNorm;
        double radius;

        Box3i iboxdone;
        iboxdone.SetNull();

        // First check the cell containing p
        if (bbox.IsInEx(p)) {
            Point3i ip;
            PToIP(p, ip);
            // Clamp
            for (int a=0;a<3;a++) { if(ip[a]<0) ip[a]=0; if(ip[a]>=siz[a]) ip[a]=siz[a]-1; }

            int ci = CellIndex(ip[0], ip[1], ip[2]);
            for (int li = gridStart[ci]; li < gridStart[ci+1]; li++) {
                int fi = links[li].faceIdx;
                if (fi < 0) continue;
                Vec3d cp;
                double d = PointTriangleDist(p,
                    mesh->verts[mesh->faces[fi].v[0]],
                    mesh->verts[mesh->faces[fi].v[1]],
                    mesh->verts[mesh->faces[fi].v[2]], cp);
                if (d < minDist) {
                    minDist = d;
                    closestPt = cp;
                    winner = fi;
                    newradius = minDist;
                }
                faceMarker[fi] = currentMark;
            }
            iboxdone = Box3i(ip, ip);
        }

        // Expanding search
        Box3i ibox(Point3i(0,0,0), Point3i(siz[0]-1, siz[1]-1, siz[2]-1));

        do {
            radius = newradius;
            // Compute box around p with radius
            Box3d boxtodo;
            boxtodo.Add(Vec3d(p.x - radius, p.y - radius, p.z - radius));
            boxtodo.Add(Vec3d(p.x + radius, p.y + radius, p.z + radius));

            Box3i iboxtodo;
            BoxToIBox(boxtodo, iboxtodo);
            iboxtodo.Intersect(ibox);

            if (!iboxtodo.IsNull()) {
                for (int ix = iboxtodo.min[0]; ix <= iboxtodo.max[0]; ix++)
                    for (int iy = iboxtodo.min[1]; iy <= iboxtodo.max[1]; iy++)
                        for (int iz = iboxtodo.min[2]; iz <= iboxtodo.max[2]; iz++) {
                            // Skip already-done cells
                            if (!iboxdone.IsNull() &&
                                ix >= iboxdone.min[0] && ix <= iboxdone.max[0] &&
                                iy >= iboxdone.min[1] && iy <= iboxdone.max[1] &&
                                iz >= iboxdone.min[2] && iz <= iboxdone.max[2])
                                continue;

                            int ci = CellIndex(ix, iy, iz);
                            for (int li = gridStart[ci]; li < gridStart[ci+1]; li++) {
                                int fi = links[li].faceIdx;
                                if (fi < 0) continue;
                                if (faceMarker[fi] == currentMark) continue;
                                faceMarker[fi] = currentMark;

                                Vec3d cp;
                                double d = PointTriangleDist(p,
                                    mesh->verts[mesh->faces[fi].v[0]],
                                    mesh->verts[mesh->faces[fi].v[1]],
                                    mesh->verts[mesh->faces[fi].v[2]], cp);
                                if (d < minDist) {
                                    minDist = d;
                                    closestPt = cp;
                                    winner = fi;
                                }
                            }
                        }
            }

            if (winner < 0) newradius = radius + voxelNorm;
            else newradius = minDist;

            iboxdone = iboxtodo;
        } while (minDist > radius);

        return winner;
    }
};

// ============================================================================
// Section 11: Mesh I/O Parsers
// ============================================================================

static std::string toLower(const std::string& s) {
    std::string r = s;
    for (char& c : r) c = (char)tolower((unsigned char)c);
    return r;
}

static std::string getExtension(const std::string& filename) {
    size_t dot = filename.rfind('.');
    if (dot == std::string::npos) return "";
    return toLower(filename.substr(dot + 1));
}

// ---- PLY Parser ----
static bool loadPLY(const std::string& filename, Mesh& m) {
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) return false;

    char line[4096];
    // Read "ply"
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); return false; }
    if (strncmp(line, "ply", 3) != 0) { fclose(fp); return false; }

    enum Format { ASCII, BIN_LE, BIN_BE } format = ASCII;
    int nverts = 0, nfaces = 0;
    bool inVertex = false, inFace = false;

    struct PropInfo {
        std::string name;
        std::string type;
        bool isList;
        std::string listCountType;
        std::string listElemType;
    };
    std::vector<PropInfo> vertProps, faceProps;
    int vertXIdx=-1, vertYIdx=-1, vertZIdx=-1;

    // Parse header
    while (fgets(line, sizeof(line), fp)) {
        std::string sline(line);
        // Remove trailing whitespace
        while (!sline.empty() && (sline.back()=='\n'||sline.back()=='\r'||sline.back()==' '))
            sline.pop_back();

        if (sline == "end_header") break;

        std::istringstream iss(sline);
        std::string token;
        iss >> token;

        if (token == "format") {
            std::string fmt;
            iss >> fmt;
            if (fmt == "ascii") format = ASCII;
            else if (fmt == "binary_little_endian") format = BIN_LE;
            else if (fmt == "binary_big_endian") format = BIN_BE;
        } else if (token == "element") {
            std::string ename;
            int ecount;
            iss >> ename >> ecount;
            if (ename == "vertex") { nverts = ecount; inVertex = true; inFace = false; }
            else if (ename == "face") { nfaces = ecount; inVertex = false; inFace = true; }
            else { inVertex = false; inFace = false; }
        } else if (token == "property") {
            PropInfo pi;
            pi.isList = false;
            std::string t;
            iss >> t;
            if (t == "list") {
                pi.isList = true;
                iss >> pi.listCountType >> pi.listElemType >> pi.name;
            } else {
                pi.type = t;
                iss >> pi.name;
            }
            if (inVertex) {
                if (pi.name == "x") vertXIdx = (int)vertProps.size();
                if (pi.name == "y") vertYIdx = (int)vertProps.size();
                if (pi.name == "z") vertZIdx = (int)vertProps.size();
                vertProps.push_back(pi);
            } else if (inFace) {
                faceProps.push_back(pi);
            }
        }
    }

    if (vertXIdx < 0 || vertYIdx < 0 || vertZIdx < 0) { fclose(fp); return false; }

    // Helper: type sizes
    auto typeSize = [](const std::string& t) -> int {
        if (t=="char"||t=="uchar"||t=="int8"||t=="uint8") return 1;
        if (t=="short"||t=="ushort"||t=="int16"||t=="uint16") return 2;
        if (t=="int"||t=="uint"||t=="int32"||t=="uint32"||t=="float"||t=="float32") return 4;
        if (t=="double"||t=="float64"||t=="int64"||t=="uint64") return 8;
        return 4;
    };

    auto swapBytes = [](char* data, int size) {
        for (int i = 0; i < size/2; i++)
            std::swap(data[i], data[size-1-i]);
    };

    auto readBinValue = [&](FILE* f, const std::string& type, bool bigEndian) -> double {
        int sz = typeSize(type);
        char buf[8] = {};
        if (fread(buf, 1, sz, f) != (size_t)sz) return 0;
        if (bigEndian) swapBytes(buf, sz);

        if (type=="float"||type=="float32") { float v; memcpy(&v,buf,4); return v; }
        if (type=="double"||type=="float64") { double v; memcpy(&v,buf,8); return v; }
        if (type=="int"||type=="int32") { int32_t v; memcpy(&v,buf,4); return v; }
        if (type=="uint"||type=="uint32") { uint32_t v; memcpy(&v,buf,4); return v; }
        if (type=="short"||type=="int16") { int16_t v; memcpy(&v,buf,2); return v; }
        if (type=="ushort"||type=="uint16") { uint16_t v; memcpy(&v,buf,2); return v; }
        if (type=="char"||type=="int8") { int8_t v; memcpy(&v,buf,1); return v; }
        if (type=="uchar"||type=="uint8") { uint8_t v; memcpy(&v,buf,1); return v; }
        return 0;
    };

    m.verts.resize(nverts);
    m.faces.resize(nfaces);

    bool bigEndian = (format == BIN_BE);

    if (format == ASCII) {
        // Read vertices
        for (int vi = 0; vi < nverts; vi++) {
            if (!fgets(line, sizeof(line), fp)) { fclose(fp); return false; }
            std::istringstream vss(line);
            std::vector<double> vals;
            double val;
            while (vss >> val) vals.push_back(val);
            if (vertXIdx < (int)vals.size()) m.verts[vi].x = vals[vertXIdx];
            if (vertYIdx < (int)vals.size()) m.verts[vi].y = vals[vertYIdx];
            if (vertZIdx < (int)vals.size()) m.verts[vi].z = vals[vertZIdx];
        }
        // Read faces
        for (int fi = 0; fi < nfaces; fi++) {
            if (!fgets(line, sizeof(line), fp)) { fclose(fp); return false; }
            std::istringstream fss(line);
            int cnt;
            fss >> cnt;
            if (cnt >= 3) {
                fss >> m.faces[fi].v[0] >> m.faces[fi].v[1] >> m.faces[fi].v[2];
            }
            // skip remaining indices for polygons > 3
        }
    } else {
        // Binary
        // Read vertices
        for (int vi = 0; vi < nverts; vi++) {
            for (int pi = 0; pi < (int)vertProps.size(); pi++) {
                if (vertProps[pi].isList) {
                    int cnt = (int)readBinValue(fp, vertProps[pi].listCountType, bigEndian);
                    for (int j = 0; j < cnt; j++)
                        readBinValue(fp, vertProps[pi].listElemType, bigEndian);
                } else {
                    double val = readBinValue(fp, vertProps[pi].type, bigEndian);
                    if (pi == vertXIdx) m.verts[vi].x = val;
                    if (pi == vertYIdx) m.verts[vi].y = val;
                    if (pi == vertZIdx) m.verts[vi].z = val;
                }
            }
        }
        // Read faces
        for (int fi = 0; fi < nfaces; fi++) {
            for (int pi = 0; pi < (int)faceProps.size(); pi++) {
                if (faceProps[pi].isList) {
                    int cnt = (int)readBinValue(fp, faceProps[pi].listCountType, bigEndian);
                    std::vector<int> indices(cnt);
                    for (int j = 0; j < cnt; j++)
                        indices[j] = (int)readBinValue(fp, faceProps[pi].listElemType, bigEndian);
                    if (pi == 0 && cnt >= 3) {
                        m.faces[fi].v[0] = indices[0];
                        m.faces[fi].v[1] = indices[1];
                        m.faces[fi].v[2] = indices[2];
                    }
                } else {
                    readBinValue(fp, faceProps[pi].type, bigEndian);
                }
            }
        }
    }

    fclose(fp);
    return true;
}

// ---- OFF Parser ----
static bool loadOFF(const std::string& filename, Mesh& m) {
    FILE* fp = fopen(filename.c_str(), "r");
    if (!fp) return false;

    char line[4096];
    // Read header
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); return false; }
    // May be "OFF" or "OFF\n" or have counts on same line
    std::string header(line);
    while (!header.empty() && (header.back()=='\n'||header.back()=='\r'))
        header.pop_back();

    int nverts = 0, nfaces = 0, nedges = 0;

    if (header == "OFF" || header == "COFF" || header == "NOFF") {
        if (fscanf(fp, "%d %d %d", &nverts, &nfaces, &nedges) < 2) { fclose(fp); return false; }
        fgets(line, sizeof(line), fp); // consume rest of line
    } else {
        // Maybe "OFF nverts nfaces nedges" on one line
        if (sscanf(line, "OFF %d %d %d", &nverts, &nfaces, &nedges) < 2 &&
            sscanf(line, "%d %d %d", &nverts, &nfaces, &nedges) < 2) {
            fclose(fp); return false;
        }
    }

    m.verts.resize(nverts);
    for (int i = 0; i < nverts; i++) {
        if (fscanf(fp, "%lf %lf %lf", &m.verts[i].x, &m.verts[i].y, &m.verts[i].z) < 3) {
            fclose(fp); return false;
        }
        fgets(line, sizeof(line), fp); // skip rest of line (optional normals/colors)
    }

    m.faces.resize(nfaces);
    for (int i = 0; i < nfaces; i++) {
        int cnt;
        if (fscanf(fp, "%d", &cnt) < 1) { fclose(fp); return false; }
        if (cnt >= 3) {
            if (fscanf(fp, "%d %d %d", &m.faces[i].v[0], &m.faces[i].v[1], &m.faces[i].v[2]) < 3) {
                fclose(fp); return false;
            }
        }
        fgets(line, sizeof(line), fp); // skip rest
    }

    fclose(fp);
    return true;
}

// ---- OBJ Parser ----
static bool loadOBJ(const std::string& filename, Mesh& m) {
    FILE* fp = fopen(filename.c_str(), "r");
    if (!fp) return false;

    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == 'v' && line[1] == ' ') {
            Vec3d v;
            if (sscanf(line + 2, "%lf %lf %lf", &v.x, &v.y, &v.z) >= 3)
                m.verts.push_back(v);
        } else if (line[0] == 'f' && line[1] == ' ') {
            // Parse face: handles "f v", "f v/vt", "f v/vt/vn", "f v//vn"
            std::istringstream iss(line + 2);
            std::string token;
            std::vector<int> indices;
            while (iss >> token) {
                int idx = atoi(token.c_str());
                if (idx < 0) idx = (int)m.verts.size() + idx + 1; // relative indexing
                indices.push_back(idx - 1); // OBJ is 1-indexed
            }
            if (indices.size() >= 3) {
                Triangle t;
                t.v[0] = indices[0]; t.v[1] = indices[1]; t.v[2] = indices[2];
                m.faces.push_back(t);
                // Triangulate quads+
                for (size_t i = 3; i < indices.size(); i++) {
                    Triangle t2;
                    t2.v[0] = indices[0]; t2.v[1] = indices[i-1]; t2.v[2] = indices[i];
                    m.faces.push_back(t2);
                }
            }
        }
    }

    fclose(fp);
    return !m.verts.empty();
}

// ---- Load mesh by extension ----
static bool LoadMesh(const std::string& filename, Mesh& m) {
    std::string ext = getExtension(filename);
    bool ok = false;
    if (ext == "ply") ok = loadPLY(filename, m);
    else if (ext == "off") ok = loadOFF(filename, m);
    else if (ext == "obj") ok = loadOBJ(filename, m);
    else {
        printf("unknown file format '%s'.\n", ext.c_str());
        return false;
    }

    if (ok) m.computeBBox();
    return ok;
}

// ============================================================================
// Section 12: Mesh Cleaning
// ============================================================================
static void CleanMesh(Mesh& m, const char* name) {
    // Remove duplicate vertices (exact match via sorting)
    int origVn = m.vn();

    // Build index map
    std::vector<int> order(m.vn());
    for (int i = 0; i < m.vn(); i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (m.verts[a].x != m.verts[b].x) return m.verts[a].x < m.verts[b].x;
        if (m.verts[a].y != m.verts[b].y) return m.verts[a].y < m.verts[b].y;
        return m.verts[a].z < m.verts[b].z;
    });

    std::vector<int> remap(m.vn());
    int newIdx = 0;
    std::vector<Vec3d> newVerts;
    for (int i = 0; i < (int)order.size(); i++) {
        if (i == 0 || m.verts[order[i]].x != m.verts[order[i-1]].x ||
            m.verts[order[i]].y != m.verts[order[i-1]].y ||
            m.verts[order[i]].z != m.verts[order[i-1]].z) {
            newVerts.push_back(m.verts[order[i]]);
            newIdx = (int)newVerts.size() - 1;
        }
        remap[order[i]] = newIdx;
    }

    int dupCount = origVn - (int)newVerts.size();

    // Remap face indices
    for (auto& f : m.faces)
        for (int j = 0; j < 3; j++)
            f.v[j] = remap[f.v[j]];

    m.verts = newVerts;

    // Remove unreferenced vertices
    std::vector<bool> used(m.vn(), false);
    for (auto& f : m.faces)
        for (int j = 0; j < 3; j++)
            used[f.v[j]] = true;

    std::vector<int> remap2(m.vn(), -1);
    std::vector<Vec3d> finalVerts;
    for (int i = 0; i < m.vn(); i++) {
        if (used[i]) {
            remap2[i] = (int)finalVerts.size();
            finalVerts.push_back(m.verts[i]);
        }
    }
    int unrefCount = m.vn() - (int)finalVerts.size();

    for (auto& f : m.faces)
        for (int j = 0; j < 3; j++)
            f.v[j] = remap2[f.v[j]];

    m.verts = finalVerts;
    m.computeBBox();

    printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n", dupCount, unrefCount, name);
}

// ============================================================================
// Section 13: PLY Error Exporter
// ============================================================================
// Color ramp: quality -> RGBA (red to yellow to green to cyan to blue)
static void QualityToColor(float q, float minQ, float maxQ, unsigned char rgba[4]) {
    if (maxQ <= minQ) { rgba[0]=rgba[1]=rgba[2]=128; rgba[3]=255; return; }
    float t = (q - minQ) / (maxQ - minQ);
    if (t < 0) t = 0; if (t > 1) t = 1;

    float r, g, b;
    if (t < 0.25f) {
        r = 0; g = t * 4.0f; b = 1;
    } else if (t < 0.5f) {
        r = 0; g = 1; b = 1.0f - (t - 0.25f) * 4.0f;
    } else if (t < 0.75f) {
        r = (t - 0.5f) * 4.0f; g = 1; b = 0;
    } else {
        r = 1; g = 1.0f - (t - 0.75f) * 4.0f; b = 0;
    }
    rgba[0] = (unsigned char)(r * 255);
    rgba[1] = (unsigned char)(g * 255);
    rgba[2] = (unsigned char)(b * 255);
    rgba[3] = 255;
}

static void SaveErrorPLY(const Mesh& m, const std::string& filename,
                         float colorMin, float colorMax) {
    FILE* fp = fopen(filename.c_str(), "w");
    if (!fp) return;

    // Find quality range if not specified
    float qmin = colorMin, qmax = colorMax;
    if (qmin == 0 && qmax == 0 && !m.vertQuality.empty()) {
        qmin = *std::min_element(m.vertQuality.begin(), m.vertQuality.end());
        qmax = *std::max_element(m.vertQuality.begin(), m.vertQuality.end());
    }

    fprintf(fp, "ply\n");
    fprintf(fp, "format ascii 1.0\n");
    fprintf(fp, "element vertex %d\n", m.vn());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "property float quality\n");
    fprintf(fp, "property uchar red\n");
    fprintf(fp, "property uchar green\n");
    fprintf(fp, "property uchar blue\n");
    fprintf(fp, "property uchar alpha\n");
    fprintf(fp, "element face %d\n", m.fn());
    fprintf(fp, "property list uchar int vertex_indices\n");
    fprintf(fp, "end_header\n");

    for (int i = 0; i < m.vn(); i++) {
        float q = (i < (int)m.vertQuality.size()) ? m.vertQuality[i] : 0;
        unsigned char rgba[4];
        QualityToColor(q, qmin, qmax, rgba);
        fprintf(fp, "%f %f %f %f %d %d %d %d\n",
            (float)m.verts[i].x, (float)m.verts[i].y, (float)m.verts[i].z,
            q, rgba[0], rgba[1], rgba[2], rgba[3]);
    }

    for (int i = 0; i < m.fn(); i++) {
        fprintf(fp, "3 %d %d %d\n", m.faces[i].v[0], m.faces[i].v[1], m.faces[i].v[2]);
    }

    fclose(fp);
}

// ============================================================================
// Section 15: SamplingFlags
// ============================================================================
struct SamplingFlags {
    enum {
        HIST                        = 0x0001,
        VERTEX_SAMPLING             = 0x0002,
        EDGE_SAMPLING               = 0x0004,
        FACE_SAMPLING               = 0x0008,
        MONTECARLO_SAMPLING         = 0x0010,
        SUBDIVISION_SAMPLING        = 0x0020,
        SIMILAR_SAMPLING            = 0x0040,
        NO_SAMPLING                 = 0x0070,
        SAVE_ERROR                  = 0x0100,
        INCLUDE_UNREFERENCED_VERTICES = 0x0200,
        USE_STATIC_GRID             = 0x0400,
    };
};

// ============================================================================
// Section 14: Sampling Engine
// ============================================================================
class Sampling {
    Mesh& S1;       // source mesh (sampled)
    Mesh& S2;       // target mesh (grid built on this)
    UniformGrid grid;

    double dist_upper_bound;
    double n_samples_per_area_unit;
    unsigned long n_samples_target;
    int Flags;

    Histogram hist;
    unsigned long n_total_samples;
    unsigned long n_total_area_samples;
    unsigned long n_total_edge_samples;
    unsigned long n_total_vertex_samples;
    double max_dist;
    double mean_dist;
    double RMS_dist;
    double volume;
    double area_S1;

    int n_samples;  // temp counter
    int n_hist_bins;
    int print_every_n_elements;

    // Track which vertices are referenced by faces
    std::vector<bool> vertReferred;

    float AddSample(const Vec3d& p) {
        double dist = dist_upper_bound;
        Vec3d closestPt;
        int fi = grid.GetClosest(p, dist_upper_bound, dist, closestPt);

        if (fi < 0 || dist >= dist_upper_bound)
            return -1.0f;

        if (dist > max_dist) max_dist = dist;
        mean_dist += dist;
        RMS_dist += dist * dist;
        n_total_samples++;

        if (Flags & SamplingFlags::HIST)
            hist.Add(fabs(dist));

        return (float)dist;
    }

    // --- Vertex sampling ---
    void VertexSampling() {
        int cnt = 0;
        printf("Vertex sampling\n");
        for (int vi = 0; vi < S1.vn(); vi++) {
            bool referred = (vi < (int)vertReferred.size()) ? vertReferred[vi] : false;
            if (referred || (Flags & SamplingFlags::INCLUDE_UNREFERENCED_VERTICES)) {
                float error = AddSample(S1.verts[vi]);
                n_total_vertex_samples++;

                if (Flags & SamplingFlags::SAVE_ERROR) {
                    if (vi >= (int)S1.vertQuality.size())
                        S1.vertQuality.resize(S1.vn(), 0);
                    S1.vertQuality[vi] = error;
                }

                if (print_every_n_elements > 0 && !(++cnt % print_every_n_elements))
                    printf("Sampling vertices %d%%\r", (100 * cnt / S1.vn()));
            }
        }
        printf("                       \r");
    }

    // --- Edge sampling ---
    void SampleEdge(const Vec3d& v0, const Vec3d& v1, int n_per_edge) {
        Vec3d e = (v1 - v0) / (double)(n_per_edge + 1);
        for (int i = 1; i <= n_per_edge; i++) {
            AddSample(v0 + e * i);
            n_total_edge_samples++;
        }
    }

    void EdgeSampling() {
        printf("Edge sampling\n");
        // Build unique edge list
        typedef std::pair<int,int> EdgePair;
        std::vector<EdgePair> edges;
        for (int fi = 0; fi < S1.fn(); fi++) {
            for (int i = 0; i < 3; i++) {
                int a = S1.faces[fi].v[i];
                int b = S1.faces[fi].v[(i+1)%3];
                if (a > b) std::swap(a, b);
                edges.push_back({a, b});
            }
        }
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

        double n_samples_per_length_unit;
        if (Flags & SamplingFlags::FACE_SAMPLING)
            n_samples_per_length_unit = sqrt(n_samples_per_area_unit);
        else
            n_samples_per_length_unit = n_samples_per_area_unit;

        double n_samples_decimal = 0.0;
        int cnt = 0;
        for (auto& e : edges) {
            double len = Distance(S1.verts[e.first], S1.verts[e.second]);
            n_samples_decimal += len * n_samples_per_length_unit;
            int ns = (int)n_samples_decimal;
            SampleEdge(S1.verts[e.first], S1.verts[e.second], ns);
            n_samples_decimal -= (double)ns;

            if (print_every_n_elements > 0 && !(++cnt % print_every_n_elements))
                printf("Sampling edge %d%%\r", (int)(100 * cnt / edges.size()));
        }
        printf("                     \r");
    }

    // --- Face sampling: Montecarlo ---
    void AddRandomSample(int fi) {
        const Vec3d& p0 = S1.verts[S1.faces[fi].v[0]];
        const Vec3d& p1 = S1.verts[S1.faces[fi].v[1]];
        const Vec3d& p2 = S1.verts[S1.faces[fi].v[2]];
        Vec3d v1 = p1 - p0;
        Vec3d v2 = p2 - p0;

        double r1 = (double)rand() / (double)RAND_MAX;
        double r2 = (double)rand() / (double)RAND_MAX;
        if (r1 + r2 > 1.0) { r1 = 1.0 - r1; r2 = 1.0 - r2; }

        AddSample(p0 + v1 * r1 + v2 * r2);
        n_total_area_samples++;
    }

    void MontecarloFaceSampling() {
        double n_samples_decimal = 0.0;
        srand((unsigned int)clock());
        for (int fi = 0; fi < S1.fn(); fi++) {
            n_samples_decimal += 0.5 * S1.faceDoubleArea(fi) * n_samples_per_area_unit;
            int ns = (int)n_samples_decimal;
            for (int i = 0; i < ns; i++)
                AddRandomSample(fi);
            n_samples_decimal -= (double)ns;
        }
    }

    // --- Face sampling: Subdivision ---
    void FaceSubdiv(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, int maxdepth) {
        if (maxdepth == 0) {
            AddSample((v0 + v1 + v2) / 3.0);
            n_total_area_samples++;
            n_samples++;
            return;
        }

        double d01 = SquaredDistance(v0, v1);
        double d12 = SquaredDistance(v1, v2);
        double d20 = SquaredDistance(v2, v0);
        int res;
        if (d01 > d12)
            res = (d01 > d20) ? 0 : 2;
        else
            res = (d12 > d20) ? 1 : 2;

        Vec3d pp;
        switch (res) {
        case 0: pp = (v0+v1)/2.0; FaceSubdiv(v0,pp,v2,maxdepth-1); FaceSubdiv(pp,v1,v2,maxdepth-1); break;
        case 1: pp = (v1+v2)/2.0; FaceSubdiv(v0,v1,pp,maxdepth-1); FaceSubdiv(v0,pp,v2,maxdepth-1); break;
        case 2: pp = (v2+v0)/2.0; FaceSubdiv(v0,v1,pp,maxdepth-1); FaceSubdiv(pp,v1,v2,maxdepth-1); break;
        }
    }

    void SubdivFaceSampling() {
        int cnt = 0;
        double n_samples_decimal = 0.0;
        printf("Subdivision face sampling\n");
        for (int fi = 0; fi < S1.fn(); fi++) {
            n_samples_decimal += 0.5 * S1.faceDoubleArea(fi) * n_samples_per_area_unit;
            int ns = (int)n_samples_decimal;
            if (ns) {
                int maxdepth = (int)(log((double)ns) / log(2.0));
                n_samples = 0;
                FaceSubdiv(S1.verts[S1.faces[fi].v[0]], S1.verts[S1.faces[fi].v[1]],
                           S1.verts[S1.faces[fi].v[2]], maxdepth);
            }
            n_samples_decimal -= (double)n_samples;

            if (print_every_n_elements > 0 && !(++cnt % print_every_n_elements))
                printf("Sampling face %d%%\r", (100 * cnt / S1.fn()));
        }
        printf("                     \r");
    }

    // --- Face sampling: Similar Triangles ---
    void SimilarTriangles(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, int n_per_edge) {
        Vec3d V1 = (v1 - v0) / (double)(n_per_edge - 1);
        Vec3d V2 = (v2 - v0) / (double)(n_per_edge - 1);

        for (int i = 1; i < n_per_edge - 1; i++)
            for (int j = 1; j < n_per_edge - 1 - i; j++) {
                AddSample(v0 + V1 * (double)i + V2 * (double)j);
                n_total_area_samples++;
                n_samples++;
            }
    }

    void SimilarFaceSampling() {
        int cnt = 0;
        double n_samples_decimal = 0.0;
        printf("Similar Triangles face sampling\n");
        for (int fi = 0; fi < S1.fn(); fi++) {
            n_samples_decimal += 0.5 * S1.faceDoubleArea(fi) * n_samples_per_area_unit;
            int ns = (int)n_samples_decimal;
            if (ns) {
                int n_per_edge = (int)((sqrt(1.0 + 8.0 * (double)ns) + 5.0) / 2.0);
                n_samples = 0;
                SimilarTriangles(S1.verts[S1.faces[fi].v[0]], S1.verts[S1.faces[fi].v[1]],
                                 S1.verts[S1.faces[fi].v[2]], n_per_edge);
            }
            n_samples_decimal -= (double)n_samples;

            if (print_every_n_elements > 0 && !(++cnt % print_every_n_elements))
                printf("Sampling face %d%%\r", (100 * cnt / S1.fn()));
        }
        printf("                     \r");
    }

public:
    Sampling(Mesh& s1, Mesh& s2) : S1(s1), S2(s2), dist_upper_bound(0),
        n_samples_per_area_unit(0), n_samples_target(0), Flags(0),
        n_total_samples(0), n_total_area_samples(0), n_total_edge_samples(0),
        n_total_vertex_samples(0), max_dist(0), mean_dist(0), RMS_dist(0),
        volume(0), n_samples(0), n_hist_bins(256), print_every_n_elements(1)
    {
        area_S1 = s1.computeArea();
        print_every_n_elements = S1.fn() / 100;
        if (print_every_n_elements <= 1) print_every_n_elements = 2;

        // Mark referred vertices
        vertReferred.assign(S1.vn(), false);
        for (int fi = 0; fi < S1.fn(); fi++)
            for (int j = 0; j < 3; j++)
                vertReferred[S1.faces[fi].v[j]] = true;
    }

    void Hausdorff() {
        // Build grid on S2
        grid.Build(S2);

        dist_upper_bound = S2.bbox.Diag();
        if (Flags & SamplingFlags::HIST)
            hist.SetRange(0.0, dist_upper_bound / 100.0, n_hist_bins);

        n_total_area_samples = n_total_edge_samples = n_total_vertex_samples = n_total_samples = 0;
        n_samples = 0;
        max_dist = -1e30;
        mean_dist = RMS_dist = 0;

        // Vertex sampling
        if (Flags & SamplingFlags::VERTEX_SAMPLING)
            VertexSampling();

        // Budget management
        if (n_samples_target > n_total_samples) {
            n_samples_target -= (unsigned long)n_total_samples;
            n_samples_per_area_unit = n_samples_target / area_S1;

            // Edge sampling
            if (Flags & SamplingFlags::EDGE_SAMPLING) {
                EdgeSampling();
                if (n_samples_target > n_total_samples)
                    n_samples_target -= (unsigned long)n_total_samples;
                else
                    n_samples_target = 0;
            }

            // Face sampling
            if ((Flags & SamplingFlags::FACE_SAMPLING) && n_samples_target > 0) {
                n_samples_per_area_unit = n_samples_target / area_S1;
                if (Flags & SamplingFlags::MONTECARLO_SAMPLING)  MontecarloFaceSampling();
                if (Flags & SamplingFlags::SUBDIVISION_SAMPLING) SubdivFaceSampling();
                if (Flags & SamplingFlags::SIMILAR_SAMPLING)     SimilarFaceSampling();
            }
        }

        // Compute per-vertex color if saving error
        if (Flags & SamplingFlags::SAVE_ERROR) {
            // vertQuality is already set during vertex sampling
            // For non-vertex-sampled vertices, quality stays at 0
        }

        // Compute final statistics
        n_samples_per_area_unit = (double)n_total_samples / area_S1;
        volume = mean_dist / n_samples_per_area_unit / 2.0;
        if (n_total_samples > 0) {
            mean_dist /= n_total_samples;
            RMS_dist = sqrt(RMS_dist / n_total_samples);
        }
    }

    void SetFlags(int flags) { Flags = flags; }
    void SetSamplesTarget(unsigned long n) {
        n_samples_target = n;
        n_samples_per_area_unit = n_samples_target / area_S1;
    }
    void SetSamplesPerAreaUnit(double n) {
        n_samples_per_area_unit = n;
        n_samples_target = (unsigned long)(n_samples_per_area_unit * area_S1);
    }

    double GetArea() const { return area_S1; }
    double GetDistMax() const { return max_dist; }
    double GetDistMean() const { return mean_dist; }
    double GetDistRMS() const { return RMS_dist; }
    unsigned long GetNSamples() const { return n_total_samples; }
    unsigned long GetNAreaSamples() const { return n_total_area_samples; }
    unsigned long GetNEdgeSamples() const { return n_total_edge_samples; }
    unsigned long GetNVertexSamples() const { return n_total_vertex_samples; }
    double GetNSamplesPerAreaUnit() const { return n_samples_per_area_unit; }
    unsigned long GetNSamplesTarget() const { return n_samples_target; }
    Histogram& GetHist() { return hist; }
};

// ============================================================================
// Section 16: main() — CLI + Orchestration
// ============================================================================
static void Usage() {
    printf("\nUsage:  "
        "metro_standalone file1 file2 [opt]\n"
        "Where opt can be:\n"
        "  -v         disable vertex sampling\n"
        "  -e         disable edge sampling\n"
        "  -f         disable face sampling\n"
        "  -u         ignore unreferred vertices\n"
        "  -sx        set the face sampling mode\n"
        "             where x can be:\n"
        "              -s0  montecarlo sampling\n"
        "              -s1  subdivision sampling\n"
        "              -s2  similar triangles sampling (Default)\n"
        "  -n#        set the required number of samples (overrides -a)\n"
        "  -a#        set the required number of samples per area unit (overrides -n)\n"
        "  -c         save a mesh with error as per-vertex colour and quality\n"
        "  -C # #     Set the min/max values used for color mapping\n"
        "  -L         Remove duplicated and unreferenced vertices before processing\n"
        "  -h         write files with histograms of error distribution\n"
        "\n"
        "Default options are to sample vertexes, edge and faces by taking \n"
        "a number of samples that is approx. 10x the face number.\n"
    );
    exit(-1);
}

static std::string SaveFileName(const std::string& filename) {
    size_t pos = filename.rfind('.');
    if (pos == std::string::npos) return filename + "_metro.ply";
    return filename.substr(0, pos) + "_metro.ply";
}

int main(int argc, char** argv) {
    Mesh S1, S2;
    float ColorMin = 0, ColorMax = 0;
    double dist1_max, dist2_max;
    unsigned long n_samples_target = 0;
    double n_samples_per_area_unit = 0;
    int flags;
    bool NumberOfSamples = false;
    bool SamplesPerAreaUnit = false;
    bool CleaningFlag = false;

    printf("-------------------------------\n"
           "   Metro Standalone V.1.00 \n"
           "   (based on VCGlib Metro)\n"
           "   release date: " __DATE__ "\n"
           "-------------------------------\n\n");

    if (argc <= 2) Usage();

    // Default flags
    flags = SamplingFlags::VERTEX_SAMPLING |
            SamplingFlags::EDGE_SAMPLING |
            SamplingFlags::FACE_SAMPLING |
            SamplingFlags::SIMILAR_SAMPLING;

    // Parse command line
    for (int i = 3; i < argc;) {
        if (argv[i][0] == '-')
            switch (argv[i][1]) {
            case 'h': flags |= SamplingFlags::HIST; break;
            case 'v': flags &= ~SamplingFlags::VERTEX_SAMPLING; break;
            case 'e': flags &= ~SamplingFlags::EDGE_SAMPLING; break;
            case 'f': flags &= ~SamplingFlags::FACE_SAMPLING; break;
            case 'u': flags |= SamplingFlags::INCLUDE_UNREFERENCED_VERTICES; break;
            case 's':
                switch (argv[i][2]) {
                case '0': flags = (flags | SamplingFlags::MONTECARLO_SAMPLING) & (~SamplingFlags::NO_SAMPLING); break;
                case '1': flags = (flags | SamplingFlags::SUBDIVISION_SAMPLING) & (~SamplingFlags::NO_SAMPLING); break;
                case '2': flags = (flags | SamplingFlags::SIMILAR_SAMPLING) & (~SamplingFlags::NO_SAMPLING); break;
                default: printf("unable to parse option '%s'\n", argv[i]); exit(0);
                }
                break;
            case 'n': NumberOfSamples = true; n_samples_target = (unsigned long)atoi(&(argv[i][2])); break;
            case 'a': SamplesPerAreaUnit = true; n_samples_per_area_unit = atof(&(argv[i][2])); break;
            case 'c': flags |= SamplingFlags::SAVE_ERROR; break;
            case 'L': CleaningFlag = true; break;
            case 'C': ColorMin = (float)atof(argv[i+1]); ColorMax = (float)atof(argv[i+2]); i += 2; break;
            default: printf("unable to parse option '%s'\n", argv[i]); exit(0);
            }
        i++;
    }

    flags |= SamplingFlags::USE_STATIC_GRID;

    // Load meshes
    printf("Loading mesh '%s'...\n", argv[1]);
    if (!LoadMesh(argv[1], S1)) {
        printf("error loading '%s'\n", argv[1]);
        return -1;
    }
    printf("read mesh `%s'\n", argv[1]);

    printf("Loading mesh '%s'...\n", argv[2]);
    if (!LoadMesh(argv[2], S2)) {
        printf("error loading '%s'\n", argv[2]);
        return -1;
    }
    printf("read mesh `%s'\n", argv[2]);

    // Cleaning
    if (CleaningFlag) {
        CleanMesh(S1, argv[1]);
        CleanMesh(S2, argv[2]);
    }

    std::string S1NewName = SaveFileName(argv[1]);
    std::string S2NewName = SaveFileName(argv[2]);

    if (!NumberOfSamples && !SamplesPerAreaUnit) {
        NumberOfSamples = true;
        n_samples_target = 10 * (unsigned long)std::max(S1.fn(), S2.fn());
    }

    // Save original bboxes for display
    Box3d tmp_bbox_M1 = S1.bbox;
    Box3d tmp_bbox_M2 = S2.bbox;

    // Compute combined bbox
    Box3d bbox;
    bbox.Add(S1.bbox);
    bbox.Add(S2.bbox);
    bbox.Offset(bbox.Diag() * 0.02);
    S1.bbox = bbox;
    S2.bbox = bbox;

    // Initialize time
    int t0 = clock();

    Sampling ForwardSampling(S1, S2);
    Sampling BackwardSampling(S2, S1);

    // Print mesh info
    printf("Mesh info:\n");
    printf(" M1: '%s'\n\tvertices  %7i\n\tfaces     %7i\n\tarea      %12.4f\n",
           argv[1], S1.vn(), S1.fn(), ForwardSampling.GetArea());
    printf("\tbbox (%7.4f %7.4f %7.4f)-(%7.4f %7.4f %7.4f)\n",
           tmp_bbox_M1.min[0], tmp_bbox_M1.min[1], tmp_bbox_M1.min[2],
           tmp_bbox_M1.max[0], tmp_bbox_M1.max[1], tmp_bbox_M1.max[2]);
    printf("\tbbox diagonal %f\n", (float)tmp_bbox_M1.Diag());
    printf(" M2: '%s'\n\tvertices  %7i\n\tfaces     %7i\n\tarea      %12.4f\n",
           argv[2], S2.vn(), S2.fn(), BackwardSampling.GetArea());
    printf("\tbbox (%7.4f %7.4f %7.4f)-(%7.4f %7.4f %7.4f)\n",
           tmp_bbox_M2.min[0], tmp_bbox_M2.min[1], tmp_bbox_M2.min[2],
           tmp_bbox_M2.max[0], tmp_bbox_M2.max[1], tmp_bbox_M2.max[2]);
    printf("\tbbox diagonal %f\n", (float)tmp_bbox_M2.Diag());

    // Forward distance
    printf("\nForward distance (M1 -> M2):\n");
    ForwardSampling.SetFlags(flags);
    if (NumberOfSamples) {
        ForwardSampling.SetSamplesTarget(n_samples_target);
        n_samples_per_area_unit = ForwardSampling.GetNSamplesPerAreaUnit();
    } else {
        ForwardSampling.SetSamplesPerAreaUnit(n_samples_per_area_unit);
        n_samples_target = ForwardSampling.GetNSamplesTarget();
    }
    printf("target # samples      : %lu\ntarget # samples/area : %f\n", n_samples_target, n_samples_per_area_unit);
    ForwardSampling.Hausdorff();
    dist1_max = ForwardSampling.GetDistMax();
    printf("\ndistances:\n  max  : %f (%f  wrt bounding box diagonal)\n",
           (float)dist1_max, (float)(dist1_max / bbox.Diag()));
    printf("  mean : %f\n", ForwardSampling.GetDistMean());
    printf("  RMS  : %f\n", ForwardSampling.GetDistRMS());
    printf("# vertex samples %9lu\n", ForwardSampling.GetNVertexSamples());
    printf("# edge samples   %9lu\n", ForwardSampling.GetNEdgeSamples());
    printf("# area samples   %9lu\n", ForwardSampling.GetNAreaSamples());
    printf("# total samples  %9lu\n", ForwardSampling.GetNSamples());
    printf("# samples per area unit: %f\n\n", ForwardSampling.GetNSamplesPerAreaUnit());

    // Backward distance
    printf("\nBackward distance (M2 -> M1):\n");
    BackwardSampling.SetFlags(flags);
    if (NumberOfSamples) {
        BackwardSampling.SetSamplesTarget(n_samples_target);
        n_samples_per_area_unit = BackwardSampling.GetNSamplesPerAreaUnit();
    } else {
        BackwardSampling.SetSamplesPerAreaUnit(n_samples_per_area_unit);
        n_samples_target = BackwardSampling.GetNSamplesTarget();
    }
    printf("target # samples      : %lu\ntarget # samples/area : %f\n", n_samples_target, n_samples_per_area_unit);
    BackwardSampling.Hausdorff();
    dist2_max = BackwardSampling.GetDistMax();
    printf("\ndistances:\n  max  : %f (%f  wrt bounding box diagonal)\n",
           (float)dist2_max, (float)(dist2_max / bbox.Diag()));
    printf("  mean : %f\n", BackwardSampling.GetDistMean());
    printf("  RMS  : %f\n", BackwardSampling.GetDistRMS());
    printf("# vertex samples %9lu\n", BackwardSampling.GetNVertexSamples());
    printf("# edge samples   %9lu\n", BackwardSampling.GetNEdgeSamples());
    printf("# area samples   %9lu\n", BackwardSampling.GetNAreaSamples());
    printf("# total samples  %9lu\n", BackwardSampling.GetNSamples());
    printf("# samples per area unit: %f\n\n", BackwardSampling.GetNSamplesPerAreaUnit());

    // Compute time info
    unsigned long elapsed_time = clock() - t0;
    unsigned long n_total_sample = ForwardSampling.GetNSamples() + BackwardSampling.GetNSamples();
    double mesh_dist_max = std::max(dist1_max, dist2_max);

    printf("\nHausdorff distance: %f (%f  wrt bounding box diagonal)\n",
           (float)mesh_dist_max, (float)(mesh_dist_max / bbox.Diag()));
    printf("  Computation time  : %d ms\n", (int)(1000.0 * elapsed_time / CLOCKS_PER_SEC));
    if (elapsed_time > 0)
        printf("  # samples/second  : %f\n\n", (float)n_total_sample / ((float)elapsed_time / CLOCKS_PER_SEC));

    // Save error files
    if (flags & SamplingFlags::SAVE_ERROR) {
        SaveErrorPLY(S1, S1NewName, ColorMin, ColorMax);
        SaveErrorPLY(S2, S2NewName, ColorMin, ColorMax);
        printf("Saved error meshes: '%s' and '%s'\n", S1NewName.c_str(), S2NewName.c_str());
    }

    // Save histograms
    if (flags & SamplingFlags::HIST) {
        ForwardSampling.GetHist().FileWrite("forward_result.csv");
        BackwardSampling.GetHist().FileWrite("backward_result.csv");
        printf("Saved histograms: 'forward_result.csv' and 'backward_result.csv'\n");
    }

    return 0;
}
