// =============================================================================
// ATTENE_IMPLICIT_POINT3D.H - 3D Implicit Point Types for Robust CSG
// =============================================================================
//
// Based on "Indirect Predicates for Geometric Constructions" by Marco Attene
// and "Interactive and Robust Mesh Booleans" by Cherchi et al.
//
// Key insight: Instead of computing explicit coordinates for intersection points
// (which introduces numerical error), we store points implicitly as the result
// of geometric constructions (e.g., "intersection of line L and plane P").
// All predicates then operate on these implicit representations.
//
// Point Types:
//   - Explicit3D: Regular (x, y, z) coordinates
//   - LPI: Line-Plane Intersection (segment intersects triangle plane)
//   - TPI: Three-Plane Intersection (three triangle planes meet)
//
// =============================================================================

#ifndef __ATTENE_IMPLICIT_POINT3D_H_
#define __ATTENE_IMPLICIT_POINT3D_H_

#include "attene_predicates.h"

// =============================================================================
// Point Type Enumeration
// =============================================================================

enum AttenePointType3D {
	ATTENE_EXPLICIT_3D = 0,
	ATTENE_LPI_3D = 1,    // Line-Plane Intersection
	ATTENE_TPI_3D = 2     // Three-Plane Intersection
};

// =============================================================================
// Forward Declarations
// =============================================================================

class AttenePoint3D;
class AtteneExplicit3D;
class AtteneLPI3D;
class AtteneTPI3D;
struct BSNumber;

// =============================================================================
// Base Class for 3D Points
// =============================================================================

class AttenePoint3D {
protected:
	AttenePointType3D type;

public:
	AttenePoint3D(AttenePointType3D t) : type(t) {}
	virtual ~AttenePoint3D() {}

	AttenePointType3D getType() const { return type; }
	bool isExplicit() const { return type == ATTENE_EXPLICIT_3D; }
	bool isLPI() const { return type == ATTENE_LPI_3D; }
	bool isTPI() const { return type == ATTENE_TPI_3D; }

	// Lambda accessors - the point is represented as (lx/d, ly/d, lz/d).
	// Interval version is the fast filter; when its signs are unreliable the
	// predicates fall back to attene_get_bs_lambda (exact, arbitrary
	// precision). Returns false if the interval representation could not
	// prove the construction non-degenerate; the exact path still works.
	virtual bool getIntervalLambda(
		attene_interval& lx, attene_interval& ly, attene_interval& lz,
		attene_interval& d) const = 0;

	// Get approximate explicit coordinates (for visualization, debugging)
	virtual bool getApproxXYZ(double& x, double& y, double& z) const = 0;

	// Get as Vec3d (may lose precision for implicit points)
	Vec3d toVec3d() const {
		double x, y, z;
		getApproxXYZ(x, y, z);
		return Vec3d(x, y, z);
	}
};

// =============================================================================
// Explicit 3D Point
// =============================================================================

class AtteneExplicit3D : public AttenePoint3D {
public:
	double x, y, z;

	AtteneExplicit3D() : AttenePoint3D(ATTENE_EXPLICIT_3D), x(0), y(0), z(0) {}
	AtteneExplicit3D(double _x, double _y, double _z)
		: AttenePoint3D(ATTENE_EXPLICIT_3D), x(_x), y(_y), z(_z) {}
	AtteneExplicit3D(const Vec3d& v)
		: AttenePoint3D(ATTENE_EXPLICIT_3D), x(v.x), y(v.y), z(v.z) {}
	AtteneExplicit3D(const Vec3f& v)
		: AttenePoint3D(ATTENE_EXPLICIT_3D), x(v.x), y(v.y), z(v.z) {}

	double X() const { return x; }
	double Y() const { return y; }
	double Z() const { return z; }

	const double* ptr() const { return &x; }

	bool getIntervalLambda(
		attene_interval& lx, attene_interval& ly, attene_interval& lz,
		attene_interval& d) const override
	{
		lx = attene_interval(x);
		ly = attene_interval(y);
		lz = attene_interval(z);
		d = attene_interval(1.0);
		return true;
	}

	bool getApproxXYZ(double& ox, double& oy, double& oz) const override {
		ox = x;
		oy = y;
		oz = z;
		return true;
	}
};

// =============================================================================
// Line-Plane Intersection (LPI)
// =============================================================================
//
// Represents the intersection of a line segment with a plane.
// The line is defined by two explicit points: P and Q
// The plane is defined by three explicit points: R, S, T
//
// The intersection point I lies on segment PQ such that:
//   I = P + t * (Q - P)
// where t is the parameter satisfying:
//   dot(I - R, normal(RST)) = 0
//
// We store the point implicitly as (lx/d, ly/d, lz/d) where:
//   d = dot(Q - P, normal(RST))
//   lx, ly, lz are computed to give I when divided by d

class AtteneLPI3D : public AttenePoint3D {
public:
	// The segment endpoints
	const AtteneExplicit3D* P;
	const AtteneExplicit3D* Q;

	// The plane-defining triangle
	const AtteneExplicit3D* R;
	const AtteneExplicit3D* S;
	const AtteneExplicit3D* T;

	// Cached interval values
	mutable attene_interval cached_lx, cached_ly, cached_lz, cached_d;
	mutable bool cache_computed;
	mutable bool cache_valid;

	AtteneLPI3D(
		const AtteneExplicit3D* _P, const AtteneExplicit3D* _Q,
		const AtteneExplicit3D* _R, const AtteneExplicit3D* _S, const AtteneExplicit3D* _T)
		: AttenePoint3D(ATTENE_LPI_3D)
		, P(_P), Q(_Q), R(_R), S(_S), T(_T)
		, cache_computed(false), cache_valid(false)
	{}

	bool getIntervalLambda(
		attene_interval& lx, attene_interval& ly, attene_interval& lz,
		attene_interval& d) const override;

	bool getApproxXYZ(double& ox, double& oy, double& oz) const override;
};

// =============================================================================
// Three-Plane Intersection (TPI)
// =============================================================================
//
// Represents the intersection of three planes.
// Each plane is defined by three explicit points.
// Plane 1: P1, Q1, R1
// Plane 2: P2, Q2, R2
// Plane 3: P3, Q3, R3
//
// The intersection point is where all three planes meet (if they do).

class AtteneTPI3D : public AttenePoint3D {
public:
	// Plane 1
	const AtteneExplicit3D* P1;
	const AtteneExplicit3D* Q1;
	const AtteneExplicit3D* R1;

	// Plane 2
	const AtteneExplicit3D* P2;
	const AtteneExplicit3D* Q2;
	const AtteneExplicit3D* R2;

	// Plane 3
	const AtteneExplicit3D* P3;
	const AtteneExplicit3D* Q3;
	const AtteneExplicit3D* R3;

	// Cached values
	mutable attene_interval cached_lx, cached_ly, cached_lz, cached_d;
	mutable bool cache_computed;
	mutable bool cache_valid;

	AtteneTPI3D(
		const AtteneExplicit3D* _P1, const AtteneExplicit3D* _Q1, const AtteneExplicit3D* _R1,
		const AtteneExplicit3D* _P2, const AtteneExplicit3D* _Q2, const AtteneExplicit3D* _R2,
		const AtteneExplicit3D* _P3, const AtteneExplicit3D* _Q3, const AtteneExplicit3D* _R3)
		: AttenePoint3D(ATTENE_TPI_3D)
		, P1(_P1), Q1(_Q1), R1(_R1)
		, P2(_P2), Q2(_Q2), R2(_R2)
		, P3(_P3), Q3(_Q3), R3(_R3)
		, cache_computed(false), cache_valid(false)
	{}

	// Helper: compute plane normal and d-coefficient from three points
	// Plane equation: n.x * x + n.y * y + n.z * z = d
	static void computePlaneInterval(
		const AtteneExplicit3D* A, const AtteneExplicit3D* B, const AtteneExplicit3D* C,
		attene_interval& nx, attene_interval& ny, attene_interval& nz, attene_interval& d);

	bool getIntervalLambda(
		attene_interval& lx, attene_interval& ly, attene_interval& lz,
		attene_interval& d) const override;

	bool getApproxXYZ(double& ox, double& oy, double& oz) const override;
};

// =============================================================================
// Exact Lambda (BSNumber)
// =============================================================================
// Computes the lambda representation (lx/d, ly/d, lz/d) of a point using
// exact arbitrary-precision arithmetic (math/bsnumber.h). This is the exact
// tier behind every indirect predicate: the fixed-size expansion arithmetic
// in attene_predicates.h overflows for products of implicit-point lambdas,
// so the exact fallback must be arbitrary precision.
void attene_get_bs_lambda(
	const AttenePoint3D* p,
	BSNumber& lx, BSNumber& ly, BSNumber& lz, BSNumber& d);

// =============================================================================
// Indirect Predicates for 3D Points
// =============================================================================

// orient3d for any combination of explicit and implicit points
// Returns > 0, < 0, or = 0 based on orientation

int attene_orient3d_indirect(
	const AttenePoint3D* a, const AttenePoint3D* b,
	const AttenePoint3D* c, const AttenePoint3D* d);

// Compare points along X axis: returns -1 if a < b, 0 if equal, 1 if a > b
int attene_lessThanOnX(const AttenePoint3D* a, const AttenePoint3D* b);
int attene_lessThanOnY(const AttenePoint3D* a, const AttenePoint3D* b);
int attene_lessThanOnZ(const AttenePoint3D* a, const AttenePoint3D* b);

// Same comparison with the coordinate selected by index (0=X, 1=Y, 2=Z).
int attene_lessThanOnAxis(const AttenePoint3D* a, const AttenePoint3D* b, int axis);

// Map a ProjectionAxis-style pair (0=XY, 1=YZ, 2=XZ) to the two 3D component
// indices it keeps. Shared by every caller that mixes projected 2D predicates
// with per-axis comparisons.
void attene_proj_axes(int axisPair, int* u, int* v);

// Exact geometric equality: true iff a and b denote the same real point.
// (Same point constructed two different ways compares equal.)
bool attene_samePoint3D(const AttenePoint3D* a, const AttenePoint3D* b);

// Projected 2D orientation of three 3D points (any mix of explicit/implicit),
// evaluated exactly on the implicit representations - NOT on rounded
// getApproxXYZ coordinates. axisPair uses the ProjectionAxis convention:
//   0 = XY (u=x, v=y), 1 = YZ (u=y, v=z), 2 = XZ (u=x, v=z)
// Returns the sign of orient2d(u_a,v_a, u_b,v_b, u_c,v_c).
int attene_orient2d_proj(
	int axisPair,
	const AttenePoint3D* a, const AttenePoint3D* b, const AttenePoint3D* c);

#endif // __ATTENE_IMPLICIT_POINT3D_H_
