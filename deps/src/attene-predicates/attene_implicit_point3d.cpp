// attene_implicit_point3d.cpp - 3D Implicit point method implementations
//
// Based on "Indirect Predicates for Geometric Constructions" by Marco Attene
// and "Interactive and Robust Mesh Booleans" by Cherchi et al.

#include "attene_implicit_point3d.h"
#include "bsnumber.h"
#include <cassert>

// Interval lambdas are computed under FE_UPWARD; the compiler must honor the
// runtime rounding mode (see attene_predicates.h).
#pragma float_control(precise, on)
#pragma fenv_access(on)

// Enclosing interval for a - b under the current (FE_UPWARD) rounding mode.
// attene_interval(a - b) wraps the ROUNDED plain-double difference in a
// degenerate interval that does not contain the true value, which makes every
// downstream "reliable" sign potentially wrong. Interval subtraction encloses.
static inline attene_interval ip3d_interval_diff(double a, double b) {
	return attene_interval(a) - attene_interval(b);
}

// =============================================================================
// AtteneLPI3D - Line-Plane Intersection
// =============================================================================

bool AtteneLPI3D::getIntervalLambda(
	attene_interval& lx, attene_interval& ly, attene_interval& lz,
	attene_interval& d) const
{
	if (cache_computed) {
		lx = cached_lx;
		ly = cached_ly;
		lz = cached_lz;
		d = cached_d;
		return cache_valid;
	}

	attene_setFPUModeToRoundUP();

	// Direction of segment: v = Q - P
	attene_interval vx = ip3d_interval_diff(Q->x, P->x);
	attene_interval vy = ip3d_interval_diff(Q->y, P->y);
	attene_interval vz = ip3d_interval_diff(Q->z, P->z);

	// Two edges of the plane triangle
	// e1 = S - R, e2 = T - R
	attene_interval e1x = ip3d_interval_diff(S->x, R->x);
	attene_interval e1y = ip3d_interval_diff(S->y, R->y);
	attene_interval e1z = ip3d_interval_diff(S->z, R->z);
	attene_interval e2x = ip3d_interval_diff(T->x, R->x);
	attene_interval e2y = ip3d_interval_diff(T->y, R->y);
	attene_interval e2z = ip3d_interval_diff(T->z, R->z);

	// Normal = e1 x e2
	attene_interval nx = e1y * e2z - e1z * e2y;
	attene_interval ny = e1z * e2x - e1x * e2z;
	attene_interval nz = e1x * e2y - e1y * e2x;

	// d = dot(v, n) = vx*nx + vy*ny + vz*nz
	attene_interval denom = vx * nx + vy * ny + vz * nz;

	// Check for parallel (line parallel to plane)
	if (denom.containsZero()) {
		attene_setFPUModeToRoundNEAR();
		cache_computed = true;
		cache_valid = false;
		return false;
	}

	// w = R - P
	attene_interval wx = ip3d_interval_diff(R->x, P->x);
	attene_interval wy = ip3d_interval_diff(R->y, P->y);
	attene_interval wz = ip3d_interval_diff(R->z, P->z);

	// t_num = dot(w, n) (numerator for parameter t = t_num / d)
	attene_interval t_num = wx * nx + wy * ny + wz * nz;

	// Intersection point I = P + t * v = P + (t_num/d) * v
	// Multiply through by d:
	// I * d = P * d + t_num * v
	// So: lx = P.x * d + t_num * vx, etc.

	cached_lx = attene_interval(P->x) * denom + t_num * vx;
	cached_ly = attene_interval(P->y) * denom + t_num * vy;
	cached_lz = attene_interval(P->z) * denom + t_num * vz;
	cached_d = denom;

	attene_setFPUModeToRoundNEAR();

	cache_computed = true;
	cache_valid = true;

	lx = cached_lx;
	ly = cached_ly;
	lz = cached_lz;
	d = cached_d;
	return true;
}

bool AtteneLPI3D::getApproxXYZ(double& ox, double& oy, double& oz) const {
	attene_interval lx, ly, lz, d;
	if (getIntervalLambda(lx, ly, lz, d)) {
		double dval = d.estimate();
		if (dval != 0) {
			ox = lx.estimate() / dval;
			oy = ly.estimate() / dval;
			oz = lz.estimate() / dval;
			return true;
		}
	}

	// Interval representation could not prove the construction non-degenerate
	// (near-parallel line/plane). Fall back to the exact lambda so callers
	// never observe uninitialized coordinates.
	BSNumber blx, bly, blz, bd;
	attene_get_bs_lambda(this, blx, bly, blz, bd);
	if (bd.sign() == 0) {
		ox = oy = oz = 0.0;
		return false;
	}
	double dval = (double)bd;
	if (dval == 0.0) dval = (bd.sign() > 0) ? DBL_MIN : -DBL_MIN;
	ox = (double)blx / dval;
	oy = (double)bly / dval;
	oz = (double)blz / dval;
	return true;
}

// =============================================================================
// AtteneTPI3D - Three-Plane Intersection
// =============================================================================

void AtteneTPI3D::computePlaneInterval(
	const AtteneExplicit3D* A, const AtteneExplicit3D* B, const AtteneExplicit3D* C,
	attene_interval& nx, attene_interval& ny, attene_interval& nz, attene_interval& d)
{
	// e1 = B - A, e2 = C - A
	attene_interval e1x = ip3d_interval_diff(B->x, A->x);
	attene_interval e1y = ip3d_interval_diff(B->y, A->y);
	attene_interval e1z = ip3d_interval_diff(B->z, A->z);
	attene_interval e2x = ip3d_interval_diff(C->x, A->x);
	attene_interval e2y = ip3d_interval_diff(C->y, A->y);
	attene_interval e2z = ip3d_interval_diff(C->z, A->z);

	// n = e1 x e2
	nx = e1y * e2z - e1z * e2y;
	ny = e1z * e2x - e1x * e2z;
	nz = e1x * e2y - e1y * e2x;

	// d = dot(n, A)
	d = nx * attene_interval(A->x) + ny * attene_interval(A->y) + nz * attene_interval(A->z);
}


bool AtteneTPI3D::getIntervalLambda(
	attene_interval& lx, attene_interval& ly, attene_interval& lz,
	attene_interval& d) const
{
	if (cache_computed) {
		lx = cached_lx;
		ly = cached_ly;
		lz = cached_lz;
		d = cached_d;
		return cache_valid;
	}

	attene_setFPUModeToRoundUP();

	// Get the three plane equations
	attene_interval n1x, n1y, n1z, d1;
	attene_interval n2x, n2y, n2z, d2;
	attene_interval n3x, n3y, n3z, d3;

	computePlaneInterval(P1, Q1, R1, n1x, n1y, n1z, d1);
	computePlaneInterval(P2, Q2, R2, n2x, n2y, n2z, d2);
	computePlaneInterval(P3, Q3, R3, n3x, n3y, n3z, d3);

	// Solve the 3x3 system using Cramer's rule
	// | n1x n1y n1z |   | x |   | d1 |
	// | n2x n2y n2z | * | y | = | d2 |
	// | n3x n3y n3z |   | z |   | d3 |

	// Determinant of coefficient matrix
	attene_interval det =
		n1x * (n2y * n3z - n2z * n3y) -
		n1y * (n2x * n3z - n2z * n3x) +
		n1z * (n2x * n3y - n2y * n3x);

	if (det.containsZero()) {
		attene_setFPUModeToRoundNEAR();
		cache_computed = true;
		cache_valid = false;
		return false;
	}

	// x = det(d, n1y, n1z; ...) / det
	cached_lx =
		d1 * (n2y * n3z - n2z * n3y) -
		n1y * (d2 * n3z - n2z * d3) +
		n1z * (d2 * n3y - n2y * d3);

	// y = det(n1x, d, n1z; ...) / det
	cached_ly =
		n1x * (d2 * n3z - n2z * d3) -
		d1 * (n2x * n3z - n2z * n3x) +
		n1z * (n2x * d3 - d2 * n3x);

	// z = det(n1x, n1y, d; ...) / det
	cached_lz =
		n1x * (n2y * d3 - d2 * n3y) -
		n1y * (n2x * d3 - d2 * n3x) +
		d1 * (n2x * n3y - n2y * n3x);

	cached_d = det;

	attene_setFPUModeToRoundNEAR();

	cache_computed = true;
	cache_valid = true;

	lx = cached_lx;
	ly = cached_ly;
	lz = cached_lz;
	d = cached_d;
	return true;
}

bool AtteneTPI3D::getApproxXYZ(double& ox, double& oy, double& oz) const {
	attene_interval lx, ly, lz, d;
	if (getIntervalLambda(lx, ly, lz, d)) {
		double dval = d.estimate();
		if (dval != 0) {
			ox = lx.estimate() / dval;
			oy = ly.estimate() / dval;
			oz = lz.estimate() / dval;
			return true;
		}
	}

	// Fall back to the exact lambda (see AtteneLPI3D::getApproxXYZ).
	BSNumber blx, bly, blz, bd;
	attene_get_bs_lambda(this, blx, bly, blz, bd);
	if (bd.sign() == 0) {
		ox = oy = oz = 0.0;
		return false;
	}
	double dval = (double)bd;
	if (dval == 0.0) dval = (bd.sign() > 0) ? DBL_MIN : -DBL_MIN;
	ox = (double)blx / dval;
	oy = (double)bly / dval;
	oz = (double)blz / dval;
	return true;
}

// =============================================================================
// Exact lambda via BSNumber (arbitrary precision, never overflows)
// =============================================================================

static void bs_lambda_explicit(const AtteneExplicit3D* p,
	BSNumber& lx, BSNumber& ly, BSNumber& lz, BSNumber& d)
{
	lx = BSNumber(p->x);
	ly = BSNumber(p->y);
	lz = BSNumber(p->z);
	d = BSNumber((int32_t)1);
}

static void bs_lambda_lpi(const AtteneLPI3D* p,
	BSNumber& lx, BSNumber& ly, BSNumber& lz, BSNumber& d)
{
	// Same construction as the interval lambda, evaluated exactly.
	BSNumber vx = BSNumber(p->Q->x) - BSNumber(p->P->x);
	BSNumber vy = BSNumber(p->Q->y) - BSNumber(p->P->y);
	BSNumber vz = BSNumber(p->Q->z) - BSNumber(p->P->z);

	BSNumber e1x = BSNumber(p->S->x) - BSNumber(p->R->x);
	BSNumber e1y = BSNumber(p->S->y) - BSNumber(p->R->y);
	BSNumber e1z = BSNumber(p->S->z) - BSNumber(p->R->z);
	BSNumber e2x = BSNumber(p->T->x) - BSNumber(p->R->x);
	BSNumber e2y = BSNumber(p->T->y) - BSNumber(p->R->y);
	BSNumber e2z = BSNumber(p->T->z) - BSNumber(p->R->z);

	BSNumber nx = e1y * e2z - e1z * e2y;
	BSNumber ny = e1z * e2x - e1x * e2z;
	BSNumber nz = e1x * e2y - e1y * e2x;

	BSNumber denom = vx * nx + vy * ny + vz * nz;

	BSNumber wx = BSNumber(p->R->x) - BSNumber(p->P->x);
	BSNumber wy = BSNumber(p->R->y) - BSNumber(p->P->y);
	BSNumber wz = BSNumber(p->R->z) - BSNumber(p->P->z);

	BSNumber t_num = wx * nx + wy * ny + wz * nz;

	lx = BSNumber(p->P->x) * denom + t_num * vx;
	ly = BSNumber(p->P->y) * denom + t_num * vy;
	lz = BSNumber(p->P->z) * denom + t_num * vz;
	d = denom;
}

static void bs_plane(const AtteneExplicit3D* A, const AtteneExplicit3D* B, const AtteneExplicit3D* C,
	BSNumber& nx, BSNumber& ny, BSNumber& nz, BSNumber& d)
{
	BSNumber e1x = BSNumber(B->x) - BSNumber(A->x);
	BSNumber e1y = BSNumber(B->y) - BSNumber(A->y);
	BSNumber e1z = BSNumber(B->z) - BSNumber(A->z);
	BSNumber e2x = BSNumber(C->x) - BSNumber(A->x);
	BSNumber e2y = BSNumber(C->y) - BSNumber(A->y);
	BSNumber e2z = BSNumber(C->z) - BSNumber(A->z);

	nx = e1y * e2z - e1z * e2y;
	ny = e1z * e2x - e1x * e2z;
	nz = e1x * e2y - e1y * e2x;

	d = nx * BSNumber(A->x) + ny * BSNumber(A->y) + nz * BSNumber(A->z);
}

static void bs_lambda_tpi(const AtteneTPI3D* p,
	BSNumber& lx, BSNumber& ly, BSNumber& lz, BSNumber& d)
{
	BSNumber n1x, n1y, n1z, d1;
	BSNumber n2x, n2y, n2z, d2;
	BSNumber n3x, n3y, n3z, d3;

	bs_plane(p->P1, p->Q1, p->R1, n1x, n1y, n1z, d1);
	bs_plane(p->P2, p->Q2, p->R2, n2x, n2y, n2z, d2);
	bs_plane(p->P3, p->Q3, p->R3, n3x, n3y, n3z, d3);

	d = n1x * (n2y * n3z - n2z * n3y) -
	    n1y * (n2x * n3z - n2z * n3x) +
	    n1z * (n2x * n3y - n2y * n3x);

	lx = d1 * (n2y * n3z - n2z * n3y) -
	     n1y * (d2 * n3z - n2z * d3) +
	     n1z * (d2 * n3y - n2y * d3);

	ly = n1x * (d2 * n3z - n2z * d3) -
	     d1 * (n2x * n3z - n2z * n3x) +
	     n1z * (n2x * d3 - d2 * n3x);

	lz = n1x * (n2y * d3 - d2 * n3y) -
	     n1y * (n2x * d3 - d2 * n3x) +
	     d1 * (n2x * n3y - n2y * n3x);
}

void attene_get_bs_lambda(
	const AttenePoint3D* p,
	BSNumber& lx, BSNumber& ly, BSNumber& lz, BSNumber& d)
{
	assert(p != NULL);
	switch (p->getType())
	{
	case ATTENE_EXPLICIT_3D:
		bs_lambda_explicit(static_cast<const AtteneExplicit3D*>(p), lx, ly, lz, d);
		break;
	case ATTENE_LPI_3D:
		bs_lambda_lpi(static_cast<const AtteneLPI3D*>(p), lx, ly, lz, d);
		break;
	case ATTENE_TPI_3D:
		bs_lambda_tpi(static_cast<const AtteneTPI3D*>(p), lx, ly, lz, d);
		break;
	}
}
