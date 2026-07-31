// attene_predicates.cpp - Expansion arithmetic and orientation predicate implementations
//
// Based on "Indirect Predicates for Geometric Constructions" by Marco Attene
// https://arxiv.org/abs/2105.09772

#include "attene_predicates.h"

// This TU evaluates interval arithmetic under FE_UPWARD and expansion
// arithmetic under FE_TONEAREST, switching modes at runtime. The compiler
// must not assume a rounding mode or move FP ops across the fesetround()
// boundaries.
#pragma float_control(precise, on)
#pragma fenv_access(on)

// =============================================================================
// Expansion Arithmetic Operations
// =============================================================================

int attene_growExpansion(int elen, const double* e, double b, double* h) {
	double Q = b;
	for (int i = 0; i < elen; i++) {
		double Qnew;
		attene_TwoSum(Q, e[i], Qnew, h[i]);
		Q = Qnew;
	}
	h[elen] = Q;
	return elen + 1;
}

int attene_growExpansionZeroElim(int elen, const double* e, double b, double* h) {
	int hindex = 0;
	double Q = b;
	for (int i = 0; i < elen; i++) {
		double Qnew, hh;
		attene_TwoSum(Q, e[i], Qnew, hh);
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}
	if (Q != 0.0 || hindex == 0) {
		h[hindex++] = Q;
	}
	return hindex;
}

int attene_fastExpansionSumZeroElim(int elen, const double* e, int flen, const double* f, double* h) {
	double Q;
	int hindex = 0;
	int eindex = 0;
	int findex = 0;

	double enow = e[0];
	double fnow = f[0];

	if ((fnow > enow) == (fnow > -enow)) {
		Q = enow;
		enow = (++eindex < elen) ? e[eindex] : 0.0;
	} else {
		Q = fnow;
		fnow = (++findex < flen) ? f[findex] : 0.0;
	}

	while (eindex < elen && findex < flen) {
		double Qnew, hh;
		if ((fnow > enow) == (fnow > -enow)) {
			attene_TwoSum(Q, enow, Qnew, hh);
			enow = (++eindex < elen) ? e[eindex] : 0.0;
		} else {
			attene_TwoSum(Q, fnow, Qnew, hh);
			fnow = (++findex < flen) ? f[findex] : 0.0;
		}
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}

	while (eindex < elen) {
		double Qnew, hh;
		attene_TwoSum(Q, enow, Qnew, hh);
		enow = (++eindex < elen) ? e[eindex] : 0.0;
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}

	while (findex < flen) {
		double Qnew, hh;
		attene_TwoSum(Q, fnow, Qnew, hh);
		fnow = (++findex < flen) ? f[findex] : 0.0;
		Q = Qnew;
		if (hh != 0.0) {
			h[hindex++] = hh;
		}
	}

	if (Q != 0.0 || hindex == 0) {
		h[hindex++] = Q;
	}
	return hindex;
}

int attene_scaleExpansionZeroElim(int elen, const double* e, double b, double* h) {
	int hindex = 0;
	double Q, sum;
	double bhi, blo;
	attene_Split(b, bhi, blo);

	double product1, product0;
	attene_TwoProd(e[0], b, product1, product0);

	double hh;
	if (product0 != 0) {
		h[hindex++] = product0;
	}
	Q = product1;

	for (int i = 1; i < elen; i++) {
		double enow = e[i];
		attene_TwoProd(enow, b, product1, product0);
		attene_TwoSum(Q, product0, sum, hh);
		if (hh != 0) {
			h[hindex++] = hh;
		}
		attene_TwoSum(product1, sum, Q, hh);
		if (hh != 0) {
			h[hindex++] = hh;
		}
	}

	if (Q != 0.0 || hindex == 0) {
		h[hindex++] = Q;
	}
	return hindex;
}

// =============================================================================
// Expansion Arithmetic Operators
// =============================================================================

attene_expansion operator+(const attene_expansion& a, const attene_expansion& b) {
	attene_expansion r;
	r.length = attene_fastExpansionSumZeroElim(a.length, a.components, b.length, b.components, r.components);
	return r;
}

attene_expansion operator-(const attene_expansion& a, const attene_expansion& b) {
	// Negate b then add
	double neg_b[ATTENE_EXPANSION_MAX_SIZE];
	for (int i = 0; i < b.length; i++) {
		neg_b[i] = -b.components[i];
	}
	attene_expansion r;
	r.length = attene_fastExpansionSumZeroElim(a.length, a.components, b.length, neg_b, r.components);
	return r;
}

attene_expansion operator*(const attene_expansion& a, double b) {
	attene_expansion r;
	r.length = attene_scaleExpansionZeroElim(a.length, a.components, b, r.components);
	return r;
}

attene_expansion operator*(double a, const attene_expansion& b) {
	return b * a;
}

attene_expansion operator*(const attene_expansion& a, const attene_expansion& b) {
	// Product of two expansions
	attene_expansion r;
	r.length = 0;

	if (a.length == 0 || b.length == 0) {
		r.components[0] = 0;
		r.length = 1;
		return r;
	}

	// Start with first term of b times all of a
	attene_expansion temp;
	temp.length = attene_scaleExpansionZeroElim(a.length, a.components, b.components[0], temp.components);
	r = temp;

	// Add each subsequent term
	for (int i = 1; i < b.length; i++) {
		temp.length = attene_scaleExpansionZeroElim(a.length, a.components, b.components[i], temp.components);
		attene_expansion sum;
		sum.length = attene_fastExpansionSumZeroElim(r.length, r.components, temp.length, temp.components, sum.components);
		r = sum;
	}

	return r;
}

// =============================================================================
// 3D Orientation Predicates
// =============================================================================

// Build a 2-component expansion holding a - b EXACTLY (TwoDiff head+tail).
// Wrapping the plain double difference (which rounds) in a 1-term expansion
// makes the "exact" expansion path lie for near-degenerate inputs.
static attene_expansion attene_expansion_diff(double a, double b) {
	attene_expansion r;
	double head, tail;
	attene_TwoDiff(a, b, head, tail);
	if (tail != 0.0) {
		r.components[0] = tail;
		r.components[1] = head;
		r.length = 2;
	} else {
		r.components[0] = head;
		r.length = 1;
	}
	return r;
}

// Enclosing interval for a - b under the current (FE_UPWARD) rounding mode.
// attene_interval(a - b) is a degenerate interval around the ROUNDED
// difference and does not contain the true value; interval subtraction does.
static inline attene_interval attene_interval_diff(double a, double b) {
	return attene_interval(a) - attene_interval(b);
}

int attene_orient3d_interval(
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz,
	double dx, double dy, double dz)
{
	attene_setFPUModeToRoundUP();

	attene_interval adx = attene_interval_diff(ax, dx);
	attene_interval ady = attene_interval_diff(ay, dy);
	attene_interval adz = attene_interval_diff(az, dz);
	attene_interval bdx = attene_interval_diff(bx, dx);
	attene_interval bdy = attene_interval_diff(by, dy);
	attene_interval bdz = attene_interval_diff(bz, dz);
	attene_interval cdx = attene_interval_diff(cx, dx);
	attene_interval cdy = attene_interval_diff(cy, dy);
	attene_interval cdz = attene_interval_diff(cz, dz);

	attene_interval m01 = adx * bdy - ady * bdx;
	attene_interval m02 = adx * bdz - adz * bdx;
	attene_interval m12 = ady * bdz - adz * bdy;

	attene_interval det = m12 * cdx - m02 * cdy + m01 * cdz;

	attene_setFPUModeToRoundNEAR();

	if (det.signIsReliable()) {
		return det.sign();
	}

	return 0;  // Uncertain - need exact computation
}

int attene_orient3d_expansion(
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz,
	double dx, double dy, double dz)
{
	attene_expansion adx = attene_expansion_diff(ax, dx);
	attene_expansion ady = attene_expansion_diff(ay, dy);
	attene_expansion adz = attene_expansion_diff(az, dz);
	attene_expansion bdx = attene_expansion_diff(bx, dx);
	attene_expansion bdy = attene_expansion_diff(by, dy);
	attene_expansion bdz = attene_expansion_diff(bz, dz);
	attene_expansion cdx = attene_expansion_diff(cx, dx);
	attene_expansion cdy = attene_expansion_diff(cy, dy);
	attene_expansion cdz = attene_expansion_diff(cz, dz);

	attene_expansion m01 = adx * bdy - ady * bdx;
	attene_expansion m02 = adx * bdz - adz * bdx;
	attene_expansion m12 = ady * bdz - adz * bdy;

	attene_expansion det = m12 * cdx - m02 * cdy + m01 * cdz;

	return det.sign();
}

int attene_orient3d(
	double ax, double ay, double az,
	double bx, double by, double bz,
	double cx, double cy, double cz,
	double dx, double dy, double dz)
{
	// Try interval arithmetic first
	int result = attene_orient3d_interval(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
	if (result != 0) {
		return result;
	}

	// Fall back to exact expansion arithmetic
	return attene_orient3d_expansion(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
}

int attene_orient3d(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d) {
	return attene_orient3d(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z);
}

// =============================================================================
// 2D Orientation Predicates
// =============================================================================

int attene_orient2d_interval(
	double ax, double ay,
	double bx, double by,
	double cx, double cy)
{
	attene_setFPUModeToRoundUP();

	attene_interval acx = attene_interval_diff(ax, cx);
	attene_interval acy = attene_interval_diff(ay, cy);
	attene_interval bcx = attene_interval_diff(bx, cx);
	attene_interval bcy = attene_interval_diff(by, cy);

	attene_interval det = acx * bcy - acy * bcx;

	attene_setFPUModeToRoundNEAR();

	if (det.signIsReliable()) {
		return det.sign();
	}

	return 0;
}

int attene_orient2d_expansion(
	double ax, double ay,
	double bx, double by,
	double cx, double cy)
{
	attene_expansion acx = attene_expansion_diff(ax, cx);
	attene_expansion acy = attene_expansion_diff(ay, cy);
	attene_expansion bcx = attene_expansion_diff(bx, cx);
	attene_expansion bcy = attene_expansion_diff(by, cy);

	attene_expansion det = acx * bcy - acy * bcx;

	return det.sign();
}

int attene_orient2d(
	double ax, double ay,
	double bx, double by,
	double cx, double cy)
{
	int result = attene_orient2d_interval(ax, ay, bx, by, cx, cy);
	if (result != 0) {
		return result;
	}
	return attene_orient2d_expansion(ax, ay, bx, by, cx, cy);
}

int attene_orient2d(const Vec2d& a, const Vec2d& b, const Vec2d& c) {
	return attene_orient2d(a.x, a.y, b.x, b.y, c.x, c.y);
}

// =============================================================================
// Comparison Predicates for Implicit Points
// =============================================================================

int attene_compare_ratios_same_denom_interval(
	const attene_interval& a, const attene_interval& b, const attene_interval& d)
{
	attene_interval diff = a - b;  // a - b

	if (!diff.signIsReliable() || !d.signIsReliable()) {
		return 0;  // Uncertain
	}

	int ds = d.sign();
	int diffs = diff.sign();

	if (ds == 0) return 0;  // Degenerate

	return (ds > 0) ? -diffs : diffs;  // Return -1 if a/d < b/d
}

int attene_compare_ratios_interval(
	const attene_interval& a, const attene_interval& da,
	const attene_interval& b, const attene_interval& db)
{
	if (!da.signIsReliable() || !db.signIsReliable()) {
		return 0;  // Uncertain
	}

	int das = da.sign();
	int dbs = db.sign();

	if (das == 0 || dbs == 0) return 0;  // Degenerate

	attene_interval lhs = a * db;
	attene_interval rhs = b * da;
	attene_interval diff = lhs - rhs;

	if (!diff.signIsReliable()) {
		return 0;
	}

	int diffs = diff.sign();

	// If signs of denominators differ, flip the comparison
	if ((das > 0) != (dbs > 0)) {
		return diffs;  // Flipped
	}
	return -diffs;  // Normal
}
