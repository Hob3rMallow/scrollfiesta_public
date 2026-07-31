// =============================================================================
// ATTENE_INDIRECT_PREDICATES3D.CPP - Indirect 3D Predicates Implementation
// =============================================================================
//
// Two-tier evaluation for predicates over implicit (lambda-represented) points:
//
//   1. Interval filter: lambdas and the predicate polynomial are evaluated in
//      interval arithmetic under FE_UPWARD. If the resulting sign is reliable,
//      return it.
//   2. Exact tier: lambdas and the polynomial are re-evaluated in BSNumber
//      arbitrary-precision arithmetic (math/bsnumber.h) and the true sign is
//      returned. Fixed-size expansion arithmetic is NOT usable here - products
//      of two implicit-point lambdas need tens of thousands of components,
//      far beyond ATTENE_EXPANSION_MAX_SIZE (the old code silently overran
//      its buffers in precisely the near-degenerate cases that matter).
//
// FPU-mode discipline: getIntervalLambda() brackets its own computation with
// FE_UPWARD/FE_TONEAREST (and cache hits touch nothing), so every predicate
// re-arms FE_UPWARD after fetching lambdas and before doing interval math,
// then restores FE_TONEAREST before returning or entering the exact tier.

#include "attene_implicit_point3d.h"
#include "bsnumber.h"

// Interval filters run under FE_UPWARD; the compiler must honor the runtime
// rounding mode (see attene_predicates.h).
#pragma float_control(precise, on)
#pragma fenv_access(on)

// =============================================================================
// Helper: 3x3 Determinant for Interval Arithmetic
// =============================================================================

static attene_interval det3x3_interval(
	const attene_interval& a00, const attene_interval& a01, const attene_interval& a02,
	const attene_interval& a10, const attene_interval& a11, const attene_interval& a12,
	const attene_interval& a20, const attene_interval& a21, const attene_interval& a22)
{
	return a00 * (a11 * a22 - a12 * a21) -
	       a01 * (a10 * a22 - a12 * a20) +
	       a02 * (a10 * a21 - a11 * a20);
}

static BSNumber det3x3_bs(
	const BSNumber& a00, const BSNumber& a01, const BSNumber& a02,
	const BSNumber& a10, const BSNumber& a11, const BSNumber& a12,
	const BSNumber& a20, const BSNumber& a21, const BSNumber& a22)
{
	return a00 * (a11 * a22 - a12 * a21) -
	       a01 * (a10 * a22 - a12 * a20) +
	       a02 * (a10 * a21 - a11 * a20);
}

static int bs_sign_i(int s) { return (s > 0) ? 1 : (s < 0) ? -1 : 0; }

// =============================================================================
// orient3d for Implicit Points
// =============================================================================
//
// We need to compute the sign of:
//   | ax-dx  ay-dy  az-dz |
//   | bx-dx  by-dy  bz-dz |
//   | cx-dx  cy-dy  cz-dz |
//
// When points are implicit (represented as lambda/denom), the differences
// become ratios with common denominators; the determinant's sign is the sign
// of the numerator determinant times the signs of the denominators.

static int orient3d_indirect_exact_bs(
	const AttenePoint3D* a, const AttenePoint3D* b,
	const AttenePoint3D* c, const AttenePoint3D* d)
{
	BSNumber lax, lay, laz, da;
	BSNumber lbx, lby, lbz, db;
	BSNumber lcx, lcy, lcz, dc;
	BSNumber ldx, ldy, ldz, dd;

	attene_get_bs_lambda(a, lax, lay, laz, da);
	attene_get_bs_lambda(b, lbx, lby, lbz, db);
	attene_get_bs_lambda(c, lcx, lcy, lcz, dc);
	attene_get_bs_lambda(d, ldx, ldy, ldz, dd);

	int sda = da.sign(), sdb = db.sign(), sdc = dc.sign(), sdd = dd.sign();
	if (sda == 0 || sdb == 0 || sdc == 0 || sdd == 0) return 0; // degenerate construction

	// Row i is (point_i - point_d) scaled by (d_i * d_d); the common positive
	// factor cancels out of the sign except for sign(da*db*dc*dd^3) =
	// sign(da)*sign(db)*sign(dc)*sign(dd).
	BSNumber r0x = lax * dd - ldx * da;
	BSNumber r0y = lay * dd - ldy * da;
	BSNumber r0z = laz * dd - ldz * da;

	BSNumber r1x = lbx * dd - ldx * db;
	BSNumber r1y = lby * dd - ldy * db;
	BSNumber r1z = lbz * dd - ldz * db;

	BSNumber r2x = lcx * dd - ldx * dc;
	BSNumber r2y = lcy * dd - ldy * dc;
	BSNumber r2z = lcz * dd - ldz * dc;

	BSNumber det = det3x3_bs(
		r0x, r0y, r0z,
		r1x, r1y, r1z,
		r2x, r2y, r2z);

	return bs_sign_i(det.sign()) * sda * sdb * sdc * sdd;
}

int attene_orient3d_indirect(
	const AttenePoint3D* a, const AttenePoint3D* b,
	const AttenePoint3D* c, const AttenePoint3D* d)
{
	// Fast path: all explicit
	if (a->isExplicit() && b->isExplicit() && c->isExplicit() && d->isExplicit()) {
		const AtteneExplicit3D* ea = static_cast<const AtteneExplicit3D*>(a);
		const AtteneExplicit3D* eb = static_cast<const AtteneExplicit3D*>(b);
		const AtteneExplicit3D* ec = static_cast<const AtteneExplicit3D*>(c);
		const AtteneExplicit3D* ed = static_cast<const AtteneExplicit3D*>(d);
		return attene_orient3d(
			ea->x, ea->y, ea->z,
			eb->x, eb->y, eb->z,
			ec->x, ec->y, ec->z,
			ed->x, ed->y, ed->z);
	}

	// Interval filter
	attene_interval lax, lay, laz, da;
	attene_interval lbx, lby, lbz, db;
	attene_interval lcx, lcy, lcz, dc;
	attene_interval ldx, ldy, ldz, dd;

	bool valid = true;
	valid &= a->getIntervalLambda(lax, lay, laz, da);
	valid &= b->getIntervalLambda(lbx, lby, lbz, db);
	valid &= c->getIntervalLambda(lcx, lcy, lcz, dc);
	valid &= d->getIntervalLambda(ldx, ldy, ldz, dd);

	if (valid) {
		attene_setFPUModeToRoundUP();

		attene_interval num_ax_dx = lax * dd - ldx * da;
		attene_interval num_ay_dy = lay * dd - ldy * da;
		attene_interval num_az_dz = laz * dd - ldz * da;
		attene_interval denom_ad = da * dd;

		attene_interval num_bx_dx = lbx * dd - ldx * db;
		attene_interval num_by_dy = lby * dd - ldy * db;
		attene_interval num_bz_dz = lbz * dd - ldz * db;
		attene_interval denom_bd = db * dd;

		attene_interval num_cx_dx = lcx * dd - ldx * dc;
		attene_interval num_cy_dy = lcy * dd - ldy * dc;
		attene_interval num_cz_dz = lcz * dd - ldz * dc;
		attene_interval denom_cd = dc * dd;

		attene_interval det_num = det3x3_interval(
			num_ax_dx, num_ay_dy, num_az_dz,
			num_bx_dx, num_by_dy, num_bz_dz,
			num_cx_dx, num_cy_dy, num_cz_dz);

		attene_interval det_denom = denom_ad * denom_bd * denom_cd;

		attene_setFPUModeToRoundNEAR();

		if (det_num.signIsReliable() && det_denom.signIsReliable()) {
			int denom_sign = det_denom.sign();
			if (denom_sign != 0) {
				return det_num.sign() * denom_sign;
			}
		}
	}

	// Exact tier
	return orient3d_indirect_exact_bs(a, b, c, d);
}

// =============================================================================
// Comparison Predicates for Sorting
// =============================================================================

// Interval comparison of la/da vs lb/db. Returns -1/0/+1, or sets *reliable
// to false when the interval information cannot decide.
static int compareLambdaOnAxis(
	const attene_interval& la, const attene_interval& da,
	const attene_interval& lb, const attene_interval& db,
	bool* reliable)
{
	*reliable = false;

	if (!da.signIsReliable() || !db.signIsReliable()) {
		return 0;
	}

	int das = da.sign();
	int dbs = db.sign();

	if (das == 0 || dbs == 0) return 0;

	attene_interval lhs = la * db;
	attene_interval rhs = lb * da;
	attene_interval diff = lhs - rhs;

	if (!diff.signIsReliable()) {
		return 0;
	}

	*reliable = true;
	int diffs = diff.sign();

	// la/da - lb/db = (la*db - lb*da) / (da*db)
	if ((das > 0) == (dbs > 0)) {
		return diffs;
	}
	return -diffs;
}

// Exact comparison of the axis-th lambda coordinate ratio (0=X, 1=Y, 2=Z).
static int compareLambdaOnAxisBS(const AttenePoint3D* a, const AttenePoint3D* b, int axis)
{
	BSNumber lax, lay, laz, da;
	BSNumber lbx, lby, lbz, db;
	attene_get_bs_lambda(a, lax, lay, laz, da);
	attene_get_bs_lambda(b, lbx, lby, lbz, db);

	int sda = da.sign(), sdb = db.sign();
	if (sda == 0 || sdb == 0) return 0; // degenerate construction

	const BSNumber& la = (axis == 0) ? lax : (axis == 1) ? lay : laz;
	const BSNumber& lb = (axis == 0) ? lbx : (axis == 1) ? lby : lbz;

	BSNumber diff = la * db - lb * da;
	return bs_sign_i(diff.sign()) * sda * sdb;
}

int attene_lessThanOnAxis(const AttenePoint3D* a, const AttenePoint3D* b, int axis)
{
	// Fast path: both explicit
	if (a->isExplicit() && b->isExplicit()) {
		const AtteneExplicit3D* ea = static_cast<const AtteneExplicit3D*>(a);
		const AtteneExplicit3D* eb = static_cast<const AtteneExplicit3D*>(b);
		double va = (axis == 0) ? ea->x : (axis == 1) ? ea->y : ea->z;
		double vb = (axis == 0) ? eb->x : (axis == 1) ? eb->y : eb->z;
		if (va < vb) return -1;
		if (va > vb) return 1;
		return 0;
	}

	attene_interval lax, lay, laz, da;
	attene_interval lbx, lby, lbz, db;

	bool valid = a->getIntervalLambda(lax, lay, laz, da);
	valid &= b->getIntervalLambda(lbx, lby, lbz, db);

	if (valid) {
		const attene_interval& la = (axis == 0) ? lax : (axis == 1) ? lay : laz;
		const attene_interval& lb = (axis == 0) ? lbx : (axis == 1) ? lby : lbz;

		attene_setFPUModeToRoundUP();
		bool reliable = false;
		int result = compareLambdaOnAxis(la, da, lb, db, &reliable);
		attene_setFPUModeToRoundNEAR();
		if (reliable) return result;
	}

	return compareLambdaOnAxisBS(a, b, axis);
}

int attene_lessThanOnX(const AttenePoint3D* a, const AttenePoint3D* b)
{
	return attene_lessThanOnAxis(a, b, 0);
}

int attene_lessThanOnY(const AttenePoint3D* a, const AttenePoint3D* b)
{
	return attene_lessThanOnAxis(a, b, 1);
}

int attene_lessThanOnZ(const AttenePoint3D* a, const AttenePoint3D* b)
{
	return attene_lessThanOnAxis(a, b, 2);
}

// =============================================================================
// Exact geometric equality
// =============================================================================

bool attene_samePoint3D(const AttenePoint3D* a, const AttenePoint3D* b)
{
	if (a == b) return true;

	if (a->isExplicit() && b->isExplicit()) {
		const AtteneExplicit3D* ea = static_cast<const AtteneExplicit3D*>(a);
		const AtteneExplicit3D* eb = static_cast<const AtteneExplicit3D*>(b);
		return ea->x == eb->x && ea->y == eb->y && ea->z == eb->z;
	}

	return attene_lessThanOnX(a, b) == 0 &&
	       attene_lessThanOnY(a, b) == 0 &&
	       attene_lessThanOnZ(a, b) == 0;
}

// =============================================================================
// Projected 2D orientation
// =============================================================================
//
// orient2d of the (u,v) projections of three 3D points, where (u,v) picks two
// of the three coordinate axes (ProjectionAxis convention: 0=XY, 1=YZ, 2=XZ).
// With A=(lau/da, lav/da) etc:
//
//   orient = (Bu-Au)(Cv-Av) - (Bv-Av)(Cu-Au)
//          = [ (lbu*da - lau*db)(lcv*da - lav*dc)
//            - (lbv*da - lav*db)(lcu*da - lau*dc) ] / (da^2 * db * dc)
//
// so sign(orient) = sign(numerator) * sign(db) * sign(dc), provided da != 0.

void attene_proj_axes(int axisPair, int* u, int* v)
{
	// 0 = XY, 1 = YZ, 2 = XZ - must match ProjectionAxis and the projection
	// mapping used by ConstrainedPolygon::getProjected2D / CoplanarProjectTo2D.
	switch (axisPair)
	{
	case 0: *u = 0; *v = 1; break;
	case 1: *u = 1; *v = 2; break;
	default: *u = 0; *v = 2; break;
	}
}

static int orient2d_proj_exact_bs(
	int u, int v,
	const AttenePoint3D* a, const AttenePoint3D* b, const AttenePoint3D* c)
{
	BSNumber la[3], da;
	BSNumber lb[3], db;
	BSNumber lc[3], dc;

	attene_get_bs_lambda(a, la[0], la[1], la[2], da);
	attene_get_bs_lambda(b, lb[0], lb[1], lb[2], db);
	attene_get_bs_lambda(c, lc[0], lc[1], lc[2], dc);

	int sda = da.sign(), sdb = db.sign(), sdc = dc.sign();
	if (sda == 0 || sdb == 0 || sdc == 0) return 0; // degenerate construction

	BSNumber n1 = lb[u] * da - la[u] * db;
	BSNumber n2 = lc[v] * da - la[v] * dc;
	BSNumber n3 = lb[v] * da - la[v] * db;
	BSNumber n4 = lc[u] * da - la[u] * dc;

	BSNumber det = n1 * n2 - n3 * n4;

	// Denominator da^2*db*dc: da^2 > 0, so only db and dc contribute sign.
	return bs_sign_i(det.sign()) * sdb * sdc;
}

int attene_orient2d_proj(
	int axisPair,
	const AttenePoint3D* a, const AttenePoint3D* b, const AttenePoint3D* c)
{
	int u, v;
	attene_proj_axes(axisPair, &u, &v);

	// Fast path: all explicit
	if (a->isExplicit() && b->isExplicit() && c->isExplicit()) {
		const AtteneExplicit3D* ea = static_cast<const AtteneExplicit3D*>(a);
		const AtteneExplicit3D* eb = static_cast<const AtteneExplicit3D*>(b);
		const AtteneExplicit3D* ec = static_cast<const AtteneExplicit3D*>(c);
		const double* pa = ea->ptr();
		const double* pb = eb->ptr();
		const double* pc = ec->ptr();
		return attene_orient2d(pa[u], pa[v], pb[u], pb[v], pc[u], pc[v]);
	}

	// Interval filter
	attene_interval la[3], da;
	attene_interval lb[3], db;
	attene_interval lc[3], dc;

	bool valid = a->getIntervalLambda(la[0], la[1], la[2], da);
	valid &= b->getIntervalLambda(lb[0], lb[1], lb[2], db);
	valid &= c->getIntervalLambda(lc[0], lc[1], lc[2], dc);

	if (valid) {
		attene_setFPUModeToRoundUP();

		attene_interval n1 = lb[u] * da - la[u] * db;
		attene_interval n2 = lc[v] * da - la[v] * dc;
		attene_interval n3 = lb[v] * da - la[v] * db;
		attene_interval n4 = lc[u] * da - la[u] * dc;

		attene_interval det = n1 * n2 - n3 * n4;
		attene_interval denom = db * dc;

		attene_setFPUModeToRoundNEAR();

		if (det.signIsReliable() && denom.signIsReliable()) {
			int denom_sign = denom.sign();
			if (denom_sign != 0) {
				return det.sign() * denom_sign;
			}
		}
	}

	// Exact tier
	return orient2d_proj_exact_bs(u, v, a, b, c);
}
