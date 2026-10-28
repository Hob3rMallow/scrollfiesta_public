/* A growth search changes its trial map many times while its starting map
 * stays fixed. Retain that face-local frame for exactly one component search.
 * The moving endpoint and the full determinant quadratic are still checked. */
typedef struct SRGrowthReference {
    double e[2], g[2], determinant;
    int positive;
} SRGrowthReference;

static void sr_growth_reference_set(SRGrowthReference *r, const double *uv,
                                     const int32_t *f)
{
    r->positive = UvGuard_positive(uv,f);
    for (size_t k = 0; k < 2; k++) {
        size_t a = 2*(size_t)f[0]+k, b = 2*(size_t)f[1]+k, c = 2*(size_t)f[2]+k;
        r->e[k] = uv[b]-uv[a]; r->g[k] = uv[c]-uv[a];
    }
    r->determinant = r->e[0]*r->g[1]-r->e[1]*r->g[0];
}

static int sr_growth_reference_safe(const SRGrowthReference *r,
                                      const double *trial, const int32_t *f)
{
    double de[2], dg[2];
    if (!r->positive) return 0;
    for (size_t k = 0; k < 2; k++) {
        size_t a = 2*(size_t)f[0]+k, b = 2*(size_t)f[1]+k, c = 2*(size_t)f[2]+k;
        de[k] = (trial[b]-trial[a])-r->e[k];
        dg[k] = (trial[c]-trial[a])-r->g[k];
    }
    double a = de[0]*dg[1]-de[1]*dg[0];
    double b = de[0]*r->g[1]-de[1]*r->g[0]+r->e[0]*dg[1]-r->e[1]*dg[0];
    double c = r->determinant;
    double minimum = fmin(c,(a+b)+c);
    if (a > 0) {
        double t = -.5*b/a;
        if (t > 0 && t < 1) minimum = fmin(minimum,(a*t+b)*t+c);
    }
    double scale = fabs(c)+fabs(b)+fabs(a);
    return isfinite(scale) && minimum > 64*DBL_EPSILON*scale;
}
