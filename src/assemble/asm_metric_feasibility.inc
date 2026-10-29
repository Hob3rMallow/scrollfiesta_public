/* Necessary original-metric condition for an injective interior vertex star.
 * Lin and Sinnamon, The Generalized Wielandt Inequality, Theorem 3.1:
 * https://www.math.uwo.ca/faculty/sinnamon/pdf/wielandt.pdf
 * For singular values in [.75,1.25], kappa <= 5/3 and a source angle a maps
 * inside [2 atan(tan(a/2)/kappa), 2 atan(kappa*tan(a/2))]. A single oriented
 * disk fan without overlap must sum to 2*pi. We certify the single closed
 * link combinatorially before using that condition. This is not a sufficient
 * feasibility test, and no source geometry or gate is changed. */
typedef struct AmfStats {
    size_t closed_stars, obstruction_stars, affected_faces;
    size_t orientation_obstructions;
    double affected_area, largest_gap;
    int io;
} AmfStats;

enum { AMF_SINGLE_TURN = 1, AMF_POSITIVE_ORIENTATION = 2 };

static void amf_bounds(double angle, double *lower, double *upper)
{
    double sn = sin(.5*angle), cs = cos(.5*angle), kappa = 1.25/.75;
    *lower = 2*atan2(sn,kappa*cs); *upper = 2*atan2(kappa*sn,cs);
}

static AmfStats amf_chart_marks(const AsmChart *ch, FILE *csv, uint8_t *marks)
{
    AmfStats out = {0}; out.io = 1; Arena_T arena = Arena_new();
    if (marks) memset(marks,0,ch->nv);
    size_t *offset = ARENA_CALLOC(arena,ch->nv+1,sizeof(size_t));
    for (size_t k = 0; k < 3*ch->nf; k++) {
        if (ch->faces[k] < 0 || (size_t)ch->faces[k] >= ch->nv) { out.io = 0; goto done; }
        offset[(size_t)ch->faces[k]+1]++;
    }
    size_t maximum = 0;
    for (size_t v = 0; v < ch->nv; v++) { if (offset[v+1] > maximum) maximum = offset[v+1]; offset[v+1] += offset[v]; }
    size_t *cursor = ARENA_ALLOC(arena,ch->nv*sizeof(size_t)); memcpy(cursor,offset,ch->nv*sizeof(size_t));
    size_t *corners = ARENA_ALLOC(arena,3*ch->nf*sizeof(size_t));
    for (size_t k = 0; k < 3*ch->nf; k++) corners[cursor[ch->faces[k]]++] = k;
    int32_t *a = ARENA_ALLOC(arena,maximum*sizeof(int32_t)), *b = ARENA_ALLOC(arena,maximum*sizeof(int32_t));
    uint8_t *visited = ARENA_ALLOC(arena,maximum), *affected = ARENA_CALLOC(arena,ch->nf,1);
    const double full = 2*acos(-1.);
    for (size_t v = 0; v < ch->nv; v++) {
        size_t n = offset[v+1]-offset[v]; if (n < 3) continue;
        for (size_t k = 0; k < n; k++) {
            size_t corner = corners[offset[v]+k], f = corner/3, local = corner%3;
            a[k] = ch->faces[3*f+(local+1)%3]; b[k] = ch->faces[3*f+(local+2)%3];
        }
        memset(visited,0,n); int32_t next = a[0]; int closed = 1;
        for (size_t k = 0; k < n; k++) {
            size_t matches = 0, chosen = 0;
            for (size_t j = 0; j < n; j++) if (a[j] == next) { matches++; chosen = j; }
            if (matches != 1 || visited[chosen] || a[chosen] == b[chosen] || (size_t)a[chosen] == v || (size_t)b[chosen] == v) { closed = 0; break; }
            visited[chosen] = 1; next = b[chosen];
        }
        if (!closed || next != a[0]) continue;
        double sum = 0, lower = 0, upper = 0;
        for (size_t k = 0; k < n; k++) {
            double x[3], y[3], cross[3], dot = 0;
            for (int d = 0; d < 3; d++) {
                x[d] = (double)ch->xyz[3*(size_t)a[k]+d]-ch->xyz[3*v+d];
                y[d] = (double)ch->xyz[3*(size_t)b[k]+d]-ch->xyz[3*v+d]; dot += x[d]*y[d];
            }
            for (int d = 0; d < 3; d++) cross[d] = x[(d+1)%3]*y[(d+2)%3]-x[(d+2)%3]*y[(d+1)%3];
            double length = hypot(hypot(cross[0],cross[1]),cross[2]);
            if (!(length > 0) || !isfinite(length) || !isfinite(dot)) { closed = 0; break; }
            double angle = atan2(length,dot), lo, hi; amf_bounds(angle,&lo,&hi);
            sum += angle; lower += lo; upper += hi;
        }
        if (!closed) continue; out.closed_stars++;
        /* The margin is deliberately much larger than accumulated roundoff. */
        double gap = fmax(lower-full,full-upper);
        if (!(gap > 1e-6 + 128*DBL_EPSILON*(double)n)) continue;
        out.obstruction_stars++; out.largest_gap = fmax(out.largest_gap,gap);
        /* Positive triangles need an integer number of complete turns.
         * A small multiply covered fan may be allowed by the material band.
         * Excluding just one turn proves non-injectivity, not necessarily
         * failure of the current admission certificate. */
        double margin = 1e-6 + 128*DBL_EPSILON*(double)n;
        double turns = fmax(1,ceil((lower-margin)/full));
        int impossible = turns*full > upper+margin;
        out.orientation_obstructions += impossible;
        if (marks) marks[v] = (uint8_t)(AMF_SINGLE_TURN | (impossible ? AMF_POSITIVE_ORIENTATION : 0));
        for (size_t k = 0; k < n; k++) affected[corners[offset[v]+k]/3] = 1;
        if (csv) {
            if (fprintf(csv,"%d,%zu,%d,%zu,%.17g,%.17g,%.17g,%.17g,%s,%.9g,%.9g,%.9g,",
                ch->id,v,ch->vid?ch->vid[v]:-1,n,sum,lower,upper,gap,upper<full?"angle_deficit":"angle_excess",
                (double)ch->xyz[3*v],(double)ch->xyz[3*v+1],(double)ch->xyz[3*v+2]) < 0) out.io = 0;
            for (size_t k = 0; k < n; k++) if (fprintf(csv,"%s%zu",k?" ":"",corners[offset[v]+k]/3) < 0) out.io = 0;
            if (fprintf(csv,"\n") < 0) out.io = 0;
        }
    }
    AsmMetricFace *metric = ARENA_ALLOC(arena,ch->nf*sizeof(AsmMetricFace));
    if (AsmMetric_prepare(ch,metric)) { out.io = 0; goto done; }
    for (size_t f = 0; f < ch->nf; f++) if (affected[f]) { out.affected_faces++; out.affected_area += metric[f].area; }
done:
    Arena_dispose(&arena); return out;
}

static AmfStats amf_chart(const AsmChart *ch, FILE *csv)
{ return amf_chart_marks(ch,csv,NULL); }

static int amf_selftest(void)
{
    int failures = 0;
    int bounds_ok = 1;
    /* Independent direct transforms, including obtuse source angles, must
     * stay inside the necessary interval used by the obstruction proof. */
    for (int a = 1; a < 200; a++) for (int direction = 0; direction < 50; direction++) {
        double angle = acos(-1.)*a/200, phi = 2*acos(-1.)*direction/50, lower, upper;
        double x[2] = {.75*cos(phi),1.25*sin(phi)};
        double y[2] = {.75*cos(phi+angle),1.25*sin(phi+angle)};
        double actual = atan2(fabs(x[0]*y[1]-x[1]*y[0]),x[0]*y[0]+x[1]*y[1]);
        amf_bounds(angle,&lower,&upper);
        if (actual < lower-1e-12 || actual > upper+1e-12) bounds_ok = 0;
    }
    fprintf(stderr,"  native metric angle bounds against 9950 direct transforms: %s\n",bounds_ok?"ok":"FAIL");
    failures += !bounds_ok;
    for (int test = 0; test < 8; test++) {
        float xyz[17*3] = {0}; int32_t faces[16*3]; size_t nf = test == 4 ? 7 : test == 5 ? 9 : test == 6 ? 16 : 8;
        for (int k = 0; k < 16; k++) {
            double angle = 2*acos(-1.)*(k%8)/8, radius = k<8?1:2;
            xyz[3*(k+1)] = (float)(radius*cos(angle)); xyz[3*(k+1)+1] = (float)(radius*sin(angle));
            xyz[3*(k+1)+2] = test == 2 ? (k%2?2.f:-2.f) : 0.f;
            faces[3*k] = 0; faces[3*k+1] = k+1; faces[3*k+2] = k/8*8+(k+1)%8+1;
        }
        if (test == 1 || test == 3 || test == 4 || test == 5 || test == 6 || test == 7) xyz[2] = 4;
        if (test == 3) for (int k = 0; k < 17*3; k++) xyz[k] += 1024;
        if (test == 5) memcpy(faces+24,faces,3*sizeof(int32_t));
        if (test == 7) for (size_t k = 0; k < nf; k++) { int32_t swap = faces[3*k+1]; faces[3*k+1] = faces[3*k+2]; faces[3*k+2] = swap; }
        AsmChart ch = {0}; ch.nv = 17; ch.nf = nf; ch.xyz = xyz; ch.faces = faces;
        AmfStats st = amf_chart(&ch,NULL); int expected = test == 1 || test == 2 || test == 3 || test == 7;
        int okay = st.io && st.obstruction_stars == (size_t)expected && st.closed_stars == (size_t)(test<4 || test==7);
        fprintf(stderr,"  native metric star bound %d: %zu closed, %zu obstructions, %s\n",test,st.closed_stars,st.obstruction_stars,okay?"ok":"FAIL");
        failures += !okay;
    }
    return failures;
}
