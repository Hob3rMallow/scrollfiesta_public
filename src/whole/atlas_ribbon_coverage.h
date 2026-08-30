/* atlas_ribbon_coverage.h -- the coverage-audit reason codes shared by the
 * producer (atlas_ribbon_fit_tool's write_coverage_audit, which rasterizes
 * ribbon_coverage_reason.tif and the JSON/CSV audit) and the viewer
 * (strip_preview --mode coverage).  Names, display priority, and palette
 * live here ONLY, so the TIF bytes, the legend CSV, and the preview colours
 * cannot drift apart.  Header-only: viewers link nothing extra.
 *
 * Display priority is used ONLY when several source layers paint the same
 * atlas pixel (and by the preview's max-rank downsample); every layer is
 * still counted separately in the JSON/CSV audit. */
#ifndef ATLAS_RIBBON_COVERAGE_H
#define ATLAS_RIBBON_COVERAGE_H

#include <stdint.h>

typedef enum {
    ATLAS_RIBBON_COVERAGE_BACKGROUND = 0,
    ATLAS_RIBBON_COVERAGE_RIBBON = 1,
    ATLAS_RIBBON_COVERAGE_PREFIT_CULL = 2,
    ATLAS_RIBBON_COVERAGE_INVALID_CHART = 3,
    ATLAS_RIBBON_COVERAGE_OUTSIDE_V = 4,
    ATLAS_RIBBON_COVERAGE_OUTSIDE_U = 5,
    ATLAS_RIBBON_COVERAGE_NO_U_SUPPORT = 6,
    ATLAS_RIBBON_COVERAGE_SUBGRID_U = 7,
    ATLAS_RIBBON_COVERAGE_U_GAP = 8,
    ATLAS_RIBBON_COVERAGE_METRIC_U = 9,
    ATLAS_RIBBON_COVERAGE_TOPOLOGY_U = 10,
    ATLAS_RIBBON_COVERAGE_CROSSING_U = 11,
    ATLAS_RIBBON_COVERAGE_VERTICAL_METRIC = 12,
    ATLAS_RIBBON_COVERAGE_GEOMETRY_MISMATCH = 13,
    ATLAS_RIBBON_COVERAGE_VERTICAL_METRIC_FILL = 14,
    ATLAS_RIBBON_COVERAGE_REASON_COUNT = 15
} AtlasRibbonCoverageReason;

static inline const char *AtlasRibbonCoverage_name(int reason)
{
    static const char *name[ATLAS_RIBBON_COVERAGE_REASON_COUNT] = {
        "background",
        "ribboned",
        "prefit_cull",
        "invalid_chart",
        "outside_v_lattice",
        "outside_u_lattice",
        "no_u_support",
        "subgrid_u_island",
        "u_gap_cut",
        "metric_u_cut",
        "topology_u_cut",
        "crossing_u_cut",
        "vertical_metric_cut",
        "geometry_mismatch",
        "vertical_metric_fill_cut"
    };
    return reason >= 0 && reason < ATLAS_RIBBON_COVERAGE_REASON_COUNT
         ? name[reason] : "invalid";
}

static inline int AtlasRibbonCoverage_priority(int reason)
{
    static const uint8_t priority[ATLAS_RIBBON_COVERAGE_REASON_COUNT] = {
        0, 1, 4, 5, 2, 2, 3, 6, 7, 9, 10, 12, 11, 13, 8
    };
    return reason >= 0 && reason < ATLAS_RIBBON_COVERAGE_REASON_COUNT
         ? priority[reason] : 0;
}

static inline void AtlasRibbonCoverage_rgb(int reason, uint8_t *rgb)
{
    switch (reason) {
    case ATLAS_RIBBON_COVERAGE_BACKGROUND:
        rgb[0] = 12;  rgb[1] = 14;  rgb[2] = 40;  break;
    case ATLAS_RIBBON_COVERAGE_RIBBON:
        rgb[0] = 155; rgb[1] = 155; rgb[2] = 155; break;
    case ATLAS_RIBBON_COVERAGE_PREFIT_CULL:
        rgb[0] = 255; rgb[1] = 80;  rgb[2] = 180; break;
    case ATLAS_RIBBON_COVERAGE_INVALID_CHART:
        rgb[0] = 30;  rgb[1] = 30;  rgb[2] = 30;  break;
    case ATLAS_RIBBON_COVERAGE_OUTSIDE_V:
        rgb[0] = 60;  rgb[1] = 100; rgb[2] = 220; break;
    case ATLAS_RIBBON_COVERAGE_OUTSIDE_U:
        rgb[0] = 40;  rgb[1] = 60;  rgb[2] = 145; break;
    case ATLAS_RIBBON_COVERAGE_NO_U_SUPPORT:
        rgb[0] = 245; rgb[1] = 220; rgb[2] = 40;  break;
    case ATLAS_RIBBON_COVERAGE_SUBGRID_U:
        rgb[0] = 235; rgb[1] = 40;  rgb[2] = 235; break;
    case ATLAS_RIBBON_COVERAGE_U_GAP:
        rgb[0] = 230; rgb[1] = 130; rgb[2] = 30;  break;
    case ATLAS_RIBBON_COVERAGE_METRIC_U:
        rgb[0] = 230; rgb[1] = 35;  rgb[2] = 35;  break;
    case ATLAS_RIBBON_COVERAGE_TOPOLOGY_U:
        rgb[0] = 135; rgb[1] = 45;  rgb[2] = 210; break;
    case ATLAS_RIBBON_COVERAGE_CROSSING_U:
        rgb[0] = 255; rgb[1] = 255; rgb[2] = 255; break;
    case ATLAS_RIBBON_COVERAGE_VERTICAL_METRIC:
        rgb[0] = 30;  rgb[1] = 220; rgb[2] = 220; break;
    case ATLAS_RIBBON_COVERAGE_GEOMETRY_MISMATCH:
        rgb[0] = 80;  rgb[1] = 255; rgb[2] = 70;  break;
    case ATLAS_RIBBON_COVERAGE_VERTICAL_METRIC_FILL:
        rgb[0] = 20;  rgb[1] = 130; rgb[2] = 140; break;
    default:
        rgb[0] = 0;   rgb[1] = 0;   rgb[2] = 0;   break;
    }
}

#endif /* ATLAS_RIBBON_COVERAGE_H */
