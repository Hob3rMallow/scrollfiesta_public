/* An alias must be connected back to its original source index by validated
 * cut obligations. Matching XYZ/IDs alone cannot excuse an arbitrary split. */
static int apr_cut_preserved(Arena_T arena, const AsmChart *old, const AsmRun *after, size_t chart)
{
    const AsmChart *next = after->charts+chart;
    if (AsmChart_registered(old) || !AsmChart_registered(next) || next->nv <= old->nv ||
        old->id != next->id || old->comp != next->comp || old->parent != next->parent || old->area3d != next->area3d ||
        !acut_source_same(old,next) || memcmp(old->xyz,next->xyz,3*old->nv*sizeof(float)) ||
        memcmp(old->vid,next->vid,old->nv*sizeof(int32_t))) return 0;
    uint8_t *member = ARENA_CALLOC(arena,after->n_charts,1); member[chart] = 1;
    AsmField *field = AsmField_new_local(arena,after,member,NULL);
    if (!field || !field->material_vertex) return 0;
    size_t offset = field->vertex_offset[chart];
    for (size_t v = 0; v < next->nv; v++) {
        size_t root = (size_t)field->material_vertex[offset+v]-offset;
        if (root >= old->nv || old->vid[root] != next->vid[v] ||
            memcmp(old->xyz+3*root,next->xyz+3*v,3*sizeof(float))) return 0;
    }
    for (size_t k = 0; k < 3*old->nf; k++)
        if (field->material_vertex[offset+(size_t)next->faces[k]] != (int32_t)(offset+(size_t)old->faces[k])) return 0;
    for (size_t i = 0; i < field->n_seams; i++) if (field->seams[i].original_cut) {
        AsmSeamMeasure seam; AsmField_seam(field,i,field->uv,&seam); if (!seam.pass) return 0;
    }
    AsmMetricFace *metric = ARENA_ALLOC(arena,next->nf*sizeof(AsmMetricFace));
    AsmMetricStats measured = {0};
    if (AsmMetric_prepare(next,metric)) return 0;
    AsmMetric_measure_encoded_uv64(next,metric,field->uv+2*offset,&measured);
    return measured.area > 0 && !measured.invalid && AsmMetric_within25(measured.minimum,measured.maximum) &&
        measured.within10 >= .95*measured.area;
}
