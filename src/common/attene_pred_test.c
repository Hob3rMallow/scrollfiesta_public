/* attene_pred_test — standalone unit test for the exact Voronoi side predicates.
 * Exercises explicit / edge-bisector / facet-bisector2 constructions and the
 * side predicate against configurations with known answers. Returns 0 on pass. */
#include "attene_predicates_wrap.h"

int main(void) {
    return AttenePred_selftest() == 0 ? 0 : 1;
}
