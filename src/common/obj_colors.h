#ifndef OBJ_COLORS_INCLUDED
#define OBJ_COLORS_INCLUDED

/* 20-color palette for per-component OBJ vertex coloring.
 * Colors from ColorBrewer qualitative palettes (Set1 + Tab20). */
#define OBJ_NUM_COLORS 20
static const float OBJ_COLORS[20][3] = {
    {0.894f, 0.102f, 0.110f},  /* #E41A1C red        */
    {0.216f, 0.494f, 0.722f},  /* #377EB8 blue       */
    {0.302f, 0.686f, 0.290f},  /* #4DAF4A green      */
    {0.596f, 0.306f, 0.639f},  /* #984EA3 purple     */
    {1.000f, 0.498f, 0.000f},  /* #FF7F00 orange     */
    {1.000f, 1.000f, 0.200f},  /* #FFFF33 yellow     */
    {0.651f, 0.337f, 0.157f},  /* #A65628 brown      */
    {0.969f, 0.506f, 0.749f},  /* #F781BF pink       */
    {0.600f, 0.600f, 0.600f},  /* #999999 gray       */
    {0.400f, 0.761f, 0.647f},  /* #66C2A5 teal       */
    {0.988f, 0.553f, 0.384f},  /* #FC8D62 salmon     */
    {0.553f, 0.627f, 0.796f},  /* #8DA0CB lavender   */
    {0.906f, 0.541f, 0.765f},  /* #E78AC3 orchid     */
    {0.651f, 0.847f, 0.329f},  /* #A6D854 lime       */
    {1.000f, 0.851f, 0.184f},  /* #FFD92F gold       */
    {0.898f, 0.769f, 0.580f},  /* #E5C494 tan        */
    {0.702f, 0.702f, 0.702f},  /* #B3B3B3 silver     */
    {0.122f, 0.471f, 0.706f},  /* #1F78B4 dark blue  */
    {0.200f, 0.627f, 0.173f},  /* #33A02C dark green  */
    {0.890f, 0.102f, 0.110f},  /* #E31A1C crimson    */
};

#endif
