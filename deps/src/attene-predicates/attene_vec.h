#ifndef ATTENE_VEC_H
#define ATTENE_VEC_H
// Minimal Vec shim for the vendored Attene predicates (replaces breadboard
// math/vector_math.h). Only the members the predicate code touches.
struct Vec2d { double x, y; Vec2d() : x(0), y(0) {} Vec2d(double _x, double _y) : x(_x), y(_y) {} };
struct Vec3d { double x, y, z; Vec3d() : x(0), y(0), z(0) {} Vec3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {} };
struct Vec3f { float x, y, z; Vec3f() : x(0), y(0), z(0) {} Vec3f(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {} };
#endif
