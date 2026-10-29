#ifndef JSON_READ_INCLUDED
#define JSON_READ_INCLUDED

#include <stddef.h>
#include "arena.h"

/* Minimal, arena-backed JSON reader.  The pipeline only needs to READ a handful
 * of manifests (report.json bbox, material_tracks/manifest.json offsets,
 * per-track geometry_report.json), so this is a compact recursive-descent parser
 * (objects / arrays / strings / numbers / true / false / null) that materializes
 * a read-only tree in the arena.  It is not a full-fidelity emitter and does not
 * preserve number formatting; strings are decoded (\" \\ \/ \b \f \n \r \t and
 * \uXXXX in the BMP).  There is no free -- dispose the arena. */

typedef enum {
    JSON_NULL = 0,
    JSON_BOOL,
    JSON_NUMBER,
    JSON_STRING,
    JSON_ARRAY,
    JSON_OBJECT
} JsonType;

typedef struct JsonValue JsonValue;

/* Parse a NUL-terminated JSON string.  On error returns NULL and, if err is
 * non-NULL, points *err at a static diagnostic.  The tree lives in `arena`. */
const JsonValue *Json_parse(Arena_T arena, const char *text, const char **err);

/* Read a whole file into the arena and parse it.  NULL on IO or parse error. */
const JsonValue *Json_parse_file(Arena_T arena, const char *path,
                                 const char **err);

JsonType Json_type(const JsonValue *v);

/* Array access. */
size_t           Json_array_len(const JsonValue *v);
const JsonValue *Json_array_get(const JsonValue *v, size_t i);

/* Object access by key, and by ordinal (for iterating unknown keys). */
const JsonValue *Json_object_get(const JsonValue *v, const char *key);
size_t           Json_object_len(const JsonValue *v);
const char      *Json_object_key(const JsonValue *v, size_t i);
const JsonValue *Json_object_value_at(const JsonValue *v, size_t i);

/* Scalar reads.  On type mismatch a number/int read returns `def`, a string read
 * returns NULL, a bool read returns `def`. */
double      Json_as_double(const JsonValue *v, double def);
long        Json_as_long(const JsonValue *v, long def);
int         Json_as_bool(const JsonValue *v, int def);
const char *Json_as_string(const JsonValue *v);

/* Convenience: object member as double/long, or `def` if absent/mismatched. */
double Json_member_double(const JsonValue *obj, const char *key, double def);
long   Json_member_long(const JsonValue *obj, const char *key, long def);

int Json_read_selftest(void);

#endif
