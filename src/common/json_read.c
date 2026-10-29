/* ============================================================================
 * json_read.c -- see json_read.h.  Compact recursive-descent JSON reader that
 * builds a read-only tree in an arena.  Containers are materialized into
 * contiguous arena arrays (O(1) index / O(members) key lookup); the only
 * transient allocation is a malloc'd growable vector used while a single
 * container is being parsed, released before the parser returns from it.
 * ==========================================================================*/

#include "json_read.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct JsonValue {
    JsonType type;
    union {
        int    boolean;      /* JSON_BOOL */
        double number;       /* JSON_NUMBER */
        const char *string;  /* JSON_STRING (NUL-terminated, arena) */
        struct {             /* JSON_ARRAY */
            const JsonValue **items;
            size_t count;
        } arr;
        struct {             /* JSON_OBJECT */
            const char      **keys;
            const JsonValue **vals;
            size_t count;
        } obj;
    } u;
};

typedef struct {
    const char *p;      /* cursor */
    const char *end;    /* one past last char */
    Arena_T     arena;
    const char *err;    /* static message on failure */
} Parser;

static const JsonValue *parse_value(Parser *ps);

/* ---- transient growable pointer vector (malloc, parse-time only) ---------- */
typedef struct {
    const void **data;
    size_t count, cap;
} Vec;

static int vec_push(Vec *v, const void *item) {
    if (v->count == v->cap) {
        size_t ncap = v->cap ? v->cap * 2 : 8;
        const void **nd = (const void **)realloc(v->data, ncap * sizeof(*nd));
        if (nd == NULL) return -1;
        v->data = nd;
        v->cap = ncap;
    }
    v->data[v->count++] = item;
    return 0;
}

static void skip_ws(Parser *ps) {
    while (ps->p < ps->end) {
        char c = *ps->p;
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            ps->p++;
        } else {
            break;
        }
    }
}

static JsonValue *new_node(Parser *ps, JsonType t) {
    JsonValue *v = ARENA_ALLOC(ps->arena, sizeof *v);
    memset(v, 0, sizeof *v);
    v->type = t;
    return v;
}

/* Encode a BMP codepoint as UTF-8 into dst; returns bytes written (1..3). */
static size_t utf8_encode(unsigned cp, char *dst) {
    if (cp < 0x80u) {
        dst[0] = (char)cp;
        return 1;
    }
    if (cp < 0x800u) {
        dst[0] = (char)(0xC0u | (cp >> 6));
        dst[1] = (char)(0x80u | (cp & 0x3Fu));
        return 2;
    }
    dst[0] = (char)(0xE0u | (cp >> 12));
    dst[1] = (char)(0x80u | ((cp >> 6) & 0x3Fu));
    dst[2] = (char)(0x80u | (cp & 0x3Fu));
    return 3;
}

static int hex4(const char *s, unsigned *out) {
    unsigned v = 0;
    for (int i = 0; i < 4; i++) {
        char c = s[i];
        unsigned d = 0;
        if (c >= '0' && c <= '9') {
            d = (unsigned)(c - '0');
        } else if (c >= 'a' && c <= 'f') {
            d = (unsigned)(c - 'a') + 10u;
        } else if (c >= 'A' && c <= 'F') {
            d = (unsigned)(c - 'A') + 10u;
        } else {
            return -1;
        }
        v = (v << 4) | d;
    }
    *out = v;
    return 0;
}

/* Parse a string starting at the opening quote; returns arena NUL-terminated
 * text (never NULL on success, "" for empty).  Advances past the close quote. */
static const char *parse_string_raw(Parser *ps) {
    if (ps->p >= ps->end || *ps->p != '"') {
        ps->err = "expected string";
        return NULL;
    }
    ps->p++;  /* opening quote */
    const char *start = ps->p;
    /* Decoded length <= raw length (escapes never expand past their source). */
    size_t raw = 0;
    const char *scan = start;
    while (scan < ps->end && *scan != '"') {
        if (*scan == '\\') scan++;  /* skip the escaped char in the length scan */
        if (scan < ps->end) scan++;
        raw++;
    }
    char *out = ARENA_ALLOC(ps->arena, raw + 1);
    size_t n = 0;
    while (ps->p < ps->end && *ps->p != '"') {
        char c = *ps->p;
        if (c != '\\') {
            out[n++] = c;
            ps->p++;
            continue;
        }
        ps->p++;  /* backslash */
        if (ps->p >= ps->end) { ps->err = "truncated escape"; return NULL; }
        char e = *ps->p++;
        switch (e) {
            case '"':  out[n++] = '"';  break;
            case '\\': out[n++] = '\\'; break;
            case '/':  out[n++] = '/';  break;
            case 'b':  out[n++] = '\b'; break;
            case 'f':  out[n++] = '\f'; break;
            case 'n':  out[n++] = '\n'; break;
            case 'r':  out[n++] = '\r'; break;
            case 't':  out[n++] = '\t'; break;
            case 'u': {
                unsigned cp = 0;
                if (ps->end - ps->p < 4 || hex4(ps->p, &cp) != 0) {
                    ps->err = "bad \\u escape";
                    return NULL;
                }
                ps->p += 4;
                if (cp >= 0xD800u && cp <= 0xDFFFu) cp = 0xFFFDu; /* no surrogates */
                n += utf8_encode(cp, out + n);
                break;
            }
            default:
                ps->err = "bad escape";
                return NULL;
        }
    }
    if (ps->p >= ps->end) { ps->err = "unterminated string"; return NULL; }
    ps->p++;  /* closing quote */
    out[n] = '\0';
    return out;
}

static const JsonValue *parse_number(Parser *ps) {
    char *endp = NULL;
    /* strtod needs a NUL-terminated region; the source buffer is NUL-terminated
     * (Json_parse requires it), so scanning past `end` cannot run away. */
    double d = strtod(ps->p, &endp);
    if (endp == ps->p) {
        ps->err = "bad number";
        return NULL;
    }
    ps->p = endp;
    JsonValue *v = new_node(ps, JSON_NUMBER);
    v->u.number = d;
    return v;
}

static const JsonValue *parse_array(Parser *ps) {
    ps->p++;  /* '[' */
    Vec vec = {0};
    skip_ws(ps);
    if (ps->p < ps->end && *ps->p == ']') {
        ps->p++;
        JsonValue *v = new_node(ps, JSON_ARRAY);
        return v;  /* empty */
    }
    for (;;) {
        const JsonValue *item = parse_value(ps);
        if (item == NULL || vec_push(&vec, item) != 0) {
            if (item != NULL) ps->err = "out of memory";
            free(vec.data);
            return NULL;
        }
        skip_ws(ps);
        if (ps->p >= ps->end) { ps->err = "unterminated array"; free(vec.data); return NULL; }
        char c = *ps->p++;
        if (c == ',') { skip_ws(ps); continue; }
        if (c == ']') break;
        ps->err = "expected ',' or ']'";
        free(vec.data);
        return NULL;
    }
    JsonValue *v = new_node(ps, JSON_ARRAY);
    v->u.arr.count = vec.count;
    v->u.arr.items = ARENA_ALLOC(ps->arena, vec.count * sizeof(*v->u.arr.items));
    memcpy(v->u.arr.items, vec.data, vec.count * sizeof(*v->u.arr.items));
    free(vec.data);
    return v;
}

static const JsonValue *parse_object(Parser *ps) {
    ps->p++;  /* '{' */
    Vec keys = {0}, vals = {0};
    skip_ws(ps);
    if (ps->p < ps->end && *ps->p == '}') {
        ps->p++;
        return new_node(ps, JSON_OBJECT);  /* empty */
    }
    for (;;) {
        skip_ws(ps);
        const char *key = parse_string_raw(ps);
        if (key == NULL) { free(keys.data); free(vals.data); return NULL; }
        skip_ws(ps);
        if (ps->p >= ps->end || *ps->p != ':') {
            ps->err = "expected ':'";
            free(keys.data); free(vals.data);
            return NULL;
        }
        ps->p++;
        const JsonValue *val = parse_value(ps);
        if (val == NULL ||
            vec_push(&keys, key) != 0 || vec_push(&vals, val) != 0) {
            if (val != NULL) ps->err = "out of memory";
            free(keys.data); free(vals.data);
            return NULL;
        }
        skip_ws(ps);
        if (ps->p >= ps->end) { ps->err = "unterminated object"; free(keys.data); free(vals.data); return NULL; }
        char c = *ps->p++;
        if (c == ',') continue;
        if (c == '}') break;
        ps->err = "expected ',' or '}'";
        free(keys.data); free(vals.data);
        return NULL;
    }
    JsonValue *v = new_node(ps, JSON_OBJECT);
    v->u.obj.count = keys.count;
    v->u.obj.keys = ARENA_ALLOC(ps->arena, keys.count * sizeof(*v->u.obj.keys));
    v->u.obj.vals = ARENA_ALLOC(ps->arena, vals.count * sizeof(*v->u.obj.vals));
    memcpy(v->u.obj.keys, keys.data, keys.count * sizeof(*v->u.obj.keys));
    memcpy(v->u.obj.vals, vals.data, vals.count * sizeof(*v->u.obj.vals));
    free(keys.data);
    free(vals.data);
    return v;
}

static int match_literal(Parser *ps, const char *lit) {
    size_t n = strlen(lit);
    if ((size_t)(ps->end - ps->p) < n || memcmp(ps->p, lit, n) != 0) return -1;
    ps->p += n;
    return 0;
}

static const JsonValue *parse_value(Parser *ps) {
    skip_ws(ps);
    if (ps->p >= ps->end) { ps->err = "unexpected end of input"; return NULL; }
    char c = *ps->p;
    switch (c) {
        case '{': return parse_object(ps);
        case '[': return parse_array(ps);
        case '"': {
            const char *s = parse_string_raw(ps);
            if (s == NULL) return NULL;
            JsonValue *v = new_node(ps, JSON_STRING);
            v->u.string = s;
            return v;
        }
        case 't': {
            if (match_literal(ps, "true") != 0) { ps->err = "bad literal"; return NULL; }
            JsonValue *v = new_node(ps, JSON_BOOL); v->u.boolean = 1; return v;
        }
        case 'f': {
            if (match_literal(ps, "false") != 0) { ps->err = "bad literal"; return NULL; }
            JsonValue *v = new_node(ps, JSON_BOOL); v->u.boolean = 0; return v;
        }
        case 'n': {
            if (match_literal(ps, "null") != 0) { ps->err = "bad literal"; return NULL; }
            return new_node(ps, JSON_NULL);
        }
        default:
            if (c == '-' || (c >= '0' && c <= '9')) return parse_number(ps);
            ps->err = "unexpected character";
            return NULL;
    }
}

const JsonValue *Json_parse(Arena_T arena, const char *text, const char **err) {
    if (err != NULL) *err = NULL;
    if (text == NULL) { if (err) *err = "null input"; return NULL; }
    Parser ps;
    ps.p = text;
    ps.end = text + strlen(text);
    ps.arena = arena;
    ps.err = NULL;
    const JsonValue *v = parse_value(&ps);
    if (v == NULL) { if (err) *err = ps.err ? ps.err : "parse error"; return NULL; }
    skip_ws(&ps);
    if (ps.p != ps.end) { if (err) *err = "trailing data after value"; return NULL; }
    return v;
}

const JsonValue *Json_parse_file(Arena_T arena, const char *path,
                                 const char **err) {
    if (err != NULL) *err = NULL;
    FILE *f = fopen(path, "rb");
    if (f == NULL) { if (err) *err = "cannot open file"; return NULL; }
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); if (err) *err = "seek failed"; return NULL; }
    long sz = ftell(f);
    if (sz < 0) { fclose(f); if (err) *err = "tell failed"; return NULL; }
    rewind(f);
    char *buf = ARENA_ALLOC(arena, (size_t)sz + 1);
    size_t rd = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    buf[rd] = '\0';
    return Json_parse(arena, buf, err);
}

JsonType Json_type(const JsonValue *v) {
    return v ? v->type : JSON_NULL;
}

size_t Json_array_len(const JsonValue *v) {
    return (v && v->type == JSON_ARRAY) ? v->u.arr.count : 0;
}

const JsonValue *Json_array_get(const JsonValue *v, size_t i) {
    if (!v || v->type != JSON_ARRAY || i >= v->u.arr.count) return NULL;
    return v->u.arr.items[i];
}

const JsonValue *Json_object_get(const JsonValue *v, const char *key) {
    if (!v || v->type != JSON_OBJECT || key == NULL) return NULL;
    for (size_t i = 0; i < v->u.obj.count; i++) {
        if (strcmp(v->u.obj.keys[i], key) == 0) return v->u.obj.vals[i];
    }
    return NULL;
}

size_t Json_object_len(const JsonValue *v) {
    return (v && v->type == JSON_OBJECT) ? v->u.obj.count : 0;
}

const char *Json_object_key(const JsonValue *v, size_t i) {
    if (!v || v->type != JSON_OBJECT || i >= v->u.obj.count) return NULL;
    return v->u.obj.keys[i];
}

const JsonValue *Json_object_value_at(const JsonValue *v, size_t i) {
    if (!v || v->type != JSON_OBJECT || i >= v->u.obj.count) return NULL;
    return v->u.obj.vals[i];
}

double Json_as_double(const JsonValue *v, double def) {
    return (v && v->type == JSON_NUMBER) ? v->u.number : def;
}

long Json_as_long(const JsonValue *v, long def) {
    if (!v || v->type != JSON_NUMBER) return def;
    double d = v->u.number;
    /* round-to-nearest, matching the pipeline's round() on placement scalars */
    return (long)(d < 0.0 ? d - 0.5 : d + 0.5);
}

int Json_as_bool(const JsonValue *v, int def) {
    return (v && v->type == JSON_BOOL) ? v->u.boolean : def;
}

const char *Json_as_string(const JsonValue *v) {
    return (v && v->type == JSON_STRING) ? v->u.string : NULL;
}

double Json_member_double(const JsonValue *obj, const char *key, double def) {
    return Json_as_double(Json_object_get(obj, key), def);
}

long Json_member_long(const JsonValue *obj, const char *key, long def) {
    return Json_as_long(Json_object_get(obj, key), def);
}

/* ============================================================================
 * Self-test
 * ==========================================================================*/
#define CK(c, m) do { if (!(c)) { fprintf(stderr, "  FAIL: %s\n", (m)); fails++; } \
                      else fprintf(stderr, "  ok: %s\n", (m)); } while (0)

int Json_read_selftest(void) {
    int fails = 0;
    Arena_T arena = Arena_new();
    fprintf(stderr, "[selftest] json_read\n");

    /* t1: scalars, nesting, arrays, escapes */
    {
        const char *txt =
            "{ \"bbox_zyx\": [4480, 3328, 2816],"
            "  \"chunk_size\": 128,"
            "  \"name\": \"z04480_y03328\","
            "  \"flag\": true, \"missing\": null,"
            "  \"neg\": -12.5, \"esc\": \"a\\\"b\\n\\u0041\" }";
        const char *err = NULL;
        const JsonValue *root = Json_parse(arena, txt, &err);
        CK(root != NULL && err == NULL, "t1 parses");
        CK(Json_type(root) == JSON_OBJECT, "t1 root is object");
        const JsonValue *bbox = Json_object_get(root, "bbox_zyx");
        CK(Json_array_len(bbox) == 3, "t1 bbox len 3");
        CK(Json_as_long(Json_array_get(bbox, 0), -1) == 4480, "t1 bbox[0]=4480");
        CK(Json_member_long(root, "chunk_size", -1) == 128, "t1 chunk_size=128");
        CK(Json_as_bool(Json_object_get(root, "flag"), 0) == 1, "t1 flag=true");
        CK(Json_type(Json_object_get(root, "missing")) == JSON_NULL, "t1 null");
        CK(Json_member_double(root, "neg", 0.0) == -12.5, "t1 neg=-12.5");
        const char *esc = Json_as_string(Json_object_get(root, "esc"));
        CK(esc != NULL && strcmp(esc, "a\"b\nA") == 0, "t1 escapes decode");
        CK(Json_object_get(root, "absent") == NULL, "t1 absent key -> NULL");
    }

    /* t2: array of objects (the manifest tracks shape) */
    {
        const char *txt =
            "{\"tracks\":["
            " {\"global_track\":7,\"members\":[{\"global_u_offset\":-64.0}]},"
            " {\"global_track\":9,\"members\":[{\"global_u_offset\":128.5}]}"
            "]}";
        const JsonValue *root = Json_parse(arena, txt, NULL);
        const JsonValue *tr = Json_object_get(root, "tracks");
        CK(Json_array_len(tr) == 2, "t2 two tracks");
        const JsonValue *t1 = Json_array_get(tr, 1);
        CK(Json_member_long(t1, "global_track", -1) == 9, "t2 track1 id=9");
        const JsonValue *m0 = Json_array_get(Json_object_get(t1, "members"), 0);
        CK(Json_member_double(m0, "global_u_offset", 0.0) == 128.5,
           "t2 track1 member0 u_offset=128.5");
    }

    /* t3: object iteration by ordinal (registration charts shape) */
    {
        const char *txt = "{\"charts\":{\"7\":{\"super_component\":0},"
                          "\"9\":{\"super_component\":0},\"20\":{\"super_component\":1}}}";
        const JsonValue *charts = Json_object_get(Json_parse(arena, txt, NULL), "charts");
        CK(Json_object_len(charts) == 3, "t3 three charts");
        long sum_super = 0;
        for (size_t i = 0; i < Json_object_len(charts); i++) {
            const char *k = Json_object_key(charts, i);
            const JsonValue *val = Json_object_value_at(charts, i);
            (void)k;
            sum_super += Json_member_long(val, "super_component", -1);
        }
        CK(sum_super == 1, "t3 super_component sum = 1");
    }

    /* t4: malformed inputs report an error, not a crash */
    {
        const char *err = NULL;
        CK(Json_parse(arena, "{\"a\":}", &err) == NULL && err != NULL, "t4a bad value");
        CK(Json_parse(arena, "[1,2", &err) == NULL, "t4b unterminated array");
        CK(Json_parse(arena, "", &err) == NULL, "t4c empty input");
        CK(Json_parse(arena, "{} junk", &err) == NULL, "t4d trailing data");
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] %s (%d failures)\n", fails == 0 ? "ALL PASS" : "FAILURES", fails);
    return fails;
}
