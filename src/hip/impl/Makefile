HIPCC ?= hipcc
AR ?= ar
NM ?= nm
READELF ?= readelf
CC ?= cc
ARCH ?= gfx1201
BUILD_DIR ?= build

CPPFLAGS ?=
CXXFLAGS ?= -O3 -std=c++17 -fPIC -Wall -Wextra -Wpedantic
HIPFLAGS ?= --offload-arch=$(ARCH)

HEADER := mls_project_clean.h
SOURCE := mls_project_hip.hip
OBJECT := $(BUILD_DIR)/mls_project_hip.o
LIBRARY := $(BUILD_DIR)/libmls_project_clean.a
RUNNER_SOURCE := mls_clean_runner.cpp
RUNNER_OBJECT := $(BUILD_DIR)/mls_clean_runner.o
RUNNER := $(BUILD_DIR)/mls_clean_runner

.PHONY: all object library runner check check-header check-object check-library check-runner check-symbols clean

all: library runner check

object: $(OBJECT)

library: $(LIBRARY)

runner: $(RUNNER)

$(BUILD_DIR):
	mkdir -p $@

$(OBJECT): $(SOURCE) $(HEADER) | $(BUILD_DIR)
	$(HIPCC) $(CPPFLAGS) $(CXXFLAGS) $(HIPFLAGS) -x hip -c $(SOURCE) -o $@

$(LIBRARY): $(OBJECT)
	$(AR) rcs $@ $<

$(RUNNER_OBJECT): $(RUNNER_SOURCE) $(HEADER) | $(BUILD_DIR)
	$(HIPCC) $(CPPFLAGS) $(CXXFLAGS) $(HIPFLAGS) -I. -c $(RUNNER_SOURCE) -o $@

$(RUNNER): $(RUNNER_OBJECT) $(LIBRARY)
	$(HIPCC) $(HIPFLAGS) $^ -o $@

check: check-header check-object check-library check-runner check-symbols

check-header: $(HEADER) | $(BUILD_DIR)
	$(CC) -std=c11 -Wall -Wextra -Werror -I. -x c -fsyntax-only $(HEADER)

check-object: $(OBJECT)
	test -s $(OBJECT)
	file $(OBJECT)
	$(READELF) -h $(OBJECT) >/dev/null

check-library: $(LIBRARY)
	test -s $(LIBRARY)
	$(AR) t $(LIBRARY) | grep -Fx 'mls_project_hip.o'

check-runner: $(RUNNER)
	test -x $(RUNNER)

check-symbols: $(OBJECT)
	$(NM) -g --defined-only $(OBJECT) | awk '$$3 == "MLS_project_verts" { found=1 } END { exit !found }'
	$(NM) -g --defined-only $(OBJECT) | awk '$$3 == "mls_clean_arena_create" { found=1 } END { exit !found }'
	$(NM) -g --defined-only $(OBJECT) | awk '$$3 == "mls_clean_arena_destroy" { found=1 } END { exit !found }'
	$(NM) -g --defined-only $(OBJECT) | awk '$$3 == "mls_clean_last_status" { found=1 } END { exit !found }'
	$(NM) -g --defined-only $(OBJECT) | awk '$$3 == "mls_clean_arena_get_info" { found=1 } END { exit !found }'
	$(READELF) -Ws $(OBJECT) | grep -Eq '[[:space:]]MLS_project_verts$$'

clean:
	rm -rf -- $(BUILD_DIR)
