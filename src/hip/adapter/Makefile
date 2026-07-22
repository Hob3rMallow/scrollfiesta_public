CXX ?= g++
HIPCC ?= /usr/bin/hipcc
OBJCOPY ?= objcopy
AR ?= ar
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
HIP_ARCH ?= gfx1201

IMPL_DIR := ../impl
BUILD_DIR := build
CLEAN_OBJ := $(IMPL_DIR)/build/mls_project_hip.o
RENAMED_OBJ := $(BUILD_DIR)/mls_project_hip_renamed.o
RENAMED_LIB := $(BUILD_DIR)/libmls_project_clean_renamed.a
ADAPTER_OBJ := $(BUILD_DIR)/mls_scrollfiesta_adapter.o
SMOKE_OBJ := $(BUILD_DIR)/adapter_smoke.o
SMOKE_BIN := $(BUILD_DIR)/adapter_smoke

.PHONY: all clean smoke
all: $(SMOKE_BIN)

$(BUILD_DIR):
	mkdir -p $@

$(RENAMED_OBJ): $(CLEAN_OBJ) | $(BUILD_DIR)
	$(OBJCOPY) --redefine-sym MLS_project_verts=mls_clean_project_verts $< $@

$(RENAMED_LIB): $(RENAMED_OBJ)
	$(AR) rcs $@ $<

$(ADAPTER_OBJ): mls_scrollfiesta_adapter.cpp mls_scrollfiesta_adapter.h $(IMPL_DIR)/mls_project_clean.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(IMPL_DIR) -c $< -o $@

$(SMOKE_OBJ): adapter_smoke.cpp mls_scrollfiesta_adapter.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SMOKE_BIN): $(ADAPTER_OBJ) $(SMOKE_OBJ) $(RENAMED_LIB)
	HSA_ENABLE_DXG_DETECTION=1 ROCPROFILER_REGISTER_ENABLED=0 $(HIPCC) --offload-arch=$(HIP_ARCH) $^ -o $@

smoke: $(SMOKE_BIN)
	HSA_ENABLE_DXG_DETECTION=1 ROCPROFILER_REGISTER_ENABLED=0 ./$(SMOKE_BIN)

clean:
	rm -rf $(BUILD_DIR)
