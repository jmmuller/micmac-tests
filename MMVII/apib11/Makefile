MMV1_DIR=../..

MMV1_LIBDIR=${MMV1_DIR}/lib
MMV2_DIR=${MMV1_DIR}/MMVII
MMV2_BINDIR=${MMV2_DIR}/bin
MMV2_INCDIR=${MMV2_DIR}/include
export MODULE_NAME = MMVII

DOC_INCLUDES=../include/MMVII_PCSens.h ../include/MMVII_memory.h ../include/MMVII_AimeTieP.h ../include/MMVII_Geom3D.h ../include/MMVII_Mappings.h ../include/MMVII_MeasuresIm.h ../include/MMVII_PoseRel.h ../include/MMVII_SysSurR.h  ../include/MMVII_nums.h


# =========== Includes & Libraries
MMV1_INCLUDES=-I ${MMV1_DIR}/include -I${MMV1_DIR}
MMV2_INCLUDES=-I ${MMV2_DIR} -I ${MMV2_INCDIR} -I ${MMV2_DIR}/ExternalInclude/eigen-3.4.0
PYBIND_INCLUDES=$(shell python3 -m pybind11 --includes) 

# =========== Compiler & Flags
ifeq ($(OS),Windows_NT)
    CXX=cl.exe
    LD=link.exe
    RM=del
    EXTRA_OBJECTS=../bin/libP2007.lib
    MMVII_BIN=../bin/MMVII.exe
    CXXFLAGS=/openmp /std:c++17 /Wall /W3 /O2
else
    CXX=clang++
    LD=clang++
    RM=rm -f
    EXTRA_OBJECTS=../bin/libP2007.a
    MMVII_BIN=../bin/MMVII
    CXXFLAGS=-fopenmp -std=c++17 -Wall -Werror -O3 -fPIC
endif

export CXXFLAGS
ALL_INC=${MMV2_INCLUDES} ${MMV1_INCLUDES}  ${PYBIND_INCLUDES}

PYTHON = python3
PIP = pip3

# ===========  Use ccache if found
CCACHE:=${shell which ccache}
ifneq (,${CCACHE})
  override CXX:=${CCACHE} ${CXX}
endif

export SRCS:=${wildcard *.cpp}
SRCS_INC:=${wildcard *.h}

# ============ Targets
export XML_MICMAC_FILES=$(wildcard ../../include/XML_MicMac/*.xml)
export XML_GEN_FILES=$(wildcard ../../include/XML_GEN/*.xml)
export XML_MMVII_LOCALPARAMETERS=$(wildcard ../MMVII-LocalParameters/*.xml)
export XML_MMVII_LOCALPARAMETERS_DEFAULT=../MMVII-LocalParameters/Default/MMVII-UserOfProfile.xml
export EXTRA_OBJECTS
export MMVII_BIN


.phony: all clean clean_python install uninstall

all: wheel

clean: clean_python
	$(RM) docstrings.h
	$(RM)r xml

clean_python:
	$(RM)R build/ dist/ $(MODULE_NAME).egg-info/ __pycache__/

docstrings.h: $(DOC_INCLUDES) Makefile
	$(RM) $@
	$(RM)r xml
	(cat Doxyfile ; echo "INPUT=$(DOC_INCLUDES)") | doxygen -
	python3 makedoc.py xml/class*xml xml/struct*xml xml/namespaceMMVII.xml > $@ || $(RM) $@
	$(RM)r xml

wheel: docstrings.h ${SRCS} ${SRCS_INC} $(MMVII_BIN) $(EXTRA_OBJECTS) $(XML_MICMAC_FILES) $(XML_GEN_FILES) $(XML_MMVII_LOCALPARAMETERS)
	$(MAKE) clean_python
	CC="$$CXX" SRCS="$$SRCS" CXXFLAGS="$$CXXFLAGS" XML_MICMAC_FILES="$$XML_MICMAC_FILES" XML_GEN_FILES="$$XML_GEN_FILES" XML_MMVII_LOCALPARAMETERS="$$XML_MMVII_LOCALPARAMETERS" XML_MMVII_LOCALPARAMETERS_DEFAULT="$$XML_MMVII_LOCALPARAMETERS_DEFAULT" EXTRA_OBJECTS="$$EXTRA_OBJECTS" MMVII_BIN="$$MMVII_BIN" $(PIP) wheel . --no-deps -w dist/


install:
	if [ ! -f dist/*.whl ]; then \
		$(MAKE) wheel; \
	fi
	$(PIP) install --force-reinstall dist/*.whl

uninstall:
	$(PIP) uninstall -y $(MODULE_NAME)
