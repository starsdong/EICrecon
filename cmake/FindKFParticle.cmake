# * Locate KFParticle library (https://github.com/alisw/KFParticle) Defines:
#
# KFPARTICLE_FOUND KFPARTICLE_INCLUDE_DIR KFPARTICLE_INCLUDE_DIRS (not cached)
# KFPARTICLE_LIBRARY KFPARTICLE_LIBRARIES (not cached) and the imported target
# KFParticle::KFParticle
#
# KFParticle headers change behaviour with the field-mode macro the library was
# compiled with. The default (non-FIXTARGET) build uses HomogeneousField, so the
# imported target propagates that definition to consumers.

find_library(
  KFPARTICLE_LIBRARY
  NAMES KFParticle
  HINTS $ENV{KFPARTICLE_ROOT}/lib ${KFPARTICLE_ROOT_DIR}/lib ${KFParticle_ROOT}/lib)

set(_kfparticle_include_hints)
if(KFPARTICLE_LIBRARY)
  # cmake-lint: disable=E1126
  file(REAL_PATH "${KFPARTICLE_LIBRARY}" _kfparticle_real_lib)
  get_filename_component(_kfparticle_lib_dir "${_kfparticle_real_lib}"
                         DIRECTORY)
  get_filename_component(_kfparticle_prefix "${_kfparticle_lib_dir}" DIRECTORY)
  list(APPEND _kfparticle_include_hints "${_kfparticle_prefix}/include")
  unset(_kfparticle_real_lib)
  unset(_kfparticle_lib_dir)
  unset(_kfparticle_prefix)
endif()

find_path(
  KFPARTICLE_INCLUDE_DIR KFParticle/KFParticle.h
  HINTS ${_kfparticle_include_hints} $ENV{KFPARTICLE_ROOT}/include
        ${KFPARTICLE_ROOT_DIR}/include ${KFParticle_ROOT}/include
  NO_DEFAULT_PATH)

unset(_kfparticle_include_hints)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KFParticle DEFAULT_MSG KFPARTICLE_INCLUDE_DIR
                                  KFPARTICLE_LIBRARY)

mark_as_advanced(KFPARTICLE_FOUND KFPARTICLE_INCLUDE_DIR KFPARTICLE_LIBRARY)

# Headers are installed flat under include/KFParticle and include each other
# without the directory prefix
set(KFPARTICLE_INCLUDE_DIRS ${KFPARTICLE_INCLUDE_DIR}
                            ${KFPARTICLE_INCLUDE_DIR}/KFParticle)
set(KFPARTICLE_LIBRARIES ${KFPARTICLE_LIBRARY})

if(KFPARTICLE_FOUND AND NOT TARGET KFParticle::KFParticle)
  add_library(KFParticle::KFParticle SHARED IMPORTED)
  set_target_properties(
    KFParticle::KFParticle
    PROPERTIES IMPORTED_LOCATION "${KFPARTICLE_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${KFPARTICLE_INCLUDE_DIRS}"
               INTERFACE_COMPILE_DEFINITIONS "HomogeneousField")
endif()
