include(BuildModules.cmake)

# Define catalogues
make_catalogue(
  NAME bbp
  SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bbp"
  OUTPUT "CAT_BBP_SOURCES"
  MECHS CaDynamics_E2 Ca_HVA Ca_LVAst Ih Im K_Pst K_Tst Nap_Et2 NaTa_t NaTs2_t SK_E2 SKv3_1
  PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"
  CXX_FLAGS_TARGET "${ARB_CXX_FLAGS_TARGET_FULL}"
  STANDALONE FALSE
  CXXIN   FALSE
  VERBOSE ${ARB_CAT_VERBOSE})

make_catalogue(
  NAME allen
  SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/allen"
  OUTPUT "CAT_ALLEN_SOURCES"
  MECHS CaDynamics Ca_HVA Ca_LVA Ih Im Im_v2 K_P K_T Kd Kv2like Kv3_1 NaTa NaTs NaV Nap SK
  PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"
  CXX_FLAGS_TARGET "${ARB_CXX_FLAGS_TARGET_FULL}"
  STANDALONE FALSE
  CXXIN   FALSE
  VERBOSE ${ARB_CAT_VERBOSE})

make_catalogue(
  NAME default
  SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/default"
  OUTPUT "CAT_DEFAULT_SOURCES"
  MECHS exp2syn expsyn expsyn_stdp hh kamt kdrmt nax nernst pas gj
  PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"
  CXX_FLAGS_TARGET "${ARB_CXX_FLAGS_TARGET_FULL}"
  STANDALONE FALSE
  CXXIN   FALSE
  VERBOSE ${ARB_CAT_VERBOSE})

make_catalogue(
        NAME smol
        SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/smol"
        OUTPUT "CAT_SMOL_SOURCES"
        CXX_FLAGS_TARGET "${ARB_CXX_FLAGS_TARGET_FULL}"
        STANDALONE FALSE
        VERBOSE ${ARB_CAT_VERBOSE}
        CXXIN   FALSE
        MECHS ca_conc cacc cah cal h k kca kdr leak na_a na_s calpid cx36 cx36temp
        PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"
    )

# Goal make it possible to add CPP catalogues with a single make_catalogue function
make_catalogue(
        NAME EMC
        SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/EMC/EMC_Mechs"  # redirected to output path ++ the emc_cat.cpp
        OUTPUT "CAT_EMC_SOURCES"                             #    see above -
        CXX_FLAGS_TARGET "${ARB_CXX_FLAGS_TARGET_FULL}"
        STANDALONE FALSE
        VERBOSE TRUE
        CXXIN   TRUE
        MECHS ou_noise gj_linear calshadow
        PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"            #   Unused NOW
)

# Join sources
set(arbor_mechanism_sources
  ${CAT_BBP_SOURCES}
  ${CAT_ALLEN_SOURCES}
  ${CAT_EMC_SOURCES}
  ${CAT_DEFAULT_SOURCES}
  ${CAT_SMOL_SOURCES}
  PARENT_SCOPE)

if(ARB_WITH_CUDA_CLANG OR ARB_WITH_HIP_CLANG)
    set_source_files_properties(${arbor_mechanism_sources} PROPERTIES LANGUAGE CXX)
    set_source_files_properties(${EMC_mechanism_sources_custom} PROPERTIES LANGUAGE CXX)
endif()
