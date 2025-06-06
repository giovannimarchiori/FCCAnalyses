# workaround for ROOT not properly exporting the VDT includes
find_package(Vdt)

message(STATUS "includes-------------------------- dataframe edm4hep: ${EDM4HEP_INCLUDE_DIRS}")
message(STATUS "includes-------------------------- dataframe podio  : ${podio_INCLUDE_DIR}")
message(STATUS "includes-------------------------- dataframe delphes: ${DELPHES_INCLUDE_DIR}")
message(STATUS "includes-------------------------- dataframe delphes EXt TrkCov: ${DELPHES_EXTERNALS_TKCOV_INCLUDE_DIR}")
message(STATUS "includes-------------------------- dataframe delphes EXt: ${DELPHES_EXTERNALS_INCLUDE_DIR}")

# sources and headers
file(GLOB sources src/*.cc)

file(GLOB headers RELATIVE ${CMAKE_CURRENT_LIST_DIR} FCCAnalyses/*.h)
list(FILTER headers EXCLUDE REGEX "LinkDef.h")

if(NOT WITH_DD4HEP)
  list(FILTER headers EXCLUDE REGEX "CaloNtupleizer.h")
  list(FILTER sources EXCLUDE REGEX "CaloNtupleizer.cc")
endif()

if(NOT WITH_ONNX)
  list(FILTER headers EXCLUDE REGEX "JetFlavourUtils.h")
  list(FILTER sources EXCLUDE REGEX "JetFlavourUtils.cc")
  list(FILTER headers EXCLUDE REGEX "WeaverUtils.h")
  list(FILTER sources EXCLUDE REGEX "WeaverUtils.cc")
endif()

if(NOT WITH_ACTS)
  list(FILTER headers EXCLUDE REGEX "VertexFitterActs.h")
  list(FILTER sources EXCLUDE REGEX "VertexFitterActs.cc")
  list(FILTER headers EXCLUDE REGEX "VertexFinderActs.h")
  list(FILTER sources EXCLUDE REGEX "VertexFinderActs.cc")
endif()

if(NOT WITH_PODIO_DATASOURCE)
  list(FILTER headers EXCLUDE REGEX "JetClusteringUtilsSource.h")
  list(FILTER sources EXCLUDE REGEX "JetClusteringUtilsSource.cc")
  list(FILTER headers EXCLUDE REGEX "LinkSource.h")
  list(FILTER headers EXCLUDE REGEX "ReconstructedParticleSource.h")
  list(FILTER sources EXCLUDE REGEX "ReconstructedParticleSource.cc")
  list(FILTER headers EXCLUDE REGEX "TrackSource.h")
  list(FILTER sources EXCLUDE REGEX "TrackSource.cc")
endif()

message(STATUS "FCCAnalyses headers:\n   ${headers}")
message(STATUS "FCCAnalyses sources:\n   ${sources}")

message(STATUS "CMAKE_CURRENT_SOURCE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR  ${CMAKE_INSTALL_INCLUDEDIR}")

add_library(FCCAnalyses SHARED ${sources})
target_include_directories(FCCAnalyses PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/addons>
                           $<INSTALL_INTERFACE:include>
                           ${VDT_INCLUDE_DIR}
                           ${DELPHES_INCLUDE_DIR}
                           ${DELPHES_EXTERNALS_INCLUDE_DIR}
                           ${DELPHES_EXTERNALS_TKCOV_INCLUDE_DIR}
                           )

message(STATUS " ====== DELPHES LIBRARY = " ${DELPHES_LIBRARY} )
message(STATUS " ====== DELPHES_EXTERNALS_TKCOV_INCLUDE_DIR = " ${DELPHES_EXTERNALS_TKCOV_INCLUDE_DIR} )

target_link_libraries(FCCAnalyses PUBLIC
                      ROOT::Physics
                      ROOT::MathCore
                      ROOT::ROOTVecOps
                      ROOT::ROOTDataFrame
                      podio::podio
                      podio::podioRootIO
                      podio::podioDataSource
                      EDM4HEP::edm4hep
                      EDM4HEP::edm4hepDict
                      EDM4HEP::utils
                      ${ADDONS_LIBRARIES}
                      ${DELPHES_LIBRARY}
                      gfortran # todo: why necessary?
)

if(WITH_DD4HEP)
  target_link_libraries(FCCAnalyses PUBLIC DD4hep::DDCore)
endif()

if(WITH_ACTS)
  target_link_libraries(FCCAnalyses PUBLIC ActsCore)
  target_compile_definitions(FCCAnalyses PRIVATE "ACTS_VERSION_MAJOR=${Acts_VERSION_MAJOR}")
endif()

set_target_properties(FCCAnalyses PROPERTIES
  PUBLIC_HEADER "${headers}")

ROOT_GENERATE_DICTIONARY(G__FCCAnalyses
                         ${headers}
                         MODULE FCCAnalyses
                         LINKDEF FCCAnalyses/LinkDef.h
)

install(TARGETS FCCAnalyses
    EXPORT FCCAnalysesTargets
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/FCCAnalyses"
    COMPONENT dev
)

install(FILES
  "${PROJECT_BINARY_DIR}/analyzers/dataframe/libFCCAnalyses.rootmap"
    DESTINATION "${INSTALL_LIB_DIR}" COMPONENT dev)

if (${ROOT_VERSION} GREATER 6)
    install(FILES
          "${PROJECT_BINARY_DIR}/analyzers/dataframe/libFCCAnalyses_rdict.pcm"
                DESTINATION "${INSTALL_LIB_DIR}" COMPONENT dev)
endif()
