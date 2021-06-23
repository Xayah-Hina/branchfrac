
add_library(Utils
        FileSystem.h
        Hashmap.h
        IndexedFaceMesh.cpp
        IndexedFaceMesh.h
        IndexedTetMesh.cpp
        IndexedTetMesh.h
        Logger.h
        OBJLoader.h
        SceneLoader.cpp
        SceneLoader.h
        StringTools.h
        SystemInfo.h
        TetGenLoader.cpp
        TetGenLoader.h
        Timing.h
        Version.h
        VolumeIntegration.cpp
        VolumeIntegration.h

        pymesh/Exception.h
        pymesh/MshLoader.h
        pymesh/MshLoader.cpp
        pymesh/MshSaver.h
        pymesh/MshSaver.cpp

        CMakeLists.txt
        )

############################################################
# GenericParameters
############################################################
include_directories(${GenericParameters_INCLUDE_DIR})
if (TARGET Ext_GenericParameters)
    add_dependencies(Utils Ext_GenericParameters)
endif ()

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

install(TARGETS Utils
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(DIRECTORY .
        DESTINATION include/Utils
        FILES_MATCHING PATTERN "*.h")




