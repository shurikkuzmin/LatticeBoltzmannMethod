CMAKE_MINIMUM_REQUIRED (VERSION 2.4)
PROJECT (FreeEnergy)
INCLUDE_DIRECTORIES("/opt/local/include/vtk")
INClUDE_DIRECTORIES("/usr/local/include/vtk-5.6")
INCLUDE_DIRECTORIES("/usr/include/vtk-5.2")
INCLUDE_DIRECTORIES("/usr/local/vtk/include/vtk-5.2")
INCLUDE_DIRECTORIES("/global/scratch/software/vtk/vtk-5.4.2/include/vtk-5.4")
INCLUDE_DIRECTORIES("/home/shurik/vtk/include/vtk-5.6")
LINK_DIRECTORIES("/opt/local/lib/vtk")
LINK_DIRECTORIES("/usr/local/lib/vtk-5.6")
LINK_DIRECTORIES("/usr/lib/vtk5-2")
LINK_DIRECTORIES("/usr/local/vtk/lib/vtk-5.2")
LINK_DIRECTORIES("/global/scratch/software/vtk/vtk-5.4.2/lib/vtk-5.4")
SET(CMAKE_CXX_COMPILER mpiicpc)
SET(CMAKE_CXX_FLAGS "-Wno-deprecated")
OPTION(MYPROJECT_DEBUG
  "Build the project using debugging code and stop for gdb"
  OFF)

OPTION(PERFORMANCE_DEBUG
  "Build the project using the code for profiling"
  OFF)

OPTION(WHOLE_PROJECT_DEBUG
  "Build the project using debugging code"
  OFF)

OPTION(HYDRODYNAMICS
    "Build the project using simple hydrodynamics"
  OFF
)

OPTION(PERFORMANCE
    "Build the project with the performance optimization"
  OFF
)
IF(WHOLE_PROJECT_DEBUG)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wno-deprecated")
ENDIF(WHOLE_PROJECT_DEBUG)

IF(PERFORMANCE_DEBUG)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
ENDIF(PERFORMANCE_DEBUG)


IF(MYPROJECT_DEBUG)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DDEBUG -Wno-deprecated")
ENDIF(MYPROJECT_DEBUG)

IF(HYDRODYNAMICS)
 SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHYDRO ")
ENDIF(HYDRODYNAMICS)

IF(PERFORMANCE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")
ENDIF(PERFORMANCE)


SET(LIBS vtkCommon vtkIO vtkFiltering)
ADD_EXECUTABLE(main.out main.cpp mpi_wrapper.cpp mpi_singleton.cpp)
TARGET_LINK_LIBRARIES(main.out ${LIBS})
