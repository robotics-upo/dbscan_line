cmake_minimum_required(VERSION 2.8.3)
project(dbscan_line_cluster)

## Compile as C++11, supported in ROS Kinetic and newer
set (CMAKE_CXX_STANDARD 14)

# find_package(Armadillo REQUIRED)
# find_package(Eigen REQUIRED)
set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3")


include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
)

add_library(dbscan_line_cluster SHARED src/dbscan.cpp src/dbscan_lines.cpp
          src/line.cpp src/detected_line.cpp
)

install( TARGETS dbscan_line_cluster DESTINATION ${INSTALL_LOCATION}lib/ PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE )

install(DIRECTORY include/dbscan_line DESTINATION ${INSTALL_LOCATION}include)
