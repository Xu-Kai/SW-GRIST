# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)

set (DIR_SRCS grist.f90                 swe_driver.f90        swe_restart_module.f90     swe_time_integration.f90
				swe_init_module.f90   swe_spatial_module.f90     swe_vars_module.f90
				swe_diagnose_module.f90   swe_inout_module.f90  swe_spatial_module_tr.f90)

link_directories(${CMAKE_BINARY_DIR}/lib)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/grist_lib/src)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/namelist)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/grid)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/hpe)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/io)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/setting)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/operators)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/drivers/phys)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/nh-dynamics)

ADD_EXECUTABLE(swe.exe ${DIR_SRCS})
SET_TARGET_PROPERTIES(swe.exe PROPERTIES LINKER_LANGUAGE Fortran)
TARGET_LINK_LIBRARIES(swe.exe grist_lib namelist_lib grid_lib hpe_lib io_lib setting_lib phys_lib utils_lib nh_dynamics_lib operators_lib)
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/run)
