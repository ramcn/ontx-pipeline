# Install script for directory: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/chakenal/ontx-pipeline/x86-bin")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithAsserts")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/pipeline/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/assembly_graph/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/modules/path_extend/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/modules/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/paired_info/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/stages/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/io/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils/kmer_mph/cmake_install.cmake")
  include("/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/modules/coverage_model/cmake_install.cmake")

endif()

