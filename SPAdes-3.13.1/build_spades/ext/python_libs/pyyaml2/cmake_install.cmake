# Install script for directory: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2

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

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "runtime")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/spades/pyyaml2" TYPE FILE FILES
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/__init__.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/composer.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/constructor.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/cyaml.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/dumper.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/emitter.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/error.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/events.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/loader.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/nodes.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/parser.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/reader.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/representer.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/resolver.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/scanner.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/serializer.py"
    "/home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/python_libs/pyyaml2/tokens.py"
    )
endif()

