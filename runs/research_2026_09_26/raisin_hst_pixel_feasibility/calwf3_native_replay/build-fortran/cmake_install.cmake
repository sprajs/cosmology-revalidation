# Install script for directory: /home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/source

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
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
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/ctegen2/cmake_install.cmake")
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos/cmake_install.cmake")
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio/cmake_install.cmake")
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib/cmake_install.cmake")
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables/cmake_install.cmake")
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
