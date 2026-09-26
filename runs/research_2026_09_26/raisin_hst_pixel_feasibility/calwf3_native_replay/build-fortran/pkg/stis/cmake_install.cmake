# Install script for directory: /home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/source/pkg/stis

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs0.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs0.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs0.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs0.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs0.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs0.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs0.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs0.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs0.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs1.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs1.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs1.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs1.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs1.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs1.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs1.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs1.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs1.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs2.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs2.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs2.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs2.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs2.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs2.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs2.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs2.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs2.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs4.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs4.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs4.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs4.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs4.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs4.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs4.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs4.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs4.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs6.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs6.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs6.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs6.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs6.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs6.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs6.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs6.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs6.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs7.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs7.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs7.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs7.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs7.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs7.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs7.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs7.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs7.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs8.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs8.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs8.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs8.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs8.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs8.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs8.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs8.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs8.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs11.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs11.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs11.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs11.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs11.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs11.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs11.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs11.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs11.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs12.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs12.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs12.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/cs12.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs12.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs12.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs12.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/cs12.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/CMakeFiles/cs12.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/lib/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/stis/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
