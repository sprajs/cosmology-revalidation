# Install script for directory: /home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/source/pkg/wfc3

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf32d.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf32d.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf32d.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/wf32d.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf32d.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf32d.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf32d.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf32d.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/CMakeFiles/wf32d.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ccd.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ccd.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ccd.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/wf3ccd.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ccd.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ccd.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ccd.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ccd.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/CMakeFiles/wf3ccd.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3cte.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3cte.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3cte.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/wf3cte.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3cte.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3cte.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3cte.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3cte.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/CMakeFiles/wf3cte.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ir.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ir.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ir.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/wf3ir.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ir.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ir.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ir.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3ir.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/CMakeFiles/wf3ir.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3rej.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3rej.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3rej.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/wf3rej.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3rej.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3rej.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3rej.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3rej.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/CMakeFiles/wf3rej.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3sum.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3sum.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3sum.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/wf3sum.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3sum.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3sum.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3sum.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/wf3sum.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/CMakeFiles/wf3sum.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calwf3.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calwf3.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calwf3.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/calwf3.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calwf3.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calwf3.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calwf3.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calwf3.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/CMakeFiles/calwf3.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/lib/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/wfc3/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
