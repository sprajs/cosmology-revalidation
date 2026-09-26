# Install script for directory: /home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/source/pkg/acs

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acs2d.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acs2d.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acs2d.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/acs2d.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acs2d.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acs2d.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acs2d.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acs2d.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/CMakeFiles/acs2d.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsccd.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsccd.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsccd.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/acsccd.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsccd.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsccd.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsccd.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsccd.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/CMakeFiles/acsccd.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acscte.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acscte.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acscte.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/acscte.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acscte.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acscte.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acscte.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acscte.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/CMakeFiles/acscte.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/acscteforwardmodel.e")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsrej.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsrej.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsrej.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/acsrej.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsrej.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsrej.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsrej.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acsrej.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/CMakeFiles/acsrej.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acssum.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acssum.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acssum.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/acssum.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acssum.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acssum.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acssum.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/acssum.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/CMakeFiles/acssum.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calacs.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calacs.e")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calacs.e"
         RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/calacs.e")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calacs.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calacs.e")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calacs.e"
         OLD_RPATH "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/lib:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/tables:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/hstio:/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/cvos:"
         NEW_RPATH "$ORIGIN/../lib:/usr/local/lib:/usr/lib:")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/calacs.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/CMakeFiles/calacs.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/lib/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/build-fortran/pkg/acs/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
