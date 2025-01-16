# This file implements the SAPPHIREPP_SETUP_TARGET_MHD macro.
#
# Usage:
#       sapphirepp_setup_target_mhd(target)
#       sapphirepp_setup_target_mhd(target config_header_dir)
#
# This appends the Sapphire++ MHD module 
# and deal.II target to the link interface 
# of the specified target,
# which in turn ensures that necessary include directories, linker flags,
# compile flags and compile definitions are set.
#
# If no "config_header_dir" argument is specified after the target, 
# the current CMAKE_CURRENT_SOURCE_DIR is used instead.
# This parameter should give the path to the "config.h" file directory,
# in order to include the "config.h" file the target.
#

macro(sapphirepp_setup_target_mhd _target)
  if(${ARGC} GREATER_EQUAL 2)
    set(_config_header_dir ${ARGV1})
  else()
    set(_config_header_dir ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  message(STATUS "Set up ${_target} for MHD module"
    " using ${_config_header_dir}/config.h"
  )

  target_sources(${_target} PRIVATE ${MHD_SOURCES})

  deal_ii_setup_target(${_target})
  target_include_directories(${_target} PUBLIC
    ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
  )
  target_include_directories(${_target} PUBLIC ${_config_header_dir})
  target_link_libraries(${_target} UtilsLib MHDLib)
endmacro()
