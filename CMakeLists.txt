##################
# Initialisation #
##################
# CMake requirements
# 3.1 is required for CMAKE_CXX_STANDARD
cmake_minimum_required(VERSION 3.1)

# set project name
cmake_policy(SET CMP0048 NEW)
project(ngshmmalign VERSION 0.1 LANGUAGES CXX)
add_definitions(-DPACKAGE_NAME="${PROJECT_NAME}" -DPACKAGE_VERSION="${PROJECT_VERSION}")

# on macOS, MACOSX_RPATH is enabled by default on more recent versions
# of CMake. Disable this behaviour, and let user enable it if need be.
cmake_policy(SET CMP0042 OLD)

# Set Release type for builds where CMAKE_BUILD_TYPE is unset
# This is usually a good default as this implictly enables
#
#   CXXFLAGS = -O3 -DNDEBUG
#
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

# Summarise found features at end
include(FeatureSummary)

################
# Dependencies #
################
# Request C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
	message("OpenMP could not be found. ${PROJECT_NAME} will still work, albeit in single-threaded mode")
endif()

# Boost
find_package(Boost 1.59.0 COMPONENTS filesystem iostreams program_options system REQUIRED)


#########
# Files #
#########
add_executable(
	ngshmmalign
		src/aligner.hpp
		src/aligner_impl.hpp
		src/debug.hpp
		src/dna_array.hpp
		src/fasta.hpp
		src/fastq.hpp
		src/hmmalign.hpp
		src/hmmalign_impl.hpp
		src/index_impl.hpp
		src/reference.hpp
		src/reference_impl.hpp
		src/utility_functions.hpp
		src/main.cpp
)

target_include_directories(ngshmmalign PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(ngshmmalign ${Boost_LIBRARIES})

install(
	TARGETS ngshmmalign
	DESTINATION ${CMAKE_INSTALL_BINDIR})

###########
# Summary #
###########
feature_summary(WHAT ALL)