project(ADBA)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(Boost 1.55 COMPONENTS system filesystem program_options REQUIRED)

include_directories(
	${Boost_INCLUDE_DIRS}
)

link_directories(
	${Boost_LIBRARY_DIRS}
)

set(COMPILER_FLAGS
	"-ggdb -Wall -L/usr/include -I/usr/include -std=c++0x -pipe -Wno-deprecated -pedantic"
)

add_executable(adba
	main.cpp
	deBruijnGraph.h
	deBruijnGraph.cpp
	Sequence.h
	Sequence.cpp
	Vertex.h
	Vertex.cpp
	UnitigGraph.h
	UnitigGraph.cpp
)

set_target_properties(adba PROPERTIES
	COMPILE_FLAGS ${COMPILER_FLAGS}
)

target_link_libraries(adba ${Boost_LIBRARIES})
