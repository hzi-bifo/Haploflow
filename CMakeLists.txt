project(Haploflow)
set(Boost_USE_STATIC_LIBS ON)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(Boost 1.54 COMPONENTS system filesystem program_options REQUIRED)

include_directories(
	${Boost_INCLUDE_DIRS}
)

link_directories(
	${Boost_LIBRARY_DIRS}
)

set(COMPILER_FLAGS
	"-ggdb -Wall -L/usr/include -I/usr/include -std=c++0x -pipe -Wno-deprecated -pedantic"
)

add_executable(haploflow
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

set_target_properties(haploflow PROPERTIES
	COMPILE_FLAGS ${COMPILER_FLAGS}
)

target_link_libraries(haploflow ${Boost_LIBRARIES})
