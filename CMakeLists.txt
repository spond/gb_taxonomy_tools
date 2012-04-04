cmake_minimum_required(VERSION 2.8)
project(TN93)

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)


set(CMAKE_CONFIGURATION_TYPES Release)

include_directories(
    src/include
)

add_executable(
    gid-taxid
    src/gid-taxid.c
    src/util/avl.c
    src/util/buffers.c
)

add_executable(
    taxonomy-reader
    src/taxonomy_reader.c
    src/util/avl.c
    src/util/buffers.c
    src/util/taxinfo.c
)


add_executable(
    taxonomy2tree
    src/taxonomy2tree.c
    src/util/avl.c
    src/util/buffers.c
    src/util/taxinfo.c
)


add_executable(
    tree2ps
    src/tree2ps.c
    src/util/avl.c
    src/util/buffers.c
)

set(DEFAULT_COMPILE_FLAGS "-O3")
set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")

install(
    TARGETS gid-taxid
	RUNTIME DESTINATION bin
)

install(
    TARGETS taxonomy-reader
	RUNTIME DESTINATION bin
)

install(
    TARGETS taxonomy2tree
	RUNTIME DESTINATION bin
)

install(
    TARGETS tree2ps
	RUNTIME DESTINATION bin
)

add_custom_target(gid2taxid 
	DEPENDS gid-taxid
)

add_custom_target(taxonomyparser
	DEPENDS taxonomy-reader
)

add_custom_target(T2T
	DEPENDS taxonomy2tree
)

add_custom_target(T2PS
	DEPENDS tree2ps
)

