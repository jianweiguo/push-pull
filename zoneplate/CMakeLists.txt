
IF(WIN32)
	# Od : no optimisation; Ot: favor code speed over code size.
	if(DO_NOT_pushpull2d_WINDOWS_RELEASE)
		SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Od /Ot")
	endif(DO_NOT_pushpull2d_WINDOWS_RELEASE)
ELSE(WIN32)
	# CGAL needs rounding math to compile well with GCC
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
ENDIF(WIN32)

AUX_SOURCE_DIRECTORY(. SUB_SOURCES)
SET(SOURCES ${SOURCES} ${SUB_SOURCES})
file(GLOB HeaderFiles *.h)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/zoneplate)

IF (WIN32)
   ADD_EXECUTABLE(zoneplate ${SOURCES}  ${HeaderFiles})
ELSE (WIN32)
   LINK_DIRECTORIES(/usr/X11R6/lib)
   ADD_EXECUTABLE(zoneplate_stk ${SOURCES}  ${HeaderFiles})
ENDIF (WIN32)
 
TARGET_LINK_LIBRARIES(zoneplate ThirdParty)
