#=============================================================================
# CMake configuration file for the STEERING_CONTROLLER demo.
# This example program works only with Irrlicht run-time visualization
#=============================================================================

IF(NOT ENABLE_MODULE_IRRLICHT)
    RETURN()
ENDIF()

MESSAGE(STATUS "...add STEERING_CONTROLLER demo")

#--------------------------------------------------------------
# List all model files for this demo

SET(MODEL_FILES
    ../hmmwv/HMMWV.h
    ../hmmwv/HMMWV.cpp
    ../hmmwv/vehicle/HMMWV_Vehicle.h
    ../hmmwv/vehicle/HMMWV_Vehicle.cpp
    ../hmmwv/vehicle/HMMWV_VehicleReduced.h
    ../hmmwv/vehicle/HMMWV_VehicleReduced.cpp
    ../hmmwv/suspension/HMMWV_DoubleWishbone.h
    ../hmmwv/suspension/HMMWV_DoubleWishbone.cpp
    ../hmmwv/suspension/HMMWV_DoubleWishboneReduced.h
    ../hmmwv/suspension/HMMWV_DoubleWishboneReduced.cpp
    ../hmmwv/steering/HMMWV_PitmanArm.h
    ../hmmwv/steering/HMMWV_PitmanArm.cpp
    ../hmmwv/steering/HMMWV_RackPinion.h
    ../hmmwv/steering/HMMWV_RackPinion.cpp
    ../hmmwv/driveline/HMMWV_Driveline2WD.h
    ../hmmwv/driveline/HMMWV_Driveline2WD.cpp
    ../hmmwv/driveline/HMMWV_Driveline4WD.h
    ../hmmwv/driveline/HMMWV_Driveline4WD.cpp
    ../hmmwv/powertrain/HMMWV_SimplePowertrain.h
    ../hmmwv/powertrain/HMMWV_SimplePowertrain.cpp
    ../hmmwv/powertrain/HMMWV_Powertrain.h
    ../hmmwv/powertrain/HMMWV_Powertrain.cpp
    ../hmmwv/brake/HMMWV_BrakeSimple.h
    ../hmmwv/brake/HMMWV_BrakeSimple.cpp
    ../hmmwv/wheel/HMMWV_Wheel.h
    ../hmmwv/wheel/HMMWV_Wheel.cpp
    ../hmmwv/tire/HMMWV_RigidTire.h
    ../hmmwv/tire/HMMWV_RigidTire.cpp
    ../hmmwv/tire/HMMWV_LugreTire.h
    ../hmmwv/tire/HMMWV_LugreTire.cpp
    ../hmmwv/tire/HMMWV_FialaTire.h
    ../hmmwv/tire/HMMWV_FialaTire.cpp
)

SET(DEMO_FILES
    demo_SteeringController.cpp
)

SOURCE_GROUP("subsystems" FILES ${MODEL_FILES})
SOURCE_GROUP("" FILES ${DEMO_FILES})

#--------------------------------------------------------------
# Additional include directories

INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/demos/vehicle")
INCLUDE_DIRECTORIES(${CH_IRRLICHTINC})


#--------------------------------------------------------------
# List of all required libraries

SET(LIBRARIES
    ChronoEngine
    ChronoEngine_irrlicht
    ChronoEngine_vehicle)

#--------------------------------------------------------------
# Create the executable

ADD_EXECUTABLE(demo_SteeringController ${DEMO_FILES} ${MODEL_FILES})
SET_TARGET_PROPERTIES(demo_SteeringController PROPERTIES 
                      COMPILE_FLAGS "${CH_CXX_FLAGS} ${CH_IRRLICHT_CXX_FLAGS}"
                      LINK_FLAGS "${LINKERFLAG_EXE}")
TARGET_LINK_LIBRARIES(demo_SteeringController ${LIBRARIES})
INSTALL(TARGETS demo_SteeringController DESTINATION bin)
