add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)

find_package(GSL)
target_include_directories(${PROJECT_NAME} PUBLIC ${GSL_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
