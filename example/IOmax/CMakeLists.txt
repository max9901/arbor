
#this actually adds the custom Catalogue
file(GLOB smolLocalMech
        "mechanisms/smolLocal/*.cpp"
    )
message(STATUS "foo include dir: ${smolmech}")
add_library(smolLocalCat mechanisms/smolLocal_catalogue.cpp ${smolLocalMech})
target_link_libraries(smolLocalCat arbor ext-random123)

add_executable(IOmax EXCLUDE_FROM_ALL IOmax.cpp)
add_dependencies(examples IOmax)
target_link_libraries(IOmax PRIVATE smolLocalCat arbor arborio arborenv arbor-sup ${json_library_name} ext-random123)
