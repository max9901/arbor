
#this actually adds the custom Catalogue
#file(GLOB LocalMech
#        "mechanisms/IOU_Mechs/*.cpp"
#    )
#add_library(LocalMech mechanisms/Local_catalogue.cpp ${LocalMech})
#target_link_libraries(LocalMech arbor ext-random123)





add_executable(IOmax EXCLUDE_FROM_ALL IOmax.cpp)
add_dependencies(examples IOmax)
#LocalMech
target_link_libraries(IOmax PRIVATE arbor arborio arborenv arbor-sup ${json_library_name} ext-random123)
