add_executable(IO_EMC EXCLUDE_FROM_ALL IO_EMC.cpp)
add_dependencies(examples IO_EMC)
#LocalMech
target_link_libraries(IO_EMC PRIVATE arbor arborio arborenv arbor-sup ${json_library_name} ext-random123)
