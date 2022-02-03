#Local_catalogue
if(ARB_WITH_GPU)
    file(GLOB Local_catalogue
            "new_default/*.cpp"
            "new_default_catalogue.cpp"
            "new_default/*.cu"
            )
else()
    file(GLOB Local_catalogue
            "new_default/*.cpp"
            "new_default_catalogue.cpp"
            )
endif()

set(new_default_sources_custom
        ${Local_catalogue}
        PARENT_SCOPE
        )