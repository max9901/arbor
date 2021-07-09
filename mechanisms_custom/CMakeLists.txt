#this actually adds the custom Catalogue

#Local_catalogue
file(GLOB Local_catalogue
        "IOU_Mechs/*.cpp"
        "IOU_catalogue.cpp"
    )

if(ARB_WITH_GPU)
    message(STATUS "custom mechs reports is has a GPU backend")
    file(GLOB Local_catalogue
            "IOU_Mechs/*.cpp"
            "IOU_catalogue.cpp"
            "IOU_Mechs/*.cu"
            )
else()
    message(STATUS "no GPU")
    file(GLOB Local_catalogue
            "IOU_Mechs/*.cpp"
            "IOU_catalogue.cpp"
            )
endif()

set(arbor_mechanism_sources_custom
        ${Local_catalogue}
        PARENT_SCOPE
   )