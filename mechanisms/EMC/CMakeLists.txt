#Local_catalogue
if(ARB_WITH_GPU)
    file(GLOB Local_catalogue
            "EMC_Mechs/*.cpp"
            "EMC_catalogue.cpp"
            "EMC_Mechs/*.cu"
            )
else()
    file(GLOB Local_catalogue
            "EMC_Mechs/*.cpp"
            "EMC_catalogue.cpp"
            )
endif()

set(EMC_mechanism_sources_custom
        ${Local_catalogue}
        PARENT_SCOPE
   )