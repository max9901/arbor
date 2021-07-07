#this actually adds the custom Catalogue

#Local_catalogue
file(GLOB Local_catalogue
        "IOU_Mechs/*.cpp"
        "IOU_catalogue.cpp"
    )

set(arbor_mechanism_sources_custom
        ${Local_catalogue}
        PARENT_SCOPE
   )