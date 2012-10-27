set(RSTKLevelSets_SRC "")

add_library(RSTKLevelSets ${RSTKLevelSets_SRC})
target_link_libraries(RSTKLevelSets ITKCommon )

#itk_module_target(RSTKLevelSets)
