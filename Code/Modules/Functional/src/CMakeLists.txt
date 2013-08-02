set(RSTKFunctional_SRC "")

add_library(RSTKFunctional ${RSTKFunctional_SRC})
target_link_libraries(RSTKFunctional ITKCommon )

#itk_module_target(RSTKFunctional)
