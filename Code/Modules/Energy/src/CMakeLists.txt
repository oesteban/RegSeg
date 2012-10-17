set(RSTKEnergy_SRC
  EnergyBase.cxx
)

add_library(RSTKEnergy ${RSTKEnergy_SRC})
target_link_libraries(RSTKEnergy ITKCommon )

#itk_module_target(RSTKEnergy)
