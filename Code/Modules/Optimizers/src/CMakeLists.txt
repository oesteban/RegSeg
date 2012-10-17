set(RSTKOptimizers_SRC
  EnergyOptimizerBase.cxx
  EnergyOptimizer.cxx
  GradientDescentEnergyOptimizer.cxx
)

add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
target_link_libraries(RSTKOptimizers
  ${RSTKEnergy_LIBRARIES}
  ${ITKOptimizers_LIBRARIES}
  )
itk_module_target(RSTKOptimizers)