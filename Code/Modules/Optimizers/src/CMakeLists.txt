set(RSTKOptimizers_SRC
  LevelSetsOptimizerBase.cxx
  LevelSetsOptimizer.cxx
  GradientDescentLevelSetsOptimizer.cxx
)

add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
target_link_libraries(RSTKOptimizers
  ${RSTKLevelSets_LIBRARIES}
  )
#itk_module_target(RSTKOptimizers)
