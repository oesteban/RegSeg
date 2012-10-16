set(DOCUMENTATION "This module contains ITK classes than encapsulate numerical optimizers")

itk_module(RSTKOptimizers
  DEPENDS
    ITKCommon
    ITKOptimizers
    ITKTransform
    ITKImageGrid
    ITKDisplacementField
#  TEST_DEPENDS
#    ITKTestKernel
#    ITKMetricsv4
  DESCRIPTION
    "${DOCUMENTATION}"
)

# ITKOptimizers dependency added to get itkCostFunction for itkSingleValuedCostFunctionv4
