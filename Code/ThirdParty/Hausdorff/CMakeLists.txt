# Copyright (c) 2011 LTSI INSERM U642
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
#     * Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#     * Neither name of LTSI, INSERM nor the names
# of any contributors may be used to endorse or promote products derived from this
# software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

PROJECT (vtkHausdorffDistancePointSetFilter)

### The following are needed to avoid warnings. See CMake doc. ######
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

if(COMMAND cmake_policy)
cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)
#####################################################################


find_package(VTK COMPONENTS vtkCommonCore vtkIOLegacy vtkFiltersSources NO_MODULE)

ADD_LIBRARY( vtkHausdorffDistance SHARED vtkHausdorffDistancePointSetFilter.cxx )
set_property(TARGET vtkHausdorffDistance APPEND PROPERTY COMPILE_DEFINITIONS "${VTK_DEFINITIONS}")

TARGET_LINK_LIBRARIES( vtkHausdorffDistance ${VTK_LIBRARIES})
set_property(TARGET vtkHausdorffDistance APPEND PROPERTY INCLUDE_DIRECTORIES ${VTK_INCLUDE_DIRS})

INSTALL(TARGETS vtkHausdorffDistance RUNTIME DESTINATION bin)
