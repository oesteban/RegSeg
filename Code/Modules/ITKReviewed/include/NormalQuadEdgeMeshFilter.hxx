/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __NormalQuadEdgeMeshFilter_hxx
#define __NormalQuadEdgeMeshFilter_hxx

#include "NormalQuadEdgeMeshFilter.h"
#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_matrix.h>

namespace rstk
{
template< typename TInputMesh, typename TOutputMesh >
NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::NormalQuadEdgeMeshFilter():
 m_IsWindingCCW(false),
 m_TotalArea(0.0)
{
  this->m_Weight = THURMER;
}

template< typename TInputMesh, typename TOutputMesh >
NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::~NormalQuadEdgeMeshFilter()
{}

template< typename TInputMesh, typename TOutputMesh >
void
NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::CheckTriangleWinding()
{
	  OutputMeshPointer  output = this->GetOutput();
	  OutputPolygonType *poly;
	  OutputPointType pt[3];

	  int j = 0;

	  for ( OutputCellsContainerConstIterator
	        cell_it = output->GetCells()->Begin();
	        cell_it != output->GetCells()->End();
	        ++cell_it )
	    {
	    poly = dynamic_cast< OutputPolygonType * >( cell_it.Value() );

	    if ( poly != ITK_NULLPTR )
	      {
	      if ( poly->GetNumberOfPoints() == 3 )
	        {
	    	  int             k(0);
	    	  OutputPointType tpt[3];
	    	  OutputQEType *edge = poly->GetEdgeRingEntry();
	    	  OutputQEType *temp = edge;
	    	  do
	    	    {
	    	    tpt[k++] = output->GetPoint( temp->GetOrigin() );
	    	    temp = temp->GetLnext();
	    	    }
	    	  while ( temp != edge );
	    	  pt[j] = TriangleType::ComputeGravityCenter(tpt[0], tpt[1], tpt[2]);
	    	  j++;
	    	  if (j==3) {
	    		  break;
	    	  }
	        }
	      }
	    }

	  double test = (pt[1][0] - pt[0][0]) * (pt[1][1] + pt[0][1]) +
			         (pt[2][0] - pt[1][0]) * (pt[2][1] + pt[1][1]) +
			         (pt[0][0] - pt[2][0]) * (pt[0][1] + pt[2][1]);
	  m_IsWindingCCW = (test < 0);
}

template< typename TInputMesh, typename TOutputMesh >
void
NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::ComputeAllFaceNormals()
{
  CheckTriangleWinding();
  OutputMeshPointer  output = this->GetOutput();
  OutputPolygonType *poly;

  for ( OutputCellsContainerConstIterator
        cell_it = output->GetCells()->Begin();
        cell_it != output->GetCells()->End();
        ++cell_it )
    {
    poly = dynamic_cast< OutputPolygonType * >( cell_it.Value() );

    if ( poly != ITK_NULLPTR )
      {
      if ( poly->GetNumberOfPoints() == 3 )
        {
        output->SetCellData( cell_it->Index(), ComputeFaceNormal(poly) );
        }
      }
    }
}

template< typename TInputMesh, typename TOutputMesh >
void
NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::ComputeAllVertexNormals()
{
  OutputMeshPointer            output = this->GetOutput();
  OutputPointsContainerPointer points = output->GetPoints();
  OutputPointIdentifier        id;

  OutputMeshType *outputMesh = this->GetOutput();

  for ( OutputPointsContainerIterator it = points->Begin();
        it != points->End();
        ++it )
    {
    id = it->Index();
    output->SetPointData( id, ComputeVertexNormal(id,outputMesh));
    }
}

template< typename TInputMesh, typename TOutputMesh >
typename NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::OutputFaceNormalType
NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::ComputeFaceNormal(OutputPolygonType *iPoly)
{
  OutputMeshPointer output = this->GetOutput();

  OutputPointType pt[3];
  int             k(0);

  OutputQEType *edge = iPoly->GetEdgeRingEntry();
  OutputQEType *temp = edge;


  do
    {
    pt[k++] = output->GetPoint( temp->GetOrigin() );
    temp = temp->GetLnext();
    }
  while ( temp != edge );

  // Test of winding, all normals should point inwards.
  OutputFaceNormalType n = TriangleType::ComputeNormal(pt[0], pt[1], pt[2]);
//  InputVectorType e1 = pt[1] - pt[0];
//  InputVectorType e2 = pt[2] - pt[0];
//  double test = dot_product(n.GetVnlVector(), vnl_cross_3d(e1.GetVnlVector(), e2.GetVnlVector()));
//
//  e1 = pt[1] - pt[2];
//  e2 = pt[0] - pt[2];
//  double test2 = dot_product(n.GetVnlVector(), vnl_cross_3d(e1.GetVnlVector(), e2.GetVnlVector()));
//
//  InputVectorType e01 = pt[1] - pt[0];
//  InputVectorType e12 = pt[2] - pt[1];
//  InputVectorType e20 = pt[0] - pt[2];
//  double test = (pt[1][0] - pt[0][0]) * (pt[1][1] + pt[0][1]) +
//		         (pt[2][0] - pt[1][0]) * (pt[2][1] + pt[1][1]) +
//		         (pt[0][0] - pt[2][0]) * (pt[0][1] + pt[2][1]);

  if (m_IsWindingCCW)
	  n*= -1.0;

  return n;
}

template< typename TInputMesh, typename TOutputMesh >
typename NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::OutputVertexNormalType
NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::ComputeVertexNormal(const OutputPointIdentifier & iId, OutputMeshType *outputMesh)
{

  OutputQEType *       edge = outputMesh->FindEdge(iId);
  OutputQEType *       temp = edge;
  OutputCellIdentifier cell_id(0);

  OutputVertexNormalType n(0.);
  OutputFaceNormalType   face_normal(0.);

  OutputPointType pt[3];

  double wi, twi;
  size_t nfaces = 0;

  do
    {
    cell_id = temp->GetLeft();
    if ( cell_id != OutputMeshType::m_NoFace )
      {
      outputMesh->GetCellData(cell_id, &face_normal);
      wi = Weight(iId, cell_id,outputMesh);
      twi= wi;
      nfaces++;
      n += face_normal * wi;
      }
    temp = temp->GetOnext();
    }
  while ( temp != edge );

  double a = twi / nfaces;
  this->m_VertexAreaContainer[iId] = a;
  this->m_TotalArea+= a;
  n.Normalize();
  return n;
}

template< typename TInputMesh, typename TOutputMesh >
typename NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >::OutputVertexNormalComponentType
NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::Weight(const OutputPointIdentifier & iPId, const OutputCellIdentifier & iCId, OutputMeshType *outputMesh)
{
  if ( m_Weight == GOURAUD )
    {
    return static_cast< OutputVertexNormalComponentType >( 1.0 );
    }
  else
    {

    OutputPolygonType *poly = dynamic_cast< OutputPolygonType * >(
      outputMesh->GetCells()->GetElement(iCId) );
    if ( poly != ITK_NULLPTR ) // this test should be removed...
      {
      // this test should be removed...
      if ( poly->GetNumberOfPoints() == 3 )
        {
        int              internal_id(0), k(0);
        OutputPointType  pt[3];
        OutputVectorType u, v;

        OutputQEType *edge = poly->GetEdgeRingEntry();
        OutputQEType *temp = edge;
        do
          {
          pt[k] = outputMesh->GetPoint( temp->GetOrigin() );
          if ( temp->GetOrigin() == iPId )
            {
            internal_id = k;
            }

          temp = temp->GetLnext();
          k++;
          }
        while ( temp != edge );

        switch ( m_Weight )
          {
          default:
          case GOURAUD:
            {
            return static_cast< OutputVertexNormalComponentType >( 1. );
            }
          case THURMER:
            {
            // this implementation may be included inside itkTriangle
            switch ( internal_id )
              {
              case 0:
                u = pt[1] - pt[0];
                v = pt[2] - pt[0];
                break;
              case 1:
                u = pt[0] - pt[1];
                v = pt[2] - pt[1];
                break;
              case 2:
                u = pt[0] - pt[2];
                v = pt[1] - pt[2];
                break;
              }
            typename OutputVectorType::RealValueType norm_u = u.GetNorm();
            if ( norm_u > vnl_math::eps )
              {
              norm_u = 1. / norm_u;
              u *= norm_u;
              }

            typename OutputVectorType::RealValueType norm_v = v.GetNorm();
            if ( norm_v > vnl_math::eps )
              {
              norm_v = 1. / norm_v;
              v *= norm_v;
              }
            return static_cast< OutputVertexNormalComponentType >(
                     std::acos(u * v) );
            }
          case AREA:
            {
            return static_cast< OutputVertexNormalComponentType >(
                     TriangleType::ComputeArea(pt[0], pt[1], pt[2]) );
            }
          }
        }
      else
        {
        std::cout << "Input should be a triangular mesh!!!" << std::endl;
        return static_cast< OutputVertexNormalComponentType >( 0. );
        }
      }
    else
      {
      return static_cast< OutputVertexNormalComponentType >( 0. );
      }
    }
}

template< typename TInputMesh, typename TOutputMesh >
void NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::GenerateData()
{
	this->CopyInputMeshToOutputMesh();
	this->m_VertexAreaContainer.SetSize(this->GetOutput()->GetNumberOfPoints());
	this->m_VertexAreaContainer.Fill(0.0);
	this->m_TotalArea = 0.0;
	this->ComputeAllFaceNormals();
	this->ComputeAllVertexNormals();
}

template< typename TInputMesh, typename TOutputMesh >
void NormalQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  std::cout << indent << "Weight: " << m_Weight << std::endl;
}
}

#endif
