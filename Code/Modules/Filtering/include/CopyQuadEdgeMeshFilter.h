/*
 * CopyQuadEdgeMeshFilter.h
 *
 *  Created on: Sep 11, 2013
 *      Author: oesteban
 */

#ifndef COPYQUADEDGEMESHFILTER_H_
#define COPYQUADEDGEMESHFILTER_H_

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>

namespace rstk
{
template<typename TInputMesh, typename TOutputMesh>
class CopyQuadEdgeMeshFilter : public itk::QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
    typedef CopyQuadEdgeMeshFilter    Self;
    typedef itk::QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh> Superclass;
    typedef itk::SmartPointer<Self>        Pointer;

    itkNewMacro(Self);
protected:
    CopyQuadEdgeMeshFilter(){}
    ~CopyQuadEdgeMeshFilter(){}

    void GenerateData() {
        this->CopyInputMeshToOutputMesh();
    }
};

}//namespace itk end

#endif /* COPYQUADEDGEMESHFILTER_H_ */
