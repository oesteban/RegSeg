/*
 * CopyQuadEdgeMeshFilter.h
 *
 *  Created on: Sep 11, 2013
 *      Author: oesteban
 */

#ifndef COPYQUADEDGEMESHFILTER_H_
#define COPYQUADEDGEMESHFILTER_H_

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>

namespace itk
{
template<typename TInputMesh, typename TOutputMesh>
class CopyQuadEdgeMeshFilter : public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
    typedef CopyQuadEdgeMeshFilter    Self;
    typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh> Superclass;
    typedef SmartPointer<Self>        Pointer;

    itkNewMacro(Self);
protected:
    void GenerateData()
    {
        this->CopyInputMeshToOutputMesh();
    }
};

}//namespace itk end

#endif /* COPYQUADEDGEMESHFILTER_H_ */
