/*
 * RadialBasisFunction.h
 *
 *  Created on: Jun 10, 2013
 *      Author: oesteban
 */

#ifndef RADIALBASISFUNCTION_H_
#define RADIALBASISFUNCTION_H_

#include "itkNumericTraits.h"
#include <cmath>

namespace rstk {
namespace RBF {


template< class TPixel, class TScalarType = double, unsigned int NDimensions = 3>
class RadialBasisFunction
{
public:
	typedef RadialBasisFunction             Self;
	typedef TPixel                          PixelType;
	typedef typename PixelType::VectorType  VectorType;

	itkStaticConstMacro( Dimension, unsigned int, NDimensions );

	RadialBasisFunction() {};
	~RadialBasisFunction(){};

	//virtual inline TScalarType operator()( const TPixel &center, const TPixel &point ) const =0;
	virtual inline TScalarType GetWeight( const TPixel &center, const TPixel &point ) const =0;
};
}
}



#endif /* RADIALBASISFUNCTION_H_ */
