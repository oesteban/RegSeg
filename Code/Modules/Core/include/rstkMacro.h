// This file is part of RegSeg
//
// Copyright 2014-2017, Oscar Esteban <code@oscaresteban.es>
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.

#ifndef RSTKMACRO_H_
#define RSTKMACRO_H_

#include <itkMacro.h>


/*
    if ( this->m_##name.size() <= id ) { \
		itkExceptionMacro( << "vector " #name "is not initialized, or position " << id << " not valid."); \
	} \
 */


#define rstkSetVectorElement(name,type) \
virtual void Set##name##Element ( size_t id, const type _arg ) { \
	if ( this->m_##name.size() <= id ) { \
		itkExceptionMacro( << "vector " #name " is not initialized, or position " << id << " not valid."); \
	} \
	if ( this->m_##name[id] != _arg ) { \
		this->m_##name[id] = _arg; \
		this->Modified(); \
	} \
}

#define rstkGetVectorElement(name,type) \
virtual type Get##name##Element ( size_t id ) { \
	return this->m_##name[id]; \
}

#define rstkGetConstVectorElement( name, type ) \
virtual type Get##name##Element( size_t id ) const { \
	return this->m_##name[id]; \
}

#define rstkFillVector(name,type) \
virtual void Fill##name ( const type _arg ) { \
	for( size_t id = 0; id < this->m_##name[id]; id++) { \
		if ( this->m_##name[id] != _arg ) { \
			this->m_##name[id] = _arg; \
			this->Modified(); \
		} \
	} \
}

#define rstkVectorMethods(name,type) \
virtual void Set##name##Element ( size_t id, const type _arg ) { \
	if ( this->m_##name.size() <= id ) { \
		itkExceptionMacro( << "vector " #name " is not initialized, or position " << id << " not valid."); \
	} \
	if ( this->m_##name[id] != _arg ) { \
		this->m_##name[id] = _arg; \
		this->Modified(); \
	} \
} \
virtual type Get##name##Element( const size_t id ) const { \
	return this->m_##name[id]; \
} \
virtual void Fill##name ( const type _arg ) { \
	for( size_t id = 0; id < this->m_##name[id]; id++) { \
		if ( this->m_##name[id] != _arg ) { \
			this->m_##name[id] = _arg; \
			this->Modified(); \
		} \
	} \
}



#define rstkGetObjectList(name, type) \
virtual type * Get##name##OfLevel ( const size_t id ) const  { \
    return this->m_##name##s[id].GetPointer();    \
}

#define rstkGetObjectListWithLast(name, type) \
virtual type * Get##name##OfLevel ( const int id ) const  { \
	if ( this->m_##name##s.empty() ) { \
		itkExceptionMacro(<< "list " #name " is empty." );\
	} \
	size_t real_id = id; \
	if( id == -1 ) { \
		real_id = this->m_##name##s.size() -1; \
	} else if ( id < -1 || id > this->m_##name##s.size() ) { \
		itkExceptionMacro(<< "attempted to access invalid position in list " #name "." );\
	} \
    return this->m_##name##s[real_id].GetPointer(); \
}


#define rstkSetObjectList(name, type) \
virtual void Set##name##OfLevel ( const size_t id, type* obj )  { \
	if ( this->m_##name##s.size() <= id ) { \
		itkExceptionMacro( << "vector " #name " is not initialized, or position " << id << " not valid."); \
	} \
	if ( this->m_##name##s[id] != _arg ) { \
		this->m_##name##s[id] = _arg; \
		this->Modified(); \
	} \
}

/** Get a smart const pointer to an object.  Creates the member
 * Get"name"() (e.g., GetPoints()). */
#define rstkGetConstObjectList(name, type) \
virtual const type * Get##name##OfLevel ( const size_t id ) const  { \
    return this->m_##name##s[id].GetPointer();    \
}
#endif /* RSTKMACRO_H_ */
