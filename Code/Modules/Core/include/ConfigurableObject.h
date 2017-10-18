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

#ifndef CONFIGURABLEOBJECT_H_
#define CONFIGURABLEOBJECT_H_

#include <map>
#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

namespace rstk {

class ConfigurableObject {
public:
	typedef bpo::variables_map  SettingsMap;
	typedef bpo::options_description SettingsDesc;

	//virtual void AddOptions( SettingsDesc& opts ) const = 0;
	void SetSettings( SettingsMap& s ) { this->m_Settings = s; }

protected:
	//virtual ConfigurableObject() {};
	virtual ~ConfigurableObject() {};


	virtual void ParseSettings() = 0;
	SettingsMap m_Settings;
private:
};


}  // namespace rstk

#endif /* CONFIGURABLEOBJECT_H_ */
