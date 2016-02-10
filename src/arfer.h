/*
==============================================================================
VCFaid,
This file contains modified functions derived from arfer and adapted to htslib
==============================================================================
Copyright (c) 2013 Adrian Tan, Erik Garrison

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
==============================================================================
*/

#ifndef ARFER_H
#define ARFER_H

namespace vcfaid
{

  template<typename TGlVector, typename TValue>
    inline void
    _estBiallelicAF(TGlVector const& glVector, TValue& hweAF0, TValue& hweAF1) {
    if (!glVector.empty()) {
      TValue numGl = glVector.size();
      TValue af0 = hweAF0;
      TValue af1 = hweAF1;
      TValue err = 1;
      while(err > 1e-20) {
	TValue AA = af0 * af0;
	TValue Aa = 2 * af0 * af1;
	TValue aa = af1 * af1;
      
	hweAF0 = 0;
	hweAF1 = 0;
	for(typename TGlVector::const_iterator itG = glVector.begin(); itG !=glVector.end();++itG) {
	  TValue gt[3];
	  gt[0] = AA * itG->at(0);
	  gt[1] = Aa * itG->at(1);
	  gt[2] = aa * itG->at(2);
	  TValue p = gt[0] + gt[1] + gt[2];
	  gt[0] /= p;
	  gt[1] /= p;
	  gt[2] /= p;
	  hweAF0 += gt[0] + 0.5 * gt[1];
	  hweAF1 += gt[2] + 0.5 * gt[1];
	}
	hweAF0 /= numGl;
	hweAF1 /= numGl;
	err = (af0-hweAF0)*(af0-hweAF0) + (af1-hweAF1)*(af1-hweAF1);
	af0 = hweAF0;
	af1 = hweAF1;
      }
    }
  }


  template<typename TGlVector, typename TValue>
  inline void
  _estBiallelicGTFreq(TGlVector const& glVector, TValue& mlehomRef, TValue& mlehet, TValue& mlehomAlt) {
    if (!glVector.empty()) {
      TValue numGl = glVector.size();
      TValue homRefPrior = mlehomRef;
      TValue hetPrior = mlehet;
      TValue homAltPrior = mlehomAlt;
      TValue err = 1;
      for(std::size_t count = 0; ((err > 1e-20) && (count<1000)); ++count) {
	mlehomRef = 0;
	mlehet = 0;
	mlehomAlt = 0;
	for(typename TGlVector::const_iterator itG = glVector.begin(); itG !=glVector.end();++itG) {
	  TValue gt[3];
	  gt[0] = homRefPrior * itG->at(0);
	  gt[1] = hetPrior * itG->at(1);
	  gt[2] = homAltPrior * itG->at(2);
	  TValue p = gt[0] + gt[1] + gt[2];
	  mlehomRef += gt[0]/p;
	  mlehet += gt[1]/p;
	  mlehomAlt += gt[2]/p;
	}
	mlehomRef /= numGl;
	mlehet /= numGl;
	mlehomAlt /= numGl;
	err = (homRefPrior-mlehomRef)*(homRefPrior-mlehomRef)+(hetPrior-mlehet)*(hetPrior-mlehet)+(homAltPrior-mlehomAlt)*(homAltPrior-mlehomAlt);
	homRefPrior = mlehomRef;
	hetPrior = mlehet;
	homAltPrior = mlehomAlt;
      }
    }
  }

}

#endif

