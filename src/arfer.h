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

  template<typename TConfig, typename TGlVector, typename TValue>
  inline void
  _estBiallelicAF(TConfig const& c, TGlVector const& glVector, TValue (&hweAF)[2]) {
    if (!glVector.empty()) {
      TValue numGl = glVector.size();
      TValue afprior[2];
      afprior[0] = 0.5;
      afprior[1] = 0.5;
      TValue gtprior[3];
      TValue gt[3];
      TValue p;
      TValue err = 1;
      for(std::size_t count = 0; ((err > c.epsilon) && (count<c.maxiter)); ++count) {
	gtprior[0] = afprior[0] * afprior[0];
	gtprior[1] = 2 * afprior[0] * afprior[1];
	gtprior[2] = afprior[1] * afprior[1];
      
	hweAF[0] = 0;
	hweAF[1] = 0;
	for(typename TGlVector::const_iterator itG = glVector.begin(); itG !=glVector.end();++itG) {
	  gt[0] = gtprior[0] * itG->at(0);
	  gt[1] = gtprior[1] * itG->at(1);
	  gt[2] = gtprior[2] * itG->at(2);
	  p = gt[0] + gt[1] + gt[2];
	  gt[0] /= p;
	  gt[1] /= p;
	  gt[2] /= p;
	  hweAF[0] += gt[0] + 0.5 * gt[1];
	  hweAF[1] += gt[2] + 0.5 * gt[1];
	}
	hweAF[0] /= numGl;
	hweAF[1] /= numGl;
	err = (afprior[0]-hweAF[0])*(afprior[0]-hweAF[0]) + (afprior[1]-hweAF[1])*(afprior[1]-hweAF[1]);
	afprior[0] = hweAF[0];
	afprior[1] = hweAF[1];
      }
    }
  }


  template<typename TConfig, typename TGlVector, typename TValue>
  inline void
  _estBiallelicGTFreq(TConfig const& c, TGlVector const& glVector, TValue (&mleGTFreq)[3]) {
    if (!glVector.empty()) {
      TValue numGl = glVector.size();
      TValue prior[3];
      prior[0] = 1.0/3.0;
      prior[1] = 1.0/3.0;
      prior[2] = 1.0/3.0;
      TValue gt[3];
      TValue p;
      TValue err = 1;
      for(std::size_t count = 0; ((err > c.epsilon) && (count<c.maxiter)); ++count) {
	mleGTFreq[0] = 0;
	mleGTFreq[1] = 0;
	mleGTFreq[2] = 0;
	for(typename TGlVector::const_iterator itG = glVector.begin(); itG !=glVector.end();++itG) {
	  gt[0] = prior[0] * itG->at(0);
	  gt[1] = prior[1] * itG->at(1);
	  gt[2] = prior[2] * itG->at(2);
	  p = gt[0] + gt[1] + gt[2];
	  mleGTFreq[0] += gt[0]/p;
	  mleGTFreq[1] += gt[1]/p;
	  mleGTFreq[2] += gt[2]/p;
	}
	mleGTFreq[0] /= numGl;
	mleGTFreq[1] /= numGl;
	mleGTFreq[2] /= numGl;
	err = (prior[0]-mleGTFreq[0])*(prior[0]-mleGTFreq[0]) + (prior[1]-mleGTFreq[1])*(prior[1]-mleGTFreq[1]) + (prior[2]-mleGTFreq[2])*(prior[2]-mleGTFreq[2]);
	prior[0] = mleGTFreq[0];
	prior[1] = mleGTFreq[1];
	prior[2] = mleGTFreq[2];
      }
    }
  }

}

#endif

