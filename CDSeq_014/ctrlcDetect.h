 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011-2013 University of Oxford
  * Version: 0.10.2
  *
  * University of Oxford means the Chancellor, Masters and Scholars of
  * the University of Oxford, having an administrative office at
  * Wellington Square, Oxford OX1 2JD, UK. 
  *
  * This file is part of Gerardus.
  *
  * This program is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
  *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details. The offer of this
  * program under the terms of the License is subject to the License
  * being interpreted in accordance with English Law and subject to any
  * action against the University of Oxford being under the jurisdiction
  * of the English Courts.
  *
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see
  * <http://www.gnu.org/licenses/>.
  */
 
#ifndef GERARDUSCOMMON_H
#define GERARDUSCOMMON_H
 
/* mex headers */
#include <mex.h>
 
/* C++ headers */
#include <iostream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
 
/* ITK headers */
//#include "itkOffset.h"
 
/*
 * utIsInterruptPending(): "undocumented MATLAB API implemented in
 * libut.so, libut.dll, and included in the import library
 * libut.lib. To use utIsInterruptPending in a mex-file, one must
 * manually declare bool utIsInterruptPending() because this function
 * is not included in any header files shipped with MATLAB. Since
 * libut.lib, by default, is not linked by mex, one must explicitly
 * tell mex to use libut.lib." -- Wotao Yin, 
 * http://www.caam.rice.edu/~wy1/links/mex_ctrl_c_trick/
 *
 */
#ifdef __cplusplus 
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif
 
/*
 * ctrlcCheckPoint(): function to check whether the user has pressed
 * Ctrl+C, and if so, terminate execution returning an error message
 * with a hyperlink to the offending function's help, and a hyperlink
 * to the line in the source code file this function was called from
 *
 * It is implemented as a C++ macro to check for the CTRL+C flag, and
 * a call to function ctrlcErrMsgTxt() inside, to throw the error. The
 * reason is that if ctrlcCheckPoint() were a function instead of a
 * macro, this would introduce a function call at every iteration of
 * the loop, which is very expensive. But then we don't want to put
 * the whole error message part inside a macro, it's bug-prone and bad
 * programming practice. And once the CTRL+C has been detected,
 * whether the error message is generated a bit faster or not is not
 * important.
 *
 * In practice, to use this function put a call like this e.g. inside
 * loops that may take for a very long time:
 *
 *    // exit if user pressed Ctrl+C
 *    ctrlcCheckPoint(__FILE__, __LINE__);
 *
 * sourceFile: full path and name of the C++ file that calls this
 *             function. This should usually be the preprocessor
 *             directive __FILE__
 *
 * lineNumber: line number where this function is called from. This
 *             should usually be the preprocessor directive __LINE__
 *
 */
inline
void ctrlcErrMsgTxt(std::string sourceFile, int lineNumber) {
 
  // run from here the following code in the Matlab side:
  //
  // >> path = mfilename('fullpath')
  //
  // this provides the full path and function name of the function
  // that called ctrlcCheckPoint()
  int nlhs = 1; // number of output arguments we expect
  mxArray *plhs[1]; // to store the output argument
  int nrhs = 1; // number of input arguments we are going to pass
  mxArray *prhs[1]; // to store the input argument we are going to pass
  prhs[0] = mxCreateString("fullpath"); // input argument to pass
  if (mexCallMATLAB(nlhs, plhs, nrhs, prhs, "mfilename")) { // run mfilename('fullpath')
    mexErrMsgTxt("ctrlcCheckPoint(): mfilename('fullpath') returned error");
  }
  if (plhs == NULL) {
    mexErrMsgTxt("ctrlcCheckPoint(): mfilename('fullpath') returned NULL array of outputs");
  }
  if (plhs[0] == NULL) {
    mexErrMsgTxt("ctrlcCheckPoint(): mfilename('fullpath') returned NULL output instead of valid path");
  }
 
  // get full path to current function, including function's name
  // (without the file extension)
  char *pathAndName = mxArrayToString(plhs[0]);
  if (pathAndName == NULL) {
    mexErrMsgTxt("ctrlcCheckPoint(): mfilename('fullpath') output cannot be converted to string");
  }
 
  // for some reason, using mexErrMsgTxt() to give this output
  // doesn't work. Instead, we have to give the output to the
  // standar error, and then call mexErrMsgTxt() to terminate
  // execution of the program
  std::cerr << "Operation terminated by user during "
	    << "<a href=\"matlab:helpUtils.errorDocCallback('"
	    << mexFunctionName()
	    << "', '" << pathAndName << ".m', " << lineNumber << ")\">"
	    << mexFunctionName()
	    << "</a> (<a href=\"matlab:opentoline('"
	    << sourceFile
	    << "'," << lineNumber << ",0)\">line " << lineNumber
	    << "</a>)"
	    << std::endl;
  mexErrMsgTxt("");
}
 
#define ctrlcCheckPoint(sourceFile, lineNumber)		\
  if (utIsInterruptPending()) {				\
    ctrlcErrMsgTxt(sourceFile, lineNumber);		\
  }
