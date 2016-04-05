/*CXXR $Id: Print.h 749 2010-01-07 18:14:35Z arr $
  *CXXR
  *CXXR This file is part of CXXR, a project to refactor the R interpreter
  *CXXR into C++.  It may consist in whole or in part of program code and
  *CXXR documentation taken from the R project itself, incorporated into
  *CXXR CXXR (and possibly MODIFIED) under the terms of the GNU General Public
  *CXXR Licence.
  *CXXR 
  *CXXR CXXR is Copyright (C) 2008-10 Andrew R. Runnalls, subject to such other
  *CXXR copyrights and copyright restrictions as may be stated below.
  *CXXR 
  *CXXR CXXR is not part of the R project, and bugs and other issues should
  *CXXR not be reported via r-bugs or other R project channels; instead refer
  *CXXR to the CXXR website.
  *CXXR */
 
/*
  *  R : A Computer Language for Statistical Data Analysis
  *  Copyright (C) 1998-2008    Robert Gentleman, Ross Ihaka
  *                             and the R Development Core Team
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU Lesser General Public License as published by
  *  the Free Software Foundation; either version 2.1 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU Lesser General Public License for more details.
  *
  *  You should have received a copy of the GNU Lesser General Public License
  *  along with this program; if not, a copy is available at
  *  http://www.r-project.org/Licenses/
*/

#ifndef R_EXT_PRINT_H_
#define R_EXT_PRINT_H_
 
#ifndef NO_C_HEADERS
# include <stdarg.h>
#endif
 
#ifdef  __cplusplus
extern "C" {
#endif
 
void Rprintf(const char *, ...);
void REprintf(const char *, ...);
void Rvprintf(const char *, va_list);
void REvprintf(const char *, va_list);
 
#ifdef  __cplusplus
}
#endif
 
#endif /* R_EXT_PRINT_H_ */
