function c = periodic_idwt(d,a,astart,at,atstart,j0)
%periodic_idwt() inverse periodic discrete wavelet transform (reconstruction)
%
%   periodic_idwt(d,a,astart,at,atstart,j0) computes the
%   inverse periodic discrete wavelet transform of d (reconstruction),
%   starting at the coarsest level j0,
%   with primal mask a starting at index astart,
%   and dual mask at starting at index atstart.
%
%   special cases:
%   periodic_idwt(d,a,astart,at,atstart) uses j0=j0_per(a,at)
%   periodic_idwt(d,a,astart,j0) uses at=a, atstart=astart
%   periodic_idwt(d,a,j0) uses at=a, atstart=astart=0
%   periodic_idwt(d,a) uses at=a, atstart=astart=0 and j0=j0_per(a)
%
%   INPUT
%   d: wavelet coefficients in V_J, row vector, length 2^J
%   a: primal refinement mask (row vector) [a_astart,...,a_{astart+L}], L\ge 2
%   astart: index of first primal mask entry
%   at: dual refinement mask (row vector) [at_atstart,...,at_{atstart+Lt}], Lt\ge 2
%   atstart: index of first dual mask entry
%   j0: coarsest level for IDWT, 0<=j0<J, default value 0
%
%   OUTPUT
%   c: generator coefficients in V_J, row vector, length 2^J
%
%   This file is part of MS-Lib, the multiscale library.
%   Copyright (c) 2011-2012 Thorsten Raasch
%
%   MS-Lib is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   MS-Lib is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MS-Lib.  If not, see <http://www.gnu.org/licenses/>.

if nargin==2,
    % orthogonal case with astart=0 and minimally possible refinement level
    j0=j0_per(a);
    at=a;
    astart=0;
    atstart=0;
end

if nargin==3,
    % orthogonal case with astart=0 and custom minimal refinement level
    j0=astart;
    at=a;    
    astart=0;
    atstart=0;    
end

if nargin==4,
    % orthogonal case with custom astart and custom j0
    j0=at;
    at=a;
    atstart=astart;
end

if nargin==5,
    % nonorthogonal case with minimally possible refinement level
    j0=j0_per(a,astart);
end

J = log2(length(d));
c = d;
for j=j0+1:J,
    c(1:pow2(j)) = c(1:pow2(j))*Mj_per(a,astart,at,atstart,j);
end
