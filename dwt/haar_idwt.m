function c = haar_idwt(d,j0)
%haar_dwt() inverse discrete Haar wavelet transform (reconstruction)
%   haar_idwt(d,j0) computes the inverse discrete Haar wavelet transform of a
%   signal d of length 2^J down to the coarsest level 0<=j0<J.
%
%   haar_idwt(d) computes the inverse discrete Haar wavelet transform of d,
%   down to the coarsest level j0=0
%
%   INPUT
%   d: wavelet coefficients in V_J, row vector, length 2^J
%   j0: coarsest level for DWT, 0<=j0<=J, default value 0
%
%   OUTPUT
%   c: generator coefficients in V_J, row vector, length 2^J

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

if nargin<2,
    j0=0;
end
    
J = log2(length(d));
c = d;
for j = j0:J-1,
    c(1:pow2(j+1)) = c(1:pow2(j+1))/sqrt(2)...
        *[kron(speye(pow2(j)), [1 1]);...
          kron(speye(pow2(j)), [1 -1])];
end
