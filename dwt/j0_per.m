function j0 = j0_per(a,at)
%j0_per() compute minimal level for discrete periodic wavelet transform
%
%   j0_per(a,at) computes the minimally possible level for the
%   decomposition steps within the discrete periodic wavelet transform
%   with masks a and at.
%   j0 is chosen such that in the corresponding transform matrices,
%   the mask entries a_k and b_k=(-1)^k at_{1-k} do not overlap.
%   In particular, j0 is chosen via 2^(j0+1)>=max(length(a),length(at)).
%   For the Haar mask N_mask(1)=[1 1], e.g., this means that j0=0.
%   Other masks will require that j0>0.
%
%   j0_per(a) computes the minimally possible level in the orthogonal case a=at.
%
%   INPUT
%   a: primal refinement mask (row vector) [a_astart,...,a_{astart+L}], L\ge 2
%   at: dual refinement mask (row vector) [at_atstart,...,at_{atstart+Lt}], Lt\ge 2
%
%   OUTPUT
%   j0: minimal level
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

% compute j0 such that 2^j0>=max(length(a),length(at)) and j0 is minimal
if nargin==1,
    j0=nextpow2(length(a))-1;
else
    j0=nextpow2(max(length(a),length(at)))-1;
end
