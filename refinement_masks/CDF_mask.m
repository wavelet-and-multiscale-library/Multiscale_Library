function [a,astart] = CDF_mask(m,mt)
%CDF_mask(m,mt) refinement mask of a dual B-spline generator
%   CDF_mask(m,mt) computes the refinement mask of a compactly supported
%   refinable function phi_tilde with approximation order mt
%   that is dual to the m-th order B-spline N_m and has minimal
%   support size m+2*mt-2, as constructed by Cohen, Daubechies and
%   Feauveau (1992). m+mt is even. The refinable function is centered
%   at 0 for even and at 1/2 for odd values of m. Therefore, the mask
%   starts at astart=-floor(m/2)-m+1.
%
%   INPUT
%   m: spline order, m>=1
%   mt: approximation order 
%
%   OUTPUT
%   a: refinement mask of phi_tilde
%   astart: start index of the refinement mask, astart=-floor(m/2)-m+1

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


if m<=0,
    error('m must be >=1')
end
if mt<=0,
    error('mt must be >=1')
end
if mod(m+mt,2)~=0,
    error('m+mt must be even')
end

uh2=[1 -2 1]; % (z-1)^2

K=(m+mt)/2;

a=[0];
uh2a=[1]; % 1
for k=0:K-1,
    zpower=[1 zeros(1,K-1-k)];
    a=polyadd(a,nchoosek(K-1+k,k)/(-4)^k*conv(zpower,uh2a));
    uh2a=conv(uh2a,uh2);
end
a=conv(a,N_mask(mt));

if nargout==2,
    astart=-floor(m/2)-mt+1;
end

