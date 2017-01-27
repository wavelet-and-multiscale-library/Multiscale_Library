function [a,astart] = N_mask(m)
%N_mask(m) refinement mask of a cardinal B-spline of order m
%   N_mask(m) computes the refinement mask of a cardinal B-spline N_m of
%   order m. The B-spline is centered at 0 for even and at 1/2 for odd
%   values of m. Therefore, the mask starts at astart=-floor(m/2).
%
%   INPUT
%   m: spline order (m=1: step function, m=2: hat function, ...)
%
%   OUTPUT
%   a: refinement mask of N_m
%   astart: start index of the refinement mask, astart=-floor(m/2)

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

a=pow2(1-m)*[1 1];
for k=2:m,
    a=conv(a,[1 1]);
end

if nargout==2,
    astart=-floor(m/2);
end
