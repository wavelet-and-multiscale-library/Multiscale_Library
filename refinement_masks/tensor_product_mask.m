function [a,astart] = tensor_product_mask(a1,astart1,a2,astart2)
%tensor_product_mask() refinement mask of tensor product between two univariate refinable functions
%   tensor_product_mask(a1,astart1,a2,astart2) takes the refinement masks of
%   two univariate refinable functions phi1, phi2 and computes the
%   2D refinement mask of
%     phi(x,y)=phi1(x)*phi2(y)
%   as a matrix with a 2D start index.
%   We use the convention that column numbers correspond to the x-direction,
%   whereas row numbers correspond to the y-direction.
%
%   INPUT
%   a1,astart1: refinement mask of phi1 and its start index
%   a2,astart2: refinement mask of phi2 and its start index
%
%   OUTPUT
%   a: refinement mask of phi
%   astart: start index of the refinement mask

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

a=a2'*a1;

if nargout==2,
    astart=[astart1 astart2];
end
