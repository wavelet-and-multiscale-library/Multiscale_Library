function [a,astart] = DM_mask(a1,astart1,a2,astart2)
%DM_mask() refinement mask of correlation function between two refinable functions
%   DM_mask(a1,astart1,a2,astart2) computes the refinement mask of
%     F(x)=int_{-infty}^infty phi1(y) phi2(y-x) dy
%   which is given, see Dahmen/Micchelli [1], by
%     a_k=sum_l(a1_l*a2_{l-k})/2
%   i.e. by discrete convolution. This refinement mask can be used, e.g.,
%   to compute the L_2 inner products of phi1 against integer translates
%   of phi2 up to machine precision.
%
%   INPUT
%   a1,astart1: refinement mask of phi1 and its start index
%   a2,astart2: refinement mask of phi2 and its start index
%
%   OUTPUT
%   a: refinement mask of F
%   astart: start index of the refinement mask, astart=astart1-aend2
%
%   REFERENCES
%   [1] W. Dahmen and C. Micchelli,
%       Using the refinement equation for evaluating integrals of wavelets,
%       SIAM J. Numer. Anal. 30(1993), 507-537

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

a=conv(a1,fliplr(a2));

if nargout==2,
    astart=astart1-(astart2+length(a2)-1);
end

