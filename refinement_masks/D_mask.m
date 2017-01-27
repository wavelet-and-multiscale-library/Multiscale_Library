function [a,astart] = D_mask(m)
%D_mask(m) refinement mask of a Daubechies scaling function
%   D_mask(m) computes the refinement mask of the orthogonal scaling
%   functions with approximation order m and minimal support [0,2*m-1],
%   as constructed by I. Daubechies [1], see also Shann/Yen [2],[3].
%   The mask starts at astart=0.
%
%   INPUT
%   m: order of Daubechies scaling function, 1<=m<=4
%
%   OUTPUT
%   a: refinement mask
%   astart: start index of the refinement mask, astart=0
%
%   REFERENCES
%   [1] I. Daubechies, Orthonormal bases of compactly supported wavelets,
%       Commun. Pure Appl. Math. 41(1988), 909-996
%   [2] W.-C. Shann: Exact solutions for Daubechies orthonormal scaling
%       coefficients, Department of Mathmematics, National Central
%       University, 1997
%   [3] W.-C. Shann/C.-C. Yen: On the exact values of orthonormal scaling
%       coefficients of lengths 8 and 10,
%       Appl. Comput. Harmon. Anal. 6(1999), 109-112

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

switch m
    case 1
        a=[1 1];
    case 2
        a=[1+sqrt(3),3+sqrt(3),3-sqrt(3),1-sqrt(3)]/4;
    case 3
        theta=sqrt(5+2*sqrt(10));
        a=[(1+sqrt(10)+theta)/16,...
           (5+sqrt(10)+3*theta)/16,...
           (5-sqrt(10)+theta)/8,...
           (5-sqrt(10)-theta)/8,...
           (5+sqrt(10)-3*theta)/16,...
           (1+sqrt(10)-theta)/16];
    case 4
        u=700+210*sqrt(15);
        v=(28+6*sqrt(35))*u^(1/3)-70*2^(1/3)+(2*u)^(2/3);
        s=(56*u^(1/3)*sqrt(v)+70*2^(1/3)*sqrt(v)-(2*u)^(2/3)*sqrt(v)...
            +12*u^(1/3)*sqrt(35*v)+336*sqrt(3*u)+48*sqrt(105*u))/(u^(1/3)*sqrt(v));
        alpha=0.5+0.5*sqrt(35)+sqrt(3*v/u^(1/3))/6+sqrt(3*s)/6;
        a=1/32*conv([1,4,6,4,1],[alpha,2-sqrt(140)+5/alpha,2+sqrt(140)-alpha,-5/alpha]);
    otherwise
        error('illegal order m');     
end

if nargout==2,
    astart=0;
end
