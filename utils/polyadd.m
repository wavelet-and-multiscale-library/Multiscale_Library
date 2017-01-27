function p=polyadd(p1,p2)
%polyadd(p1,p2) adds two polynomials of possibly different degree
%   polyadd(p1,p2) computes the coefficients of the polynomial p1+p2,
%   where p1,p2 are polynomials of possibly different degree.
%
%   INPUT
%   p1,p2: coefficients of the polynomials p1, p2
%   
%   OUTPUT
%   p: coefficients of p1+p2

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

if length(p1)<length(p2),
    p=p2+[zeros(1,length(p2)-length(p1)) p1];
else
    p=p1+[zeros(1,length(p1)-length(p2)) p2];
end
