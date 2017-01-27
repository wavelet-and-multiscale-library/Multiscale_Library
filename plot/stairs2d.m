function stairs2d(x,y,C)
%stairs2d() bivariate stairstep plot
%   stairs2d(x,y,C) yields a bivariate stairstep graph of the elements of
%   the m-by-n matrix C, where the coordinates of the stairsteps are given
%   by the vectors x and y of length n+1 and m+1, respectively.
%
%   stairs2d(C) yields a bivariate stairstep graph of the elements of
%   the m-by-n matrix C, where the stairsteps have integer coordinates
%   and start at (0,0).
%
%   INPUT
%   C: m-by-n matrix of z values of a 2D function
%      (row number <-> y coordinate, column number <-> x coordinate)
%
%   OUTPUT
%   no output

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

if nargin==1,
    C=x;
end
m = size(C,1);
n = size(C,2);

if nargin==1,
    x=0:n;
    y=0:m;
end

assert(length(x)==n+1,'length of input vector x must be equal to the number of columns+1');
assert(length(y)==m+1,'length of input vector y must be equal to the number of rows+1');
    
% create 3D staircase plot from single 3D patches
cmin=min(min(C));
clf;
for k1=1:n,
    for k2=1:m,
        xk=x(k1:k1+1);
        yk=y(k2:k2+1);
        zk=[cmin C(k2,k1)];
        FV.vertices= ...
         [xk([1 2 2 1 1 2 2 1]); ...
          yk([1 1 2 2 1 1 2 2]); ...
          zk([1 1 1 1 2 2 2 2])]';
        FV.faces=[1 2 3 4; 5 6 7 8; 1 2 6 5; 4 3 7 8; 1 5 8 4; 2 6 7 3];
        FV.facecolor=[1 1 1];
        patch(FV)
    end
end
axis tight
axis equal
xlabel('x');
ylabel('y');
