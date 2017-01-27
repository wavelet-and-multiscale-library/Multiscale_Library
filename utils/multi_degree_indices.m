function I = multi_degree_indices(d,k)
%multi_degree_indices() compute all multi-indices with a given degree
%   multi_degree_indices(d,k) yields a matrix which contains
%   in its columns all multi-indices nu in dimension d whose
%   degree |k| is less or equal to k
%
%   INPUT
%   d: dimension of the multi-indices
%   k: degree bound of the multi-indices (sum of all entries)
%
%   OUTPUT
%   I: matrix of multi-indices (columnwise)

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

if d<=0,
    error('dimension of multi-indices must be positive')
end

% We iterate over the dimension to obtain all multi-indices below a certain degree.
% The only multi-indices in dimension 1 with degree <=k are 0:k.
I=[0:k];

for dd=2:d,
    % multi-indices in dimension dd of degree kk<=k decompose into a
    % multi-index of degree kkk<=kk and a last entry kk-kkk
    J=[];
    for jj=1:size(I,2),
        kkk=sum(I(:,jj));
        J=[J [kron(ones(1,k-kkk+1),I(:,jj)); 0:k-kkk]];
    end
    I=J;
end
