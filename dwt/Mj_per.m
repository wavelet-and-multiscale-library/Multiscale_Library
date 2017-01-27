function M = Mj_per(a,astart,at,atstart,j)
%Mj_per() matrix representation of periodic discrete wavelet transform (decomposition)
%   
%   Mj_per(a,astart,at,atstart,j) computes the (sparse) matrix M_j representing
%   one transformation step of the periodic discrete wavelet transform
%   (decomposition) as a basis transformation within V_j;
%   the inverse transform uses M_j' with reversed roles of a and at.
%   We assume that supp(a)={astart,...,astart+N}, supp(at)={atstart,...,atstart+Nt}
%
%   Mj_per(a,j) computes the matrix M_j in the orthogonal case a=at
%   and astart=atstart=0
%
%   INPUT
%   a: primal refinement mask (row vector) [a_astart,...,a_{astart+L}], L\ge 2
%   astart: index of first primal mask entry
%   at: dual refinement mask (row vector) [at_atstart,...,at_{atstart+Lt}], Lt\ge 2
%   atstart: index of first dual mask entry
%   j: level in which to transform
%
%   OUTPUT
%   M: transformation matrix (sparse format)
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

if nargin==2,
    j=astart;
    at=a;
    astart=0;
    atstart=0;
end

N=size(a,2)-1;
Nt=size(at,2)-1;
M=sparse(pow2(j),pow2(j));

% abort if level is too small for M to hold the mask
if (max(N,Nt)+1>pow2(j)),
    error('j too small!');
end

% setup first row of M
abuf=[a,sparse(1,pow2(j)-N-1)];
abuf=[abuf abuf abuf]; % provides three copies of the mask
% the first row of M is a subvector of abuf...
M(1,:)=abuf(pow2(j)+1-astart:pow2(j+1)-astart);
% ..., the other rows are circulant copies thereof
for row=2:pow2(j-1),
    M(row,:)=circshift(M(row-1,:),[0 2]);
end

% setup (2^(j-1)+1)-th row of M, using
%   b_k = (-1)^k at_{1-k}
%       = -(-1)^{1-k} at_{1-k}
bbuf=[fliplr((-1).^[atstart+1:atstart+Nt+1].*at),sparse(1,pow2(j)-Nt-1)];
bbuf=[bbuf bbuf bbuf]; % provides three copies of the wavelet expansion coefficients
% the (2^(j-1)+1)-th row of M is a subvector of bbuf...
offset=pow2(j)+Nt+atstart;
M(pow2(j-1)+1,:)=bbuf(offset:offset+pow2(j)-1);
% ..., the other rows are circulant copies thereof
for row=pow2(j-1)+2:pow2(j),
    M(row,:)=circshift(M(row-1,:),[0 2]);
end

M=M./sqrt(2);
