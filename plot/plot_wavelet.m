function [x,y] = plot_wavelet(a,astart,at,atstart,j)
%plot_wavelet() evaluate (bi)orthogonal wavelet on 2^{-j}Z
%   plot_wavelet(a,astart,j) evaluates an orthogonal wavelet psi
%   with refinement mask a starting at astart on a dyadic grid 2^{-j}Z
%
%   plot_wavelet(a,astart,at,atstart,j) evaluates a biorthogonal wavelet psi
%   with primal refinement mask a starting at astart and dual refinement
%   mask at starting at atstart on a dyadic grid 2^{-j}Z.
%
%   INPUT
%   a: primal refinement mask (row vector) [a_astart,...,a_{astart+L}], L\ge 1
%   at: dual refinement mask (row vector) [at_atstart,...,at_{atstart+Lt}], Lt\ge 1
%   j: desired dyadic resolution
%
%   OUTPUT
%   x: (astart+1-atstart-Lt)/2:2^{-j}:(astart+L+1-atstart)/2
%   y: phi(x)

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

if nargin==3,
    % orthogonal case
    j=at; % since j should be the last argument
    at=a;
    atstart=astart;
end

L=size(a,2)-1;
if L<=0,
   error('mask a is too short')
end
aend=astart+L;

Lt=size(at,2)-1;
if Lt<=0,
    error('mask at is too short')
end
atend=atstart+Lt;

b=zeros(size(at));
for k=1-Lt:1,
    b(k+Lt)=(-1)^k*at(2-k);
end
bstart=1-atend;
bend=1-atstart;

if L==1,
    % special case of Haar-type (piecewise constant) wavelets, i.e., CDF(1,*)
    % psi is piecewise constant w.r.t. Z/2
    if nargout==0,
        stairs([(astart+bstart):(aend+bend)]/2,[b 0]);
    else
        x=(astart+bstart)/2:pow2(-j):(aend+bend)/2;
        y=[kron(bk,ones(1,pow2(j-1))) 0];
    end
else
    % first we evaluate phi at resolution j-1
    [xphi,yphi]=plot_refinable_function(a,astart,j-1);

    % psi is a linear combination of phi(2.-k),
    % k running from bstart to bend, i.e., psi has support length
    %   (aend+bend-astart-bstart)/2 = (L+Lt)/2
    y=zeros(1,pow2(j-1)*(L+Lt)+1);
    for k=1-Lt:1,
        y=y+(-1)^k*at(2-k)*[zeros(1,(k-1+Lt)*pow2(j-1)) yphi zeros(1,(1-k)*pow2(j-1))];
    end

    if nargout==0,
        plot((astart+bstart)/2:pow2(-j):(aend+bend)/2,y);
    else
        x=(astart+bstart)/2:pow2(-j):(aend+bend)/2;
    end
end
