function [x,y] = plot_refinable_function(a,astart,j,mu)
%plot_refinable_function() evaluate (derivative of) 1d refinable function on 2^{-j}Z
%   plot_refinable_function(a) evaluates the refinable function phi
%   with refinement mask a and support [0,length(a)+1] on the integers,
%   up to machine precision
%
%   plot_refinable_function(a,astart) does the same, using a start index
%   astart of the mask and support [astart,astart+length(a)+1] of phi.
%
%   plot_refinable_function(a,astart,j) evaluates phi on the dyadic grid
%   2^{-j}Z.
%
%   plot_refinable_function(a,astart,j,mu) evaluates the mu-th derivative
%   of phi on the dyadic grid 2^{-j}Z
%
%   INPUT
%   a: refinement mask (row vector) [a_astart,...,a_{astart+L}], L>=1
%   astart: start index of mask, integer
%   j: desired dyadic resolution
%   mu: derivative order, mu>=0
%
%   OUTPUT
%   x: astart:2^{-j}:astart+L
%   y: phi^{(mu)}(x)

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

L=size(a,2)-1;
if L<=0,
   error('mask a is too short')
end

if nargin==1,
    astart=0;
    j=0;
    mu=0;
end

if nargin==2,
    j=0;
    mu=0;
end

if nargin==3,
    mu=0;
end

if L==1,
    % exclude Haar case, we assume that then mu=0
    if nargout==0,
        stairs([astart astart+1],[1 0]);
    else,
        x=astart+[0:pow2(j)]/pow2(j);
        y=[ones(1,pow2(j)) 0];
    end
else,
    % first compute the values of phi^{(mu)} on the integers
    A=zeros(L-1);
    for k=1:L-1,
      % A_{k,l}=a_{2k-l}
      for l=max((2*k-L),1):min((2*k),L-1),
        A(k,l)=a(2*k-l+1);
      end
    end; % loop could be sped up
    % solve extended eigenvalue problem
    V=rot90(vander(1:max([L-1,mu+1])));
    y=[0 ([pow2(mu)*A-eye(L-1);V(1:mu+1,1:L-1)]\[zeros(L+mu-1,1);(-1)^mu*factorial(mu)])' 0];

    % compute values of phi^{(mu)} on 2^{-j}Z inductively, using refinement relation
    for j1=1:j,
        % embed old values at odd positions
        y=[y(1) y(2:end)*kron(speye(size(y,2)-1),[0 1])];
        % compute values at odd points 2^{-j}(2k+1) via refinement relation
        for k=0:(pow2(j1-1)*L-1),
            result=0;
            for m=max(0,ceil(pow2(1-j1)*(2*k+1)-L)):min(L,floor(pow2(1-j1)*(2*k+1))),
                result=result+pow2(mu)*a(m+1)*y(2*(2*k+2-pow2(j1-1)*m)-1);
            end
            y(2*k+2)=result;
        end
    end

    if nargout==0,
        plot(astart+[0:L*pow2(j)]/pow2(j),y);
    else,
        x=astart+[0:L*pow2(j)]/pow2(j);
    end
end
