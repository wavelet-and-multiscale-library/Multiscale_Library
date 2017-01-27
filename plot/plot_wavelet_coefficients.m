function plot_wavelet_coefficients(d,j,j0)
%plot_wavelet_coefficients() visualize 1D wavelet coefficient distribution
%   plot_wavelet_coefficients() visualizes a 1D wavelet coefficient array d
%   by a 2-dimensional plot. A coefficient d_{j,k} is visualized by a
%   filled rectangle at (j,2^{-j}k) with width 2^{-j} and height 1. The
%   absolute value of d_{j,k} determines the fill colour in a logarithmic scale.
%
%   INPUT
%   d: wavelet coefficients in V_j, row vector, length approximately 2^j
%   j: finest level of wavelet coefficients
%   j0: coarsest level, 0<=j0<=j, default value 0

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

if nargin<3, j0=0; end
if nargin<2, j=log2(length(d)); end

% deduce number of generator coefficients on coarsest scale
dimVj0=pow2(j0)+length(d)-pow2(j);

% compute maximal absolute value
maxabs=max(abs(d));

% draw the support boxes
i=1;
for k=0:dimVj0-1, % generator coefficients
    patch([k/dimVj0;(k+1)/dimVj0;(k+1)/dimVj0;k/dimVj0],...
        [j0-0.5;j0-0.5;j0+0.5;j0+0.5],...
        log10(abs(d(i))/maxabs));
    i=i+1;
end
for l=j0:j-1, % wavelet coefficients
    for k=0:pow2(l)-1,
        patch([pow2(-l)*k;pow2(-l)*(k+1);pow2(-l)*(k+1);pow2(-l)*k],...
            [l+0.5;l+0.5;l+1.5;l+1.5],...
            log10(abs(d(i))/maxabs));
        i=i+1;
    end;
end

box on;
set(gca,'Layer','top');
set(gca,'YLim',[j0-0.5 j+0.5]);
xlabel '2^{-j}k';
ylabel 'level j (*: generators)';
set(gca,'XTick',[]);
set(gca,'YTick',[j0:j]);
set(gca,'YTickLabel',{[int2str(j0),'*'],j0:j-1});
set(gca,'CLim',[-15 0]);
colorbar;
