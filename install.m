% Setup script for MS-Lib
%
% This file is part of MS-Lib, the multiscale library.
% Copyright (c) 2011-2012 Thorsten Raasch
%
% MS-Lib is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% MS-Lib is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with MS-Lib.  If not, see <http://www.gnu.org/licenses/>.

% add some subdirectories to the Matlab search path
mfilepath=fileparts(which(mfilename));
addpath(fullfile(mfilepath,'dwt'))
addpath(fullfile(mfilepath,'plot'))
addpath(fullfile(mfilepath,'refinement_masks'))
addpath(fullfile(mfilepath,'unit_tests'))
addpath(fullfile(mfilepath,'utils'))
clear mfilepath;
