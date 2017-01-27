function y = multi_factorial(alpha)
%multi_factorial() factorial of a multi-index
%   multi_factorial(alpha) yields the factorial of the multi-index alpha,
%   being defined as
%
%       alpha! = alpha_1! * ... * alpha_d!
%
%   INPUT
%   alpha: multi-index of dimension d
%
%   OUTPUT
%   y: factorial of alpha

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

y=prod(factorial(alpha));
