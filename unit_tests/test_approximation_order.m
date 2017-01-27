function result = test_approximation_order()
%test_approximation_order() test approximation orders of builtin refinement masks
%   test_approximation_order() is one of the unit tests in MS-Lib.
%   We check whether the builtin univariate refinement masks have
%   the correct approximation order.
%   
%   OUTPUT
%   result: true or false

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

disp('MS-Lib unit test: check approximation orders of builtin refinement masks')

% check approximation orders of m-th order B-splines N_m
disp('- checking approximation orders of m-th order B-splines...')
result=true;
for m=1:10,
    a=N_mask(m);
    k=compute_approximation_order(a);
    result=result&(k==m); % N_m has approximation order m
end
if result,
    disp('  ... ok!');
else
    disp('  ... ERROR!');
end

% check approximation orders of dual B-spline generators, CDF(m,mt)
disp('- checking approximation orders of dual B-spline generators, CDF(m,mt)...')
result=true;
for m=1:4,
    for mt=m:2:m+4,
        a=CDF_mask(m,mt);
        k=compute_approximation_order(a);
        result=result&(k==mt); % CDF(m,mt) scaling functions have order mt
    end
end
if result,
    disp('  ... ok!');
else
    disp('  ... ERROR!');
end

% check approximation orders of m-th order Daubechies scaling functions
disp('- checking approximation orders of m-th order Daubechies scaling functions...')
result=true;
for m=1:4,
    a=D_mask(m);
    k=compute_approximation_order(a);
    result=result&(k==m); % Daubechies scaling functions have order m
end
if result,
    disp('  ... ok!');
else
    disp('  ... ERROR!');
end


function order=compute_approximation_order(a)
% computes approximation order of a mask, by checking whether the
% Strang-Fix conditions
%   p(0)=1, p^{(mu)}(pi)=0, 1<=mu<=k-1
% hold for some k, with
%   p(x)=sum_k (a_k*exp(-i*k*x))/2.
% Note that
%   p^{(mu)}(x)=sum_k (a_k*(-i)^mu*k^mu*exp(-i*k*x))/2
order=0;
if abs(polyval(a/2,1)-1)<1e-15,
    order=1;
    while true,
        if abs(polyval(a/2.*[1:length(a)].^order.*(-i)^order,-1))<1e-14,
            order=order+1;
        else
            return
        end
    end
end
