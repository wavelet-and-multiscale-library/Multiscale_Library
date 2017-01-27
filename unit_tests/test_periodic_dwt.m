function result = test_periodic_dwt()
%test_periodic_dwt() test periodic discrete wavelet transform
%   test_periodic_dwt() is one of the unit tests in MS-Lib.
%   We check whether the builtin routines for the discrete transform
%   routines work properly for different periodic wavelet bases.
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

more off

disp('MS-Lib unit test: check periodic discrete wavelet transform routines')

% check discrete Haar wavelet transform
disp('- checking discrete Haar wavelet transform...')
result=true;
for J=8:10,
    for j0=0:3,
        c=randn(1,pow2(J));
        d=haar_dwt(c,j0);
        c0=haar_idwt(d,j0);        
        result=result&(norm(c-c0,'inf')/norm(c,'inf')<1e-14);
        d=haar_idwt(c,j0);
        c0=haar_dwt(d,j0);
        result=result&(norm(c-c0,'inf')/norm(c,'inf')<1e-14);
    end
end
if result,
    disp('  ... ok!');
else
    disp('  ... ERROR!');
end

% check some periodic discrete wavelet transform matrices M for (bi)orthogonality
disp('- checking some periodic discrete wavelet transform matrices M_j for orthogonality...')
result=true;
for m=1:4,
    a=D_mask(m);
    for j=4:6,
        M=Mj_per(a,j);
        result=result&(norm(speye(pow2(j))-M*M','inf')<1e-14);
    end
end
if result,
    disp('  ... ok!');
else
    disp('  ... ERROR!');
end
disp('- checking some periodic discrete wavelet transform matrices M_j for biorthogonality...')
result=true;
for m=1:4,
    [a,astart]=N_mask(m);
    for mt=m:2:m+4,
        [at,atstart]=CDF_mask(m,mt);        
        for j=5:6,
            M=Mj_per(a,astart,at,atstart,j);
            Mt=Mj_per(at,atstart,a,astart,j);            
            result=result&(norm(speye(pow2(j))-M*Mt','inf')<1e-14);
        end
    end
end
if result,
    disp('  ... ok!');
else
    disp('  ... ERROR!');
end

% check some periodic discrete wavelet transforms
disp('- checking some periodic discrete wavelet transforms...')
result=true;
for m=1:4,
    a=D_mask(m);    
    for J=6:8,
        for j0=j0_per(a):j0_per(a)+2,
            c=randn(1,pow2(J));
            d=periodic_dwt(c,a,j0);
            c0=periodic_idwt(d,a,j0);        
            result=result&(norm(c-c0,'inf')/norm(c,'inf')<1e-14);
            d=periodic_idwt(c,a,j0);
            c0=periodic_dwt(d,a,j0);
            result=result&(norm(c-c0,'inf')/norm(c,'inf')<1e-14);
        end
    end
end
if result,
    disp('  ... ok!');
else
    disp('  ... ERROR!');
end
