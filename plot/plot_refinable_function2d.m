function [X,Y,Z] = plot_refinable_function2d(A,Astart,j,mu)
%plot_refinable_function2d() evaluate (partial derivative of) 2d refinable function on 2^{-j}Z^2
%   plot_refinable_function(A) evaluates the bivariate refinable function phi
%   with refinement mask A and support [0,length(a)+1]
%   on the integers, up to machine precision
%
%   INPUT
%   A: refinement mask ((M+1)-by-(N+1) matrix), M,N>=1
%   Astart: start indices of the mask, integer row vector of length 2
%   j: desired dyadic resolution
%   mu: derivative order, row vector, mu>=(0,0)
%
%   OUTPUT
%   X,Y: x,y coordinates of the Z values, ready for surf()
%   Z: point values of mu-th derivative of phi

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

M=size(A,1)-1;
N=size(A,2)-1;
if min(M,N)<=0,
   error('mask a is too short')
end

if nargin==1,
    Astart=[0 0];
    j=0;
    mu=[0 0];
end

if nargin==2,
    j=0;
    mu=[0 0];
end

if nargin==3,
    mu=[0 0];
end

if min(M,N)==1,
    % exclude Haar case, we assume that then mu=0
    if nargout==0,
        stairs2d(Astart(1)+[0:1],Astart(2)+[0:1],1)
    else
        X=kron(ones(pow2(j)+1,1),Astart(1)+[0:pow2(-j):1]);
        Y=kron(Astart(2)+[0:pow2(-j):1]',ones(1,pow2(j)+1));
        Z=[ones(pow2(j)) zeros(pow2(j),1); zeros(1,pow2(j)) 0];
    end
else
    % First compute the values Z(alpha) of the mu-th derivative of phi
    % on the integers alpha in a support cube [1:N-1]x[1:M-1].
    % Z(alpha) is uniquely determined by the following two sets of
    % equations, see Dahmen/Micchelli 1992:
    %
    % a) eigenvalue conditions ((M-1)*(N-1) equations)
    %    2^{-|mu|}*Z(alpha) = sum_beta A_{2*alpha-beta}*Z(beta)
    %
    % b) orthogonality conditions (#{nu:|nu|<=|mu|} equations)
    %    sum_alpha (-alpha)^nu Z(alpha) = mu!*delta(mu,nu),|nu|<=|mu|
    %
    % We stack these two sets of conditions as a linear condition for
    %    z=reshape(Z',(M-1)*(N-1),1), with Z=reshape(z,N-1,M-1)',
    % so that    
    %    Z(m,n)=z_{(m-1)*(N-1)+n}, 1<=m<=M-1, 1<=n<=N-1

    z=zeros((M-1)*(N-1),1);
    
    B=zeros((M-1)*(N-1)); % start with a square matrix holding the refinement coeffs
    % B_{alpha,beta} = A_{2*alpha-beta}-2^{-|mu|}delta(alpha,beta)
    for alpham=1:M-1,
        for alphan=1:N-1,
            % alpha=[alphan alpham]
            for betam=1:M-1,
                for betan=1:N-1,
                    % beta=[betan betam]
                    % 2*alpha-beta=[2*alphan-betan 2*alpham-betam]
                    if(2*alphan-betan>=0 & 2*alphan-betan<=N & 2*alpham-betam>=0 & 2*alpham-betam<=M),
                        B((alpham-1)*(N-1)+alphan,(betam-1)*(N-1)+betan)=...
                        A(2*alpham-betam+1,2*alphan-betan+1);
                    end
                    if(alpham==betam & alphan==betan),
                        B((alpham-1)*(N-1)+alphan,(betam-1)*(N-1)+betan)=...
                        B((alpham-1)*(N-1)+alphan,(betam-1)*(N-1)+betan)-pow2(-sum(mu));
                    end
                end
            end
        end
    end
    % right-hand side r starts with zeros
    r=zeros((M-1)*(N-1),1);
    % stack B and r with orthogonality conditions
    nus=multi_degree_indices(2,sum(mu));
    for jj=1:size(nus,2),
        v=zeros(1,(M-1)*(N-1));
        for betam=1:M-1,
            for betan=1:N-1,
                % beta=[betan betam]
                v((betam-1)*(N-1)+betan)=(-betan)^nus(1,jj)*(-betam)^nus(2,jj); % =(-beta)^nu
            end
        end
        B=[B;v];
        if nus(:,jj)==mu',
            r=[r;multi_factorial(mu)];
        else
            r=[r;0];
        end
    end
    z=B\r;
    
    if nargout==0,
    else
        X=kron(ones(M*pow2(j)+1,1),Astart(1)+[0:pow2(-j):N]);
        Y=kron(Astart(2)+[0:pow2(-j):M]',ones(1,N*pow2(j)+1));
        Z=[zeros(1,N+1); zeros(M-1,1) reshape(z,N-1,M-1)' zeros(M-1,1); zeros(1,N+1)];
    end
end
