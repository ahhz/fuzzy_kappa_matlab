function [fk, sim, P, E] = fuzzy_kappa(A, B, M, Mask, f)
% FUZZY_KAPPA Applies the Fuzzy Kappa method (Hagen-Zanker, 2009)
% [fk, sim, P, E] = fuzzy_kappa(A, B, M, Mask, f)
%
% =======================================================================
% Copyright 2022
% Author: Alex Hagen-Zanker
% University of Surrey
%
% Distributed under the MIT Licence (http://opensource.org/licenses/MIT)
% =======================================================================
%
% OUTPUTS
% fk is the fuzzy kappa statistic
% sim is the similarity raster layer 
% P is the mean agreement
% E is the expected agreement.
%
% INPUTS
% A is first categorical raster layer
% B is second categorical raster layer
% Mask is binary raster layer to delineate the study area 
% M is category similarity matrix
% f is a distance decay function
%
% PRECONDITIONS
% Classes in A and B are coded as integers, consecutively numbered and 
% starting at 0. It is allowed for classes to be fully absent in the maps.
% Values in M must be in the range [0,1]; The number of rows must be
% equal to the number of classes in A; The number of columns must be equal 
% to the number of classes in B.
% The argument f must be function handle for a unary function with inputs 
% in the range [0,inf] with f(0) = 1 and f(inf) = 0, and f(x) <= f(y) for 
% all x > y
% 
% POSTCONDITIONS
% P in range [0,1]
% E in range [0,1]
% fk in range [-1, 1]
% sim of the same dimension as A, B, Mask; Values inside study area in 
% range [0,1]; Values outside study area NaN
% Corner case: The study area is empty: P = E = fk = NaN, sim = NaN matrix 
% Corner case: Both maps are uniform and same class: P = E = 1, fk = NaN
% 
% REFERENCE
% Hagen‐Zanker, A., 2009. An improved Fuzzy Kappa statistic that accounts 
% for spatial autocorrelation. International Journal of Geographical 
% Information Science, 23(1), pp.61-73.

A = A + 1;
B = B + 1;
Mask = logical(Mask);
[m,n] = size(M);
r = nnz(Mask);

if r == 0
    P = NaN; E = NaN; fk = NaN;
    sim = NaN(size(Mask));
    return
end
dA = arrayfun(@(i)f(bwdist(A==i & Mask)), 1:m,'UniformOutput',false);
dB = arrayfun(@(j)f(bwdist(B==j & Mask)), 1:n,'UniformOutput',false);

cellmax = @(c)max(cat(3,c{:}),[],3);
muA = arrayfun(@(j)cellmax(cellfun(@times, dA, mat2cell(M(:,j) ,...
    ones(m,1))','UniformOutput',false)), 1:n,'UniformOutput',false);
muB = arrayfun(@(i)cellmax(cellfun(@times, dB, mat2cell(M(i,:)',...
    ones(n,1))','UniformOutput',false)), 1:m,'UniformOutput',false);

cellsum = @(c)sum(cat(3,c{:}),3);
simA = cellsum(arrayfun(@(j) (B==j) .* muA{j}, 1:n,'UniformOutput',false));
simB = cellsum(arrayfun(@(i) (A==i) .* muB{i}, 1:m,'UniformOutput',false));

sim = min(simA,simB);
sim(~Mask) = nan;

E = 0;
for i = 1:m
    for j = 1:n
        xA = muA{j}(A==i & Mask);
        xB = muB{i}(B==j & Mask);
        x = flip(unique([xA; xB]));
        cA = histcounts(-xA, [-x; inf]);
        cB = histcounts(-xB, [-x; inf]);
        p = diff([0, cumsum(cA) .* cumsum(cB)]);
        E = E +  p * x; 
    end
end
E = E / r^2;
P = mean(sim(Mask),'all');
if(E == 1)
    fk = NaN;
else
    fk =  (P - E) / (1 - E);
end