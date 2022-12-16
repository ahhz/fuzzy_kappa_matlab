% =======================================================================
% Copyright 2022
% Author: Alex Hagen-Zanker
% University of Surrey
%
% Distributed under the MIT Licence (http://opensource.org/licenses/MIT)
% =======================================================================
%
function  [fks, sim, P, E] = fuzzy_kappa_simulation(A_before, A_after, B_before, B_after, Mask, M, f)
% FUZZY_KAPPA_SIMULATION Applies the Fuzzy Kappa Simulationmethod (Van
% Vliet et al, 2013)

Mask = logical(Mask);

[m,n] = size(M);
A = A_after * m + A_before;
B = B_after * n + B_before;

Msim = eye(m*m, n*n);

for ia = 1:m
    for ja = 1:m
        for ib = 1:n
            for jb = 1:n
                index_a = (ia-1) * m + ja;
                index_b = (ib-1) * n + jb;
                both_change = ia ~= ja && ib ~= jb;
                both_persist = ia == ja && ib == jb;
                if both_change || both_persist
                    Msim(index_a,index_b) = M(ja,jb);
                end
            end
        end
    end
end

[fks, sim, P, E] = fuzzy_kappa(A, B, Msim, Mask, f);

