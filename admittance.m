function Y = admittance(nfrom, nto, r, x, b)
% ADMITTANCE Calculate the admittance matrix for an AC power network.
%
%   Y = ADMITTANCE(NFROM, NTO, R, X, B) computes the nodal admittance matrix
%   for a power system network using the incidence matrix approach.
%
%   Inputs:
%     nfrom - Mx1 vector of branch starting nodes
%     nto   - Mx1 vector of branch ending nodes  
%     r     - Mx1 vector of branch resistances (p.u.)
%     x     - Mx1 vector of branch reactances (p.u.)
%     b     - Mx1 vector of branch shunt susceptances (p.u.)
%
%   Output:
%     Y     - NxN nodal admittance matrix (p.u.)
%
%   The function follows these steps:
%     1. Determines the number of nodes (N) from the connectivity data
%     2. Calculates branch admittances from series R and X
%     3. Builds the incidence matrix A
%     4. Constructs the branch admittance matrix Yb
%     5. Computes Y = A * Yb * A' + Yshunt
%     6. Removes the reference node (node with highest number) to make Y invertible
%
%   Example:
%     Y = admittance(nfrom, nto, r, x, b);

    % Determine total number of nodes
    all_nodes = unique([nfrom; nto]);
    N = max(all_nodes);
    M = length(nfrom);
    
    % Validate input dimensions
    if length(nto) ~= M || length(r) ~= M || length(x) ~= M || length(b) ~= M
        error('All input vectors must have the same length');
    end
    
    % Calculate series admittance for each branch: y_series = 1/(r + jx)
    y_series = 1./(r + 1i*x);
    
    % Initialize incidence matrix A (N x M)
    A = zeros(N, M);
    
    % Build incidence matrix
    for m = 1:M
        i = nfrom(m);
        j = nto(m);
        A(i, m) = 1;      % Current leaves node i
        A(j, m) = -1;     % Current enters node j
    end
    
    % Create branch admittance matrix Yb (M x M diagonal)
    Yb = diag(y_series);
    
    % Calculate preliminary admittance matrix
    Y_prelim = A * Yb * A';
    
    % Add shunt admittances (half at each end of the branch)
    Y_shunt = zeros(N, N);
    for m = 1:M
        i = nfrom(m);
        j = nto(m);
        shunt_admittance = 1i * b(m) / 2;  % Half susceptance at each end
        
        Y_shunt(i, i) = Y_shunt(i, i) + shunt_admittance;
        Y_shunt(j, j) = Y_shunt(j, j) + shunt_admittance;
    end
    
    % Final admittance matrix
    Y_full = Y_prelim + Y_shunt;
    
    % Remove reference node (highest numbered node) to make matrix invertible
    % In power systems, we typically use node N as reference (V_N = 0)
    ref_node = N;
    non_ref_nodes = setdiff(1:N, ref_node);
    Y = Y_full(non_ref_nodes, non_ref_nodes);
end
