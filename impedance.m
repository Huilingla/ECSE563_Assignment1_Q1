function Z = impedance(nfrom, nto, r, x, b)
% IMPEDANCE Calculate the impedance matrix for an AC network
% 
% Inputs:
%   nfrom - vector of starting nodes for each branch
%   nto   - vector of ending nodes for each branch  
%   r     - vector of line resistances (pu)
%   x     - vector of line reactances (pu)
%   b     - vector of line shunt susceptances (pu)
%
% Output:
%   Z     - nodal impedance matrix

    % Calculate admittance matrix using previous function
    Y = admittance(nfrom, nto, r, x, b);
    
    % Number of buses
    nbus = size(Y, 1);
    
    % Create identity matrix for solving YZ = I
    I_matrix = eye(nbus);
    
    % Solve for impedance matrix using linsolve()
    Z = linsolve(Y, I_matrix);
end
