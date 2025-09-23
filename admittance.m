function Y = admittance(nfrom, nto, r, x, b)
% ADMITTANCE Calculate the admittance matrix for an AC network
% 
% Inputs:
%   nfrom - vector of starting nodes for each branch
%   nto   - vector of ending nodes for each branch  
%   r     - vector of line resistances (pu)
%   x     - vector of line reactances (pu)
%   b     - vector of line shunt susceptances (pu)
%
% Output:
%   Y     - nodal admittance matrix

    % Determine number of nodes (find maximum node number)
    nbus = max([nfrom; nto]);
    
    % Initialize admittance matrix
    Y = zeros(nbus, nbus) + 1i*zeros(nbus, nbus);
    
    % Process each branch
    for k = 1:length(nfrom)
        % Calculate series admittance
        z = r(k) + 1i*x(k);  % series impedance
        y_series = 1/z;       % series admittance
        
        % Get node indices
        i = nfrom(k);
        j = nto(k);
        
        % Add series admittance to off-diagonals
        Y(i,j) = Y(i,j) - y_series;
        Y(j,i) = Y(j,i) - y_series;
        
        % Add series admittance to diagonals
        Y(i,i) = Y(i,i) + y_series;
        Y(j,j) = Y(j,j) + y_series;
        
        % Add shunt admittance (half at each end)
        y_shunt = 1i*b(k)/2;
        Y(i,i) = Y(i,i) + y_shunt;
        Y(j,j) = Y(j,j) + y_shunt;
    end
end
