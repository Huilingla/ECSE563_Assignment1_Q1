% SOLVE_IEEE9 Solve the IEEE 9-bus system using nodal analysis
%
% This script:
%   1. Loads the IEEE 9-bus test system data
%   2. Computes the admittance matrix using the admittance function
%   3. Solves for node voltages using the internal current sources
%   4. Displays results in polar coordinates

% Load the IEEE 9-bus system data
ieee9_A1;

% Calculate admittance matrix
Y = admittance(nfrom, nto, r, x, b);

% Determine system dimensions
all_nodes = unique([nfrom; nto]);
N_total = max(all_nodes);
N = N_total - 1; % After removing reference node

% Extract current injections for non-reference nodes
% The reference node is node N_total (highest numbered node)
ref_node = N_total;
non_ref_nodes = setdiff(1:N_total, ref_node);

% Internal current sources (given in the data file)
% We need to extract only the non-reference nodes
I_non_ref = Iint(non_ref_nodes);

% Solve for node voltages: Y * V = I
% Use linsolve() instead of matrix inversion
V_non_ref = linsolve(Y, I_non_ref);

% Create complete voltage vector with reference node voltage = 0
V_complete = zeros(N_total, 1);
V_complete(non_ref_nodes) = V_non_ref;
V_complete(ref_node) = 0; % Reference node voltage

% Convert to polar coordinates
V_mag = abs(V_complete);
V_angle_deg = angle(V_complete) * 180/pi;

% Display results
fprintf('IEEE 9-Bus System Node Voltages\n');
fprintf('================================\n\n');
fprintf('Node\tMagnitude (p.u.)\tAngle (degrees)\n');
fprintf('----\t---------------\t-------------\n');

for i = 1:N_total
    fprintf('%2d\t%12.4f\t%12.2f\n', i, V_mag(i), V_angle_deg(i));
end

% Additional validation: Calculate power mismatch
I_calculated = Y * V_non_ref;
mismatch = norm(I_calculated - I_non_ref);
fprintf('\nCurrent injection mismatch norm: %e\n', mismatch);
