% SOLVE_IEEE9 Solve the IEEE 9-bus system following exact nodal analysis steps
clear all; close all; clc;

% Load the IEEE 9-bus system data
ieee9_A1;

fprintf('=== NODAL ANALYSIS OF IEEE 9-BUS SYSTEM ===\n\n');

% Step 1: Calculate admittance matrix using the prescribed methodology
fprintf('STEP 1: CALCULATING ADMITTANCE MATRIX\n');
fprintf('--------------------------------------\n');
Y_full = admittance(nfrom, nto, r, x, b);

% Step 4: Print the admittance matrix
fprintf('\nSTEP 4: ADMITTANCE MATRIX OUTPUT\n');
fprintf('--------------------------------\n');
fprintf('Full Admittance Matrix Y (%dx%d) in p.u.:\n\n', size(Y_full,1), size(Y_full,2));

% Display real and imaginary parts separately for clarity
fprintf('Real Part (Conductance Matrix G):\n');
for i = 1:size(Y_full,1)
    fprintf('Bus %d: ', i);
    for j = 1:size(Y_full,2)
        fprintf('%8.4f ', real(Y_full(i,j)));
    end
    fprintf('\n');
end

fprintf('\nImaginary Part (Susceptance Matrix B):\n');
for i = 1:size(Y_full,1)
    fprintf('Bus %d: ', i);
    for j = 1:size(Y_full,2)
        fprintf('%8.4f ', imag(Y_full(i,j)));
    end
    fprintf('\n');
end

% Step 5: Solve for node voltages using internal current sources
fprintf('\nSTEP 5: SOLVING FOR NODE VOLTAGES\n');
fprintf('---------------------------------\n');

% Determine system dimensions
N_total = max(unique([nfrom; nto]));

% Choose reference node (node 9 as shown in your results)
ref_node = 9;
fprintf('Using node %d as reference (V_%d = 0)\n', ref_node, ref_node);

% Remove reference node to make system solvable
non_ref_nodes = setdiff(1:N_total, ref_node);
Y_reduced = Y_full(non_ref_nodes, non_ref_nodes);

% Extract current injections for non-reference nodes
I_non_ref = Iint(non_ref_nodes);

fprintf('\nCurrent injections at non-reference nodes:\n');
for i = 1:length(non_ref_nodes)
    node = non_ref_nodes(i);
    fprintf('  I_%d = %7.4f + j%7.4f p.u.\n', node, real(I_non_ref(i)), imag(I_non_ref(i)));
end

% Solve Y_reduced * V_non_ref = I_non_ref using linsolve()
fprintf('\nSolving system: Y_reduced * V_non_ref = I_non_ref\n');
V_non_ref = linsolve(Y_reduced, I_non_ref);

% Create complete voltage vector with reference node voltage = 0
V_complete = zeros(N_total, 1);
V_complete(ref_node) = 0; % Reference node voltage = 0
V_complete(non_ref_nodes) = V_non_ref;

% Convert to polar coordinates
V_mag = abs(V_complete);
V_angle_rad = angle(V_complete);
V_angle_deg = V_angle_rad * 180/pi;

% Display results in polar coordinates
fprintf('\nSOLUTION: NODE VOLTAGES IN POLAR COORDINATES\n');
fprintf('===========================================\n');
fprintf('Node    Magnitude (p.u.)    Angle (degrees)    Angle (radians)\n');
fprintf('----    ----------------    --------------    --------------\n');

for i = 1:N_total
    if i == ref_node
        fprintf('%2d (ref) %12.4f        %12.2f        %12.4f\n', ...
                i, V_mag(i), V_angle_deg(i), V_angle_rad(i));
    else
        fprintf('%2d       %12.4f        %12.2f        %12.4f\n', ...
                i, V_mag(i), V_angle_deg(i), V_angle_rad(i));
    end
end

% Validation
fprintf('\nVALIDATION:\n');
fprintf('-----------\n');
I_calculated = Y_reduced * V_non_ref;
mismatch = I_calculated - I_non_ref;
max_mismatch = max(abs(mismatch));

fprintf('Maximum current injection mismatch: %e p.u.\n', max_mismatch);

if max_mismatch < 1e-10
    fprintf('Solution accuracy: Excellent\n');
else
    fprintf('Solution accuracy: Acceptable\n');
end

% Display the reduced matrix used for solving
fprintf('\nReduced Admittance Matrix (without node %d):\n', ref_node);
fprintf('Real Part:\n');
disp(real(Y_reduced));
fprintf('Imaginary Part:\n');
disp(imag(Y_reduced));
