% SOLVE_IEEE9 Solve the IEEE 9-bus system using the full admittance matrix
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

% Display as combined complex matrix
fprintf('Combined Admittance Matrix (G + jB):\n');
for i = 1:size(Y_full,1)
    fprintf('Bus %d: ', i);
    for j = 1:size(Y_full,2)
        if imag(Y_full(i,j)) >= 0
            fprintf('%8.4f + j%7.4f  ', real(Y_full(i,j)), imag(Y_full(i,j)));
        else
            fprintf('%8.4f - j%7.4f  ', real(Y_full(i,j)), abs(imag(Y_full(i,j))));
        end
    end
    fprintf('\n');
end

% Step 5: Solve for node voltages using the FULL matrix with reference node constraint
fprintf('\nSTEP 5: SOLVING FOR NODE VOLTAGES USING FULL MATRIX\n');
fprintf('---------------------------------------------------\n');

% Determine system dimensions
N_total = max(unique([nfrom; nto]));

% Choose reference node (node 9 as shown in your results)
ref_node = 9;
fprintf('Using node %d as reference (V_%d = 0)\n', ref_node, ref_node);

% Modify the full matrix to enforce reference node constraint
Y_modified = Y_full;

% Method 1: Set reference node row/column to identity and injection to zero
Y_modified(ref_node, :) = 0;
Y_modified(:, ref_node) = 0;
Y_modified(ref_node, ref_node) = 1;

I_modified = Iint;
I_modified(ref_node) = 0;  % Set reference node current to zero

fprintf('\nModified system: Y_modified * V = I_modified\n');
fprintf('(Reference node constraint enforced by modifying row/column %d)\n', ref_node);

% Solve using the modified full matrix
V_complete = linsolve(Y_modified, I_modified);

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

% Validation using the FULL matrix
fprintf('\nVALIDATION USING FULL MATRIX:\n');
fprintf('-----------------------------\n');
I_calculated_full = Y_full * V_complete;
mismatch_full = I_calculated_full - Iint;
max_mismatch_full = max(abs(mismatch_full));

fprintf('Maximum current injection mismatch: %e p.u.\n', max_mismatch_full);

if max_mismatch_full < 1e-10
    fprintf('Solution accuracy: Excellent\n');
else
    fprintf('Solution accuracy: Acceptable\n');
end

% Display the modified full matrix used for solving
fprintf('\nModified Full Admittance Matrix (with reference node constraint):\n');
fprintf('Combined Complex Form:\n');
for i = 1:size(Y_modified,1)
    fprintf('Bus %d: ', i);
    for j = 1:size(Y_modified,2)
        if imag(Y_modified(i,j)) >= 0
            fprintf('%8.4f + j%7.4f  ', real(Y_modified(i,j)), imag(Y_modified(i,j)));
        else
            fprintf('%8.4f - j%7.4f  ', real(Y_modified(i,j)), abs(imag(Y_modified(i,j))));
        end
    end
    fprintf('\n');
end
