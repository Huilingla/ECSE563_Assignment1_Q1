% Validation script for IEEE 9-bus test system
clear; clc;

% Load test system data
ieee9_A1;

fprintf('IEEE 9-Bus System Validation\n');
fprintf('============================\n\n');

% Calculate admittance matrix
Y = admittance(nfrom, nto, r, x, b);

% Calculate impedance matrix  
Z = impedance(nfrom, nto, r, x, b);

% Calculate node voltages using linsolve()
V = linsolve(Y, Iint);

% Display results in polar coordinates
fprintf('Admittance Matrix Y (pu):\n');
disp(Y);

fprintf('\nImpedance Matrix Z (pu):\n');
disp(Z);

fprintf('\nNode Voltages (polar coordinates):\n');
fprintf('Bus\tMagnitude (pu)\tAngle (degrees)\n');
for k = 1:length(V)
    mag = abs(V(k));
    ang = angle(V(k)) * 180/pi;
    fprintf('%d\t%.4f\t\t%.2f\n', k, mag, ang);
end

% Validation checks
fprintf('\nValidation Results:\n');
fprintf('||Y*Z - I|| = %.2e (should be near zero)\n', norm(Y*Z - eye(9)));
fprintf('Maximum voltage magnitude: %.4f pu\n', max(abs(V)));
fprintf('Minimum voltage magnitude: %.4f pu\n', min(abs(V)));
