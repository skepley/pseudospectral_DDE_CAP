%WRIGHT_EQN_ODE_VALIDATION - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       WRIGHT_EQN_ODE_VALIDATION description
%
%   Output:
%       WRIGHT_EQN_ODE_VALIDATION output
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 30-Oct-2023; 

clearvars
close all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))

tic
% set up parameters for the pseudospectral approximation
n = 10;
numCoef = 25;
alpha = 2;
alpha_enclosure = intval(alpha);

% set up linearization and get eigenvalue/eigenvector data
p = zeros(1, n+1);
cheb_diff_matrix = cheb(n, -1, 0);
D = cheb_diff_matrix(2:end, 2:end);
D_enclosure = intval(D);
row1 = [zeros(1, n), -alpha];
A = cat(1, row1, cheb_diff_matrix(2:end, :));
[vec, val] = eig(A);


% get vector of ordered eigenvalues and make sure exactly 2 are unstable
val = diag(val);
unstable_idx = find(real(val) == max(real(val)), 1);
lambda_1 = val(unstable_idx);
lambda_1 = real(lambda_1) + 1i*abs(imag(lambda_1)); % force the index lambda_1 to refer to the eigenvalue with positive imaginary part
lambda_2 = conj(lambda_1);
lambda = [lambda_1, lambda_2];
% lambda_enclosure = midrad(lambda, 0);
lambda_1_enclosure = PSA_eigenvalue_enclosure(alpha_enclosure, lambda_1, D_enclosure);
lambda_enclosure = [lambda_1_enclosure, conj(lambda_1_enclosure)];
if ~isequal(sum(real(val) > 0), 2)
    error('The number of unstable PSA eigenvalues is something other than 2')
end
disp('PSA has exactly 2 (numerical) unstable eigenvalues')

% Rigorously verify that each numerical eigenvalue which has negative real part is actually in the left half plane
for eigenvalue = val
    if real(eigenvalue) < 0
        eigenvalue_enclosure = PSA_eigenvalue_enclosure(alpha_enclosure, eigenvalue, D_enclosure);
        if sup(real(eigenvalue_enclosure)) >= 0
            error('Validation of stable eigenvalues failed')
        end
    end
    disp('Validation of stable eigenvalues passed')
end

% Choose M and epsilon as required to satisfy Lemma 18
M = 1000; % cutoff for computation in the Z1 bound (see paper)
epsilon = real(lambda_enclosure(1))*M - norm(D_enclosure, 1);  % choose the largest epsilon which satisfies Proposition 16
if sup(epsilon) < inf(alpha_enclosure*norm(D_enclosure, 1)/(M*real(lambda_enclosure(1))))
    error('Chosen M and epsilon do not satisfy Lemma 20')
end

% get scaled eigenvector for lambda_1
ODE_eigenvector_norm_scaling = .5;  % adjust scaling of PSA eigenvectors
pos_eig_idx = find(val == lambda_1);
xi_1 = vec(:, pos_eig_idx);
xi_1 = ODE_eigenvector_norm_scaling*xi_1;  % scale eigenvector
ODE_eigenvector_scaling = real(xi_1(1));
ODE_eigenvector_scaling_enclosure = intval(ODE_eigenvector_scaling);

a = wright_solve_homological_diagonal(alpha, lambda, A, xi_1, numCoef, n);
a_gt = a(1);
b_gt = a(2:end);
x_hat = intval(a_gt);  % An interval numerical approximation where our proofs are centered

F_PSA = wright_psa_zero_map_full(alpha, ODE_eigenvector_scaling, lambda, D, a_gt);
F_PSA_enclosure = wright_psa_zero_map_full(alpha_enclosure, ODE_eigenvector_scaling_enclosure, lambda_enclosure, D_enclosure, x_hat);
disp(norm(F_PSA))
disp(norm(F_PSA_enclosure))

disp(norm(D, inf)/real(lambda_1))

%% Test some code for evaluating DFn(x) NONRIGOROUSLY
% 
x = a_gt;
% build operator evaluations
% first matrix representation of v |---> Rn(v)
allOnes = Scalar(ones(size(x.Coefficient)), {'Taylor', 'Taylor'});
Rn_PSA_values = Rn_PSA(lambda, D, allOnes);
Rn_PSA_matrix = diag(reshape(Rn_PSA_values.Coefficient, [], 1));

% second is a matrix representation of the map v |---> Rn(x)*v
Rnx = Rn_PSA(lambda, D, x);  % apply R to x to obtain R(x)
left_times_Rx = Rnx.leftmultiplicationoperator();  

% third is a matrix representation of the map v |---> x*R(v)
left_times_x = x.leftmultiplicationoperator();
x_times_Rnv = left_times_x.Matrix*Rn_PSA_matrix;

% last is a matrix representation of the map v_beta | ---> <lambda, beta>v_beta
lambeta_values = lambeta_map(lambda, allOnes);
lambeta_matrix = diag(reshape(lambeta_values.Coefficient, [], 1));

% construct the derivative correct for order |beta| >= 2
An_dagger_bar = lambeta_matrix + alpha * Rn_PSA_matrix + alpha * (left_times_Rx.Matrix + x_times_Rnv);

% input derivative of F at orders |beta| <= 1 manually
An_dagger_bar(1, :) = [1, zeros(1, numCoef^2 - 1)];
An_dagger_bar(2, :) = [0, 1, zeros(1, numCoef^2 - 2)];
An_dagger_bar(1 + numCoef, :) = [zeros(1, numCoef), 1, zeros(1, numCoef^2 - numCoef - 1)];
An_bar = inv(An_dagger_bar);

%% Code for computing Y0 bound
% compute low order terms of An*Fn(x_hat) i.e. An o Pi^N o Fn o Pi^N
An_bar_enclosure = intval(inv(An_dagger_bar)); % upper left block of An_bar
Fx_hat_trunc = Scalar(F_PSA_enclosure.Coefficient(1:numCoef, 1:numCoef), {'Taylor', 'Taylor'});
AFx_lower_order = An_bar_enclosure * reshape(Fx_hat_trunc.Coefficient, [], 1);

% compute high order terms of AF(x_hat)
allOnes_full_intval = intval(Scalar(ones([2*numCoef-1, 2*numCoef-1]), {'Taylor', 'Taylor'}));
dot_terms = lambeta_map(lambda_enclosure, allOnes_full_intval);
R_PSA_terms = Rn_PSA(lambda_enclosure, D_enclosure, allOnes_full_intval);
higher_order_coefficients = Scalar(1./(dot_terms.Coefficient + alpha_enclosure*R_PSA_terms.Coefficient), {'Taylor', 'Taylor'});  % overflow terms of A
AFx_full = Scalar(higher_order_coefficients.Coefficient.*F_PSA_enclosure.Coefficient,{'Taylor', 'Taylor'});
AFx_high_order = AFx_full;
AFx_high_order.Coefficient(1:numCoef, 1:numCoef) = intval(zeros(numCoef, numCoef));  % zero out the finite block which is already computed above

% add them up
Y0 = AFx_high_order.norm() + sum(abs(AFx_lower_order))

%% Compute the Z0 bound
An_dagger_bar_enclosure = intval(An_dagger_bar);
Z0 = norm(eye(size(An_dagger_bar_enclosure, 1)) - An_dagger_bar_enclosure*An_bar_enclosure, 1)

%% Compute the Z1 bound
% We construct the Z1 bound as a sum of the form Z11 + alpha*Z12*Z13 

% We begin with the first term in the Z1 bound denoted by Z11 := ||An(Pi^N o DF(x) o Pi^N - AnDagger)||
% First we evaluate a rigorous enclosure of DFn(xHat) via some sub-operator evaluations

% first matrix representation of v |---> Rn(v)
allOnes_enclosure = intval(Scalar(ones(size(x_hat.Coefficient)), {'Taylor', 'Taylor'}));
Rn_PSA_values = Rn_PSA(lambda_enclosure, D_enclosure, allOnes_enclosure);
Rn_PSA_matrix = diag(reshape(Rn_PSA_values.Coefficient, [], 1));

% second is a matrix representation of the map v |---> Rn(x)*v
Rnx = Rn_PSA(lambda_enclosure, D_enclosure, x_hat);  % apply Rn to x to obtain Rn(x)
left_times_Rx = Rnx.leftmultiplicationoperator();  

% third is a matrix representation of the map v |---> x*Rn(v)
left_times_x = x_hat.leftmultiplicationoperator();
x_times_Rnv = left_times_x.Matrix*Rn_PSA_matrix;

% last is a matrix representation of the map v_beta | ---> <lambda, beta>v_beta
% this line is identical to the DDE case
lambeta_values = lambeta_map(lambda_enclosure, allOnes_enclosure);
lambeta_matrix = diag(reshape(lambeta_values.Coefficient, [], 1));

% construct the derivative correct for order |beta| >= 2
DFn_x_enclosure = lambeta_matrix + alpha * Rn_PSA_matrix + alpha * (left_times_Rx.Matrix + x_times_Rnv);

% input derivative of F at orders |beta| <= 1 manually
DFn_x_enclosure(1, :) = [1, zeros(1, numCoef^2 - 1)];
DFn_x_enclosure(2, :) = [0, 1, zeros(1, numCoef^2 - 2)];
DFn_x_enclosure(1 + numCoef, :) = [zeros(1, numCoef), 1, zeros(1, numCoef^2 - numCoef - 1)];
% DFn_x_enclosure is an interval enclosure of DFn(xHat);

Z11 = norm(An_bar_enclosure*(DFn_x_enclosure - An_dagger_bar_enclosure),1)

% Next is Z12 
% RIGOROUS VERSION:
allOnes_to_order_M_enclosure = intval(Scalar(ones(M,M), {'Taylor', 'Taylor'}));
Rn_to_order_M_values = Rn_PSA(lambda_enclosure, D_enclosure, allOnes_to_order_M_enclosure);
lambeta__to_order_M_values = lambeta_map(lambda_enclosure, allOnes_to_order_M_enclosure);

% NONRIGOROUS VERSION
% allOnes_to_order_M_enclosure = Scalar(ones(M,M), {'Taylor', 'Taylor'});
% Rn_to_order_M_values = Rn_PSA(lambda, D, allOnes_to_order_M_enclosure);
% lambeta_to_order_M_values = lambeta_map(lambda, allOnes_to_order_M_enclosure);

%% 
Z12_sum = 1./(alpha_enclosure*Rn_to_order_M_values.Coefficient + lambeta__to_order_M_values.Coefficient);
% % Z12_sum.Coefficient(1:numCoef, 1:numCoef) = intval(0);
Z12_sum(1:numCoef, 1:numCoef) = 0;

Z12_max1 = max(abs(Z12_sum(:)));
Z12_max2 = 1./(real(lambda_enclosure(1)*M - alpha_enclosure*norm(D_enclosure, 1)/epsilon)); 
Z12 = max(Z12_max1, Z12_max2);

% Finally Z13
Z13_max1 = max(abs(Rn_to_order_M_values.Coefficient(:)));
Z13_max2 = norm(D_enclosure, 1)/epsilon;
Z13_max = max(Z13_max1, Z13_max2);
Z13 = Rnx.norm() + x_hat.norm()*Z13_max;

Z1 = Z11 + alpha_enclosure*Z12*Z13
% Z1 = (2*alpha_enclosure*x_hat.norm())/(real(lambda_enclosure(1))*(numCoef-1)) + norm(A_bar_enclosure*(DF_x_enclosure - A_dagger_bar_enclosure),1)


% Z2 bound
Z21 = max(Z12, norm(An_bar_enclosure, 1)); 
Z22 = Z13_max;
Z2 = 2*alpha_enclosure*Z21*Z22

% solve radii polynomial
r = solveradiipoly([Z2, (Z0 + Z1 -1), Y0])
toc
