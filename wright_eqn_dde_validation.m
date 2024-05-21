%WRIGHT_EQN_DDE_VALIDATION - Compute local unstable manifold for the pseudospectral ODE associated to Wright's equation
%
%   Description:
%       WRIGHT_EQN_DDE_VALIDATION assumes that the script "wright_eqn_ode_validation.m" has already been run
%       successfully.
%
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: wright_eqn_ode_validation.m
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 09-Jun-2022;

clearvars
close all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
load ode_data_final



% rename eigenvalues from the PSA to avoid naming conflicts
lambda_PSA = mid(lambda);
lambda_PSA_enclosure = lambda_enclosure;


%% Get an interval enclosure of the DDE eigenvalues

T_eig_DDE = @(z)z - (z + alpha_enclosure.*exp(-z))./(1 - alpha_enclosure.*exp(-z));
DT_eig_DDE = @(z)(z + alpha_enclosure*exp(-z))*(alpha_enclosure*exp(-z))/(1 - alpha_enclosure*exp(-z)).^2;
lambda_DDE = mid(lambda_PSA_enclosure(1));
for niter = 1:10
    lambda_DDE = T_eig_DDE(lambda_DDE);
end
lambda_DDE = mid(lambda_DDE);  % coerce back into a double
z_hat_enclosure = midrad(lambda_DDE, 0);

% Y, Z bounds
Y = sup(norm(T_eig_DDE(z_hat_enclosure) - z_hat_enclosure));
test_ball = midrad(lambda_DDE, 1e-3);
Z = sup(norm(DT_eig_DDE(test_ball)));
rEigenvalue = sup(intval(Y)/intval(1 - Z));
lambda_enclosure = midrad([lambda_DDE, conj(lambda_DDE)], [rEigenvalue, rEigenvalue]);



%% This will be replaced by the implementation from wright_eqn_ode_validation
% a = wright_solve_homological_diagonal(alpha, lambda_DDE, A, xi_1, numCoef, n);
% a = wright_solve_homological_diagonal(alpha, mid(lambda_enclosure), A, xi_1, numCoef, n);
% 
% a_gt = a(1);
% b_gt = a(2:end);
% x_hat = intval(a_gt);  % An interval numerical approximation where our proofs are centered
% 
% 
DDE_eigenvector_scaling = real(xi_1(1));
DDE_eigenvector_scaling_enclosure = intval(real(xi_1(1)));
F_DDE_enclosure = wright_DDE_zero_map_full(alpha_enclosure, DDE_eigenvector_scaling_enclosure, lambda_enclosure, x_hat);
% % ab = wright_DDE_zero_map_full(alpha, ODE_eigenvector_scaling, lambda_PSA, a_gt);
%  
% disp(norm(F_DDE_enclosure))


%% Construct a floating point approximation for DF(x_bar)
% 
x = a_gt;
% build operator evaluations
% first matrix representation of v |---> R(v)
allOnes = Scalar(ones(size(x.Coefficient)), {'Taylor', 'Taylor'});
R_DDE_values = R_DDE(lambda_PSA, allOnes);
R_DDE_matrix = diag(reshape(R_DDE_values.Coefficient, [], 1));

% second is a matrix representation of the map v |---> R(x)*v
Rx = R_DDE(lambda_PSA, x);  % apply R to x to obtain R(x)
left_times_Rx = Rx.leftmultiplicationoperator();  

% third is a matrix representation of the map v |---> x*R(v)
left_times_x = x.leftmultiplicationoperator();
x_times_Rv = left_times_x.Matrix*R_DDE_matrix;

% last is a matrix representation of the map v_beta | ---> <lambda, beta>v_beta
lambeta_values = lambeta_map(lambda_PSA, allOnes);
lambeta_matrix = diag(reshape(lambeta_values.Coefficient, [], 1));

% construct the derivative correct for order |beta| >= 2
A_dagger_bar = lambeta_matrix + alpha * R_DDE_matrix + alpha * (left_times_Rx.Matrix + x_times_Rv);

% input derivative of F at orders |beta| <= 1 manually
A_dagger_bar(1, :) = [1, zeros(1, numCoef^2 - 1)];
A_dagger_bar(2, :) = [0, 1, zeros(1, numCoef^2 - 2)];
A_dagger_bar(1 + numCoef, :) = [zeros(1, numCoef), 1, zeros(1, numCoef^2 - numCoef - 1)];

%% Code for computing Y0 bound
% compute low order terms of AF(x_hat) i.e. A o Pi^N o F o Pi^N
A_bar_enclosure = intval(inv(A_dagger_bar)); % upper left block of A
Fx_hat_trunc = Scalar(F_DDE_enclosure.Coefficient(1:numCoef, 1:numCoef), {'Taylor', 'Taylor'});
AFx_lower_order = A_bar_enclosure * reshape(Fx_hat_trunc.Coefficient, [], 1);

% compute high order terms of AF(x_hat)
allOnes_full_intval = intval(Scalar(ones([2*numCoef-1, 2*numCoef-1]), {'Taylor', 'Taylor'}));
dot_terms = lambeta_map(lambda_enclosure, allOnes_full_intval);
R_DDE_terms = R_DDE(lambda_enclosure, allOnes_full_intval);
higher_order_coefficients = Scalar(1./(dot_terms.Coefficient + alpha_enclosure*R_DDE_terms.Coefficient), {'Taylor', 'Taylor'});  % overflow terms of A
AFx_full = Scalar(higher_order_coefficients.Coefficient.*F_DDE_enclosure.Coefficient,{'Taylor', 'Taylor'});
AFx_high_order = AFx_full;
AFx_high_order.Coefficient(1:numCoef, 1:numCoef) = intval(zeros(numCoef, numCoef));

% add them up
Y0 = AFx_high_order.norm() + sum(abs(AFx_lower_order))


%% Code for computing Z0
A_dagger_bar_enclosure = intval(A_dagger_bar);
Z0 = norm(eye(size(A_dagger_bar_enclosure, 1)) - A_dagger_bar_enclosure*A_bar_enclosure, 1)


%% Code for computing Z1

% First we evaluate a rigorous enclosure of DF(xHat) 
% Build some sub-operator evaluations
% first matrix representation of v |---> R(v)
allOnes_intval = intval(Scalar(ones(size(x_hat.Coefficient)), {'Taylor', 'Taylor'}));
R_DDE_values = R_DDE(lambda_enclosure, allOnes_intval);
R_DDE_matrix = diag(reshape(R_DDE_values.Coefficient, [], 1));

% second is a matrix representation of the map v |---> R(x)*v
Rx = R_DDE(lambda_enclosure, x_hat);  % apply R to x to obtain R(x)
left_times_Rx = Rx.leftmultiplicationoperator();  

% third is a matrix representation of the map v |---> x*R(v)
left_times_x = x_hat.leftmultiplicationoperator();
x_times_Rv = left_times_x.Matrix*R_DDE_matrix;

% last is a matrix representation of the map v | ---> <lambda, beta>v_beta
lambeta_values = lambeta_map(lambda_enclosure, allOnes_intval);
lambeta_matrix = diag(reshape(lambeta_values.Coefficient, [], 1));

% construct the derivative correct for order |beta| >= 2
DF_x_enclosure = lambeta_matrix + alpha * R_DDE_matrix + alpha * (left_times_Rx.Matrix + x_times_Rv);

% input derivative of F at orders |beta| <= 1 manually
DF_x_enclosure(1, :) = [1, zeros(1, numCoef^2 - 1)];
DF_x_enclosure(2, :) = [0, 1, zeros(1, numCoef^2 - 2)];
DF_x_enclosure(1 + numCoef, :) = [zeros(1, numCoef), 1, zeros(1, numCoef^2 - numCoef - 1)];
% DF_x_enclosure is an interval enclosure of DF(xHat);

% Z1 = (2*alpha_enclosure*x_hat.norm())/(real(lambda_enclosure(1))*(numCoef-1)) + norm(A_bar_enclosure*(DF_x_enclosure - A_dagger_bar_enclosure),1)

if real(lambda_enclosure(1))*(numCoef-1) - alpha_enclosure*exp(-real(lambda_enclosure(1))*(numCoef-1)) < 0
    error('Lemma 29 estimate is not satisfied')
end


Z1 = (2*alpha_enclosure*x_hat.norm())/(real(lambda_enclosure(1))*(numCoef-1) - alpha_enclosure*exp(-real(lambda_enclosure(1))*(numCoef-1)))...
    + norm(A_bar_enclosure*(DF_x_enclosure - A_dagger_bar_enclosure),1)


%% Code for computing Z2
Z2 = 2*alpha_enclosure*max(1./(real(lambda_enclosure(1))*(numCoef-1) - alpha_enclosure*exp(-real(lambda_enclosure(1))*(numCoef-1))), norm(A_bar_enclosure, 1))

%% Solve radii polynomial
r_DDE = solveradiipoly([Z2, (Z0 + Z1 -1), Y0])


% disp(size(DF_matrix, 1))
% disp(rank(DF_matrix))
return
% %% lets look at it
% close all
% figure;
% hold on
% [R, T] = meshgrid(linspace(0, 1, 100), linspace(0, 2*pi, 101));
% z = R.*exp(1i*T);
% ZZ = reshape(z, [], 1);
% ZZc = conj(ZZ);
% y1 = a(1).eval([ZZ, ZZc]);
% y2 = a(floor(end/2)).eval([ZZ, ZZc]);
% y3 = a(end).eval([ZZ, ZZc]);
% % y4 = a(4).eval([ZZ, ZZc]);
% scatter3(real(y1), real(y2), real(y3), 5, 'filled');


% % make a test orbit near equilibrium and check tht it stays on the manifold
% r = 0.9;
% theta = 0;
% z0 = r*exp(1i*theta);
% y_init = real(cell2mat(a.eval([z0, conj(z0)])));
% f = @(t,y)wright_eqn_ode(t,y,alpha);
% scatter3(y_init(1), y_init(floor((n+1)/2)), y_init(n+1), 'g')
% plot_orbit(f, [0, -2], y_init, [1,floor((n+1)/2),n+1],'OdeOptions', odeset('RelTol',1e-13,'AbsTol',1e-13), 'PlotOptions', {'r', 'LineWidth', 3})
% plot_orbit(f, [0, 2], y_init, [1,floor((n+1)/2),n+1],'OdeOptions', odeset('RelTol',1e-13,'AbsTol',1e-13), 'PlotOptions', {'g', 'LineWidth', 3})
% view([-24,22])





% % plot coefficients of the manifolds
% cheb_nodes = 0.5 * (cos(linspace(0, pi, n+1)) - 1);
% figure
% hold on
% for idx_sum = 2:numCoef
%     for jj = 1:(idx_sum+1)/2
%         c = coef_idx(idx_sum - jj, jj);
%         if norm(c) > 0
%             plot(cheb_nodes, log(abs(c)), 'LineWidth', 3);
%             disp([idx_sum - jj,jj])
%         end
%     end
% end
% dealfig()
% legend
