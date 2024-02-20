function lambda_rigorous = PSA_eigenvalue_enclosure(alpha, lambda, D)
% PSA_EIGENVALUE_ENCLOSURE - Compute a rigorous enclosure of a numerical eigenvalue of the PSA problem associated with
% Wright's equation
%    
% Inputs:
% 	alpha - (intval) Interval enclosure of parameter
%   lambda - (double) A numerical approximation of an eigenvalue
%   D - (intval) Interval enclosure of differentiation matrix
%
% Outputs:
% 	lambda_rigorous - If successful then this is an interval enclosure centered at lambda which is guranteed to contain
% 	exactly 1 true eigenvalue of the linearized equation.
%
% Subfunctions: none
% Classes required: none
% Other m-files required: none
% MAT-files required: none

% Author: Shane Kepley
% email: s.kepley@vu.nl
% Date: 11-Dec-2023; 

% some helper functions
n = size(D, 1);
D_one = D*ones(n, 1);
res = @(z)D - z*eye(n);  % define the residual function
proj_n = @(z)z(n);  % define the projection onto the last coordinate
% test_ball = midrad(lambda(1), 1e-10);

% define Newton function for the eigenvalue problem in terms of Delta_n and its deriatives
Delta_n = @(z)z + alpha*proj_n(res(z)\D_one);
D_Delta_n = @(z)intval('1') + alpha*proj_n(((res(z)^2)\D_one));
D2_Delta_n = @(z)intval('2')*alpha*proj_n((res(z)^3)\D_one);
T_eig_ODE = @(z)z - Delta_n(z)./D_Delta_n(z);
DT_eig_ODE = @(z)(Delta_n(z).*D2_Delta_n(z))./((D_Delta_n(z)).^2);

% Rigorously enclose the eigenvalues of the PSA 
lambda_enclosure = midrad(lambda, 0);  % coerce lambda to an intval

% Y, Z bounds
test_ball = midrad(lambda, 1e-7);
Y = sup(norm(T_eig_ODE(lambda_enclosure) - lambda_enclosure));
Z = sup(norm(DT_eig_ODE(test_ball)));
rBound = sup(intval(Y)/intval(1 - Z));

if rBound < 1e-7
    lambda_rigorous = midrad(lambda, rBound);
else
    error('Rigorous eigenvalue enclosure failed')
end


end % end PSA_eigenvalue_enclosure

