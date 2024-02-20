function Fn_x = wright_psa_zero_map_full(alpha, ODE_eigenvector_scaling, lambda, D, x)
%WRIGHT_PSA_ZERO_MAP_FULL - Evaluate the scalar zero finding map for the pseudo spectral approximation (ODE) of Wright's equation
%
%   Syntax:
%       Fn_x = WRIGHT_PSA_ZERO_MAP_FULL (alpha, ODE_eigenvector_scaling, lambda, D, x)
%
%   Inputs:
%       alpha - The parameter in Wrights equation
%       lambda - The vector of eigenvalues for the PS approximation
%       D - The differentiation matrix
%       x - a sequence of coefficients
%
%   Outputs:
%       T(x) - Another sequence of coefficients
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 27-Feb-2023;

N = x.Truncation(1) - 1;
Rnx_n = Rn_PSA(lambda, D, x);
% Rx_Full = blowup(lambda, D, x);  % evaluate blowup function
% Rx_Last = Rx_Full(end);  % get last coefficient
% Rnx_n = Scalar(squeeze(Rx_Last.Coefficient(:, :, end)), {'Taylor', 'Taylor'});
% switch x.NumericalClass
%     case 'double'
%         dotTerms = Scalar.zeros({'Taylor', 'Taylor'}, size(x.Coefficient));
%     case 'intval'
%         dotTerms = Scalar(intval(zeros(size(x.Coefficient))), {'Taylor', 'Taylor'});
% end
% 
% for beta_1 = 0:N
%     for beta_2 = 0:N
%         dotTerms.Coefficient(beta_1 + 1, beta_2 + 1) = x.Coefficient(beta_1 + 1, beta_2 + 1)*my_dot(lambda, [beta_1, beta_2]);
%     end
% end
dotTerms = lambeta_map(lambda, x);
Fn_x = Rnx_n*alpha + mtimes(x, Rnx_n, 'Full')*alpha + dotTerms; % Compute Fn(xHat) including higher order terms
Fn_x.Coefficient(1, 1) = x.Coefficient(1, 1);  % Equilibrium is at 0
Fn_x.Coefficient(1, 2) = x.Coefficient(1, 2) - ODE_eigenvector_scaling;  % eigenvector scalings in both directions
Fn_x.Coefficient(2, 1) = x.Coefficient(2, 1) - ODE_eigenvector_scaling;
end % end wright_psa_zero_map_full

