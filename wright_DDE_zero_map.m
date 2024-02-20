function F_x = wright_DDE_zero_map(alpha, DDE_eigenvector_scaling, lambda, x)
%WRIGHT_DDE_ZERO_MAP - Evaluate the DDE zero finding problem for Wright's equation including truncation.
%
%   Syntax:
%       output = WRIGHT_DDE_ZERO_MAP(input)
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

error('Do not call this anuymore. Use wright_dde_zero_map_full')

N = x.Truncation(1) - 1;
R_x = R_DDE(lambda, x);
dotTerms = lambeta_map(lambda, x);
F_x = alpha*R_x + dotTerms + alpha*(x*R_x);
F_x.Coefficient(1, 1) = x.Coefficient(1, 1);  % Equilibrium is at 0
F_x.Coefficient(1, 2) = x.Coefficient(1, 2) - DDE_eigenvector_scaling;
F_x.Coefficient(2, 1) = x.Coefficient(2, 1) - DDE_eigenvector_scaling;
end % end wright_DDE_zero_map

