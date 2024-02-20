function a = wright_solve_homological_diagonal(alpha, lambda, A, xi_1, numCoef, n)
%WRIGHT_SOLVE_HOMOLOGICAL_DIAGONAL - Computes a Taylor parameterization for the unstable manifold of the pseudospectral
% approximation of Wright's equation by solving the homological equations for (anti)diagonal data shape. This means compute
% coefficients satisfying |beta| <= N (see wright_solve_homological_rectangular for the alternative).

% Inputs:
%    n - An integer specifying the number of basis functions to use in the
%        pseudospectral approximation.
%    numCoef - An integer specifying the number of coefficients to use in
%              the Taylor parameterization.
%    alpha - A real number specifying the value of the parameter alpha in
%            Wright's equation.
%    lambda - A vector of eigenvalues associated with the unstable manifold.
%    xi_1 - The scaled eigenvector for the pseudospectral approximation associated with lambda(1).
%
% Outputs:
%    a - An array of Scalars of length n+1, where each entry is a Taylor polynomial with
%        numCoef coefficients. The Taylor parameterization represents the
%        unstable manifold of the pseudospectral approximation of Wright's
%        equation.
%
% Other m-files required: homological_rhs.m, homological_linear_map.m
% Subfunctions: none
% MAT-files required: none
%
% See also: homological_rhs, homological_linear_map

% Author: Shane Kepley
% email: s.kepley@vu.nl
% Date: 08-May-2023;

xi_2 = conj(xi_1);  % get eigenvector for lambda_2

% Initialize the first element of the array 'a' with zeros, where each entry
% of 'a' is a Taylor polynomial with 'numCoef' coefficients.
a(1) = Scalar.zeros({'Taylor', 'Taylor'}, [numCoef, numCoef]);
for idx = 2:(n+1)
    a(idx) = Scalar.zeros({'Taylor', 'Taylor'}, [numCoef, numCoef]);
end

% set up the constant and first order data (i.e. equilibrium and eigenvector data)
for idx = 1:n+1
    a(idx).Coefficient(1, 1) = 0; % equilibrium is at zero
    a(idx).Coefficient(1, 2) = xi_2(idx);
    a(idx).Coefficient(2, 1) = xi_1(idx);
end

% Update the coefficients of the 'a' polynomials for orders higher than 1 by recursively solving the homological equations.
for order = 2:numCoef-1
    hom_rhs = homological_rhs(a(1), a(end), order, alpha);
    
    for beta_1 = 0:order
        multiIndex = [beta_1, order - beta_1];
        A_idx = homological_linear_map(A, lambda, multiIndex);
        rhs_idx = cat(1, hom_rhs(1 + multiIndex(2)), zeros(n, 1));
        p_multiIndex = A_idx\rhs_idx;
        
        % update Taylor coefficients
        for idx = 1:n+1
            a(idx).Coefficient(1 + multiIndex(1), 1 + multiIndex(2)) = p_multiIndex(idx);
        end
    end
end
end % end wright_solve_homological_diagonal

