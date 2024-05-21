function R_u = blowup(lambda, D, u)
%BLOWUP - Evaluate the blowup operator for the pseudospectral ODE associated with some DDE.
%
%   Syntax:
%       Ru = BLOWUP(u, lambda, D) evaluates the function R_n(u)_beta = u_beta(D - <lambda, beta>I)^{-1}*D*ONE)
%       where ONE is the vector of all ones in R^n.
%
%   Inputs:
%       u - A Scalar
%       lambda - the vector of eigenvalues
%       D - Chebyshev differentiation matrix of order n
%
%   Outputs:
%       Ru - A scalar quantity
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 13-Feb-2023;



n = size(D, 1);
maxOrder = size(u.Coefficient, 1)-1;

switch u.NumericalClass
    
    case 'double'
        D_one = D*ones(n,1);
        Ru_coeffs = zeros([size(u.Coefficient), n]);
        for order = 1:maxOrder
            for beta_2 = 0:order
                beta = [order - beta_2, beta_2];
                Ru_coeffs(1 + beta(1), 1 + beta(2), :) = u.Coefficient(1 + beta(1), 1 + beta(2))*...
                    ((D - sum(lambda.*beta)*eye(size(D))) \ D_one);
            end
        end
        
        R_u = Scalar(squeeze(Ru_coeffs(:, :, 1)), {'Taylor', 'Taylor'});
        for j = 2:n
            R_u(j) = Scalar(squeeze(Ru_coeffs(:, :, j)), {'Taylor', 'Taylor'});
        end
    case 'intval'
        if ~isequal(class(D), 'intval')
            error('D must be an interval enclosure of a Chebyshev differentiation matrix')
        end
        D_one = D*intval(ones(n,1));
        Ru_coeffs = intval(zeros([size(u.Coefficient), n]));
        Id = intval(eye(size(D)));
        for order = 1:maxOrder
            for beta_2 = 0:order
                beta = [order - beta_2, beta_2];
                Ru_coeffs(1 + beta(1), 1 + beta(2), :) = u.Coefficient(1 + beta(1), 1 + beta(2))*...
                    ((D - sum(lambda.*beta)*Id) \ D_one);
            end
        end
        
        R_u = Scalar(squeeze(Ru_coeffs(:, :, 1)), {'Taylor', 'Taylor'});
        for j = 2:n
            R_u(j) = Scalar(squeeze(Ru_coeffs(:, :, j)), {'Taylor', 'Taylor'});
        end
end
end % end blowup

