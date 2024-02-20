function output = R_DDE(lambda, x)
%DELTA - Evaluate the resolvent-like function for Wrights equation. The output is what is currently denoted by R in the paper, not delta.
%
%   Syntax:
%       output = DELTA(lambda, beta)
%
%   Inputs:
%      lambda - A vector of eigenvalues (or interval enclosures of eigenvalues)
%       x - A Scalar sequence (possibly interval coeficients)
%
%   Outputs:
%       output - Evaluation of R(x). If x and lambda are interval enclosures then this returns a rigorous interval
%         enclosure of R(x)
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 14-Nov-2022;

if isa(x, 'Scalar')
    order = x.Truncation(1);

    switch x.NumericalClass
        case 'double'
            scaleFactor = zeros(order, order);
        case 'intval'
            scaleFactor = intval(zeros(order, order));
    end
    
    for row_idx = 0:order-1
        for col_idx = 0:order-1
            beta = [row_idx, col_idx];
            scaleFactor(1 + row_idx, 1 + col_idx) = R_DDE(lambda, beta);
        end
    end
    output = Scalar(scaleFactor.*x.Coefficient, x.Basis);
    
else
    output = exp(-my_dot(lambda, x));
end
end % end delta

