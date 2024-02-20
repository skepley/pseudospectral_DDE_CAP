function q_idx = homological_rhs(a_0, a_n, order, alpha)
%HOMOLOGICAL_RHS - One line description of what the function or script performs (H1 line)
%
%   HOMOLOGICAL_RHS() - A more detailed description of the function
%
%   Syntax:
%       output = HOMOLOGICAL_RHS(input1, input2)
%       [output1, output2] = HOMOLOGICAL_RHS(input1, input2, input3)
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 09-Jun-2022; 

nonLinearTerm = a_0*a_n; 
q_idx = zeros(1, order + 1);  % initialize q_beta for all |beta| = N
for k = 1:(order + 1)
    q_idx(k) = nonLinearTerm.Coefficient(order + 2 - k, k);
end
q_idx = alpha*q_idx;
end % end homological_rhs

