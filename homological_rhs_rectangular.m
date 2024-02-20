function q_idx = homological_rhs_rectangular(a_0, a_n, multiIndex, alpha)
%HOMOLOGICAL_RHS_RECTANGULAR - One line description of what the function or script performs (H1 line)
%
%   Syntax:
%       output = HOMOLOGICAL_RHS_RECTANGULAR(input)
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
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
%   Date: 13-May-2023; 

nonLinearTerm = a_0*a_n; 
q_idx = nonLinearTerm.Coefficient(1 + multiIndex(1), 1 + multiIndex(2));
q_idx = alpha*q_idx;

end % end homological_rhs_rectangular

