function linearMap = homological_linear_map(Df_p, lambda, multiIndex)
%HOMOLOGICAL_LINEAR_MAP - One line description of what the function or script performs (H1 line)
%
%   HOMOLOGICAL_LINEAR_MAP() - A more detailed description of the function
%
%   Syntax:
%       output = HOMOLOGICAL_LINEAR_MAP(input1, input2)
%       [output1, output2] = HOMOLOGICAL_LINEAR_MAP(input1, input2, input3)
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

linearMap = Df_p - sum(lambda.*multiIndex)*eye(size(Df_p));

end % end homological_linear_map

