function F_x = zero_map(alpha, lambda, x)
%ZERO_MAP - One line description of what the function or script performs (H1 line)
%
%   Syntax:
%       output = ZERO_MAP(input)
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
%   Date: 14-Nov-2022; 

F_x = delta(alpha, lambda, x) + alpha*(L(lambda, x)*x);
end % end zero_map

