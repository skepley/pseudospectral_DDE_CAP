function lambeta = my_dot(lambda, beta)
%MY_DOT - Compute the dot product. This is rigorous if lambda and beta are intervals
%
%   Syntax:
%       output = MY_DOT(input)
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
%   Date: 13-Feb-2023;


lambeta = sum(lambda.*beta);


end % end my_dot

