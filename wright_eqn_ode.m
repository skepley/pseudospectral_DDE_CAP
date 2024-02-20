function f_y = wright_eqn_ode(t, y, alpha)
%WRIGHT_EQN_ODE - One line description of what the function or script performs (H1 line)
%
%   WRIGHT_EQN_ODE() - A more detailed description of the function
%
%   Syntax:
%       output = WRIGHT_EQN_ODE(input1, input2)
%       [output1, output2] = WRIGHT_EQN_ODE(input1, input2, input3)
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

D = cheb(length(y)-1, -1, 0);
D = D(2:end, :);
f_y = cat(1, -alpha*y(end) - alpha*y(1)*y(end), D*y);
end % end wright_eqn_ode

