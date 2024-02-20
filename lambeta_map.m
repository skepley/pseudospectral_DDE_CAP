function dotTerms = lambeta_map(lambda, x)
%LAMBETA_MAP - Evaluate the linear map v_beta |---> <lambda, beta>v_beta
%
%   Syntax:
%       output = LAMBETA_MAP(input)
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
%   Date: 08-May-2023; 

N = x.Truncation(1);
dotTerms = Scalar.zeros({'Taylor', 'Taylor'}, size(x.Coefficient));

if isequal(x.NumericalClass, 'intval')
    dotTerms = intval(dotTerms);
end

for beta_1 = 0:N-1
    for beta_2 = 0:N-1
        dotTerms.Coefficient(beta_1 + 1, beta_2 + 1) = x.Coefficient(beta_1 + 1, beta_2 + 1)*my_dot(lambda, [beta_1, beta_2]);
    end
end
end % end lambeta_map

