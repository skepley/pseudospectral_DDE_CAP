function Rnx = Rn_PSA(lambda, D, x)
%RN_PSA - Evaluate Rn(x) for Wright's equation which is the blowup map followed by projection on the last coordinate
%    
% Inputs:
%       input - Description
%
% Outputs:
%       output: - Description
%
% Subfunctions: none
% Classes required: none
% Other m-files required: none
% MAT-files required: none

% Author: Shane Kepley
% email: s.kepley@vu.nl
% Date: 30-Oct-2023; 

Rnx_Full = blowup(lambda, D, x);  % evaluate blowup function
Rnx = Rnx_Full(end);  % get last coefficient
end % end Rn_PSA

