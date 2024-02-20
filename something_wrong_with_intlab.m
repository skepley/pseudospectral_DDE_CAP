%SOMETHING_WRONG_WITH_INTLAB - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       SOMETHING_WRONG_WITH_INTLAB description
%
%   Output:
%       SOMETHING_WRONG_WITH_INTLAB output
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 16-Oct-2023; 

T_eig = @(z)z - (z + alpha_enclosure.*exp(-z))./(1 - alpha_enclosure.*exp(-z));
% T_eig = @(z)-alpha_enclosure.*exp(-z).*(1+z)./(1 - alpha_enclosure*exp(-z));
DT_eig = @(z)(z + alpha_enclosure*exp(-z))*(alpha_enclosure*exp(-z))/(1 - alpha_enclosure*exp(-z)).^2;
z_hat = T_eig(lambda(1));
z_hat_ball = midrad(mid(z_hat), 1e-6);
x1 = inf(real(z_hat_ball));
x2 = sup(real(z_hat_ball));
y1 = inf(imag(z_hat_ball));
y2 = sup(imag(z_hat_ball));
t = linspace(0, 1, 100);
south = [(1-t).*x1 + t.*x2; y1*ones(size(t))];
north = [(1-t).*x1 + t.*x2; y2*ones(size(t))];
east = [x1*ones(size(t)); (1-t).*y1 + t.*y2];
west = [x2*ones(size(t)); (1-t).*y1 + t.*y2];
Ts = mid(T_eig(south(1,:) + 1i*south(2,:)));
Tn = mid(T_eig(north(1,:) + 1i*north(2,:)));
Te = mid(T_eig(east(1,:) + 1i*east(2,:)));
Tw = mid(T_eig(west(1,:) + 1i*west(2,:)));

close all
T_z = T_eig(z_hat_ball);
abs(DT_eig(z_hat_ball))
in(T_eig(z_hat_ball),z_hat_ball)
figure
hold on
plotintval(z_hat_ball);
pause(1)
plotintval(T_z);
pause(1)
plotintval(T_eig(T_z));

% test exponential
zc = midrad(1e-15 + 1e-15*1i, 5);
figure
hold on
plotintval(zc);
plotintval(zc/2);