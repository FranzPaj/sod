function [dg_dt] = dy_dt_ecef(g,GM,we)

% disentangle state vector and state transition matrix
y = g(1:6);
Phi = reshape(g(7:end),6,6);

% norm of position vector
r = sqrt(sum(y(1:3).^2));

% first part of ODE
a00_ecef = -y(1:3)*GM/r^3;

W = zeros(3,3);
W(2,1) = -we;
W(1,2) = +we;

a = a00_ecef - W*W*y(1:3) + 2*W*y(4:6);

dy_dt = [y(4:6); a];


% construct matrix df/dy for state transition matrix --> C00 effect
F00 = zeros(6,6);
F00(1:3,4:6) = eye(3);
F00(4:6,1:3) = GM/r^5 * (3*y(1:3)*y(1:3)' - eye(3)*r^2);

% part that comes from the ECEF
FW = zeros(6,6);
FW(4:6,1:3) = - W*W;
FW(4:6,4:6) = 2*W;

% F = F00 + F20 + FW;
F = F00 + FW;

% differential equation for state transition matrix
dPhi_dt = F*Phi;

% output needs to be a vector
dPhi_dt = reshape(dPhi_dt,36,1);

% combine state vector and state transition matrix
dg_dt = [dy_dt; dPhi_dt];

