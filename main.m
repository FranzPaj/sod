clc, clear, close all

%%
% Useful constants and parameters, distances in km, times in sec,
% angles in rad

C = 299792.458;  % Speed of light - km/s
WE = 7.292115e-5;  % Rotation of the Earth - rad/s
R_E = 6378;  % Earth radius - km
GM = 3.986e5;  % Gravitational parametre Earth - km^3/s^2

% Read the GNSS data

load('data/t.txt') % epochs - s

load('data/CA_range.txt') % pseudorange observations from CA code - km

load('data/PRN_ID.txt') % PRN ID of tracked GPS satellites

load('data/rx_gps.txt') % GPS satellite positions (transmitters) - km
load('data/ry_gps.txt')
load('data/rz_gps.txt')

load('data/vx_gps.txt') % GPS satellite velocities (transmitters) - km/s
load('data/vy_gps.txt')
load('data/vz_gps.txt')

load('data/rx.txt') % precise positions (receivers) - km
load('data/ry.txt')
load('data/rz.txt')

load('data/vx.txt') % precise velocities (receivers) - km/s
load('data/vy.txt')
load('data/vz.txt')

% Define useful intermediate params
num_epoch = size(t,1);

%% SECTION A
%

% Define the covariance matrix
sigma_pos = 1.e-3;  % Standard deviation of initial position - km
sigma_vel = 0.5e-3;  % Standard deviation of initial position - km
rho_rv = 0.5;  % Correlation of pos and vel along same coor at t1
v_diag = [sigma_pos^2 sigma_pos^2 sigma_pos^2 sigma_vel^2 sigma_vel^2 sigma_vel^2];
P0 = diag(v_diag);  
corr_term = sigma_pos * sigma_vel * rho_rv;  % Correlation term
row = [1 2 3 4 5 6];
col = [4 5 6 1 2 3];
sz = [6 6];
ind = sub2ind(sz,row,col);
P0(ind) = corr_term;

% Define the observations covariance matrix
% Find visible satellites
epoch = 1;
vis = find(PRN_ID(epoch,:));
num_vis = size(vis,2);  % Number of observations
sigma_obs = 1.5e-3;  % Standard deviation of pseudorange observations - km
rho_obs = 0.;  % Correlation of observations
idx_obs = eye(num_vis,num_vis);
R0 = idx_obs.*sigma_obs^2;  % Observation covariance matrix
% R0 = idx_obs.*sigma_obs^2 + (1-idx_obs).*(sigma_obs^2*rho_obs);  % Observation covariance matrix

%% SECTION B
%

% Define inital condition
% Initiale state vector
y0 = [rx(1) ry(1) rz(1) vx(1) vy(1) vz(1)]';
% Initial transition matrix
Phi0 = eye(6);
% Reshape to single coloumn
Phi0 = reshape(Phi0,36,1);
g_0 = [y0; Phi0];

% Integration
% Define timespan - epoch 1 and 2
tspan = [t(1) t(2)];
options = odeset('AbsTol',1e-08);
[tb,gb] = ode45(@(t,g) dy_dt_ecef(g,GM,WE), tspan, g_0, options);

% Extract results
gf = gb(end,:);
y_b = gf(1:6);
Phi_b = reshape(gf(7:end),6,6);

%% SECTION C
%

% Define initial conditions
% Initiale state vector
xhat = [rx(1) ry(1) rz(1) vx(1) vy(1) vz(1)]';
% Initial transition matrix
Phi = eye(6);
% Already defined covariance matrices
Phat = P0;
R = R0;
% Initialise result matrices
xhat_c = zeros(num_epoch,6);
xhat_c(1,:) = xhat;
Phat_c = zeros(num_epoch,6,6);
Phat_c(1,:,:) = Phat;

for epoch=2:num_epoch
% for epoch=2:3

    % Integration
    % Reshape to single coloumn
    Phi = reshape(Phi,36,1);
    g0 = [xhat; Phi];
    % Define timespan - epoch 1 and 2
    tspan = [t(epoch-1) t(epoch)];
    options = odeset('AbsTol',1e-12);
    [t_temp,g] = ode45(@(t,g) dy_dt_ecef(g,GM,WE), tspan, g0, options);
    % Extract results
    gf = g(end,:);
    xbar = gf(1:6)';
    Phi = reshape(gf(7:end),6,6);
    % Update covariance matrix
    P = Phi * Phat * Phi';

    % Update with observations
    % Find the visible GPS sats
    vis = find(PRN_ID(epoch,:));
    num_vis = size(vis,2);
    zbar = CA_range(epoch,vis)';
    rx_vis = rx_gps(epoch,vis);
    ry_vis = ry_gps(epoch,vis);
    rz_vis = rz_gps(epoch,vis);
    % Create opbservations covariance matrix
    idx_obs = eye(num_vis,num_vis);
    R = idx_obs.*sigma_obs^2;  % Observation covariance matrix
    % Calculate design matrix
    H = zeros(num_vis,6);        
    rho = sqrt((xbar(1)-rx_vis).^2 + (xbar(2)-ry_vis).^2 + (xbar(3)-rz_vis).^2);
    H(:,1) = (xbar(1) - rx_vis)./rho;
    H(:,2) = (xbar(2) - ry_vis)./rho;
    H(:,3) = (xbar(3) - rz_vis)./rho;
    % Calculate residual
    dzbar = zbar - rho';
    % Calculate Kalman gain matrix
    K = (P * H') / (H * P * H' + R);
    % Update state vector and covariance matrix
    xhat = xbar + K * dzbar;
    Phat = (eye(6) - K * H) * P;
    % Save state vector and covariance matrix
    xhat_c(epoch,:) = xhat;
    Phat_c(epoch,:,:) = Phat;
 
    % Initialize parametres for next integration
    Phi = eye(6);

end

% Save individual values of state vector
format long
x_c1 = xhat_c(10,:)
x_c2 = xhat_c(30,:)
x_c3 = xhat_c(50,:)
% Convert time to hours from t0
t_ax = (t - t(1)) / (60);
% Calculate residuals - m and m/s
resr_c = sqrt((xhat_c(:,1) - rx).^2 + (xhat_c(:,2) - ry).^2 + (xhat_c(:,3) - rz).^2) * 1e3;
resv_c = sqrt((xhat_c(:,4) - vx).^2 + (xhat_c(:,5) - vy).^2 + (xhat_c(:,6) - vz).^2) * 1e3;
% Create plots
figure(1)
subplot(1,2,1)
plot(t_ax,resr_c,'-b','LineWidth',1.5)
xlabel('time [min]');ylabel('\deltar [m]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])
ylim([0 inf])
subplot(1,2,2)
plot(t_ax,resv_c,'-r','LineWidth',1.5)
xlabel('time [min]');ylabel('\deltav [m/s]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])
ylim([0 inf])

%% SECTION D
%

% Calculate std devs - m and m/s
sigmar_c = sqrt(Phat_c(:,1,1) + Phat_c(:,2,2) + Phat_c(:,3,3)) * 1e3;
sigmav_c = sqrt(Phat_c(:,4,4) + Phat_c(:,5,5) + Phat_c(:,6,6)) * 1e3;
% Create plots
figure(2)
subplot(1,2,1)
plot(t_ax,sigmar_c,'-b','LineWidth',1.5)
xlabel('time [min]');ylabel('\sigma_r [m]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])
ylim([0 inf])
subplot(1,2,2)
plot(t_ax,sigmav_c,'-r','LineWidth',1.5)
xlabel('time [min]');ylabel('\sigma_v [m/s]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])
ylim([0 inf])

%% SECTION E
%

% Define initial conditions
% Initiale state vector
xhat = [rx(1) ry(1) rz(1) vx(1) vy(1) vz(1)]';
% Initial transition matrix
Phi = eye(6);
% Already defined covariance matrices
Phat = P0;
R = R0;
% Define process noise matrix
sigma_a = var(vecnorm(diff(xhat_c(:,1:3)),2,2))  % position noise - km
sigma_b = var(vecnorm(diff(xhat_c(:,4:6)),2,2))  % position noise - km/s
sigma_a = 0.0046;
% sigma_b = 0.01e-3;
Q = zeros(6,6);
row = [1 2 3];
col = [1 2 3];
sz = [6 6];
ind = sub2ind(sz,row,col);
Q(ind) = sigma_a^2;
row = [4 5 6];
col = [4 5 6];
sz = [6 6];
ind = sub2ind(sz,row,col);
Q(ind) = sigma_b^2;
% Initialise result matrices
xhat_e = zeros(num_epoch,6);
xhat_e(1,:) = xhat;
Phat_e = zeros(num_epoch,6,6);
Phat_e(1,:,:) = Phat;
e_rms = zeros(num_epoch,1);
e_rms(1) = NaN;

for epoch=2:num_epoch

    % Integration
    % Reshape to single coloumn
    Phi = reshape(Phi,36,1);
    g0 = [xhat; Phi];
    % Define timespan - epoch 1 and 2
    tspan = [t(epoch-1) t(epoch)];
    options = odeset('AbsTol',1e-08);
    [t_temp,g] = ode45(@(t,g) dy_dt_ecef(g,GM,WE), tspan, g0, options);
    % Extract results
    gf = g(end,:);
    xbar = gf(1:6)';
    Phi = reshape(gf(7:end),6,6);
    % Update covariance matrix
    P = Phi * Phat * Phi' + Q;

    % Update with observations
    % Find the visible GPS sats
    vis = find(PRN_ID(epoch,:));
    num_vis = size(vis,2);
    zbar = CA_range(epoch,vis)';
    rx_vis = rx_gps(epoch,vis);
    ry_vis = ry_gps(epoch,vis);
    rz_vis = rz_gps(epoch,vis);
    % Create opbservations covariance matrix
    idx_obs = eye(num_vis,num_vis);
    R = idx_obs.*sigma_obs^2;  % Observation covariance matrix
    % Calculate design matrix
    H = zeros(num_vis,6);        
    rho = sqrt((xbar(1)-rx_vis).^2 + (xbar(2)-ry_vis).^2 + (xbar(3)-rz_vis).^2);
    H(:,1) = (xbar(1) - rx_vis)./rho;
    H(:,2) = (xbar(2) - ry_vis)./rho;
    H(:,3) = (xbar(3) - rz_vis)./rho;
    % Calculate residual
    dzbar = zbar - rho';
    % Calculate Kalman gain matrix
    K = (P * H') / (H * P * H' + R);
    % Update state vector and covariance matrix
    xhat = xbar + K * dzbar;
    Phat = (eye(6) - K * H) * P;
    % Save state vector and covariance matrix
    xhat_e(epoch,:) = xhat;
    Phat_e(epoch,:,:) = Phat;
    % Calculate observation residual rms
    dxhat = K * dzbar;
    ebar = dzbar - H * dxhat;
    ebar_rms = rms(ebar);
    e_rms(epoch) = ebar_rms;
 
    % Initialize parametres for next integration
    Phi = eye(6);

end

% Calculate residuals - m and m/s
resr_e = sqrt((xhat_e(:,1) - rx).^2 + (xhat_e(:,2) - ry).^2 + (xhat_e(:,3) - rz).^2) * 1e3;
resv_e = sqrt((xhat_e(:,4) - vx).^2 + (xhat_e(:,5) - vy).^2 + (xhat_e(:,6) - vz).^2) * 1e3;
% Create plots
figure(3)
subplot(1,2,1)
plot(t_ax,resr_e,'-b','LineWidth',1.5)
xlabel('time [min]');ylabel('\deltar [m]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])
subplot(1,2,2)
plot(t_ax,resv_e,'-r','LineWidth',1.5)
xlabel('time [min]');ylabel('\deltav [m/s]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])

% Calculate std devs - m and m/s
sigmar_e = sqrt(Phat_e(:,1,1) + Phat_e(:,2,2) + Phat_e(:,3,3)) * 1e3;
sigmav_e = sqrt(Phat_e(:,4,4) + Phat_e(:,5,5) + Phat_e(:,6,6)) * 1e3;
% Create plots
figure(4)
subplot(1,2,1)
plot(t_ax,sigmar_e,'-b','LineWidth',1.5)
xlabel('time [min]');ylabel('\sigma_r [m]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])
subplot(1,2,2)
plot(t_ax,sigmav_e,'-r','LineWidth',1.5)
xlabel('time [min]');ylabel('\sigma_v [m/s]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])

%% SECTION F
%

% Convert observation residual to m
e_rms = e_rms * 1e3;
% Plot the observation residual RMS
figure(5)
subplot(1,2,1)
plot(t_ax,e_rms,'-k','LineWidth',1.5)
xlabel('time [min]');ylabel('e_{RMS} [m]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])
ylim([0 4])
subplot(1,2,2)
plot(t_ax,resr_e,'-b','LineWidth',1.5)
xlabel('time [min]');ylabel('\deltar [m]')
set(gca,'FontSize',14)
xlim([t_ax(1) t_ax(end)])
ylim([0 4])


