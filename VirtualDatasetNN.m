clc
clear
close all

%% DATASET 1 (Creating dataset using analytic solution)
% Define constants
G = [9.81,20];                % acceleration due to gravity [m/s^2]
L = [0.1,0.2];           % length of pendulum [m]
M = [0.1];              % mass of pendulum [kg]
C = [0.1];     % damping coefficient [kg*m^2/s]
THETA0 = [75]./180*pi;    % initial angle [radians]
OMEGA0 = [0];                     % initial angular velocity [radians/s]
samples = 1;

% Add Guassian noise
mu = 0;
sigma = 0;

% Define time specificities
dt = 0.01; % time step
t_end = 100; % [s]
tspan = [0:dt:t_end]; % time interval [s]

% Allocate memory for the data vectors X and Y
% X = zeros(2*length(tspan),length(C)*length(M)*length(g)*length(L)*length(THETA0)*length(OMEGA0)*samples);
X = zeros(length(tspan),length(C)*length(M)*length(G)*length(L)*length(THETA0)*length(OMEGA0)*samples);
Y = zeros(5,length(C)*length(M)*length(G)*length(L)*length(THETA0)*length(OMEGA0)*samples);

% PCA: L, M, C, theta0, omega0, g, H1, H2, H3, T1, T2
Y_PCA = zeros(11,length(C)*length(M)*length(G)*length(L)*length(THETA0)*length(OMEGA0)*samples);

tic
sample = 1;
for j_G = 1:length(G)
    g = G(j_G);
for j_OMEGA0 = 1:length(OMEGA0)
    omega0 = OMEGA0(j_OMEGA0);
for j_THETA0 = 1:length(THETA0)
    theta0 = THETA0(j_THETA0);
for j_C = 1:length(C)
    c = C(j_C);
for j_M = 1:length(M)
    m = M(j_M);
for j_L = 1:length(L)
    l = L(j_L);
for j_sample = 1:samples
    % Define initial conditions
    theta_init = [theta0; omega0];
    
    % Define differential equation
    f = @(t, theta) [theta(2); -g/l*sin(theta(1))-c/m*theta(2)];
    
    % Solve differential equation
    [t, theta] = ode45(f, tspan, theta_init);
    theta = theta(:,1); % take only the angles not the angular velocities
    theta = theta + mu + sigma.*randn(length(theta),1);     % add gaussian noise to theta
%     x = l * sin(theta);
%     y = l * cos(theta);
% 
%     % Create a row vector with alternating elements of x and y
%     xy = zeros(length(x) + length(y),1);
%     xy(1:2:end) = x;
%     xy(2:2:end) = y;
% 
%     X(:,sample) = xy;
    X(:,sample) = theta;
    Y(1,sample) = l;
    Y(2,sample) = m;
    Y(3,sample) = c;
    Y(4,sample) = theta0;
    Y(5,sample) = omega0;

    Y_PCA(1,sample) = l;
    Y_PCA(2,sample) = m;
    Y_PCA(3,sample) = c;
    Y_PCA(4,sample) = theta0;
    Y_PCA(5,sample) = omega0;
    Y_PCA(6,sample) = g;
    [pks, locs] = findpeaks(X(:,sample));
    Y_PCA(7,sample) = pks(1);
    Y_PCA(8,sample) = pks(2);
    Y_PCA(9,sample) = pks(3);
    Y_PCA(10,sample) = (locs(2)-locs(1))*dt; % Time between first and second peak
    Y_PCA(11,sample) = (locs(3)-locs(2))*dt;
    
    sample = sample + 1;
end
end
end
end
end
end
end
toc

X_L = X;
Y_L = Y;

% Extracting the data for 1 sample/observation
observation = 4;

figure
plot(tspan,X(:,observation)/pi*180)
xlabel('t (s)')
ylabel('theta (deg)')
title(['Single pendulum with friction (observation: ',num2str(observation),')'])
grid on

