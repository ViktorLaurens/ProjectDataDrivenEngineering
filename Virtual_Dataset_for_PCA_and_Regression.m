clc
clear
close all

%% DATASET 1 (Creating dataset using analytic solution)
% Define constants
% Gravity constants 
% G = [9.81;     % Earth
%      0.41;     % Moon
%      26.0;];   % Sun

G = [9.81];

% length of pendulum [m]
% L = [0.25; 0.5; 0.75]; 
L = [0.5]; 

% mass of pendulum [kg]
% M = [0.057;    % Tennis ball
%      0.42;     % Soccer ball (standard size 5)
%      0.17];    % Basketball
M = [0.057];
               
C = [0.05];
% C = [0.05;   % Damping coefficient 1 [kg*m^2/s]
%      0.1;    % Damping coefficient 2 [kg*m^2/s]
%      0.15];  % Damping coefficient 5 [kg*m^2/s]


% THETA0 = [75; 45; 30; 15; 5]./180*pi;     % initial angle [radians]
% OMEGA0 = [0; 0.05; 0.1];              % initial angular velocity [radians/s]
THETA0 = [45]./180*pi;
OMEGA0 = [0.02]; 
samples = 300;

% Add Guassian noise
mu = 0;
sigma = 0;
noise_level_g = 0;  % 1% 

% Define time specificities
dt = 0.01; % time step
t_end = 100; % [s]
tspan = [0:dt:t_end]; % time interval [s]

% Allocate memory for the data vectors X and Y
% X = zeros(2*length(tspan),length(C)*length(M)*length(g)*length(L)*length(THETA0)*length(OMEGA0)*samples);
X = zeros(length(tspan),length(C)*length(M)*length(G)*length(L)*length(THETA0)*length(OMEGA0)*samples);
Y = zeros(5,length(C)*length(M)*length(G)*length(L)*length(THETA0)*length(OMEGA0)*samples);

% PCA: L, M, C, theta0, omega0, g, H1, H2, H3, T1, T2
Y_PCA = zeros(16,length(C)*length(M)*length(G)*length(L)*length(THETA0)*length(OMEGA0)*samples);

tic
sample = 1;
% for j_G = 1:length(G)
%     g = G(j_G) + noise_level_rel * abs(G(j_G)) * randn();
% for j_OMEGA0 = 1:length(OMEGA0)
%     omega0 = OMEGA0(j_OMEGA0) + noise_level_rel * abs(OMEGA0(j_OMEGA0)) * randn();
% for j_THETA0 = 1:length(THETA0)
%     theta0 = THETA0(j_THETA0) + noise_level_rel * abs(THETA0(j_THETA0)) * randn();
% for j_C = 1:length(C)
%     c = C(j_C) + noise_level_rel * abs(C(j_C)) * randn();
% for j_M = 1:length(M)
%     m = M(j_M) + noise_level_rel * abs(M(j_M)) * randn();
% for j_L = 1:length(L)
%     l = L(j_L) + noise_level_rel * abs(G(j_G)) * randn();
% for j_sample = 1:samples
for j_sample = 1:samples
for j_G = 1:length(G)
    g = (0.3 + 3*rand())*G(j_G);
for j_OMEGA0 = 1:length(OMEGA0)
    omega0 = randn()*OMEGA0(j_OMEGA0);
for j_THETA0 = 1:length(THETA0)
    theta0 = (1 + 0.5*rand())*THETA0(j_THETA0);
for j_C = 1:length(C)
    c = (1 + 2*rand())*C(j_C);
for j_M = 1:length(M)
    m = (1 + 4*rand())*M(j_M);
for j_L = 1:length(L)
    l = (1 + rand())*L(j_L);
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
    Y_PCA(4,sample) = g;
    Y_PCA(5,sample) = c/m;
    Y_PCA(6,sample) = g/l;
    Y_PCA(7,sample) = theta0;
    Y_PCA(8,sample) = omega0;
    
    [pks, locs] = findpeaks(X(:,sample));
    Y_PCA(9,sample) = pks(1);
    Y_PCA(10,sample) = pks(2);
    Y_PCA(11,sample) = pks(3);
    Y_PCA(12,sample) = (locs(2)-locs(1))*dt; % Time between first and second peak
    Y_PCA(13,sample) = (locs(3)-locs(2))*dt;
    
    % Compute logarithmic decrement
    delta = log(pks(1)/pks(2));
    Y_PCA(14,sample) = delta;
    
    % Define the limit for peak detection
    peak_limit = 5*pi/180; % Adjust the value as needed
    
    % Find the time and number of oscillations to reach the peak limit
    peak_indices = find(pks < peak_limit);
    upperpeaks = find(pks > peak_limit);
    if ~isempty(peak_indices)
        T = locs(peak_indices(1)) * dt; % Time to reach the peak limit
        N = length(upperpeaks); % Number of oscillations to reach the peak limit
    else
        T = NaN;
        N = NaN;
    end
    
    Y_PCA(15,sample) = T;
    Y_PCA(16,sample) = N;
    
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
Y_PCA = real(Y_PCA)



%% Plotting
observations = [1, 2, 3]; % Specify the indices of the observations

figure
xlabel('$t$ [s]', 'Interpreter', 'latex')
ylabel('$\theta$ [deg]', 'Interpreter', 'latex')
hold on

for i = 1:numel(observations)
    observation = observations(i);
    plot(tspan, X(:,observation)/pi*180)
end

legend(sprintf('Motion %i: $m=%.2f$, $L=%.2f$, $\\theta_0=%.2f$, $\\omega_0=%.2f$', 1, Y_PCA(2,observations(1)), Y_PCA(1,observations(1)), Y_PCA(4,observations(1)), Y_PCA(5,observations(1))),...
       sprintf('Motion %i: $m=%.2f$, $L=%.2f$, $\\theta_0=%.2f$, $\\omega_0=%.2f$', 2, Y_PCA(2,observations(2)), Y_PCA(1,observations(2)), Y_PCA(4,observations(2)), Y_PCA(5,observations(2))),...
       sprintf('Motion %i: $m=%.2f$, $L=%.2f$, $\\theta_0=%.2f$, $\\omega_0=%.2f$', 3, Y_PCA(2,observations(3)), Y_PCA(1,observations(3)), Y_PCA(4,observations(3)), Y_PCA(5,observations(3))), 'Interpreter', 'latex');
title('Single pendulum with friction (multiple observations)')
grid on



