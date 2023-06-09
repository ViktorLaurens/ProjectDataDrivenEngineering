clear; close all; clc
%% Initilization
value       = 3;                    % the number of initial conditions

% Initial angular velocity of the ball
theta_dot0  = (-pi/2 + (pi/2+pi/2)*rand(1,value))*0.01;

% Initial Position of the ball
theta0      = -pi/2 + (pi/2+pi/2)*rand(1,value);

% Initial Condition Vector
IC          = [theta0;theta_dot0];

% Earth, Mercury, Venus, Moon
g           = [9.8100 3.7278 8.8290 1.6677];


% Masses different balls
% Tennis, baseball, volleyball, rugby, football, bowling
m           = [59.4 149 270 420 500 7260]/1000;

L           = linspace(0.5,1.5,3);    % Different length of the rope

b           = linspace(0.01,0.2,4);   % Different damping coefficients

RESULT      = [];                   % Empty varibale for saving the informations

ball_name   = ["tennis" "baseball" "volleyball" "rugby" "football" "bowling"];
header      = {'X_coordinate [m]','Y_Coordinate[m]','Position [rad]',...
               'Angular Velocity','Time [s]',...
               'Initial Position [rad]','Initial Velocity[rad/s]',...
               'Damping Coefficient [Ns/m]','g [m/s^2]',...
               'Rope Lenght [m]','Ball Mass [kg]'};
% ***********************************************************************%
%% Dataset Creation 
baseFileName = 'Pendulum1.csv';

for ii = 1:length (m)                       % different masses 
    for jj = 1:value                        % different initial position
        for kk = 1:value                    % different initial angluar velocity
            for ll = 1:length (g)           % different g
                for mm = 1:length(L)  % Different length of the rope
                    for nn = 1:length(b)    % Different damping

                        t_span     = linspace(0,5,500);
                        [t,result] = ode45(@(t,theta)odefun(t,theta,b(nn),g(ll),L(mm),m(ii)),t_span,IC(:,jj));
                        
                        outlierIndices = round(linspace(1, 500, 100));
                        outlierValues = rand(size(outlierIndices));
                        result(outlierIndices, 1) = outlierValues;

                        x_position = L(mm)*cos(result(:,1));
                        y_position = L(mm)*sin(result(:,1));
                        
                        % Store the results in RESULT matrix
                        C = [x_position,y_position,t,result(:,1),...
                             result(:,2),theta0(jj)*ones(length(x_position),1),...
                             theta_dot0(kk)*ones(length(x_position),1),...
                             b(nn)*ones(length(x_position),1),...
                             g(ll)*ones(length(x_position),1),...
                             L(mm)*ones(length(x_position),1),...
                             m(ii)*ones(length(x_position),1)];
%                         *ones(length(x_position),1)
                        RESULT = [RESULT;C];


                    end                     % End of b
                end                         % End of L
            end                             % end of g
        end                                 % end of theta_dot_0
    end                                     % end of theta_0
end                                         % end of m
%%
% Write the results to an Excel file
% writecell([header; num2cell(RESULT)], baseFileName, 'Sheet', 1);
writecell([header; num2cell(RESULT)], baseFileName);
%%
% plot(RESULT(1:250,3),RESULT(1:250,4))
% plot(RESULT(1:250,3),RESULT(1:250,4))
% hold on
% plot(RESULT(1251:1500,3),RESULT(1251:1500,4))
% plot(RESULT(2251:2500,3),RESULT(2251:2500,4))
% plot(RESULT(3251:3500,3),RESULT(3251:3500,4))
% plot(RESULT(4251:4500,3),RESULT(4251:4500,4))