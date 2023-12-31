%% Two-stage Robust Rostering Problem
% L. Zhao and B. Zeng, "An exact algorithm for two-stage robust optimization
% with mixed integer recourse problems," submitted, available on Optimization
% Online. org, 2012.

% In the rostering problem, an organization such as a call center or a 
% clinic needs to allocate its staff members to shifts to meet service 
% demands with a minimized operation cost and satisfy governmental or 
% industrial regulations/restrictions. Sometimes, to deal with demand 
% surges, overtime from its regular staff or part-time staff (or agency 
% staff) will be called.
clc; clear
%% Rostering problem struct definition
rostcase.T = 21; % Total number of shifts
rostcase.I = 12; % Full time staff index
rostcase.J = 3; % Part time staff index. Working hours are to be determined
rostcase.N = 8; % Work length of full time staff

% cost parameter
rostcase.c = 10 * rand(rostcase.I, rostcase.T) + 5;
rostcase.f = 10 * rand(rostcase.J, rostcase.T) + 20;
rostcase.h = 4 * rand(rostcase.J, rostcase.T) + 4;
rostcase.M = 10 * rand(1, rostcase.T) + 40;

% shift bounds
rostcase.l = randi([4, 8], rostcase.I, 1);
rostcase.u = randi([8, 14], rostcase.I, 1);
rostcase.a = randi([2, 4], rostcase.J, 1);
rostcase.b = randi([4, 6], rostcase.J, 1);
dt = randi([30, 80], 1, rostcase.T);

% Choose to save log file
time_now = string(datetime("now", 'Format', "yyyy-MM-dd-HH-mm-ss"));
f = fopen(string("./log/"+time_now+".txt"), 'w+');
fclose(f);
rostcase.logName = fopen(string("./log/"+time_now+".txt"), 'a+');

%% Deterministic Rostering Problem
[~, xMIP, yMIP] = rosteringMIP(rostcase, dt);

%% Two-stage Robust Rostering Problem
dt = 50 + 30 * rand(1, rostcase.T);
gamma = 3;
[objRO, xRO, yRO] = rosteringRO(rostcase, dt, gamma);

