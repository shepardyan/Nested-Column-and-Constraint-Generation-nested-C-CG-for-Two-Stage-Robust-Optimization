%% Deterministic Rostering Problem
% service demand
T = 21;
d = 50 * rand(1, T) + 25;  % d(t) is the known service demand
% staff
I = 12; % Full time staff index
J = 3; % Part time staff index. Working hours are to be determined
N = 8; % Work length of full time staff

% cost parameter
c = 10 * rand(I, T) + 5;
f = 10 * rand(J, T) + 20;
h = 4 * rand(J, T) + 4;
M = 10 * rand(1, T) + 40;

% shift bounds
l = randi([4, 8], I, 1);
u = randi([8, 14], I, 1);
a = randi([2, 4], J, 1);
b = randi([4, 6], J, 1);

%% Variable definition
x = binvar(I, T, 'full');
y = binvar(J, T, 'full');
z = sdpvar(J, T, 'full');
w = sdpvar(1, T, 'full');

%% Constraints
constrs = [];
% regular staff cannot work through any three consecutive shifts
for t = 1:T-2
    constrs = [constrs, x(:, t) + x(:, t + 1) + x(:, t + 2) <= 2];
end
% similar requirement on part-time staff
for t = 1:T-1
    constrs = [constrs, y(:, t) + y(:, t + 1) <= 1];
end
% shift bounds
constrs = [constrs, l <= sum(x, 2), sum(x, 2) <= u];
constrs = [constrs, a <= sum(y, 2), sum(y, 2) <= b];
% link binary activation and continuous working-hour decisions for part-time staff
constrs = [constrs, z <= N.*y];
% coverage of demands
constrs = [constrs, N * sum(x, 1) + sum(z, 1) + w >= d];
% variable bounds
constrs = [constrs, w >= 0, z >= 0];
%% Objective
obj = sum(c.*x, 'all') + sum(f.*y + h.*z, 'all') + sum(M.*w);

%% Optimization
optimize(constrs, obj);
