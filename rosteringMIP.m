function [obj, x, y] = rosteringMIP(rc, d)
% Deterministic Rostering Problem
[T, I, J, N, c, f, h, M, l, u, a, b] = dealRosteringCase(rc);

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
obj = value(obj);
x = value(x);
y = value(y);
end