function [obj, x] = MP(rc, cuts)
% Outer-level master problem implementation
[T, I, J, N, c, f, h, M, l, u, a, b] = dealRosteringCase(rc);

%% Define constants
na = size(cuts, 2);

%% Variable definition
x = binvar(I, T, 'full');
y = binvar(J, T, na, 'full');
z = sdpvar(J, T, na, 'full');
w = sdpvar(1, T, na, 'full');
eta = sdpvar(1, 1);

%% Objective
obj = sum(c.*x, 'all') + eta;

%% Constraints
constrs = [];
% regular staff cannot work through any three consecutive shifts
for t = 1:T-2
    constrs = [constrs, x(:, t) + x(:, t + 1) + x(:, t + 2) <= 2];
end
% shift bounds
constrs = [constrs, l <= sum(x, 2), sum(x, 2) <= u];


for i = 1:na
    constrs = [constrs, eta >= sum(f.*y(:, :, i) + h.*z(:, :, i), 'all') + sum(M.* w(:, :, i))];
    % similar requirement on part-time staff
    for t = 1:T-1
        constrs = [constrs, y(:, t, i) + y(:, t + 1, i) <= 1];
    end
    constrs = [constrs, a <= sum(y(:, :, i), 2), sum(y(:, :, i), 2) <= b];
    % link binary activation and continuous working-hour decisions for part-time staff
    constrs = [constrs, z(:, :, i) <= N.*y(:, :, i)];
    % coverage of demands
    constrs = [constrs, N * sum(x, 1) + sum(z(:, :, i), 1) + w(:, :, i) >= cuts{i}];
    % variable bounds
    constrs = [constrs, w(:, :, i) >= 0, z(:, :, i) >= 0];
end

optimize(constrs, obj);

%% Return values
obj = value(obj);
x = value(x);

end