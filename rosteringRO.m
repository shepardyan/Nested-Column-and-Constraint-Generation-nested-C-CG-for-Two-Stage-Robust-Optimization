function [obj, x, y] = rosteringRO(rc, dt, varargin)
% Outer-level master problem implementation
[T, I, J, N, c, f, h, M, l, u, a, b] = dealRosteringCase(rc);
%% Start info
if isfield(rc, 'logName') && ~isempty(rc.logName)
    fprintf(rc.logName, "Nested C&CG starts\n");
else
    fprintf("Nested C&CG starts");
end

%% C&CG variables
LB  = -inf;
UB  =  inf;
delta = 1e-3;

%% Find a initial cut
cuts = {dt};

%% Outer-level C&CG
k = 0;
tic;
while UB - LB > delta
    k = k + 1;
    % Solve outer problem
    [objOu, x] = MP(rc, cuts);
    % Update lower bound
    LB = objOu;
    % Solve inner problem
    if nargin == 3
        [objIn, d] = InCCG(rc, x, dt, varargin{1});
    elseif nargin == 5
        [objIn, d] = InCCG(rc, x, dt, varargin{1}, varargin{2}, varargin{3});
    else
        error("Param wrong!")
    end
    % Update cut set
    cuts{k + 1} = d;
    % Update upper bound
    UB = min([UB, objIn + sum(c.*x, 'all')]);
    if isfield(rc, 'logName') && ~isempty(rc.logName)
        fprintf(rc.logName, "Outer Iteration %2d, bound is %10.2f. UB is %10.2f, LB is %10.2f\n", k, UB - LB, UB, LB);
    else
        fprintf("Outer Iteration %2d, bound is %10.2f. UB is %10.2f, LB is %10.2f\n", k, UB - LB, UB, LB);
    end
end
%% Return values
obj = UB;
x = value(x);
[~, y] = InSP(rc, x, cuts{k});

elapsed_time = toc;
if isfield(rc, 'logName') && ~isempty(rc.logName)
    fprintf(rc.logName, "Nested C&CG: Time elapsed %.2f seconds\n", elapsed_time);
else
    fprintf("Nested C&CG: Time elapsed %.2f seconds\n", elapsed_time);
end
end