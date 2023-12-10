function [obj, x, y] = rosteringRO(rc, dt, varargin)
% Outer-level master problem implementation
[T, I, J, N, c, f, h, M, l, u, a, b] = dealRosteringCase(rc);

%% C&CG variables
LB  = -inf;
UB  =  inf;
delta = 1e-4;

%% Find a initial cut
cuts = {dt};

%% Outer-level C&CG
k = 0;
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
    fprintf("Outer Iteration %d, bound is %6.2f. UB is %6.2f, LB is %6.2f\n", k, UB - LB, UB, LB);
end
%% Return values
obj = UB;
x = value(x);
[~, y] = InSP(rc, x, cuts{k});
end