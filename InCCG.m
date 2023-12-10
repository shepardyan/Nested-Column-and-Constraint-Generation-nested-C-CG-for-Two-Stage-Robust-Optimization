function [objIn, d] = InCCG(rc, x, dt, varargin)
% Inner-level master problem implementation
[T, I, J, N, c, f, h, M, l, u, a, b] = dealRosteringCase(rc);

%% C&CG variables
inLB  = -inf;
inUB  =  inf;
delta = 1e-4;

%% Initialize inner-level scenario cuts array
y  = {zeros(J, T)};


%% Upper-level variables and constraints
g = binvar(1, T, 'full');
theta = sdpvar(1, 1);
constrsMP = [theta >= 0];
% Define uncertainty set
xi = 0.05 * dt;
d = dt + xi .* g;
if nargin == 4  % cardinality based uncertainty set (gamma)
    constrsMP = [constrsMP, sum(g) <= varargin{1}];
elseif nargin == 6  % cardinality based uncertainty set (T1, rho1, rho2)
    constrsMP = [constrsMP, sum(g(1:varargin{1}+2)) <= varargin{2}, ...
        sum(g(varargin{1}:end)) <= varargin{3}];
else
    error('parameters of inner-level master problem are wrong!')
end
objMP = -theta;

%% Initialize inner-level C&CG iteration procedure
k = 0;
constrsSP = {};
InCCG_constrs = [];
while inUB - inLB > delta
    k = k + 1;
    z{k} = sdpvar(J, T, 'full');
    w{k} = sdpvar(1, T, 'full');
    constrsMP = [constrsMP, theta <= sum(f.*y{k} + h.*z{k}, 'all') + sum(M.*w{k})];

    % Generate subproblem constraints
    constrsSP{k} = [];
    constrsSP{k} = [constrsSP{k}, z{k} <= N.* y{k}, ...
        N * sum(x, 1) + sum(z{k}, 1) + w{k} >= d, ...
        w{k} >= 0, z{k} >= 0];
    [KKTSystem, ~] = kkt(constrsSP{k}, sum(h.*z{k}, 'all') + sum(M.*w{k}), g);
    InCCG_constrs = [InCCG_constrs, KKTSystem, constrsSP{k}];
    % Optimize master problem
    optimize([InCCG_constrs, constrsMP], objMP);
    
    % Update upper bound
    inUB = -value(objMP);

    % Optimize subproblem and update lower bound
    [objSP, ySP] = InSP(rc, x, value(d)); 
    inLB = max([inLB, objSP]);

    % Update cuts
    y{k + 1} = ySP; 

    fprintf("  Inner Iteration %d, bound is %6.2f. UB is %6.2f, LB is %6.2f\n", k, inUB - inLB, inUB, inLB);

end

%% Return values
objIn = inLB;
d = value(d);
end