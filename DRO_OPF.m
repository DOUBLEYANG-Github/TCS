%% Note: This file solves the dual problem of the original DRO problem
%  Developed by Yang Yang at 25/8/2021.
%  The code is based on the IEEE-39 bus system.
%  Instead simply using the RSOME package, we provide here a direct reformulation
%  program for the original DRO problem.
%  See file DRO_reformulation.pdf for the mathmatical deduction.

%  As the variance upper bound 'sigma' increases, the overall reseve cost
%  will also increase. This is because more reserve is required to compensate for
%  the larger residue error. For instance, when sigma = 10 *
%  residue_variance, the overall reserve cost will increase by 10.

clear all; close all;
%% Case Parameter
IEEE39_Parameter;
sigma = 10 * gamma2 * residue_variance;  
%% The equivalence of the DRO model
%%Part I: construct the first stage constraints
disp('start first stage model');

% define first stage variables
Pg          = sdpvar(gen_num,1,'full');                     % generator output 
Pgs         = sdpvar(gen_num,3,'full');                     % generator output at each linearization block
PFR         = sdpvar(gen_num,1,'full');                     % PFR of each gernatiing unit 
pf          = sdpvar(line_num,1,'full');                    % power flow at each branch
theta       = sdpvar(node_num,1,'full');                    % angle at each node

% first layer constraints
constr = [];

% generation constraints
for i = 1:gen_num 
    for j = 1:3
        constr = constr + [0 <= Pgs(i,j)<= Pgs_max(i)];
    end
    constr = constr + [Pg(i)==sum(Pgs(i,:))+Pg_min(i)];
end

% power flow constr
constr = constr + [-l_max <= pf <= l_max];
for l = 1:line_num
     constr = constr + [pf(l)== line_map(l,:)*theta*admittance(l)]; % matrix multiplication    
end

% energy balance constr
for s=1:node_num
     constr = constr + [pf'*line_map(:,s) == Pg' * gen_map(:,s) - demand(s)];
end

% reserve contr
constr = constr + [0 <= PFR <= PFR_max];
constr = constr + [PFR + Pg <= Pg_max];
constr = constr + [ones(1,gen_num) * PFR >= min_gen_PFR];

% angle constr
constr = constr + [theta(slack_bus)==0];
constr = constr + [-pi<= theta <=pi];


%%Part II: construct the matrix for the second stage problem
% x = [PFR_g]            % size: 10 * 1
% z = [DSR_s]            % size: 2

x_num = gen_num;         % 10 generators
z_num = 2;               % 2 DSR
s_num = 1;               % one uncertain parameter

% Q(x,w) = min d * z, need to formulate d
d_sec = C_s;

% Gz <= h - Ex - Ms, need to formulate G, h , E, M, this will be done by enumarating constraints
% Constraint 1: - DSR_s <= 0 
G1 = -eye(z_num); h1 = zeros(z_num,1); E1 = zeros(z_num,x_num); M1 = zeros(z_num,1);

% Constraint 2: DSR_s <= DSR_max
G2 = eye(z_num); h2 = DSR_max'; E2 = zeros(z_num,x_num); M2 = zeros(z_num,1);

% Constraint 3: - \sum DSR_s <= - min_total_reserve - (-\sum_g PFR_g) 
G3 = -ones(1,z_num); h3 = - min_total_reserve; E3 = -ones(1,x_num); M3 = 0;

% Constraint 4:  coefficient * DSR <= rho * FD_max - intercept - Coefficient * PFRg  - s 
G4 = slope(11:12)'; h4 = rho * FD_max - intercept; E4 = slope(1:10)'; M4 = 1;

% Construct matrices
G = [G1;G2;G3;G4];
h = [h1;h2;h3;h4];      
E = [E1;E2;E3;E4];
M = [M1;M2;M3;M4];

constr_num = size(h,1);

% define second stage decision variables
disp('start building second stage model');

% scalar form
alpha_d  = sdpvar(1,1,'full');
beta_d   = sdpvar(s_num,1,'full');
delta_d  = sdpvar(s_num,1,'full');
lambda_d = sdpvar(s_num,1,'full');
eps_d    = sdpvar(s_num,1,'full');
kappa_d  = sdpvar(s_num,1,'full');
pai_d    = sdpvar(s_num,1,'full');
omega_d    = sdpvar(s_num,1,'full');
eta_d    = sdpvar(s_num,1,'full');

% matrix form
delta_dm  = sdpvar(s_num,constr_num,'full');
eps_dm    = sdpvar(s_num,constr_num,'full');
kappa_dm  = sdpvar(s_num,constr_num,'full');
pai_dm    = sdpvar(s_num,constr_num,'full');
lambda_dm = sdpvar(s_num,constr_num,'full');
eta_dm    = sdpvar(s_num,constr_num,'full');

% define affine parameters
z_0       = sdpvar(z_num,1,'full');
z_xi      = sdpvar(z_num,1,'full');
z_v       = sdpvar(z_num,1,'full');

% second stage constraints
constr1 = [];
constr1 = constr1 + [lambda_d >= 0];
constr1 = constr1 + [delta_d >= 0];
constr1 = constr1 + [eps_d >= 0];
constr1 = constr1 + [omega_d >= 0];

% alpha - d^T z_0 >= - xi_{min}^T delta + xi_{max}^T eps - 1^T kappa +
% nu_{max}^T * lambda
nu_max = 20;
constr1 = constr1 + [alpha_d - d_sec' * z_0 >= - xi_min'*delta_d + xi_max'* eps_d - kappa_d + pai_d + nu_max *lambda_d];
% eps - delta + 2 eta = z_xi^T d - beta
constr1 = constr1 + [eps_d - delta_d - 2*eta_d == z_xi'* d_sec - beta_d];
% -kappa - pai + pho = z_v^T d - lambda
constr1 = constr1 + [-kappa_d - pai_d + lambda_d == z_v'* d_sec - omega_d];
% socp constraint ||eta_d, kappa_d||2 <= pai_d
constr1 = constr1 + [cone([eta_d,kappa_d],pai_d)];

% dual constraint for each secondary stage constraint
constr1 = constr1 + [delta_dm >= 0];
constr1 = constr1 + [lambda_dm >= 0];
constr1 = constr1 + [eps_dm >= 0];

for j = 1: constr_num
    % (h - Ex - G z_0)_k >= - xi_{min}^T delta_k + xi_{max}^T eps_k - 1^T
    % kappa_k + 1^T pai_k + nu_{max}^T * pho_k
    constr1 = constr1 + [h(j)-E(j,:)* PFR - G(j,:)* z_0 >= - xi_min * delta_dm(j) + xi_max * eps_dm(j) - kappa_dm(j) + pai_dm(j) + nu_max * lambda_dm(j)];
    constr1 = constr1 + [eps_dm(j) - delta_dm(j) - 2*eta_dm(j) == G(j,:)*z_xi + M(j)];
    constr1 = constr1 + [-kappa_dm(j) - pai_dm(j) + lambda_dm(j) == G(j,:)*z_v];
    % socp constraint ||eta_dm(k),kappa_dm(k)|| <= pai_dm(k);
    constr1 = constr1 + [cone([eta_dm(j),kappa_dm(j)],pai_dm(j))];
end

% objective function
totalcost1 = 0;  % energy cost
totalcost2 = 0;  % gen reserve cost
totalcost3 = 0;  % expcted DSR reserve cost
for i= 1:gen_num
    for m = 1:3
        totalcost1 = totalcost1 + Pgs(i,m) * k(i,m);    
    end 
end
totalcost2 =  C_g' * PFR; 
totalcost3 =  alpha_d + omega_d'*sigma;
totalcost_all = totalcost1+totalcost2+totalcost3;
% optimization model setting

constr_all = constr + constr1;
ops = sdpsettings('solver','cplex','verbose',0,'debug',1); 
disp('send to solver')
result = optimize(constr_all,totalcost_all,ops);
% disp('solution status:')
% result.info
%disp('solvertime:')
%result.solvertime
%disp('gen cost:')
%value(totalcost1)
disp('gen reserve cost:')
value(totalcost2)
disp('expected reserve cost:')
value(totalcost2 + totalcost3)
z_0 = value(z_0);
%calculate the expected reserve
z_expect = value(z_0)+value(z_v)*sigma;

disp('exptected DSR reserve:')
DSR_expected = sum(z_expect)
gen_reserve = sum(value(PFR))

