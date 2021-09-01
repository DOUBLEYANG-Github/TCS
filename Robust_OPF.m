%% Note: This file solves the robust OPF problem
%  Developed by Yang Yang at 25/8/2021.
%  The code is based on the IEEE-39 bus system.
%  The Yalmip package should be installed before running the code.
%  The regression residue is taken as the worst one 'xi_max', i.e., requires more reserve
%  to compensate for the uncertain residue.
%  For the deterministic OPF model, simply set the 'xi_max' as zero.

clear all; close all;
%% Case Parameter
IEEE39_Parameter;

%xi_max = 0;                                                 % change to deterministic model

%% The robust OPF model
% define first stage variables
Pg          = sdpvar(gen_num,1,'full');                     % generator output 
Pgs         = sdpvar(gen_num,3,'full');                     % generator output at each linearization block
PFR         = sdpvar(gen_num,1,'full');                     % PFR of each gernatiing unit 
pf          = sdpvar(line_num,1,'full');                    % power flow at each branch
theta       = sdpvar(node_num,1,'full');                    % angle at each node
DSR         = sdpvar(DSR_num,1,'full');  
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

% generator reserve contr
constr = constr + [0 <= PFR <= PFR_max'];
constr = constr + [PFR + Pg <= Pg_max];
constr = constr + [ones(1,gen_num) * PFR >= min_gen_PFR];

% DSR reserve constr
constr = constr + [0 <= DSR <= DSR_max'];

% frequency security constr
constr = constr + [ones(1,gen_num) * PFR + ones(1,DSR_num) * DSR >= min_total_reserve];
constr = constr + [slope' * [PFR;DSR] + intercept + xi_max <= rho * FD_max];

% angle constr
constr = constr + [theta(slack_bus)==0];
constr = constr + [-pi<= theta <=pi];

% objective function 
totalcost1 = 0;
for i= 1:gen_num
    for m = 1:3
        totalcost1 = totalcost1 + Pgs(i,m) * k(i,m);    % energy cost
    end 
end
totalcost2 =  C_g' * PFR;                               % gen reserve cost
totalcost3 =  C_s' * DSR;                               % DSR reserve cost

totalcost_all = totalcost1+totalcost2+totalcost3;
% optimization model setting
ops = sdpsettings('solver','cplex','verbose',0,'debug',1); 
disp('send to solver')
result = optimize(constr,totalcost_all,ops);
%disp('solution status:')
%result.info
%disp('solvertime:')
%result.solvertime
disp('gen cost:')
value(totalcost1)
disp('gen reserve cost:')
value(totalcost2)
disp('reserve cost:')
value(totalcost2 + totalcost3)