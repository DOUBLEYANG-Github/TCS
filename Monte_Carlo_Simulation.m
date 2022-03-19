%% Note: This program calls the SIMULINK function to generate data and solves the linear regression model

%  In the IEEE 39-bus system, there are 6 generators. We assume all the
%  Inverter Air Condiioner (IAC) in the system are aggregated as one
%  demand-side resource (DSR) to provide primary frequency reserve to the system.

%  The output limit of the 10 generators are: [1040,646,725,652,508,687,580,564,865,1100] MW. In total, 7367 MW. 

%  The overall load demand is 6,254 MW. We consider a maximum power imbalance 'dP' that is 10% of the total load, i.e.,
%  625 MW.

%  Meanwhile, the maximum reserve capacity of each generator is set as 10%,
%  i.e. [100,60,70,65,50,70,60,55,85,110] MW, in total 1500 MW. The maximum
%  reserve provided by a DSR is assumed to be 100 MW. We assume the there 
%  are two DSRs in the system. 

%  The load damping parameter "D" is set as 2% of the total load demand, i.e. 125 MW/Hz.

%  Parameters for the IACs are chosen according to a realistic system in Haining,
%  China, reported by reference: X. Zhuang et. al, "Data-Driven Reserve Allocation With Frequency Security 
%  Constraint Considering Inverter Air Conditioners," in IEEE Access, vol. 7, pp. 120014-120022, 2019

%  Note that due to the random number generator, the simulation results of each execution may
%  differ !!!  


clear all; close all;

%% System Parameter for primary frequency regulation process
dP = -300;    % 5% * 6254 MW
D  = 125;     % 2% * 6254 MW
H  = 300;     % RoCoF = dP/2H >= 0.5 Hz/s => H = 1000 MW*Hz/s

%% Gen parameters
gen_num = 10;

R   = 0.002* ones(gen_num);                              % governor speed regulation (Hz/MW)
TG  = [0.3;0.4;0.5;0.3;0.4;0.5;0.5;0.3;0.4;0.5]*40;      % time constant of generator (s)
TCH = 0.6 * ones(10);                                    % steam chest reaction time (s)
TRH = 5;                                                 % reheat turbine time constant (s)
FHP = 0.8;                                               % high-pressure power fraction of reheat turbine

PFR_max = [100,60,70,65,50,70,60,55,85,110];             % MW

%% IAC parameters
Ta  = 5;                                                 % temperature difference between room and ambient air (K)
AC  = 0.01;                                              % IAC control constant

A   = 80;                                                % speed regulator, similar to 1/R of generators (kW/Hz)  
Tc  = 0.27;                                              % compressor time constant of IAC (s)
N   = 6000/1000;                                         % 6000 IACs, divide by 1000 to change KW to MW

DSR_max = 60;                                            % MW
DSR_num = 1;       

%% Battery parameters
TB = 0.5;
B_max = DSR_max;
B_num = 1;

%% simulation
options = simset('SrcWorkspace','current');
sim_time = 60;

sim_num = 11000;
Input = zeros(sim_num,12);
Integration = zeros(sim_num,1);

i = 1;
while i <= sim_num
    Re = PFR_max .* rand(1,10);    % upper bound of Gen reserve
    Pa = DSR_max * rand(1,1);      % upper bound of IAC reserve
    Pb = B_max * rand(1,1);        % upper bound of battery reserve
    % We nned to guarantee adequate reserve, so that the quasi-steady state frequency
    % can satisfy the requirement. delta_f_qss = (dP + reserve)/D >= -0.5,
    % i.e.,  reserve >= -dP*0.8;
    % Meanwhile, we need to leave some reserve margin for the secondary frequency
    % regulation process, i.e., reserve <= -dP*0.95.
    
    reserve = sum(Re)+Pa+Pb;
    if reserve >= -dP*0.8 &&  reserve <= -dP*0.95       
       [T,X,Y] =sim('Simulink_IEEE_39_two_DSR',sim_time,options);
       Frequency = Y(:,1);
       Input(i,:) = [Re,Pa,Pb];
       Integration(i) = -mean(Frequency)*60; 
       i = i+1
    end
    
end

disp('simulation is finished!')

%% Linear regression
dataset = [Input, Integration];

disp('linear regression coefficients and Hypothesis testing based on T-statistics:')
mdl = fitlm(Input,Integration);
mdl.Coefficients
Coefficient = mdl.Coefficients.Estimate;
%save('Coefficient.mat','Coefficient');

if mdl.coefTest <= 0.01 && max(mdl.Coefficients.pValue) <=0.01
    disp('The linear regression model is statistically significant!')
elseif mdl.coefTest >= 0.05 || max(mdl.Coefficients.pValue) >=0.05
    disp('Warning: The linear regression model may not be statistically significant!')
end
%scatter3(Re_input(:,1),Input(:,2),Integration,'.')

disp('average integral frequency deviation is:')
mean(abs(Integration))

disp('largest estimation error is:')
max_regression_residue = max(abs(mdl.Residuals.Raw))
%save('max_regression_residue.mat','max_regression_residue');

disp('the mean square error is:')
residue_variance = var(mdl.Residuals.Raw)
%save('residue_variance.mat','residue_variance');

%csvwrite('simulation_large_10000.csv',dataset)

%% 10-fold Cross Validation to check the effectiveness of the linear regression

n = length(dataset);
k = 10;   
cv = cvpartition(n,'KFold',k);                % create an object for cvpartition class

max_rel_estimation_error = zeros(k,1);
mean_rel_estimation_error = zeros(k,1);

for i = 1:k
    test_classes = dataset(cv.test(i),:);     % cv.test(i) returns the label of test set     
    train_classes = dataset(~cv.test(i),:);
    
    X_train = train_classes(:,1:gen_num+DSR_num+B_num);
    Y_train = train_classes(:,gen_num+DSR_num+B_num+1);
    X_test  = test_classes(:,1:gen_num+DSR_num+B_num);
    Y_test  = test_classes(:,gen_num+DSR_num+B_num+1);
    
    mdl = fitlm(X_train,Y_train);            % fit the linear regression model
    Y_pred = predict(mdl,X_test);
    
    % the maximum estimation error and the mean estimation error
    max_rel_estimation_error(i) = max(abs(Y_test - Y_pred)./abs(Y_test)); 
    mean_rel_estimation_error(i) = mean(abs(Y_test - Y_pred)./abs(Y_test)); 
end 

figure;
plot(max_rel_estimation_error)
figure;
plot(mean_rel_estimation_error)
