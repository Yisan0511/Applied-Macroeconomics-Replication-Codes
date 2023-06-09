clear;
%% 0. Preparation
%% 0.1. Adding Paths
addpath(genpath("VAR-Toolbox-main/v3dot0"));
addpath(genpath("Kilian Toolbox"));

%% 0.2. Loading Data
load Data/kiliandata.txt

% Please use "readtable('Data/')" to replace loading for reading the new data
% The new dataset should be stored as a obs*3(long) format.

%% 1. Estimating VAR Model by OLS
det = 1;
nlags = 24;
Y = kiliandata;
[VAR, VARopt] = VARmodel(Y, nlags, det);

%% 2. Identification
%% 2.1. Recursive VAR
%% 2.1.1. Cholesky Decompostion
Chol_B0inv = chol(VAR.sigma, 'lower'); 
disp(Chol_B0inv);
%% 2.1.2. Cholesky IRF
h = 16;

Chol_IRF = irfvar(VAR.Fcomp,Chol_B0inv,nlags,3,h-1);

for i = 1:9
    Chol_IRF(i,:) = cumsum(Chol_IRF(i,:));
end

% 2.1.3. Cholesky IRF Plotting (Not finished) -> NO NEED TO DO

% 2.1.4. Cholesky FEVD (Not finished) -> NO NEED TO DO

% 2.1.5. Cholesky FEVD Plotting (Not finished) -> NO NEED TO DO

%% 2.2 Sign Restriction (Agnostic)
%% 2.2.1. Agnostirc Sign Restriction Identification
% ** from Ambrogio Cesa-Bianchi's toolbox: https://github.com/ambropo/VAR-Toolbox
VARopt.vnames = {'prod','rea', 'rpo'};
VARopt.nsteps = 16;
VARopt.firstdate = 1; % Need to revise
VARopt.frequency = 'm';
VARopt.snames = {'s','ad','osd'};

Sign = [ -1, 1, 1; -1, 1, -1; 1, 1, 1];

% Setting the parameters
VARopt.ndraws = 3000;
VARopt.sr_hor = 1;
VARopt.pctg = 68; 


SRout = SR(VAR,Sign,VARopt);

%% 2.2.2. Sign Restriction (Agnostic) IRF and IRF Plots
figure;
subplot(3,1,1)
plot(squeeze(SRout.IRall(:,3,1,:))); hold on
plot(zeros(VARopt.nsteps),'--k','LineWidth',0.5); hold on
xlim([1 VARopt.nsteps]);
title('Response of Real Oil Price to \epsilon^{os}')

subplot(3,1,2)
plot(squeeze(SRout.IRall(:,3,2,:))); hold on
plot(zeros(VARopt.nsteps),'--k','LineWidth',0.5); hold on
xlim([1 VARopt.nsteps]);
title('Response of Real Oil Price to \epsilon^{ad}')

subplot(3,1,3)
plot(squeeze(SRout.IRall(:,3,3,:))); hold on
plot(zeros(VARopt.nsteps),'--k','LineWidth',0.5); hold on
xlim([1 VARopt.nsteps]);
title('Response of Real Oil Price to \epsilon^{osd}')
%% 2.3. Sign Restriction with Elasticity Bounds
%% 2.3.1. Sign Restirction (Elasticity-Bounded) Filtration

Bounded_SRout = struct();
Bounded_SRout.IRall = zeros(16,3,3,0);

k = 1;

for i = 1:size(SRout.Ball,3)
    Imp_Mat = SRout.Ball(:,:,i);
    a13 = Imp_Mat(1,3); a33 = Imp_Mat(3,3); a12 = Imp_Mat(1,2); a32 = Imp_Mat(3,2); a23 = Imp_Mat(2,3);
    if (a13/a33) < 0.0258 && (a12/a32) < 0.0258 && a23 > -1.5
        Bounded_SRout.IRall(:,:,:,k) = SRout.IRall(:,:,:,i);
        k = k+1;
        disp(i);
    end
end

%% 2.3.2. Sign Restriction (Elasticity-Bounded) IRF and IRF Plots
figure;
subplot(3,1,1)
plot(squeeze(Bounded_SRout.IRall(:,3,1,:))); hold on
plot(zeros(VARopt.nsteps),'--k','LineWidth',0.5); hold on
xlim([1 VARopt.nsteps]);
title('Response of Real Oil Price to \epsilon^{os}')

subplot(3,1,2)
plot(squeeze(Bounded_SRout.IRall(:,3,2,:))); hold on
plot(zeros(VARopt.nsteps),'--k','LineWidth',0.5); hold on
xlim([1 VARopt.nsteps]);
title('Response of Real Oil Price to \epsilon^{ad}')

subplot(3,1,3)
plot(squeeze(Bounded_SRout.IRall(:,3,3,:))); hold on
plot(zeros(VARopt.nsteps),'--k','LineWidth',0.5); hold on
xlim([1 VARopt.nsteps]);
title('Response of Real Oil Price to \epsilon^{osd}')

%% 2.4 Max-Share Identification - Maximize: Production Change / Supply Shock
addpath(genpath("BVAR_-master"));

h = 12; % Horizon chosen for FEVD
options.hor=h;
options.priors.name='Minnesota';
BVAR = bvar(Y, nlags, options);

%%
counter=0;
bvar_size = size(BVAR.Phi_draws,3); % Number of total Bayesian draws

figure;

i = 1;
j = 1;

counter=0;
bvar_size = size(BVAR.Phi_draws,3);
jump = 100;
irfs=nan(3,h,3,(bvar_size/jump));

for j=1:jump:size(BVAR.Phi_draws,3)
	Phi_draws=squeeze(BVAR.Phi_draws(:,:,j));
	Sigma_draws=squeeze(BVAR.Sigma_draws(:,:,j));
	Qbar=max_fevd(v,h,s,Phi_draws, Sigma_draws,1000);
	
    counter=counter+1;
    disp(j)
	
    irfs(:,:,:,counter)=iresponse(Phi_draws, Sigma_draws, h, Qbar);
end

options.varnames = {'Production Change', 'Real Economic Activities', 'Real Oil Price'};
options.saveas_dir = './Plots';
options.saveas_strng = 'Max - Production Change - Supply Shock';
options.shocksnames = {'Supply Shock','AD Shock','OSD Shock'};
options.conf_sig_2 = 0.95;

plot_all_irfs_(irfs,options);

%% 2.5 Max-Share Identification - Maximize: Real Oil Price / Oil-Specific Demand Shock

figure;

v = 3;
s = 3;

counter=0;
bvar_size = size(BVAR.Phi_draws,3);
jump = 100;
irfs=nan(3,h,3,(bvar_size/jump));

for j=1:jump:size(BVAR.Phi_draws,3)
	Phi_draws=squeeze(BVAR.Phi_draws(:,:,j));
	Sigma_draws=squeeze(BVAR.Sigma_draws(:,:,j));
	Qbar=max_fevd(v,h,s,Phi_draws, Sigma_draws,1000);
	
    counter=counter+1;
	disp(j)

    irfs(:,:,:,counter)=iresponse(Phi_draws, Sigma_draws, h, Qbar);
end

options.varnames = {'Production Change', 'Real Economic Activities', 'Real Oil Price'};
options.saveas_dir = './Plots';
options.saveas_strng = 'Max - Real Price of Oil - OSD Shock';
options.shocksnames = {'Supply Shock','AD Shock','OSD Shock'};
options.conf_sig_2 = 0.95;
plot_all_irfs_(irfs,options);