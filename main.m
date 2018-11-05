%========== Simulation for a nonlocal advection-diffusion model  ==========
%---------- used for NOFLUX boundary condition of MT-NF phase segregation
%
%   Step 1: load parameter values (saved in "Data" folder)
%   Step 2: prelocate some terms for calculating integral terms
%   Step 3: for each fixed value of rho got in Step 1, solve corresponding
%       Davidenko differential equation using predictor-corrector method;
%   Step 4: save data and plot the outcome.
%
%   by Xige Yang
%==========================================================================

clear; close all; clc;
addpath('src');

% run params.m first if you haven't yet!
load Data/params.mat
p.r1 = p.r2*1.5; 
alpha = 1/4/10/5;
%p.eta1 = 0;
%p.eta2 = 0;

%%%--- change p.flag to 1 if you want to see reversable
input_msg = ['Do you want to see washout phenomena?\n',...
                '1 -- Yes\n'...
                '0 -- No\n'];
p.flag = input(input_msg);  

k1_fold = 0;
k2_fold = -10;
p.k1 = p.k1*k1_fold; 
p.k2 = p.k2*k2_fold;

% define domain and spatial mesh
p.xmin = 0;  % Initial Length
p.xmax = 5; % Final Length L=1 micron
p.L = p.xmax-p.xmin;
p.dx = 0.02; % Step in x
p.x = [p.xmin:p.dx:p.xmax]'; 
p.Nx = length(p.x);

% initial conditions
init = getInit(p.x,p);

% solve semi-discretized system 
p.T = 50;   % computational time
p.dt = 0.02; % step size for integration
p.tspan = [0:p.dt:p.T];

%% Pre-difine integral indexes
%cal_integral_coeff(p);
p.ind_range1 = floor(p.r1/p.dx);
p.ind_range2 = floor(p.r2/p.dx);

fprintf('%% Simulating for k1 = %.3f, k2 = %.3f...\n', p.k1, p.k2);
tic
[t, n]=ode15s(@(t,n)noflux_rhs(t,n,p),p.tspan,init);
toc

u = n(:,1:p.Nx); 
v = n(:,p.Nx+1:end);

up = u(:,2:end);      um = u(:,1:end-1);
vp = v(:,2:end);      vm = v(:,1:end-1);
energy_u = sum(((up-um)/p.dx).*(up-um),2);
energy_v = sum(((vp-vm)/p.dx).*(vp-vm),2);

%% Save data and draw graph
if p.flag == 0  % no reverse seggregation (wash out)
    filename = sprintf('Data/noflux_L_%d_k1_%d_k2_%d.mat',p.L,k1_fold,k2_fold);
else
    if ~exist('Data/reverse')
        mkdir('Data/reverse')
    end
    filename = sprintf('Data/reverse/noflux_L_%d_k1_%d_k2_%d.mat',p.L,k1_fold,k2_fold);
end
save(filename,'u','v','k1_fold','k2_fold','p','energy_u','energy_v');

img_noflux;
