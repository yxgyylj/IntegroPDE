%% plot outcomes from the main file

close all;
set(0,'Defaultaxesfontsize',14)
%%%!!! run 'main.m' first !

% default setting
if ~exist('k1_fold') || ~exist('k2_fold')
    k1_fold = 1;
    k2_fold = 1;
end

% load data
if p.flag == 0  % no reverse seggregation (wash out)
    filename = sprintf('Data/noflux_L_%d_k1_%d_k2_%d.mat',p.L,k1_fold,k2_fold);
else
    filename = sprintf('Data/reverse/noflux_L_%d_k1_%d_k2_%d.mat',p.L,k1_fold,k2_fold);
end
load(filename);

step = floor(1/p.dt); 
ind = 1:step:step*p.tspan(end);
ind_sur = 1:10:step*p.tspan(end);

%% Plots
if p.flag == 0  % no reverse seggregation (wash out)
    prefix_img = sprintf('Imgs/noflux_L_%d_k1_%d_k2_%d',p.L, k1_fold,k2_fold);
else
    prefix_img = sprintf('Imgs/reverse/noflux_L_%d_k1_%d_k2_%d',p.L, k1_fold,k2_fold);
end
if ~exist(prefix_img)
    mkdir(prefix_img);
end

% contour plot
tmp = figure('position',[200,300,300,400]);
subplot(2,1,1)
contourf(p.x,p.tspan(ind),u(ind,:));
colormap(hot)
colorbar
ylabel('t in hours')
title('Microtubule Density')

subplot(2,1,2)
contourf(p.x,p.tspan(ind),v(ind,:));
colormap(hot)
colorbar
xlabel(['x in 200nm'])
ylabel('t in hours')
title('Neurofilament Density')
savefig(sprintf('%s/contour_plot.fig',prefix_img));
saveas(tmp,sprintf('%s/contour_plot_L_%d.jpg',prefix_img,p.L));

% surface plot
tmp = figure('position',[200,300,300,400]);
subplot(2,1,1);
waterfall(p.x,p.tspan(ind_sur),u(ind_sur,:));
shading interp;
xlabel('x')
ylabel('t')
zlabel('MT')
title('MT Density Evolution')

subplot(2,1,2);
waterfall(p.x,p.tspan(ind_sur),v(ind_sur,:));
shading interp;
xlabel('x')
ylabel('t in hours')
zlabel('NF')
title('NF Density Evolution')
savefig(sprintf('%s/surface_plot.fig',prefix_img));
saveas(tmp,sprintf('%s/surface_plot_L_%d.jpg',prefix_img,p.L));

% Total energy plot
tmp = figure('position',[200,300,300,400]);
subplot(2,1,1);
plot(p.tspan(2:end),energy_u(2:end),'linewidth',2);
ylim([-.1*max(energy_u(2:end)), 1.1*max(energy_u(2:end))]);
xlabel('t')
ylabel('MT energy')
xlim([0,p.T]);
title('\int_ (\nabla u(x,t))^2 dx')

subplot(2,1,2);
plot(p.tspan(2:end),energy_v(2:end),'linewidth',2);
ylim([-.1*max(energy_v(2:end)), 1.1*max(energy_v(2:end))]);
xlabel('t')
ylabel('NF energy')
xlim([0,p.T]);
title('\int_ (\nabla v(x,t))^2 dx')
savefig(sprintf('%s/energy_plot.fig',prefix_img));
saveas(tmp,sprintf('%s/energy_plot_L_%d.jpg',prefix_img,p.L));