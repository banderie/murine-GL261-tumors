function [err,voxelnums]=gbm_minimizerindypts(par, timepoint, mouse)

%function [err,tumor_bg,savetumor_bg,brain,voxelnums]=gbm_minimizerindypts(par, timepoint, mouse)

%% Model parameters
p.initial_condition_type=1;
p.manual_initial_condition=2+mouse;
p.timept_ic=timepoint;
if (timepoint==1 || timepoint==3)
    p.tspan=[0.84,1];
else
    p.tspan=[0.88,1];
end
%% Tumor parameters
p.D = par(1);                            % Diffusion rate. 0.3; % 0.7;
p.D_sigma = 0;                      % Diffusion rate stddev 2
p.rho = par(2);                         % Growth rate. 24; % 50;
p.rho_sigma = 0;                    % Growth rate stddev. 4
p.CC = 1;                           % Carrying capacity is 0% - 100% (0-1 in code).
p.num_vox = 4;

%% Vizualization parameters
p.threshold = 0.16;                 % percent carrying capacity achieved (0 - 1) (do not use 0 or 1)
p.pause_time = 0;                   % pause time between frames    
p.play_type = 'full';               % last = only final tumor volume. full = full tumor movie.
p.fontsize = 25;                    % fontzise on graphs
p.zslice = 19; 
p.yslice = 39; 
p.xslice = 60; 

%% Choose computational geometry files
% G2: M1, M2
% G3: M1, M2, M3

p.tp(1) = 11.0/25.0;
p.tp(2) = 15.0/25.0;
p.tp(3) = 18.0/25.0;
p.tp(4) = 22.0/25.0;
p.tp(5) = 25.0/25.0;

p.group = 3;
p.mouse = mouse;
p.res = 0.5;        % Only use 0.5 as the higher resolutions are interpolated
p.timept = 5;       % 1,2,3,4,5 Choose which timepoint to use as simulation geometry

%% Ensemble iterations
n = 1;

%% Ensemble
% Time has the fomat T(n, # time steps); and U has the format U(n, # time
% steps, x, y, z)

for i = 1:n
%    tic;
    [T(i,:), U(i,:,:,:,:)] = run_simulation(p);
%    disp(['Loop ',num2str(i),' completed']);
%    t(i) = toc;
end

%keyboard

%% Calculate "average" tumor from ensemble
ensemble_tumor = squeeze(U(1,:,:,:,:));

% %keyboard

%% Vizualize results
p.timept=timepoint+1;
[filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
[~,~,m.brain,m.tumor,~,~,~,~] = mask(filenames);

brain_bg = m.brain;
brain_bg = scaledata(brain_bg, 0, 1);

brain_bg(ensemble_tumor(end,:,:,:) > p.threshold) = 1;

tumor_bg = squeeze(ensemble_tumor(end,:,:,:));

tumor_bg(tumor_bg < p.threshold) = 0;

plot_tumor = brain_bg + tumor_bg;

m.tumor(m.tumor>0)=1;
savetumor_bg=tumor_bg;
tumor_bg(tumor_bg>0)=1;
brain=m.brain;
figuremaker(tumor_bg,brain);
% figure
% %colormap(map)
% % surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[13 13; 13 13], 'CData',m.brain(:,:,13),'FaceColor','texturemap')
% % colormap(gray)
% % hold on
% surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[14 14; 14 14], 'CData',m.brain(:,:,14),'FaceColor','texturemap')
% colormap(gray)
% hold on
% surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[15 15; 15 15], 'CData',m.brain(:,:,15),'FaceColor','texturemap')
% colormap(gray)
% hold on
% surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[16 16; 16 16], 'CData',m.brain(:,:,16),'FaceColor','texturemap')
% colormap(gray)
% hold on
% surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[17 17; 17 17], 'CData',m.brain(:,:,17),'FaceColor','texturemap')
% colormap(gray)
% surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[17 17; 17 17], 'CData',m.brain(:,:,17),'FaceColor','texturemap')
% colormap(jet)
% surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[18 18; 18 18], 'CData',m.brain(:,:,18),'FaceColor','texturemap')
% colormap(gray)
% hold on
% surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[19 19; 19 19], 'CData',m.brain(:,:,19),'FaceColor','texturemap')
% colormap(gray)
% hold on
% surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[20 20; 20 20], 'CData',m.brain(:,:,20),'FaceColor','texturemap')
% colormap(gray)
% % hold on
% % surface('Xdata',[0 1;0 1],'Ydata',[0 0; 1 1],'Zdata',[21 21; 21 21], 'CData',m.brain(:,:,21),'FaceColor','texturemap')
% % colormap(gray)
% %colormap(vertcat(gray(250),jet(256)))
% %surface(plot_tumor(:,:,13))
% brain=m.brain;
% view(3)
% 
% figure
% imagesc(squeeze(plot_tumor(:,:,p.zslice)));
% caxis([0 2.03])
% colormap(vertcat(gray(250),jet(256)));
% figure
% surf(tumor_bg(:,:,14:20))

%% Error Function

counts=find(m.tumor==1 & tumor_bg==1);
intersection=length(counts);
counts=find(m.tumor==1 | tumor_bg==1);
union=length(counts);
vistumor=find(tumor_bg==1);
voxelnums=length(vistumor)
err=1-intersection/union;

par
