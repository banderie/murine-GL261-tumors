function [err,voxelnums]=gbm_minimizereachtimepts(par,timepoint,mouse)
if (par(1)<0 || par(2)<0)
    err=4;
    return
end
%% Model parameters
if timepoint==1
    p.initial_condition_type = 1;
    p.manual_initial_condition = 2+mouse;
    p.timept_ic = 1;    
    p.tspan=[0.84,1];
    endtp=2;
elseif timepoint==2
    p.initial_condition_type = 2;
    p.tspan=[0.88,1];
    endtp=3;
    p.timept_ic = 2;    

elseif timepoint==3
    p.initial_condition_type = 2;
    p.tspan=[0.84,1];
    endtp=4;
    p.timept_ic = 3;    
    if (mouse==3)
        p.initial_condition_type=1;
        p.manual_initial_condition=5;
    end
else 
    p.initial_condition_type = 2;
    p.tspan=[0.88,1];
    endtp=5;
    p.timept_ic = 4;    

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
p.zslice = 18; 
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
tumorsy=squeeze(ensemble_tumor);
max(max(max(tumorsy(end,:,:,:))));
title=sprintf('S1G%dM%d_timepoint%d.mat',p.group,p.mouse,endtp);
save(title,'tumorsy')

% %keyboard

%% Vizualize results

    p.timept=endtp;
    [filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
    [~,~,m.brain,m.tumor,~,~,~,~] = mask(filenames);

brain_bg = m.brain;
brain_bg = scaledata(brain_bg, 0, 1);

brain_bg(ensemble_tumor(end,:,:,:) > p.threshold) = 1;

tumor_bg = squeeze(ensemble_tumor(end,:,:,:));
plot_tumorfirst=tumor_bg+brain_bg;

tumor_bg(tumor_bg < p.threshold) = 0;

plot_tumor = brain_bg + tumor_bg;

m.tumor(m.tumor>0)=1;
tumor_bg(tumor_bg>0)=1;
% brain=m.brain;
% figuremaker(tumor_bg,brain);
% figure
% imagesc(squeeze(plot_tumorfirst(:,:,p.zslice)));
% caxis([0 2.03])
% colormap(vertcat(gray(250),jet(256)));
% 
% figure
% imagesc(squeeze(plot_tumor(:,:,p.zslice)));
% caxis([0 2.03])
% colormap(vertcat(gray(250),jet(256)));
% 
% figure
% imagesc(squeeze(brain_bg(:,:,p.zslice)));
% caxis([0 2.03])
% colormap(vertcat(gray(250),jet(256)));

%% Error Function

counts=find(m.tumor==1 & tumor_bg==1);
intersection=length(counts);
counts=find(m.tumor==1 | tumor_bg==1);
union=length(counts);

err=1-intersection/union;
vistumor=find(tumor_bg==1);
voxelnums=length(vistumor)
par;
