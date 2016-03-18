function err=gbm_minimizer(par,mouse)
%function [err,tumor_bg,savetumor_bg,brain,voxelnums] = gbm_minimizer(par,mouse)
%function [err,tumor_bg]=gbm_minimizer(par)

%% Model parameters
% Time parameters
start = 0.44;                          % Days 18-25 --> start = 0.72;
finish = 1.0;                       % Finish time should always be one, so that time is scaled between 0-1.
% p.tspan = [start,finish];             % tspan is used by ODE45.
p.tspan = [0.44,0.6,0.72,0.88,1]; % linspace(start,finish, 4);  % fixed steps are required for ensemble simulations
%p.tspan=[0.72,0.84,1];
p.t_end = p.tspan(end);               % Left over parameter for when finish was not equal to 1.

%% Tumor parameters
p.D = par(1);                            % Diffusion rate. 0.3; % 0.7;
p.D_sigma = 0;                      % Diffusion rate stddev 2
p.rho = par(2);                         % Growth rate. 24; % 50;
p.rho_sigma = 0;                    % Growth rate stddev. 4
p.CC = 1;                           % Carrying capacity is 0% - 100% (0-1 in code).
p.num_vox = 4;


if ((p.D<0) | (p.rho<0))
    err=4;
    return
end
%% Initial condition parameters
p.initial_condition_type = 1;
% Manual = 0
% Automatic = 1

% Manual
p.snr = 30; % Signal to noise ratio for logrand distribution for tumor cell distribution along needle lesion.
% for manual initial conditions select one of the following (corresponding
% to the geometries chosen below);
% G2_M1 = 1;
% G2_M2 = 2;
% G3_M1 = 3;
% G3_M2 = 4;
% G3_M3 = 5;
p.manual_initial_condition = 2+mouse;

% Automatic
p.timept_ic = 1;    % Choose which scan to use as initial condition (may require adjusting the timespan (tspan))
if(mouse==3)
    p.timept_ic=3;
    p.tspan = [0.72,0.88,1];
end
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
% 
% %% Ensemble iterations
 n = 1;
% 
%% Display all parameter choices
% disp('Review the following summary before proceeding:')
% disp('Computational geometry:')
% disp(['     Group: ', num2str(p.group)])
% disp(['     Mouse: ', num2str(p.mouse)])
% disp(['     Scan: ', num2str(p.timept)])
% disp(['     Z resolution: ', num2str(p.res), 'mm'])
switch p.initial_condition_type
    case 0 
        switch p.manual_initial_condition
            case 1
                disp('Manual initial condition: G2_M1')
            case 2
                disp('Manual initial condition: G2_M2')
            case 3
                disp('Manual initial condition: G3_M1')
            case 4
                disp('Manual initial condition: G3_M2')
            case 5
                disp('Manual initial condition: G3_M3')
        end
    case 1
%        disp(['Automatic initial condition (scan): ', num2str(p.timept_ic)])
end
% disp(['Time range: ', num2str(start), ' - ', num2str(finish)])
% disp(['Signal to noise ratio for lesion tail: ', num2str(p.snr)])
% disp('Growth and diffusion parameters')
% disp(['     Diffusion rate: ', num2str(p.D)])
% disp(['     Diffusion stddev: ', num2str(p.D_sigma)])
% disp(['     Growth rate: ', num2str(p.rho)])
% disp(['     Growth stddev: ', num2str(p.rho_sigma)])
% disp(['Number of iterations to be performed: ', num2str(n)])
% disp('Type "return" to continue')

% %keyboard
% 
% 
%% Ensemble
% Time has the fomat T(n, # time steps); and U has the format U(n, # time
% steps, x, y, z)

for i = 1:n
%    tic;
    [T(i,:), U(i,:,:,:,:)] = run_simulation(p);
%    disp(['Loop ',num2str(i),' completed']);
%    t(i) = toc;
end
% 
% %keyboard
% 
%% Calculate "average" tumor from ensemble
tic
ensemble_tumor = U(1,:,:,:,:);
if n > 1
    for i = 2:n
        ensemble_tumor = ensemble_tumor + U(i,:,:,:,:);
    end
end

ensemble_tumor = squeeze(ensemble_tumor./n);

toc
% 
% % %keyboard

%% Vizualize results

% p.timept=1;
% [filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
% disp('filenames'),filenames.tumor
%[~,~,m.brain1,m.tumor1,~,~,~,~] = mask(filenames);
%  figure;
%  imshow(initbrain(:,:,p.zslice)+inittumor(:,:,p.zslice),[])
%  title('Brain+Tumor')
% 
% figure;
% imshow(initbrain(:,:,p.zslice),[])
% title('Brain')
% 
% figure;
% imshow(inittumor(:,:,p.zslice),[])
% title('Tumor')
% 
% figure;
% imshow(inittumor(:,:,p.zslice))
% title('Unbracket')

%err=1;
%brain_bg=1;
% 
% 
if (mouse==3)
p.timept=5;
[filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
[~,~,m.brain,m.tumor,~,~,~,~] = mask(filenames);
p.timept=4;
[filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
[~,~,m.brain4,m.tumor4,~,~,~,~] = mask(filenames);

brain_bg = m.brain;
brain_bg = scaledata(brain_bg, 0, 1);
brain_bg(ensemble_tumor(end,:,:,:) > p.threshold) = 1;

brain_bg4 = m.brain4;
brain_bg4 = scaledata(brain_bg4, 0, 1);
brain_bg4(ensemble_tumor(end-1,:,:,:) > p.threshold) = 1;

 tumor_bg = squeeze(ensemble_tumor(end,:,:,:));
 tumor_bg(tumor_bg < p.threshold) = 0;
tumor_bg4 = squeeze(ensemble_tumor(end-1,:,:,:));
tumor_bg4(tumor_bg4 < p.threshold) = 0;
%tumor_bg1 = squeeze(ensemble_tumor(end-4,:,:,:));
%tumor_bg1(tumor_bg1 < p.threshold) = 0;
% 
 plot_tumor = brain_bg + tumor_bg;
 plot_tumor4=brain_bg4 + tumor_bg4;
% plot_tumor1=brain_bg1 + tumor_bg1
% 
 m.tumor(m.tumor>0)=1;
 m.tumor4(m.tumor4>0)=1;
%  m.tumor1(m.tumor1>0)=1;

 tumor_bg(tumor_bg>0)=1;
 tumor_bg4(tumor_bg4>0)=1;
%  tumor_bg1(tumor_bg1>0)=1;
 counts=find(m.tumor==1 & tumor_bg==1);
 intersection=length(counts);
 counts=find(m.tumor==1 | tumor_bg==1);
 union=length(counts);
 counts=find(m.tumor4==1 & tumor_bg4==1);
 intersection4=length(counts);
 counts=find(m.tumor4==1 | tumor_bg4==1);
 union4=length(counts);
 
   err1=1-intersection/union
  err4=1-intersection4/union4
%  err5=1-intersection1/union1
  err=err1+err4%+err5
  par
 brain=m.brain;
 brain4=m.brain4;
 figuremaker(tumor_bg,brain);
 figuremaker(tumor_bg4, brain4);
    vistumor=find(tumor_bg4==1);
    voxelnums(1)=length(vistumor);
    vistumor=find(tumor_bg==1);
    voxelnums(2)=length(vistumor);
savetumor_bg=tumor_bg;
else
    p.timept=5;
    [filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
    [~,~,m.brain,m.tumor,~,~,~,~] = mask(filenames);
    p.timept=4;
[filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
[~,~,m.brain4,m.tumor4,~,~,~,~] = mask(filenames);
p.timept=2;
[filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
[~,~,m.brain2,m.tumor2,~,~,~,~] = mask(filenames);
p.timept=3;
[filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);
[~,~,m.brain3,m.tumor3,~,~,~,~] = mask(filenames);
brain_bg = m.brain;
brain_bg = scaledata(brain_bg, 0, 1);

brain_bg(ensemble_tumor(end,:,:,:) > p.threshold) = 1;

brain_bg2 = m.brain2;
brain_bg2 = scaledata(brain_bg2, 0, 1);

brain_bg2(ensemble_tumor(end-3,:,:,:) > p.threshold) = 1;

brain_bg3 = m.brain3;
brain_bg3 = scaledata(brain_bg3, 0, 1);

brain_bg3(ensemble_tumor(end-2,:,:,:) > p.threshold) = 1;

brain_bg4 = m.brain4;
brain_bg4 = scaledata(brain_bg4, 0, 1);

brain_bg4(ensemble_tumor(end-1,:,:,:) > p.threshold) = 1;
 tumor_bg = squeeze(ensemble_tumor(end,:,:,:));
 tumor_bg(tumor_bg < p.threshold) = 0;
tumor_bg2 = squeeze(ensemble_tumor(end-3,:,:,:));
tumor_bg2(tumor_bg2 < p.threshold) = 0;
tumor_bg3 = squeeze(ensemble_tumor(end-2,:,:,:));
tumor_bg3(tumor_bg3 < p.threshold) = 0;
tumor_bg4 = squeeze(ensemble_tumor(end-1,:,:,:));
tumor_bg4(tumor_bg4 < p.threshold) = 0;
%tumor_bg1 = squeeze(ensemble_tumor(end-4,:,:,:));
%tumor_bg1(tumor_bg1 < p.threshold) = 0;
% 
 plot_tumor = brain_bg + tumor_bg;
plot_tumor2=brain_bg2 + tumor_bg2;
plot_tumor3=brain_bg3 + tumor_bg3;
 plot_tumor4=brain_bg4 + tumor_bg4;
% plot_tumor1=brain_bg1 + tumor_bg1
% 
 m.tumor(m.tumor>0)=1;
 savetumor_bg=tumor_bg;
m.tumor2(m.tumor2>0)=1;
m.tumor3(m.tumor3>0)=1;
 m.tumor4(m.tumor4>0)=1;
%  m.tumor1(m.tumor1>0)=1;

 tumor_bg(tumor_bg>0)=1;
tumor_bg2(tumor_bg2>0)=1;
tumor_bg3(tumor_bg3>0)=1;
 tumor_bg4(tumor_bg4>0)=1;
%  tumor_bg1(tumor_bg1>0)=1;
 counts=find(m.tumor==1 & tumor_bg==1);
 intersection=length(counts);
 counts=find(m.tumor==1 | tumor_bg==1);
 union=length(counts);
counts=find(m.tumor2==1 & tumor_bg2==1);
intersection2=length(counts);
counts=find(m.tumor2==1 | tumor_bg2==1);
union2=length(counts);
counts=find(m.tumor3==1 & tumor_bg3==1);
intersection3=length(counts);
counts=find(m.tumor3==1 | tumor_bg3==1);
union3=length(counts);
 counts=find(m.tumor4==1 & tumor_bg4==1);
 intersection4=length(counts);
 counts=find(m.tumor4==1 | tumor_bg4==1);
 union4=length(counts);
   err1=1-intersection/union;
  err2=1-intersection2/union2;
  err3=1-intersection3/union3;
  err4=1-intersection4/union4;
%  err5=1-intersection1/union1
  err=err1+err2+err3+err4%+err5
  par
%  brain=m.brain;
%  brain4=m.brain4;
%  brain3=m.brain3;
%  brain2=m.brain2;
%  figuremaker(tumor_bg,brain);
%  figuremaker(tumor_bg4, brain4);
%  figuremaker(tumor_bg3, brain3);
%  figuremaker(tumor_bg2, brain2);
%      vistumor=find(tumor_bg2==1);
%     voxelnums(1)=length(vistumor);
%     vistumor=find(tumor_bg3==1);
%     voxelnums(2)=length(vistumor);
%     vistumor=find(tumor_bg4==1);
%     voxelnums(3)=length(vistumor);
%     vistumor=find(tumor_bg==1);
%     voxelnums(4)=length(vistumor);
end
% 

% 

%brain_bg1 = m.brain1;
%brain_bg1 = scaledata(brain_bg1, 0, 1);

%brain_bg1(ensemble_tumor(end-4,:,:,:) > p.threshold) = 1;

 
% 
% % figure;
% % imshow(abs(m.tumor(:,:,18)-tumor_bg(:,:,18)))
% % figure;
% % imshow(abs(m.tumor2(:,:,18)-tumor_bg2(:,:,18)))
% % figure;
% % imshow(abs(m.tumor3(:,:,18)-tumor_bg3(:,:,18)))
% % figure;
% % imshow(abs(m.tumor4(:,:,18)-tumor_bg4(:,:,18)))
% 
% figure
% imagesc(squeeze(plot_tumor(:,:,p.zslice)));
% caxis([0 2.03])
% colormap(vertcat(gray(250),jet(256)));
% title('Time point 5')
% % 
% figure
% imagesc(squeeze(plot_tumor2(:,:,p.zslice)));
% caxis([0 2.03])
% colormap(vertcat(gray(250),jet(256)));
% title('Time point 2')
% 
% figure
% imagesc(squeeze(plot_tumor3(:,:,p.zslice)));
% caxis([0 2.03])
% colormap(vertcat(gray(250),jet(256)));
% title('Time point 3') 
% % 
% figure
% imagesc(squeeze(plot_tumor4(:,:,p.zslice)));
% caxis([0 2.03])
% colormap(vertcat(gray(250),jet(256)));
% title('Time point 4') 
% % 
% % figure;
% % imshow((m.tumor(:,:,18)));
% % hold on
% % plot(tumor_bg(:,:,18));
% % figure;
% % imshow(abs(m.tumor2(:,:,18)-tumor_bg2(:,:,18)))
% % figure;
% % imshow(abs(m.tumor3(:,:,18)-tumor_bg3(:,:,18)))
% % figure;
% % imshow(abs(m.tumor4(:,:,18)-tumor_bg4(:,:,18)))
% 
% %% Error Function
% 

%  counts=find(m.tumor1==1 & tumor_bg1==1);
% intersection1=length(counts);
% counts=find(m.tumor1==1 | tumor_bg1==1);
% union1=length(counts);
% 
% err=4-intersection/union-intersection2/union2-intersection3/union3-intersection4/union4
end
