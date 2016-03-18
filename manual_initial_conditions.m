function [xslice, yslice, zslice, u0] = manual_initial_conditions(p)
%% load_initial_conditions
% The needle used was 26 gauge, or 0.4636mm. The needle was withdrawn 0.4mm
% before cells were injected. Hence, we will have the initial bolus take up
% a 0.4mm x 0.4mm x 0.4mm space. This corresponds to approximately one
% voxel in the Z axis, and four voxels in each the X and Y axes. 
% To ensure a consistent number of cells (or consistant total tumor cell
% density) is applied for each geometry, the initial condition densities
% are normalized to a chosen value (currently the full density of four
% voxels) before being returned.

switch p.manual_initial_condition
    case 1
        %% S1_G2_M1 initial conditions
        % Prescan
        zslice = 20; 
        yslice = 39; 
        xslice = 60; 
        u0 = zeros(p.x,p.y,p.z);
        u0(xslice-2:xslice+1,yslice-2:yslice+1,zslice) = 1; % make cell deposite less wide and
        R = logrand_IC(p.snr,80-60 + 1, p.fontsize);
        u0(60:80,yslice,zslice) = R; 
        % u0(82,yslice,zslice) = 0.5; % deposite a superficial cell clump 
        u0 = reshape(u0, numel(u0),1); 
        disp('S1_G2_M1 initial condition generated');
    case 2
        %% S1_G2_M2 initial conditions
        % Prescan
        zslice = 20;
        yslice = 39;
        xslice = 63;
        u0 = zeros(p.x,p.y,p.z);
        u0(xslice-2:xslice+1,yslice-2:yslice+1,zslice) = 1; % make cell deposite less wide and longer.
        R = logrand_IC(p.snr,80-xslice + 1, p.fontsize);
        u0(xslice:80,yslice,zslice) = R; 
        %u0(82,yslice,zslice) = 0.5; % deposite a superficial cell clump
        u0 = reshape(u0, numel(u0),1);
        disp('S1_G2_M2 initial condition generated');
    case 3
        %% S1_G3_M1 initial conditions
        % S1_G3_M1_Prescan
        zslice = 20; % ESTIMATE! NEED COMBINED TUMOR DATA TO KNOW FOR SURE.
        yslice = 33;
        xslice = 38;
        u0 = zeros(p.x,p.y,p.z);
        u0(xslice-2:xslice+1,yslice:yslice+1,zslice) = 1; % make cell deposite less wide and longer.
        R = logrand_IC(p.snr,65-xslice + 1, p.fontsize);
        u0(xslice:65,yslice,zslice) = R; 
        %u0(65,yslice,zslice) = 0.5; % deposite a superficial cell clump
        u0 = reshape(u0, numel(u0),1);
        disp('S1_G3_M1 initial condition generated');
    case 4
        %% S1_G3_M2 initial conditions
        % S1_G3_M2_Prescan
        xslice = 60;
        yslice = 33;
        zslice = 18;
        u0 = zeros(p.x,p.y,p.z);
        u0(xslice-2:xslice+1,yslice-2:yslice+1,zslice) = 1; % make cell deposite less wide and longer.
        R = logrand_IC(p.snr,80-60 + 1, p.fontsize);
        u0(xslice:80,yslice,zslice) = R./1000; 
        %u0(82,yslice,zslice) = 0.5; % deposite a superficial cell clump
        u0 = reshape(u0, numel(u0),1);
        fprintf('S1_G3_M2 initial condition generated');
    case 5
        %% S1_G3_M3 initial conditions
        % Prescan
        zslice = 20;
        yslice = 33;
        xslice = 60;
        u0 = zeros(p.x,p.y,p.z);
        u0(xslice-2:xslice+1,yslice-2:yslice+1,zslice) = 1; % make cell deposite less wide and longer.
        R = logrand_IC(p.snr,80-xslice + 1, p.fontsize);
        u0(xslice:80,yslice,zslice) = R; 
        %u0(82,yslice,zslice) = 0.5; % deposite a superficial cell clump
        u0 = reshape(u0, numel(u0),1);
        disp('S1_G3_M3 initial condition generated');
end

%% OLD STUFF
%{
%% S1_G3_M2_Scan_1
%{
xslice = 60;
yslice = 33;
zslice = 18;
q.zslice = zslice;
q.filename_brain     = '/Users/banderies/Google Drive/School/ASU/Research/Mouse Data/Preprocessed JDP_S1_G3_M2/3_T2w Registration/00/001 - 20140407 195722 3_T2w Registration 00/Preprocessed MATLAB/003 01/Interp2ImageGrid/Grid from 003 01 - 3 - 20140414 205836 3_T2w Registration 00/Interpolated XYZ Resolution 0.100000 0.100000 0.500000 - Transformed Intensitp.y Brain_grap.yvalues.mat';
q.filename_tumor     = '/Users/banderies/Google Drive/School/ASU/Research/Mouse Data/Preprocessed JDP_S1_G3_M2/3_T2w Registration/00/001 - 20140407 195722 3_T2w Registration 00/Preprocessed MATLAB/003 01/Interp2ImageGrid/Grid from 003 01 - 3 - 20140414 205836 3_T2w Registration 00/Interpolated XYZ Resolution 0.100000 0.100000 0.500000 - Transformed Tumor_grap.yvalues.mat';
q.filename_ventricle = '/Users/banderies/Google Drive/School/ASU/Research/Mouse Data/Preprocessed JDP_S1_G3_M2/3_T2w Registration/00/001 - 20140407 195722 3_T2w Registration 00/Preprocessed MATLAB/003 01/Interp2ImageGrid/Grid from 003 01 - 3 - 20140414 205836 3_T2w Registration 00/Interpolated XYZ Resolution 0.100000 0.100000 0.500000 - Transformed Ventricle_grap.yvalues.mat';
u0 = initial_conditions(q);
u0 = reshape(u0, numel(u0),1);    
%}

%% Ep.xperimental initial conditions
% For testing for leaks
%{
zslice = 3;
yslice = 33;
xslice = 60;
u0 = zeros(p.x,p.y,p.z);
u0(45:50,40,zslice) = 1;
u0 = reshape(u0, numel(u0),1);
fprintf('Initial condition generated \n\n');
%}
%}

%% Normalize for consistant number of cells
u0 = p.num_vox*u0./(sum(u0));




end

