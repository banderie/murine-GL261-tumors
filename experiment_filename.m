function filestring = experiment_filename(group,mouse,res,timept,tissue)
% filestring = experiment_filename(group,mouse,res,timept,tissue)
%
% This function outputs the file name of the experimental data requested
% based on the group and mouse number, resolution of data, time point, and
% tissue type.
%
% Inputs:
%   group  = group number (options: 2, 3)
%   mouse  = mouse number (options: 1, 2, 3)
%   res    = resolution (options: 0.1, 0.3, 0.5)
%   timept = time point (options: 1, 2, 3, 4, 5)
%   tissue = string, what type of brain tissue (options: 'brain',
%            'tumor inferior', 'tumor superior', 'tumor all','ventricles')
%
% Output:
%   filestring = .mat file string

s1 = '/Users/rutter/Documents/UniformMouseData/Animal Data/';  %%% only line that would need to be changed
                              %%% depending on where UniformMouseData
                              %%% folder is stored

s2 = 'UniformMouseData/Animal Data/';

%%% G2
if group==2
    s3 = 'Animal S1 G2 M12/Preprocessed Data - 1/';
    s8 = 'Grid from 003 01 - 003 - 20140310 161941 3_T2w Registration 00/';
    
    %%% M1
    if mouse==1
        s4 = 'Preprocessed JDP_S1_G2_M1/';
        if timept==1
            s6 = '001 - 20140303 180105 3_T2w Registration 00/';
        elseif timept==2
            s6 = '002 - 20140307 202945 3_T2w Registration 00/';
        elseif timept==3
            s6 = '003 - 20140310 161941 3_T2w Registration 00/';
        elseif timept==4
            s6 = '004 - 20140314 194204 3_T2w Registration 00/';
        elseif timept==5
            s6 = '005 - 20140317 151941 3_T2w Registration 00/';
        end
    
    %%% M2
    elseif mouse==2
        s4 = 'Preprocessed JDP_S1_G2_M2/';
        if timept==1
            s6 = '001 - 20140303 203408 3_T2w Registration 00/';
        elseif timept==2
            s6 = '002 - 20140307 225053 3_T2w Registration 00/';
        elseif timept==3
            s6 = '003 - 20140310 185841 3_T2w Registration 00/';
        elseif timept==4
            s6 = '004 - 20140314 223517 3_T2w Registration 00/';
        elseif timept==5
            s6 = '005 - 20140317 172757 3_T2w Registration 00/';
        end
    end
    
    
%%% G3
elseif group==3
    s3 = 'Animal S1 G3 M123/Preprocessed Data - 1/';
    
    %%% M1
    if mouse==1
        s4 = 'Preprocessed JDP_S1_G3_M1/';
        s8 = 'Grid from 003 01 - 003 - 20140414 181037 3_T2w Registration 00/';
        if timept==1
            s6 = '001 - 20140407 165342 3_T2w Registration 00/';
        elseif timept==2
            s6 = '002 - 20140411 153622 3_T2w Registration 00/';
        elseif timept==3
            s6 = '003 - 20140414 181037 3_T2w Registration 00/';
        elseif timept==4
            s6 = '004 - 20140418 194316 3_T2w Registration 00/';
        elseif timept==5
            s6 = '005 - 20140421 165214 3_T2w Registration 00/';
        end
        
    %%% M2
    elseif mouse==2
        s4 = 'Preprocessed JDP_S1_G3_M2/';
        s8 = 'Grid from 003 01 - 3 - 20140414 205836 3_T2w Registration 00/';
        if timept==1
            s6 = '001 - 20140407 195722 3_T2w Registration 00/';
        elseif timept==2
            s6 = '002 - 20140411 174006 3_T2w Registration 00/';
        elseif timept==3
            s6 = '003 - 20140414 205836 3_T2w Registration 00/';
        elseif timept==4
            s6 = '004 - 20140418 221213 3_T2w Registration 00/';
        elseif timept==5
            s6 = '005 - 20140421 190226 3_T2w Registration 00/';
        end
        
    %%% M3
    elseif mouse==3
        s4 = 'Preprocessed JDP_S1_G3_M3/';
        s8 = 'Grid from 003 01 - 003 - 20140414 235724 3_T2w Registration 00/';
        if timept==1
            s6 = '001 - 20140407 221252 3_T2w Registration 00/';
        elseif timept==2
            s6 = '002 - 20140411 194109 3_T2w Registration 00/';
        elseif timept==3
            s6 = '003 - 20140414 235724 3_T2w Registration 00/';
        elseif timept==4
            s6 = '004 - 20140419 002843 3_T2w Registration 00/';
        elseif timept==5
            s6 = '005 - 20140421 211515 3_T2w Registration 00/';
        end
    end
end

s5 = '3_T2w Registration/00/'; % NOTE: S1G3M2 has other folders besides 00:
                               % 01,02,03,04,05,06    (but this is the
                               % only one with multiple folders)

s7 = 'Preprocessed MATLAB/003 01/Interp2ImageGrid/';

if res==0.1
    s9temp = 'Interpolated XYZ Resolution 0.100000 0.100000 0.100000 - Transformed';
elseif res==0.3
    s9temp = 'Interpolated XYZ Resolution 0.100000 0.100000 0.300000 - Transformed';
elseif res==0.5
    s9temp = 'Interpolated XYZ Resolution 0.100000 0.100000 0.500000 - Transformed';
end

if strcmp(tissue,'brain')==1
    s9 = strcat(s9temp,' Intensity Brain_grayvalues.mat');
elseif strcmp(tissue,'tumor inferior')==1
    if group==3 && mouse==1 && timept==3
        s9 = strcat(s9temp,' Tumor inferior_grayvalues.mat');
    else
        s9 = strcat(s9temp,' Tumor Inferior_grayvalues.mat');
    end
    if group==2 && mouse==2
        disp('ERROR: There are no Tumor Inferior files for G2 M2');
    elseif group==3 && mouse==2
        disp('ERROR: There are no Tumor Inferior files for G3 M2');
    end
elseif strcmp(tissue,'tumor superior')==1
    s9 = strcat(s9temp,' Tumor Superior_grayvalues.mat');
    if group==2 && mouse==1
        disp('ERROR: There are no Tumor Superior files for G2 M1');
    elseif group==2 && mouse==2
        disp('ERROR: There are no Tumor Superior files for G2 M2');
    elseif group==3 && mouse==2
        disp('ERROR: There are no Tumor Superior files for G3 M2');
    end
elseif strcmp(tissue,'tumor all')==1
    if (group==2 && mouse==2) || ...
            (group==3 && mouse==1 && (timept==4 || timept==5)) || ...
            (group==3 && mouse==2) || (group==3 && mouse==3)
        s9 = strcat(s9temp,' Tumor_grayvalues.mat');
    else
        s9 = strcat(s9temp,' Tumors_grayvalues.mat');
    end
elseif strcmp(tissue,'ventricles')==1
    if (group==3 && mouse==2)
        s9 = strcat(s9temp,' Ventricle_grayvalues.mat');
    else
        s9 = strcat(s9temp,' Ventricles_grayvalues.mat');
    end
end

% filestring = strcat(s1,s2,s3,s4,s5,s6,s7,s8,s9);
filestring = strcat(s1,s3,s4,s5,s6,s7,s8,s9);