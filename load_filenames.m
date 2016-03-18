function [filename_brain, filename_tumor, filename_ventricle, filename_initial_condition] = load_filenames(group,mouse,res,timept,timept_ic)
%% load_filenames

filename_brain = experiment_filename(group,mouse,res,timept,'brain');
filename_tumor = experiment_filename(group,mouse,res,timept,'tumor all');
filename_ventricle = experiment_filename(group,mouse,res,timept,'ventricles');

filename_initial_condition = experiment_filename(group,mouse,res,timept_ic,'tumor all');

end

