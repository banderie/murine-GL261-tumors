function [T,U] = run_simulation(p)
%% Run simulation
% The 'p' object contains all parameters (single numbers, no vectors,
% matricies or tensors). The 'm' object contains all vectors, matricies
% and tensors. The 'filenames' object contains all file paths

%% Load filenames
[filenames.brain,filenames.tumor,filenames.ventricle,filenames.initial_conditions] = load_filenames(p.group,p.mouse,p.res,p.timept,p.timept_ic);

%% Generate mask and load geometries (i.e. brain, tumor, ventricle and boundary conditions)
[m.mask,m.boundary,m.brain,m.tumor,m.ventricle,p.x,p.y,p.z] = mask(filenames);

%% Load initial conditions
switch p.initial_condition_type
    case 0
        [p.xslice, p.yslice, p.zslice, u0] = manual_initial_conditions(p);
        u0(u0>0)=0.5;
    case 1
        u0 = load_initial_conditions(filenames,p);
        u0(u0>0)=0.5;
    case 2
        title=sprintf('S1G%dM%d_timepoint%d.mat',p.group,p.mouse,p.timept_ic);
        load(title)
        u0 = squeeze(tumorsy(end,:,:,:));
end

%% Generate Stochastic Parameters
%Normal (anything < 0 set to 0.1)
m.stochastic_D = normrnd(p.D,p.D_sigma,p.x,p.y,p.z);
m.stochastic_rho = normrnd(p.rho,p.rho_sigma,p.x,p.y,p.z);
m.stochastic_D(m.stochastic_D < 0) = 0.1;
m.stochastic_rho(m.stochastic_rho < 0) = 0.1;

%% Run Model
OPTIONS = odeset('RelTol', 10^-6, 'AbsTol', 10^-9);
[T, U] = ode45(@(t,u) discretization(t, u, p, m), p.tspan, u0, OPTIONS);
U = reshape(U,length(T),p.x,p.y,p.z);

end