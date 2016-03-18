function uprime = discretization(t, u, p, m)

uprimexf = zeros(p.x-1,p.y,p.z); % Preallocate 6 matrices to store the forward 
uprimexb = zeros(p.x,p.y,p.z); % and backward differences in each dimension.
uprimeyf = zeros(p.x,p.y-1,p.z);
uprimeyb = zeros(p.x,p.y,p.z);
uprimezf = zeros(p.x,p.y,p.z-1);
uprimezb = zeros(p.x,p.y,p.z);

uprimex = zeros(p.x,p.y,p.z); % Preallocate 3 matrices to store the concatenation
uprimey = zeros(p.x,p.y,p.z); % of the finite differences in each dimension.
uprimez = zeros(p.x,p.y,p.z);

u = reshape(u, p.x, p.y, p.z);

deltaX = 1; % voxels are in mm, but using 0.1 makes problem stiff, so use 1 and 5 (10ths of mm). Need to adjust 
deltaY = 1; % other parameters to be consistent.
deltaZ = p.res*10;

%% Calculate forward and backward differences
for k = 1:p.z
    for j = 1:p.y
        uprimexf(:,j,k) = u(2:p.x,j,k) - u(1:p.x-1,j,k);
        uprimexb(2:p.x,j,k) = u(1:p.x-1,j,k) - u(2:p.x,j,k);
        if j ~= p.y
            uprimeyf(:,j,k) = u(:,j+1,k) - u(:,j,k);
        end
        if j ~= 1
            uprimeyb(:,j,k) = u(:,j-1,k) - u(:,j,k);
        end
        if k ~= p.z
            uprimezf(:,j,k) = u(:,j,k+1) - u(:,j,k);
        end
        if k ~= 1
            uprimezb(:,j,k) = u(:,j,k-1) - u(:,j,k);
        end
    end
end

%% Apply boundary conditions
%{
% x - direction
% both
uprimexf(m.boundary(:,:,:,1) == 1) = 0;
uprimexb(m.boundary(:,:,:,1) == 1) = 0;
% forward
uprimexf(m.boundary(:,:,:,1) == 2) = 0;
% backwards
uprimexb(m.boundary(:,:,:,1) == 3) = 0;

% y - direction
% both
uprimeyf(m.boundary(:,:,:,2) == 1) = 0;
uprimeyb(m.boundary(:,:,:,2) == 1) = 0;
% forward
uprimeyf(m.boundary(:,:,:,2) == 2) = 0;
% backwards
uprimeyb(m.boundary(:,:,:,2) == 3) = 0;

% z - direction
% both
uprimezf(m.boundary(:,:,:,3) == 1) = 0;
uprimezb(m.boundary(:,:,:,3) == 1) = 0;
% forward
uprimezf(m.boundary(:,:,:,3) == 2) = 0;
% backwards
uprimezb(m.boundary(:,:,:,3) == 3) = 0;
%}

%
for k = 1:p.z
    for j = 1:p.y
        for i = 1:p.x
            if m.boundary(i,j,k,1) ~= 0
                if m.boundary(i,j,k,1) == 1 % X-both
                    uprimexf(i,j,k) = 0;
                    uprimexb(i,j,k) = 0;
                elseif m.boundary(i,j,k,1) == 2 % X-forwards
                    uprimexf(i,j,k) = 0;
                elseif m.boundary(i,j,k,1) == 3 % X-backwards
                    uprimexb(i,j,k) = 0;
                end
            end
            if m.boundary(i,j,k,2) ~= 0
                if m.boundary(i,j,k,2) == 1 % Y-both
                    uprimeyf(i,j,k) = 0;
                    uprimeyb(i,j,k) = 0;
                elseif m.boundary(i,j,k,2) == 2 % Y-forwards
                    uprimeyf(i,j,k) = 0;
                elseif m.boundary(i,j,k,2) == 3 % Y-backwards
                    uprimeyb(i,j,k) = 0;
                end
            end
            if m.boundary(i,j,k,3) ~= 0
                if m.boundary(i,j,k,3) == 1 % Z-both
                    uprimezf(i,j,k) = 0;
                    uprimezb(i,j,k) = 0;
                elseif m.boundary(i,j,k,3) == 2 % Z-forwards
                    uprimezf(i,j,k) = 0;
                elseif m.boundary(i,j,k,3) == 3 % Z-backwards
                    uprimezb(i,j,k) = 0;
                end
            end
        end
    end
end
%}

%% Construct center differences from forward and backward differences
uprimex(1:p.x-1,:,:) = uprimexf;
uprimex = uprimex + uprimexb;
uprimey(:,1:p.y-1,:) = uprimeyf;
uprimey = uprimey + uprimeyb;
uprimez(:,:,1:p.z-1) = uprimezf;
uprimez = uprimez + uprimezb;

%% Calculate uprime
uprime = (m.stochastic_D(:,:,:,1)).*(uprimex/(deltaX.^2) + ... 
    uprimey/(deltaY.^2) + ...
    uprimez/(deltaZ.^2)) + ...
    (m.stochastic_rho(:,:,:,1)).*u.*(p.CC - u);

%% Remove anything outside brain geometry
uprime(m.mask == 0) = 0;

%% Reshape for ODE45
uprime = reshape(uprime, p.x*p.y*p.z, 1);

%% Print percent completion
%fprintf('\n%3.2f%%',100*t/p.t_end);

end

