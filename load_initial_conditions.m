function u0 = load_initial_conditions(filenames,q)
%INITIAL CONDITIONS creates initial conditions for tumor
%   Creates initial conditions based on observed tumors so that the
%   simulation can be started at a later scan.

[raw,tumor,ventricle,x,y,z] = loadData(filenames.brain,filenames.initial_conditions,filenames.ventricle);

%
for k = 1:z
    for j = 1:y
        for i = 1:x
            if ~isnan(ventricle(i,j,k))
                raw(i,j,k) = NaN;
            end
        end
    end
end
%}

tumor = scaledata(tumor,0.84,1);
tumor(isnan(tumor)) = 0;
% 
% figure
% subplot(1,2,1);
% imagesc(tumor(:,:,q.zslice));
% caxis([0,1]);
% subplot(1,2,2)
% 
% imagesc(raw(:,:,q.zslice))
% colormap(bone(128));


%fprintf('\nInitial condition generated from tumor. \nShown at specified zslice: %2.0f \nType "return" to continue.\n', q.zslice);
%disp(['Initial tumor volume = ', num2str(sum(sum(sum(tumor > 0))))]);
%u0(~isnan(raw)) = tumor(~isnan(raw));
u0 = tumor;


end

