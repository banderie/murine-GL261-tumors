function [mask, boundary, brain, tumor, ventricle, x, y, z] = mask(filenames)
% Generates a mask for a slice of brain data given a threshold and tags
% boundaries with boudary information.

[brain,tumor,ventricle,x,y,z] = loadData(filenames.brain, filenames.tumor, filenames.ventricle);

mask = zeros(x,y,z);
boundary = zeros(x,y,z,3);
for k = 1:z
    for j = 1:y
        for i = 1:x
            if ~isnan(ventricle(i,j,k))
                brain(i,j,k) = NaN;
            end
        end
    end
end

for k = 1:z
    for j = 1:y
        for i = 1:x
            if ~isnan(brain(i,j,k))
                mask(i,j,k) = 4;
                if isnan(brain(i-1,j,k)) || isnan(brain(i+1,j,k)) || ...
                        isnan(brain(i,j-1,k)) || isnan(brain(i,j+1,k)) || ...
                        isnan(brain(i,j,k-1)) || isnan(brain(i,j,k+1))
                    mask(i,j,k) = 2;
                end
            end
        end
    end
end

for k = 2:z-1
    for j = 2:y-1
        for i = 2:x-1
            if mask(i,j,k) == 2 %On a boundary
                % x-direction
                if mask(i+1,j,k) == 0 && mask(i-1,j,k) == 0 % both
                    boundary(i,j,k,1) = 1;
                elseif mask(i+1,j,k) == 0 % forwards
                    boundary(i,j,k,1) = 2;
                elseif mask(i-1,j,k) == 0 % backwards
                    boundary(i,j,k,1) = 3;
                end
                
                % y-direction
                if mask(i,j+1,k) == 0 && mask(i,j-1,k) == 0 % both
                    boundary(i,j,k,2) = 1;
                elseif mask(i,j+1,k) == 0 % forwards
                    boundary(i,j,k,2) = 2;
                elseif mask(i,j-1,k) == 0 % backwards
                    boundary(i,j,k,2) = 3;
                end
                
                % z-direction
                if mask(i,j,k+1) == 0 && mask(i,j,k-1) == 0 % both
                    boundary(i,j,k,3) = 1;
                elseif mask(i,j,k+1) == 0 % forwards
                    boundary(i,j,k,3) = 2;
                elseif mask(i,j,k-1) == 0 % backwards
                    boundary(i,j,k,3) = 3;
                end
            end
        end
    end
end

ventricle(isnan(ventricle)) = 0;
tumor(isnan(tumor)) = 0;
% disp('Mask generated')

end