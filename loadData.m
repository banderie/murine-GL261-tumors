function [B,T,V,x,y,z] = loadData(filename_brain,filename_tumor,filename_ventricle)

Bs = load(filename_brain);
Ts = load(filename_tumor);
Vs = load(filename_ventricle);
[x,y,z] = size(Bs.Volume);
[xt,yt,zt] = size(Ts.Volume);

if x ~= xt || y ~= yt || z ~= zt
    error('\n\nDimensions of brain and tumor data do not match.\n\n'); 
end
clear xt; clear yt; clear zt;

maxB = max(max(max(Bs.Volume)));
minB = min(min(min(Bs.Volume)));
maxT = max(max(max(Ts.Volume)));
minT = min(min(min(Ts.Volume)));
maxV = max(max(max(Vs.Volume)));
minV = min(min(min(Vs.Volume)));

x = x+4;
y = y+4;
z = z+4;

%fprintf('\nMax BD = %10.10f. \nMin BD = %10.10f.\n', maxB, minB);
%fprintf('\nMax TD is %4.4f. \nMin TD = %10.10f \n', maxT, minT);
%fprintf('\nMax VD is %4.4f. \nMin VD = %10.10f \n\n', maxV, minV);
%fprintf('X: %4.0f\nY: %4.0f\nZ: %4.0f\n\n',x,y,z);

B(1:x,1:y,1:z) = NaN;
T(1:x,1:y,1:z) = NaN;
V(1:x,1:y,1:z) = NaN;
B(3:end-2,3:end-2,3:end-2) = Bs.Volume;
T(3:end-2,3:end-2,3:end-2) = Ts.Volume;
V(3:end-2,3:end-2,3:end-2) = Vs.Volume;



%T=Ts.Volume;


end