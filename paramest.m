%% Driver for optimization of Mouse GBM Tumor
%S3G2

for mm=1:1
initialguess=[2.5,3; 2.5, 6; 2.5, 9; 23.75, 3; 23.75, 6; 23.75, 9; 50, 3; 50, 6; 50, 9];
err=zeros(9,1);
optimizedgbmhyp1=zeros(9,2);

for k =1:9
[optimizedgbmhyp1(k,:),err(k)]=fminsearch(@(par) gbm_minimizer(par,mm),initialguess(k,:)) 
[mm,k]
end
title=sprintf('S1G3M%d_errorshyp1IC70.mat',mm);
save(title,'err');
title=sprintf('S1G3M%d_optimizedhyp1IC70.mat',mm);
save(title,'optimizedgbmhyp1');
end

 30 D 0.2485
 70 D 0.1702
 30 rho 0.2836
 70 rho 0.1014
 8%
 
 5%