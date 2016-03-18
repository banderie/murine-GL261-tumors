% param estimator for the full mouse case
% literature estimated values of D range from 2-3.5 mm^2/hr, rho approx
% 0.01 .hr. I convert to D and rho for this simulatioin with 600*rho,
% 36000*D
for mm=1:2

initialguess=[2.5,3; 2.5, 6; 2.5, 9; 23.75, 3; 23.75, 6; 23.75, 9; 50, 3; 50, 6; 50, 9];
err1=zeros(9,1);
% optimizedgbmtp1=zeros(9,2);
% err2=zeros(9,1);
% optimizedgbmtp2=zeros(9,2);
% err3=zeros(9,1);
% optimizedgbmtp3=zeros(9,2);
% err4=zeros(9,1);
% optimizedgbmtp4=zeros(9,2);

% first time point
for k=1:9
[optimizedgbmtp1(k,:),err1(k)]=fminsearch(@(par) gbm_minimizerhyp2(par,0,mm),initialguess(k,:));
[mm,k,1]
title=sprintf('S1G3M%d_errorstp1hyp2.mat',mm);
save(title,'err1');
title=sprintf('S1G3M%d_optimizedtp1hyp2.mat',mm);
save(title,'optimizedgbmtp1');
end
[blah,j]=(min(err1));
err1final=gbm_minimizereachtimepts(optimizedgbmtp1(j,:),1,mm); % so that it saves for tp 2!
title=sprintf('S1G3M%d_errorstp1hyp2.mat',mm);
save(title,'err1');
title=sprintf('S1G3M%d_optimizedtp1hyp2.mat',mm);
save(title,'optimizedgbmtp1');


%Totalerr=err1final+err2final+err3final+err4final
 err1final
%  err2final
%  err3final
%  err4final
 optimizedgbmtp1
%  optimizedgbmtp2
%  optimizedgbmtp3
%  optimizedgbmtp4
end