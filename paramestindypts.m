% param estimator for the full mouse case
% literature estimated values of D range from 2-3.5 mm^2/hr, rho approx
% 0.01 .hr. I convert to D and rho for this simulatioin with 600*rho,
% 36000*D
for mm=3:3

initialguess=[2.5,3; 2.5, 6; 2.5, 9; 23.75, 3; 23.75, 6; 23.75, 9; 50, 3; 50, 6; 50, 9];
err1=zeros(9,1);
optimizedgbmtp1=zeros(9,2);
err2=zeros(9,1);
optimizedgbmtp2=zeros(9,2);
err3=zeros(9,1);
optimizedgbmtp3=zeros(9,2);
err4=zeros(9,1);
optimizedgbmtp4=zeros(9,2);

% % first time point
% for k=1:9
% [optimizedgbmtp1(k,:),err1(k)]=fminsearch(@(par) gbm_minimizerindypts(par,1,mm),initialguess(k,:));
% [mm,k,1]
% end
% title=sprintf('S1G3M%d_errorstp1hyp3a.mat',mm);
% save(title,'err1');
% title=sprintf('S1G3M%d_optimizedtp1hyp3a.mat',mm);
% save(title,'optimizedgbmtp1');
% 
% %initialguess=optimizedgbmtp2;
% 
% % second time point
% for k =1:9
% [optimizedgbmtp2(k,:),err2(k)]=fminsearch(@(par) gbm_minimizerindypts(par,2,mm),initialguess(k,:)); 
% [mm,k,2]
% end
% title=sprintf('S1G3M%d_errorstp2hyp3a.mat',mm);
% save(title,'err2');
% title=sprintf('S1G3M%d_optimizedtp2hyp3a.mat',mm);
% save(title,'optimizedgbmtp2');
% %initialguess=optimizedgmbtp3;

% third time point
for k=1:9
[optimizedgbmtp3(k,:),err3(k)]=fminsearch(@(par) gbm_minimizerindypts(par,3,mm),initialguess(k,:)); 
[mm,k,3]
end
title=sprintf('S1G3M%d_errorstp3hyp3a.mat',mm);
save(title,'err3');
title=sprintf('S1G3M%d_optimizedtp3hyp3a.mat',mm);
save(title,'optimizedgbmtp3');
%initialguess=optimizedgbmtp4;

% fourth time point
for k =1:9
[optimizedgbmtp4(k,:),err4(k)]=fminsearch(@(par) gbm_minimizerindypts(par,4,mm),initialguess(k,:)); 
[mm,k,4]
end
title=sprintf('S1G3M%d_errorstp4hyp3a.mat',mm);
save(title,'err4');
title=sprintf('S1G3M%d_optimizedtp4hyp3a.mat',mm);
save(title,'optimizedgbmtp4');
%Totalerr=err1final+err2final+err3final+err4final
% err1final
% err2final
% err3final
% err4final
% optimizedgbmtp1
% optimizedgbmtp2
% optimizedgbmtp3
% optimizedgbmtp4
end