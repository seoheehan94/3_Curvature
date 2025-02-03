%% prfSampleModel_curvature

addpath(genpath('/home/hanseohe/Documents/GitHub/stimulusVignetting'))
addpath(genpath('/home/hanseohe/Documents/GitHub/nsdOtopy'))
% 
% for sub = 5:8
%     tic;
%     fprintf('%s. %d. %d ...\n','prfSampleModel',sub,1);
%     prfSampleModel_length(sub,1);
%     toc;
% 
%     tic;
%     fprintf('%s. %d. %d ...\n','prfSampleModel',sub,2);
%     prfSampleModel_length(sub,2);
%     toc;
% 
%     tic;
%     fprintf('%s. %d. %d ...\n','prfSampleModel',sub,3);
%     prfSampleModel_length(sub,3);
%     toc;
% 
%     tic;
%     fprintf('%s. %d. %d ...\n','prfSampleModel',sub,4);
%     prfSampleModel_length(sub,4);
%     toc;
% end

%% regressPrfSplit
for sub = 2:8
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 1);
regressPrfSplit_curvMLV(sub, 1);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 2);
regressPrfSplit_curvMLV(sub, 2);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 3);
regressPrfSplit_curvMLV(sub, 3);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 4);
regressPrfSplit_curvMLV(sub, 4);

fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 5);
regressPrfSplit_curvMLV(sub, 5);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 6);
regressPrfSplit_curvMLV(sub, 6);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 7);
regressPrfSplit_curvMLV(sub, 7);

end

%% getVoxPref
% for sub = 1:8
%     getVoxPref_curvature(sub,7)
% end