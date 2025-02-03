clear all;

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/';

filterCurv_lh = squeeze(load_mgh([filedir, 'Curvature/surfaceData/curvBrainbyROI_groupmean_lh_fsaverage.mgh']));
filterCurv_rh = squeeze(load_mgh([filedir, 'Curvature/surfaceData/curvBrainbyROI_groupmean_rh_fsaverage.mgh']));
MLVCurv_lh = squeeze(load_mgh([filedir, 'Curvature_MLV/surfaceData/curvMLVBrainbyROI_groupmean_lh_fsaverage.mgh']));
MLVCurv_rh = squeeze(load_mgh([filedir, 'Curvature_MLV/surfaceData/curvMLVBrainbyROI_groupmean_rh_fsaverage.mgh']));

pearsonTotal = [];
spearmanTotal = [];

for cur = 1:7
    curFilter = filterCurv_lh(:,cur);
    curMLV = MLVCurv_lh(:,cur);
    validIdx = ~isnan(curFilter) & ~isnan(curMLV);
    curFilter = curFilter(validIdx);
    curMLV = curMLV(validIdx);

    [pearsonTotal(cur,1), pearsonTotal(cur,2)] = corr(curFilter, curMLV);


    [spearmanTotal(cur,1), spearmanTotal(cur,2)] = corr(curFilter, curMLV, 'Type', 'Spearman');

end

saveName = [filedir, 'Curvature_MLV/compareFilterMLVbrain.mat'];
save(saveName, "spearmanTotal","pearsonTotal");
%% Is my data normally distributed?
% Example for visual inspection
cur = 2; % Choose which curvature data to inspect
curFilter = filterCurv_lh(:, cur);
curMLV = MLVCurv_lh(:, cur);

% Remove NaNs
validIdx = ~isnan(curFilter) & ~isnan(curMLV);
curFilter = curFilter(validIdx);
curMLV = curMLV(validIdx);

% Histogram
figure;
subplot(2, 2, 1);
histogram(curFilter, 30); % Adjust the number of bins as needed
title('Histogram of curFilter');
xlabel('Values');
ylabel('Frequency');

subplot(2, 2, 2);
histogram(curMLV, 30);
title('Histogram of curMLV');
xlabel('Values');
ylabel('Frequency');

% Q-Q Plot
subplot(2, 2, 3);
qqplot(curFilter);
title('Q-Q Plot of curFilter');

subplot(2, 2, 4);
qqplot(curMLV);
title('Q-Q Plot of curMLV');

% Example using the Lilliefors test (an adaptation of the Kolmogorov-Smirnov test for normality)
[hFilter, pFilter] = lillietest(curFilter);
[hMLV, pMLV] = lillietest(curMLV);

if hFilter == 0
    fprintf('curFilter data is normally distributed (Lilliefors test, p = %.4f).\n', pFilter);
else
    fprintf('curFilter data is not normally distributed (Lilliefors test, p = %.4f).\n', pFilter);
end

if hMLV == 0
    fprintf('curMLV data is normally distributed (Lilliefors test, p = %.4f).\n', pMLV);
else
    fprintf('curMLV data is not normally distributed (Lilliefors test, p = %.4f).\n', pMLV);
end
