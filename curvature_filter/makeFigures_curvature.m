% makeFigures_length.m
clear all;
cd('/home/hanseohe/Documents/GitHub/nsdOtopy/curvature');

prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/prfsample_Curv/';  
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};

%% 1.  Coef voxel preference %%

for isub = 1:1
    load(fullfile([prffolder, 'voxCurvCoef_sub', num2str(isub), '.mat']));
 
    subName = ['sub', num2str(isub)];
    voxPref.(subName)=maxCoefCurv;
    % %% check consistency between splits
    % for visualRegion = 1:7
    %     thisfield = ['V', num2str(visualRegion)];
    % 
    %     figure;
    %     subplot(1,2,1)
    %     bar(mean_split1.(combinedRoiNames{visualRegion}),1);
    %     title('split1');
    % 
    %     subplot(1,2,2)
    %     bar(mean_split2.(combinedRoiNames{visualRegion}),1);
    %     title('split2');
    % 
    %     totalTitle = ['sub', num2str(isub), ' V', num2str(visualRegion)];
    %     sgtitle(totalTitle);
    % 
    %     % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/split_sub', num2str(isub), 'v', num2str(visualRegion), '.png']);
    % end
    % 
    %% get average for each ROI
    figure;
    subplot(4,2,1)
    bar(mean_all.V1,1);
    title('v1');

    subplot(4,2,2)
    bar(mean_all.V2,1);
    title('v2');

    subplot(4,2,3)
    bar(mean_all.V3,1);
    title('v3');

    subplot(4,2,4)
    bar(mean_all.hV4,1);
    title('v4');

    subplot(4,2,5)
    bar(mean_all.OPA,1);
    title('OPA');

    subplot(4,2,6)
    bar(mean_all.PPA,1);
    title('PPA');

    subplot(4,2,7)
    bar(mean_all.RSC,1);
    title('RSC');

    totalTitle = ['sub', num2str(isub), ' average for ROI'];
    sgtitle(totalTitle);
    % 
    % % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/all_sub', num2str(isub), '.png']);

    %% preferred length from coef
    figure;
    subplot(4,2,1)
    histogram(maxCoefCurv.V1);
    title('v1');

    subplot(4,2,2)
    histogram(maxCoefCurv.V2);
    title('v2');

    subplot(4,2,3)
    histogram(maxCoefCurv.V3);
    title('v3');

    subplot(4,2,4)
    histogram(maxCoefCurv.hV4);
    title('v4');

    subplot(4,2,5)
    histogram(maxCoefCurv.OPA);
    title('OPA');

    subplot(4,2,6)
    histogram(maxCoefCurv.PPA);
    title('PPA');

    subplot(4,2,7)
    histogram(maxCoefCurv.RSC);
    title('RSC');


    totalTitle = ['sub', num2str(isub), ' number of preferred curvature bin'];
    sgtitle(totalTitle);

    saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/curv_hist/prefCurv_sub', num2str(isub), '.png']);
end

%% combine all participants - preferred length from coef
allvoxPref={};
for visualRegion = 1:7
    currRoi = [];
    thisfield = combinedRoiNames{visualRegion};
    for isub = 1:8
        subName = ['sub', num2str(isub)];

        currRoi = [currRoi, voxPref.(subName).(thisfield)];

    end
    allvoxPref{visualRegion} = currRoi;
end

figure;
subplot(4,2,1)
histogram(allvoxPref{1});
title('v1');

subplot(4,2,2)
histogram(allvoxPref{2});
title('v2');

subplot(4,2,3)
histogram(allvoxPref{3});
title('v3');

subplot(4,2,4)
histogram(allvoxPref{4});
title('v4');

subplot(4,2,5)
histogram(allvoxPref{5});
title('OPA');

subplot(4,2,6)
histogram(allvoxPref{6});
title('PPA');

subplot(4,2,7)
histogram(allvoxPref{7});
title('RSC');


totalTitle = 'all sub number of preferred curvature bin';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/curv_hist/prefCurv_allsub.png');

