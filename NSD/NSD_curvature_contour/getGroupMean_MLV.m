% getGroupMean.m
% For the group map, we computed the weighted mean across subjects for each vertex,
% weighted by the full model R2 values.
clear all;
nSub = 8;
hemis = {'lh','rh'};
modelList = { ...
    % 'curvMLVBrain_max', ...
    % 'curvMLVBrain_wmean', ...
    'curvMLVBrain_max_posR2', ...
    % 'curvMLVBrain_wmean_posR2'
    };

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/surfaceData/';

[~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/surfaceData/sub1/curvMLVBrain_max_sub1_lh_fsaverage.mgh');

for imodel = 1:numel(modelList)
    modelName = modelList{imodel};
   
    %% load volumes
    clear vol_lh vol_rh R2_lh R2_rh

    for isub = 1:nSub
        vol_lh(:,isub) = load_mgh([filedir,'sub',num2str(isub),'/', ...
            modelName,'_sub',num2str(isub),'_lh_fsaverage.mgh']);

        vol_rh(:,isub) = load_mgh([filedir,'sub',num2str(isub),'/', ...
            modelName,'_sub',num2str(isub),'_rh_fsaverage.mgh']);

        % R2 (same for all models unless you have model-specific R2)
        R2_lh(:,isub) = load_mgh([filedir,'sub',num2str(isub),'/', ...
            'curvMLVBrain_R2_sub',num2str(isub),'_lh_fsaverage.mgh']);

        R2_rh(:,isub) = load_mgh([filedir,'sub',num2str(isub),'/', ...
            'curvMLVBrain_R2_sub',num2str(isub),'_rh_fsaverage.mgh']);
    end


    %% weighted mean
    R2_lh_pos = R2_lh;
    % R2_lh_pos = R2_lh_pos - min(R2_lh_pos, [], 2);
    R2_rh_pos = R2_rh;
    % R2_rh_pos = R2_rh_pos - min(R2_rh_pos, [], 2);
    R2_lh_pos(R2_lh_pos <= 0) = NaN;
    R2_rh_pos(R2_rh_pos <= 0) = NaN;


    prefCurv_lh = nansum(vol_lh .* R2_lh_pos, 2) ./ nansum(R2_lh_pos, 2);
    prefCurv_rh = nansum(vol_rh .* R2_rh_pos, 2) ./ nansum(R2_rh_pos, 2);



    %% save
    save_mgh(prefCurv_lh, ...
        [filedir, modelName, '_groupmean_lh_fsaverage.mgh'], M, mr);

    save_mgh(prefCurv_rh, ...
        [filedir, modelName, '_groupmean_rh_fsaverage.mgh'], M, mr);


    % %% Top-50% R² voxel survival
    % 
    % meanR2_lh = nanmean(R2_lh, 2);
    % meanR2_rh = nanmean(R2_rh, 2);
    % 
    % thr_lh = prctile(meanR2_lh, 50);
    % thr_rh = prctile(meanR2_rh, 50);
    % 
    % mask_lh = meanR2_lh >= thr_lh;
    % mask_rh = meanR2_rh >= thr_rh;
    % 
    % prefCurv_lh_top50 = prefCurv_lh;
    % prefCurv_rh_top50 = prefCurv_rh;
    % 
    % prefCurv_lh_top50(~mask_lh) = NaN;
    % prefCurv_rh_top50(~mask_rh) = NaN;
    % 
    % 
    % save_mgh(prefCurv_lh_top50, ...
    %     [filedir, modelName, '_groupmean_lh_fsaverage_top50R2.mgh'], M, mr);
    % 
    % save_mgh(prefCurv_rh_top50, ...
    %     [filedir, modelName, '_groupmean_rh_fsaverage_top50R2.mgh'], M, mr);
end

%%
clear all;
nSub = 8;
hemis = {'lh','rh'};
modelList = { ...
    % 'curvMLVBrainbyROI_max', ...
    % 'curvMLVBrainbyROI_wmean', ...
    'curvMLVBrainbyROI_max_posR2', ...
    % 'curvMLVBrainbyROI_wmean_posR2'
    };

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/surfaceData/';

[~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/surfaceData/sub1/curvMLVBrain_max_sub1_lh_fsaverage.mgh');

for imodel = 1:numel(modelList)
    modelName = modelList{imodel};
   
    %% load volumes
    clear vol_lh vol_rh R2_lh R2_rh

    for isub = 1:nSub
        vol_lh(:,:,isub) = squeeze(load_mgh([filedir,'sub',num2str(isub),'/', ...
            modelName,'_sub',num2str(isub),'_lh_fsaverage.mgh']));

        vol_rh(:,:,isub) = squeeze(load_mgh([filedir,'sub',num2str(isub),'/', ...
            modelName,'_sub',num2str(isub),'_rh_fsaverage.mgh']));

        % R2 (same for all models unless you have model-specific R2)
        R2_lh(:,:,isub) = squeeze(load_mgh([filedir,'sub',num2str(isub),'/', ...
            'curvMLVBrain_R2_sub',num2str(isub),'_lh_fsaverage.mgh']));

        R2_rh(:,:,isub) = squeeze(load_mgh([filedir,'sub',num2str(isub),'/', ...
            'curvMLVBrain_R2_sub',num2str(isub),'_rh_fsaverage.mgh']));
    end


    %% weighted mean
    R2_lh_pos = R2_lh;
    % R2_lh_pos = R2_lh_pos - min(R2_lh_pos, [], 3);
    R2_rh_pos = R2_rh;
    % R2_rh_pos = R2_rh_pos - min(R2_rh_pos, [], 3);
    R2_lh_pos(R2_lh_pos <= 0) = NaN;
    R2_rh_pos(R2_rh_pos <= 0) = NaN;


    prefCurv_lh = nansum(vol_lh .* R2_lh_pos, 2) ./ nansum(R2_lh_pos, 2);
    prefCurv_rh = nansum(vol_rh .* R2_rh_pos, 2) ./ nansum(R2_rh_pos, 2);



    %% save
    save_mgh(prefCurv_lh, ...
        [filedir, modelName, '_groupmean_lh_fsaverage.mgh'], M, mr);

    save_mgh(prefCurv_rh, ...
        [filedir, modelName, '_groupmean_rh_fsaverage.mgh'], M, mr);


    % %% Top-50% R² voxel survival
    % 
    % meanR2_lh = nanmean(R2_lh, 2);
    % meanR2_rh = nanmean(R2_rh, 2);
    % 
    % thr_lh = prctile(meanR2_lh, 50);
    % thr_rh = prctile(meanR2_rh, 50);
    % 
    % mask_lh = meanR2_lh >= thr_lh;
    % mask_rh = meanR2_rh >= thr_rh;
    % 
    % prefCurv_lh_top50 = prefCurv_lh;
    % prefCurv_rh_top50 = prefCurv_rh;
    % 
    % prefCurv_lh_top50(~mask_lh) = NaN;
    % prefCurv_rh_top50(~mask_rh) = NaN;
    % 
    % 
    % save_mgh(prefCurv_lh_top50, ...
    %     [filedir, modelName, '_groupmean_lh_fsaverage_top50R2.mgh'], M, mr);
    % 
    % save_mgh(prefCurv_rh_top50, ...
    %     [filedir, modelName, '_groupmean_rh_fsaverage_top50R2.mgh'], M, mr);
end
