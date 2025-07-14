% getGroupMean.m
% For the group map, we computed the circular mean across subjects for each vertex in V1–V4, 
% weighted by the full model R2 values.
clear all;

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/surfaceData/';

[~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/surfaceData/sub1/curvMLVBrain_sub1_lh_fsaverage.mgh');
%save_mgh(vol, 'orig.stripped.mgz', M,mr);
%result = MRIread('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/angleBrain_sub1_lh_fsaverage.mgh');

%% get all subjects volume in radians
% what to do with -1 and 0?
for isub=1:8
    vol_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/curvMLVBrain_sub', num2str(isub), '_lh_fsaverage.mgh']);
    vol_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/curvMLVBrain_sub', num2str(isub), '_rh_fsaverage.mgh']);
end

vol_lh(vol_lh==0)=NaN;
vol_lh(vol_lh==-1)=NaN;
vol_rh(vol_rh==0)=NaN;
vol_rh(vol_rh==-1)=NaN;

%% weighted circ mean
prefCurv_lh = [];
prefCurv_rh = [];
for ivox=1:size(vol_lh,1)
    prefCurv_lh(ivox) = mean(vol_lh(ivox,:)', 'omitnan');
end
for ivox=1:size(vol_rh,1)
    prefCurv_rh(ivox) = mean(vol_rh(ivox,:)', 'omitnan');
end


prefCurv_lh = prefCurv_lh';
prefCurv_rh = prefCurv_rh';

%% save
fileName_lh = [filedir,'curvMLVBrain__groupmean_lh_fsaverage.mgh'];
fileName_rh = [filedir, 'curvMLVBrain__groupmean_rh_fsaverage.mgh'];
save_mgh(prefCurv_lh, fileName_lh, M,mr);
save_mgh(prefCurv_rh, fileName_rh, M,mr);


%%
clear all;

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/surfaceData/';

[mghfile,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/surfaceData/sub1/curvMLVBrainbyROI_sub1_lh_fsaverage.mgh');
%save_mgh(vol, 'orig.stripped.mgz', M,mr);
%result = MRIread('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/angleBrain_sub1_lh_fsaverage.mgh');

%% get all subjects volume in radians
% what to do with -1 and 0?
for isub=1:8
    vol_lh(:,:,isub) = squeeze(load_mgh([filedir, 'sub', num2str(isub), '/curvMLVBrainbyROI_sub', num2str(isub), '_lh_fsaverage.mgh']));
    vol_rh(:,:,isub) = squeeze(load_mgh([filedir, 'sub', num2str(isub), '/curvMLVBrainbyROI_sub', num2str(isub), '_rh_fsaverage.mgh']));
end

vol_lh(vol_lh==0)=NaN;
vol_lh(vol_lh==-1)=NaN;
vol_rh(vol_rh==0)=NaN;
vol_rh(vol_rh==-1)=NaN;

%% weighted circ mean
prefCurv_lh = [];
prefCurv_rh = [];
for iroi=1:size(vol_lh,2)
    for ivox=1:size(vol_lh,1)
        prefCurv_lh(ivox,iroi) = mean(squeeze(vol_lh(ivox,iroi,:)), 'omitnan');
    end
end
for iroi=1:size(vol_rh,2)
    for ivox=1:size(vol_rh,1)
        prefCurv_rh(ivox,iroi) = mean(squeeze(vol_rh(ivox,iroi,:)), 'omitnan');
    end
end

prefCurv_lh = reshape(prefCurv_lh, size(mghfile));
prefCurv_rh = reshape(prefCurv_rh, size(mghfile));

%% save
fileName_lh = [filedir,'curvMLVBrainbyROI_groupmean_lh_fsaverage.mgh'];
fileName_rh = [filedir, 'curvMLVBrainbyROI_groupmean_rh_fsaverage.mgh'];
save_mgh(prefCurv_lh, fileName_lh, M,mr);
save_mgh(prefCurv_rh, fileName_rh, M,mr);