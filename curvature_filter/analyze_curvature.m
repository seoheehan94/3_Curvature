%analyze_curvature.m

%   uses files created by: regressPrfSplit_curvate.m
%   creates files used by:
clear all;
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/prfsample_Curv/';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};


%% 1.  Coef voxel preference %%
for isub = 1:8
    cd('/home/hanseohe/Documents/GitHub/nsdOtopy/Curvature');

    fprintf('%d ...\n',isub);
    clearvars -except isub roiNames combinedRoiNames prffolder
    %% set up

    bandpass = 1; bandMin = 1; bandMax = 1;
    bandpassStr = '';
    if bandpass
        bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
    end

    mean_split1 = struct;
    mean_split2 = struct;
    mean_all = struct;
    allCurvCoef = struct;
    maxCoefCurv = struct;

    %% load file
    for visualRegion = 1:7
        thisfield = combinedRoiNames{visualRegion};
        load(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
            'nsd', 'numLevels', 'numCurvs','rois','nvox','roiPrf','nsplits');

        % get mean coef
        if visualRegion == 4 || visualRegion == 5 || visualRegion == 6 || visualRegion == 7

            mean_split1.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(1,:,1:8),2));
            mean_split2.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(2,:,1:8),2));

            % average splits
            nsd.voxCurvCoef{1}(nsplits+1,:,:) = mean(nsd.voxCurvCoef{1},1);
            mean_all.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(3,:,1:8),2));
            allCurvCoef.(thisfield) = nsd.voxCurvCoef{1}(:,:,1:8);

        else %for other regions, combine ventral and dorsal
            % combine ventral and dorsal
            oldNsd = nsd;
            nsd.voxCurvCoef{1} = [];
            nsd.voxCurvCoef{2} = [];
            for iroi=1:length(rois)
                nsd.voxCurvCoef{1} = cat(2,nsd.voxCurvCoef{1},oldNsd.voxCurvCoef{iroi});
            end

            mean_split1.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(1,:,1:8),2));
            mean_split2.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(2,:,1:8),2));

            % average splits
            nsd.voxCurvCoef{1}(nsplits+1,:,:) = mean(nsd.voxCurvCoef{1},1);
            mean_all.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(3,:,1:8),2));
            allCurvCoef.(thisfield) = nsd.voxCurvCoef{1}(:,:,1:8);
        end

        % find preferred curvature for each voxel: find the bin with max coef for
        % each voxel
        [~, nvox, ~] = size(allCurvCoef.(thisfield));
        for ivox = 1: nvox
            [~, icurv] = max(allCurvCoef.(thisfield)(3,ivox,:));
            maxCoefCurv.(thisfield)(1,ivox) = icurv;
        end
    end

    saveName = [prffolder, 'voxCurvCoef_sub', num2str(isub), '.mat'];
    save(saveName, 'allCurvCoef', 'mean_all', 'mean_split1', 'mean_split2','maxCoefCurv');

end

%% make a brain volume
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoicd ../..Names = {'V1','V2','V3','hV4','OPA','PPA','RSC'};
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/prfsample_Curv/';

for isub = 1:8
    clearvars -except isub roiNames combinedRoiNames prffolder

    % load(fullfile([prffolder, 'voxCurvCoef_sub', num2str(isub), '.mat']));
    % 
    % % save all ROIs to create overlay
    % roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    % visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
    % visRoiData = niftiread(visualRoisFile);
    % placesRoisFile = fullfile(roifolder,'roi/floc-places.nii.gz'); %OPA, PPA, RSC
    % placeRoiData = niftiread(placesRoisFile);
    % 
    % allRoiData = visRoiData;
    % allRoiData(placeRoiData == 1) = 8;
    % allRoiData(placeRoiData == 2) = 9;
    % allRoiData(placeRoiData == 3) = 10;
    % 
    % ourBrain = allRoiData;
    % ourBrain(ourBrain == 2) = 1;
    % ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
    % ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
    % ourBrain(ourBrain == 7) = 4;
    % ourBrain(ourBrain == 8) = 5;
    % ourBrain(ourBrain == 9) = 6;
    % ourBrain(ourBrain == 10) = 7;
    % 
    % % make a brain volume
    % newBrain = ourBrain;
    % newBrain(newBrain > 0) = 0;
    % for visualRegion = 1:7
    %     curOurBrain = ourBrain;
    %     % if visualRegion == 2
    %     %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
    %         % elseif visualRegion == 3
    %         %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
    %     % end
    %     thisfield = combinedRoiNames{visualRegion};
    %     newBrain(curOurBrain == visualRegion) = maxCoefCurv.(thisfield)(1,:);
    % end
    % newBrain(newBrain < 0) = -1;
    % 
    % for visualRegion = 1:7
    %     curOurBrain = ourBrain;
    %     % if visualRegion == 2
    %     %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
    %     % elseif visualRegion == 3
    %     %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
    %     % end
    %     curNewBrain = curOurBrain;
    %     curNewBrain(curOurBrain ~= visualRegion) = -1;
    %     thisfield = combinedRoiNames{visualRegion};
    % 
    %     curNewBrain(curOurBrain == visualRegion) = maxCoefCurv.(thisfield)(1,:);
    % 
    %     newBrainbyROI(:,:,:,visualRegion) =curNewBrain;
    % end
    % 
    % saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/brainVolume/curvBrain_sub', num2str(isub), '.mat'];
    % save(saveName, 'newBrain');
    % 
    % saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/brainVolume/curvBrainbyROI_sub', num2str(isub), '.mat'];
    % save(saveName, 'newBrainbyROI');
    % 
    % 
    % 
    % %% save afni file
    % % size(newBrain)
    cd('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/brainVolume');
    % %3dinfo betas_session01.nii.gz
    % %3dcalc -a betas_session01.nii.gz[1] -expr a -prefix oneBeta
    % % command = ['3dinfo betas_session01_sub', num2str(isub), '.nii.gz'];
    % % system(command);
    % % command = ['3dcalc -a betas_session01_sub', num2str(isub), '.nii.gz[1] -expr a -prefix oneBeta_sub', num2str(isub)];
    % % system(command);
    % 
    % currOneBeta = ['oneBeta_sub', num2str(isub), '+orig'];
    % [err,V,Info] = BrikLoad(currOneBeta);
    % 
    % Info.RootName = ['curvBrain_sub', num2str(isub), '+orig'];
    % opt.Prefix = ['curvBrain_sub', num2str(isub)];
    % WriteBrik(newBrain,Info,opt);
    % Info.RootName = ['curvBrainbyROI_sub', num2str(isub), '+orig'];
    % opt.Prefix = ['curvBrainbyROI_sub', num2str(isub)];
    % WriteBrik(newBrainbyROI,Info,opt);

 % %% save nifti
 %  load(['curvBrain_sub', num2str(isub), '.mat']);
 % 
 %  niftiwrite(newBrain,['curvBrain', '_sub', num2str(isub),'.nii']);
 %  info_old = niftiinfo(['betas_session01_sub', num2str(isub),'.nii.gz']);
 %  info_new = niftiinfo(['curvBrain', '_sub', num2str(isub),'.nii']);
 % 
 %  info_new.PixelDimensions = [1.8, 1.8, 1.8];
 %  info_new.TransformName = info_old.TransformName;
 %  info_new.SpatialDimension = info_old.SpatialDimension;
 %  info_new.Transform = info_old.Transform;
 %  info_new.Qfactor = info_old.Qfactor;
 %  info_new.AuxiliaryFile = info_old.AuxiliaryFile;
 %  info_new.raw.pixdim = info_old.raw.pixdim;
 %  info_new.raw.aux_file = info_old.raw.aux_file;
 %  info_new.raw.sform_code = info_old.raw.sform_code;
 %  info_new.raw.srow_x = info_old.raw.srow_x;
 %  info_new.raw.srow_y = info_old.raw.srow_y;
 %  info_new.raw.srow_z = info_old.raw.srow_z;
 % 
 %  niftiwrite(newBrain,['curvBrain', '_sub', num2str(isub),'.nii'],info_new);

  %% save nifti
  load(['curvBrainbyROI_sub', num2str(isub), '.mat']);

  niftiwrite(newBrainbyROI,['curvBrainbyROI', '_sub', num2str(isub),'.nii']);
  info_old = niftiinfo(['betas_session01_sub', num2str(isub),'.nii.gz']);
  info_new = niftiinfo(['curvBrainbyROI', '_sub', num2str(isub),'.nii']);

  % info_new.ImageSize = size(newBrainbyROI);
  info_new.PixelDimensions = info_old.PixelDimensions;
  info_new.TransformName = info_old.TransformName;
  info_new.SpatialDimension = info_old.SpatialDimension;
  info_new.Transform = info_old.Transform;
  info_new.Qfactor = info_old.Qfactor;
  info_new.AuxiliaryFile = info_old.AuxiliaryFile;
  info_new.raw.pixdim = info_old.raw.pixdim;
  info_new.raw.aux_file = info_old.raw.aux_file;
  info_new.raw.sform_code = info_old.raw.sform_code;
  info_new.raw.srow_x = info_old.raw.srow_x;
  info_new.raw.srow_y = info_old.raw.srow_y;
  info_new.raw.srow_z = info_old.raw.srow_z;

  niftiwrite(newBrainbyROI,['curvBrainbyROI', '_sub', num2str(isub),'.nii'],info_new);


end