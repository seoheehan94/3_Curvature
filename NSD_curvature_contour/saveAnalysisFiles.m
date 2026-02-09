
clear all;
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/prfsample_CurvMLV/';
save_brainfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_MLV/brainVolume';
save_R2folder = '/home/hanseohe/Documents/GitHub/3_Curvature/NSD_curvature_contour/R2_analysis/';


%% 1.  get mean R2 %%
totalR2CurvSplit = [];
for isub = 1:8
    fprintf('isub:%d. ...\n',isub);
    load([prffolder 'voxCurvCoef_sub' num2str(isub) '.mat']);
    % total values of R
    totalR2CurvSplit{end+1} = roiNsdCurvR2;

end

curRoiR2CurvSplit = [];
curV1R2CurvSplit  = [];
curV2R2CurvSplit  = [];
curV3R2CurvSplit  = [];
curhV4R2CurvSplit  = [];
curPPAR2CurvSplit = [];
curOPAR2CurvSplit = [];
curRSCR2CurvSplit = [];

for j = 1:numel(totalR2CurvSplit)

    % loop over all ROIs for the "all-roi" vector
    for r = 1:numel(combinedRoiNames)
        roiName = combinedRoiNames{r};
        curRoiR2CurvSplit = [curRoiR2CurvSplit, ...
            totalR2CurvSplit{j}.(roiName)];
    end

    % keep ROI-specific vectors if you still want them
    curV1R2CurvSplit  = [curV1R2CurvSplit,  totalR2CurvSplit{j}.V1];
    curV2R2CurvSplit  = [curV2R2CurvSplit,  totalR2CurvSplit{j}.V2];
    curV3R2CurvSplit  = [curV3R2CurvSplit,  totalR2CurvSplit{j}.V3];
    curhV4R2CurvSplit  = [curhV4R2CurvSplit,  totalR2CurvSplit{j}.hV4];
    curPPAR2CurvSplit = [curPPAR2CurvSplit, totalR2CurvSplit{j}.PPA];
    curOPAR2CurvSplit = [curOPAR2CurvSplit, totalR2CurvSplit{j}.OPA];
    curRSCR2CurvSplit = [curRSCR2CurvSplit, totalR2CurvSplit{j}.RSC];
end
%
writematrix(curRoiR2CurvSplit', ...
    fullfile(save_R2folder, ['allroiR2_MLV','.csv']));

for r = 1:numel(combinedRoiNames)
    roiName = combinedRoiNames{r};
    dataVar = eval(['cur' roiName 'R2CurvSplit']);

    writematrix(dataVar', ...
        fullfile(save_R2folder, [roiName 'R2_MLV', '.csv']));
end

% mean R² voxels
meanR2Curv = struct;
meanR2Curv.all = mean(curRoiR2CurvSplit, 'omitnan');

for r = 1:numel(combinedRoiNames)
    roiName = combinedRoiNames{r};
    dataVar = eval(['cur' roiName 'R2CurvSplit']);
    meanR2Curv.(roiName) = mean(dataVar, 'omitnan');
end

% number of positive R² voxels
nPosR2Curv = struct;

validAll = ~isnan(curRoiR2CurvSplit);
nPosR2Curv.all = sum(curRoiR2CurvSplit(validAll) > 0);
nPosR2Curv.all_total = sum(validAll);

for r = 1:numel(combinedRoiNames)
    roiName = combinedRoiNames{r};
    dataVar = eval(['cur' roiName 'R2CurvSplit']);

    valid = ~isnan(dataVar);

    nPosR2Curv.(roiName) = sum(dataVar(valid) > 0);
    nPosR2Curv.([roiName '_total']) = sum(valid);
end

% proportion of positive R² voxels
propPosR2Curv = struct;

propPosR2Curv.all = ...
    nPosR2Curv.all / nPosR2Curv.all_total;

for r = 1:numel(combinedRoiNames)
    roiName = combinedRoiNames{r};
    propPosR2Curv.(roiName) = ...
        nPosR2Curv.(roiName) / nPosR2Curv.([roiName '_total']);
end

%mean R² of positive voxels only
meanPosR2Curv = struct;

maskAll = curRoiR2CurvSplit > 0 & ~isnan(curRoiR2CurvSplit);
meanPosR2Curv.all = mean(curRoiR2CurvSplit(maskAll), 'omitnan');

for r = 1:numel(combinedRoiNames)
    roiName = combinedRoiNames{r};
    dataVar = eval(['cur' roiName 'R2CurvSplit']);

    mask = dataVar > 0 & ~isnan(dataVar);

    if any(mask)
        meanPosR2Curv.(roiName) = mean(dataVar(mask), 'omitnan');
    else
        meanPosR2Curv.(roiName) = NaN;
    end
end

%% make a brain volume

for isub = 1:8
    fprintf('isub:%d. ...\n',isub);
    clearvars -except isub roiNames combinedRoiNames prffolder save_brainfolder

    load(fullfile([prffolder, 'voxCurvCoef_sub', num2str(isub), '.mat']));

    % save all ROIs to create overlay
    roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
    visRoiData = niftiread(visualRoisFile);
    placesRoisFile = fullfile(roifolder,'roi/floc-places.nii.gz'); %OPA, PPA, RSC
    placeRoiData = niftiread(placesRoisFile);

    allRoiData = visRoiData;
    allRoiData(placeRoiData == 1) = 8;
    allRoiData(placeRoiData == 2) = 9;
    allRoiData(placeRoiData == 3) = 10;

    ourBrain = allRoiData;
    ourBrain(ourBrain == 2) = 1;
    ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
    ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
    ourBrain(ourBrain == 7) = 4;
    ourBrain(ourBrain == 8) = 5;
    ourBrain(ourBrain == 9) = 6;
    ourBrain(ourBrain == 10) = 7;

    % maxcoef
    newBrain_max = ourBrain;
    newBrain_max(newBrain_max <= 0) = NaN;
    newBrain_max(newBrain_max > 0) = 0;
    for visualRegion = 1:7
        curOurBrain = ourBrain;
        % if visualRegion == 2
        %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
        % elseif visualRegion == 3
        %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
        % end
        thisfield = combinedRoiNames{visualRegion};
        newBrain_max(curOurBrain == visualRegion) = maxCoefCurv.(thisfield)(1,:);
    end

    for visualRegion = 1:7
        curOurBrain = ourBrain;
        % if visualRegion == 2
        %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
        % elseif visualRegion == 3
        %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
        % end
        curNewBrain = curOurBrain;
        curNewBrain(curOurBrain ~= visualRegion) = NaN;
        thisfield = combinedRoiNames{visualRegion};

        curNewBrain(curOurBrain == visualRegion) = maxCoefCurv.(thisfield)(1,:);

        newBrainbyROI_max(:,:,:,visualRegion) =curNewBrain;
    end


    % wmean
    newBrain_wmean = ourBrain;
    newBrain_wmean(newBrain_wmean <= 0) = NaN;
    newBrain_wmean(newBrain_wmean > 0) = 0;
    for visualRegion = 1:7
        curOurBrain = ourBrain;
        % if visualRegion == 2
        %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
        % elseif visualRegion == 3
        %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
        % end
        thisfield = combinedRoiNames{visualRegion};
        newBrain_wmean(curOurBrain == visualRegion) = wmeanCoefCurv.(thisfield)(1,:);
    end

    for visualRegion = 1:7
        curOurBrain = ourBrain;
        % if visualRegion == 2
        %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
        % elseif visualRegion == 3
        %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
        % end
        curNewBrain = curOurBrain;
        curNewBrain(curOurBrain ~= visualRegion) = NaN;
        thisfield = combinedRoiNames{visualRegion};

        curNewBrain(curOurBrain == visualRegion) = wmeanCoefCurv.(thisfield)(1,:);

        newBrainbyROI_wmean(:,:,:,visualRegion) =curNewBrain;
    end


    % R2
    r2Brain = ourBrain;
    r2Brain(ourBrain <= 0) = NaN;
    r2Brain(ourBrain > 0) = 0;
    for visualRegion = 1:7
        thisfield = combinedRoiNames{visualRegion};
        r2Brain(ourBrain == visualRegion) = roiNsdCurvR2.(thisfield)(1,:);
    end

    % only R2 >0 voxels
    r2Mask = r2Brain > 0;
    newBrain_max_R2pos        = newBrain_max;
    newBrain_max_R2pos(~r2Mask) = NaN;

    newBrain_wmean_R2pos        = newBrain_wmean;
    newBrain_wmean_R2pos(~r2Mask) = NaN;

    newBrainbyROI_max_R2pos = newBrainbyROI_max;
    for v = 1:size(newBrainbyROI_max,4)
        tmp = newBrainbyROI_max(:,:,:,v);
        tmp(~r2Mask) = NaN;
        newBrainbyROI_max_R2pos(:,:,:,v) = tmp;
    end

    newBrainbyROI_wmean_R2pos = newBrainbyROI_wmean;
    for v = 1:size(newBrainbyROI_wmean,4)
        tmp = newBrainbyROI_wmean(:,:,:,v);
        tmp(~r2Mask) = NaN;
        newBrainbyROI_wmean_R2pos(:,:,:,v) = tmp;
    end





    % save nifti
    % reference NIfTI (orientation / spatial metadata)
    refFile = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/brainVolume/betas_session01_sub', ...
        num2str(isub), '.nii.gz'];
    info_ref = niftiinfo(refFile);

    % helper: copy spatial fields from reference
    copySpatialInfo = @(info) setfield(info, 'PixelDimensions', [1.8 1.8 1.8], ...
        'TransformName', info_ref.TransformName, ...
        'SpatialDimension', info_ref.SpatialDimension, ...
        'Transform', info_ref.Transform, ...
        'Qfactor', info_ref.Qfactor, ...
        'AuxiliaryFile', info_ref.AuxiliaryFile, ...
        'raw', info_ref.raw);

    % outputs to write: {data, filename}
    outputs = {
        newBrain_max,          'curvMLVBrain_max',        false;
        newBrainbyROI_max,     'curvMLVBrainbyROI_max',   true;
        newBrain_wmean,        'curvMLVBrain_wmean',      false;
        newBrainbyROI_wmean,   'curvMLVBrainbyROI_wmean', true;
        r2Brain,               'curvMLVBrain_R2',         false;

        % R2 > 0 masked versions
        newBrain_max_R2pos,           'newBrain_max',                      false;
        newBrainbyROI_max_R2pos,      'newBrainbyROI_max',                 true;
        newBrain_wmean_R2pos,         'newBrain_wmean',                    false;
        newBrainbyROI_wmean_R2pos,    'newBrainbyROI_wmean',               true
        };

    for i = 1:size(outputs,1)
        data = outputs{i,1};
        fname = fullfile(save_brainfolder, ...
            [outputs{i,2}, '_sub', num2str(isub), '.nii']);
        useRef = outputs{i,3};

        % write once to create header
        niftiwrite(data, fname);

        % read header and overwrite spatial metadata
        info_new = niftiinfo(fname);

        % pixel dimensions (conditional)
        if useRef
            info_new.PixelDimensions = info_ref.PixelDimensions;
        else
            info_new.PixelDimensions = [1.8 1.8 1.8];
        end
        info_new.TransformName     = info_ref.TransformName;
        info_new.SpatialDimension  = info_ref.SpatialDimension;
        info_new.Transform         = info_ref.Transform;
        info_new.Qfactor           = info_ref.Qfactor;
        info_new.AuxiliaryFile     = info_ref.AuxiliaryFile;
        info_new.raw.pixdim        = info_ref.raw.pixdim;
        info_new.raw.aux_file      = info_ref.raw.aux_file;
        info_new.raw.sform_code    = info_ref.raw.sform_code;
        info_new.raw.srow_x        = info_ref.raw.srow_x;
        info_new.raw.srow_y        = info_ref.raw.srow_y;
        info_new.raw.srow_z        = info_ref.raw.srow_z;

        % rewrite with corrected header
        niftiwrite(data, fname, info_new);
    end

end