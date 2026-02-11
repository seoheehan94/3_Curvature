
clear all;
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_filter/prfsample_Curv/';
save_brainfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_filter/brainVolume';
save_analysisfolder = '/home/hanseohe/Documents/GitHub/3_Curvature/NSD/analysis/';


%% 1.  get R2 %%

totalR2CurvSplit = cell(1, 8);
totalMaxCurv     = cell(1, 8);

for isub = 1:8
    fprintf('isub:%d...\n', isub);
    load([prffolder 'voxCurvCoef_sub' num2str(isub) '.mat']); % should load roiNsdCurvR2
    totalR2CurvSplit{isub} = roiNsdCurvR2;
    totalMaxCurv{isub} = maxCoefCurv;
end

% Preallocate containers for long-format data
subj_all = [];
roi_all  = {};
R2_all   = [];
curv_all = [];

for isub = 1:8
    roiStruct = totalR2CurvSplit{isub};
    roiStructCurv = totalMaxCurv{isub};

    for r = 1:numel(combinedRoiNames)
        roiName = combinedRoiNames{r};
        R2_vals = roiStruct.(roiName);  % vector of R2 for this subj & ROI
        curv_vals = roiStructCurv.(roiName);
        
        nvox = numel(R2_vals);
        nvoxCurv = numel(curv_vals);

        if nvox ~= nvoxCurv
            warning('Voxel count mismatch: subj %d, ROI %s: nR2=%d, nCurv=%d', ...
                isub, roiName, nvox, nvoxCurv);
        end


        subj_all = [subj_all;  repmat(isub, nvox, 1)];
        roi_all  = [roi_all;   repmat({roiName}, nvox, 1)];
        R2_all   = [R2_all;    R2_vals(:)];
        curv_all   = [curv_all;    curv_vals(:)];
    end
end

% Build table: subj, ROI, R2, curv
T_all = table(subj_all, roi_all, R2_all, curv_all,...
             'VariableNames', {'subj', 'ROI', 'R2', 'curvPref'});

% Save one "all ROI" file for MLV
writetable(T_all, fullfile(save_analysisfolder, 'allROI_filter.csv'));




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

        curNewBrain(curOurBrain == visualRegion) = roiNsdCurvR2.(thisfield)(1,:);

        r2BrainbyROI(:,:,:,visualRegion) =curNewBrain;
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
        'AuxiliaryFile', info_ref.AuxiliaryFile);

    % outputs to write: {data, filename}
    outputs = {
        newBrain_max,          'curvfilterBrain_max',        false;
        newBrainbyROI_max,     'curvfilterBrainbyROI_max',   true;
        newBrain_wmean,        'curvfilterBrain_wmean',      false;
        newBrainbyROI_wmean,   'curvfilterBrainbyROI_wmean', true;

        r2Brain,               'curvfilterBrain_R2',         false;
        r2BrainbyROI,               'curvfilterBrain_R2byROI',         true;

        % R2 > 0 masked versions
        newBrain_max_R2pos,           'curvfilterBrain_max_posR2',                      false;
        newBrainbyROI_max_R2pos,      'curvfilterBrainbyROI_max_posR2',                 true;
        newBrain_wmean_R2pos,         'curvfilterBrain_wmean_posR2',                    false;
        newBrainbyROI_wmean_R2pos,    'curvfilterBrainbyROI_wmean_posR2',               true
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