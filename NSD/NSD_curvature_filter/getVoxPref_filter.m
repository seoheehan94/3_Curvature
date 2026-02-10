%getVoxPref.m

%   uses files created by: regressPrfSplit_curvate.m
%   creates files used by:
clear all;
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature_filter/prfsample_Curv/';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};


%% 1.  Coef voxel preference %%
for isub = 1:8
    
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
    roiNsdCurvR2 = struct;
    wmeanCoefCurv = struct;



    %% load file
    for visualRegion = 1:7
        fprintf('%d. %d...\n',isub, visualRegion);
        thisfield = combinedRoiNames{visualRegion};
        
        load(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
            'nsd', 'numLevels', 'numCurvs','rois','nvox','roiPrf','nsplits');
        nsd.r2curvSplit = nsd.r2lenSplit;

        % get mean coef/R2
        if visualRegion == 4 || visualRegion == 5 || visualRegion == 6 || visualRegion == 7

            mean_split1.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(1,:,1:8),2));
            mean_split2.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(2,:,1:8),2));

            % average splits
            nsd.voxCurvCoef{1}(nsplits+1,:,:) = mean(nsd.voxCurvCoef{1},1);
            mean_all.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(3,:,1:8),2));
            allCurvCoef.(thisfield) = nsd.voxCurvCoef{1}(:,:,1:8);
            roiNsdCurvR2.(thisfield) = mean(nsd.r2curvSplit{1},1);

        else %for other regions, combine ventral and dorsal
            % combine ventral and dorsal
            oldNsd = nsd;
            nsd.voxCurvCoef{1} = [];
            nsd.voxCurvCoef{2} = [];
            nsd.r2curvSplit{1} = [];
            nsd.r2curvSplit{2} = [];
            for iroi=1:length(rois)
                nsd.voxCurvCoef{1} = cat(2,nsd.voxCurvCoef{1},oldNsd.voxCurvCoef{iroi});
                nsd.r2curvSplit{1} = cat(2,nsd.r2curvSplit{1},oldNsd.r2curvSplit{iroi});
            end

            mean_split1.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(1,:,1:8),2));
            mean_split2.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(2,:,1:8),2));

            % average splits
            nsd.voxCurvCoef{1}(nsplits+1,:,:) = mean(nsd.voxCurvCoef{1},1);
            mean_all.(thisfield) = squeeze(mean(nsd.voxCurvCoef{1}(3,:,1:8),2));
            allCurvCoef.(thisfield) = nsd.voxCurvCoef{1}(:,:,1:8);
            roiNsdCurvR2.(thisfield) = mean(nsd.r2curvSplit{1},1);
        end

        % max coef curvature preference per voxel
        [~, maxCoefCurv.(thisfield)] = max(allCurvCoef.(thisfield)(3,:,:), [], 3);
        %[~, nvox, ~] = size(allCurvCoef.(thisfield));
        % for ivox = 1: nvox
        %     [~, icurv] = max(allCurvCoef.(thisfield)(3,ivox,:));
        %     maxCoefCurv.(thisfield)(1,ivox) = icurv;
        % 
        % end


        % weighted mean curvature preference per voxel
        values = reshape(1:8, 1, 1, 8);  % matches 3rd dimension
        coef = allCurvCoef.(thisfield)(3,:,1:8);  % averaged split only
        coefShift = coef - min(coef, [], 3);
        num = sum(coefShift .* values, 3);
        den = sum(coefShift, 3);
        wmeanCoefCurv.(thisfield) = num ./ den;
        % handle flat or zero-weight voxels
        wmeanCoefCurv.(thisfield)(den == 0) = NaN;


    end

    saveName = [prffolder, 'voxCurvCoef_sub', num2str(isub), '.mat'];
    save(saveName, 'allCurvCoef', 'roiNsdCurvR2', 'mean_all', 'mean_split1', 'mean_split2','maxCoefCurv','wmeanCoefCurv');

end





