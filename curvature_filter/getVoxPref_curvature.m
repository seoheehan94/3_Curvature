% getVoxPref_curvature.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: getVoxPref(1,1)
%   by: zvi roth
%   date: 7/29/2022
%   purpose: extract preferred orientation from regression weights, and sum
%   across partitions
%   uses files created by: regressPrfSplit.m
%   creates files used by: fig$$.m


function getVoxPref_curvature(isub,numregions)

%uses data from regressPrfSplit.m

%chooses preferred orientation/levels by averaging across filters
%mrQuit
close all
% clear all
global interpSz;
global backgroundSz;
global degPerPix;
global prefAnalysis;

prefAnalysis = 1;

toSavePdf = 0;

numCurvs = 8;
figFolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/figures/'];
nperms=1000;
prffolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/prfsample_Len/'];

interpSz= 714;
backgroundSz= 1024;
bandpass = 1; bandMin = 1; bandMax = 1;
bandpassStr = '';
if bandpass
    bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
end

imgScaling = 0.5;
interpSz= 714*imgScaling;
backgroundSz= 1024*imgScaling;
degPerPix = 8.4/(714*imgScaling);

gratings = load(['gratings_length.mat'],'numCurvs','numLevels',...
    'sumLenEnergy','modelLenEnergy');
for iregion=1:numregions
    visualRegion = iregion;%V1,V2,V3,V4,OPA,PPA,RSC
    load(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
        'nsd',...
        'numLevels', 'numCurvs','rois','nvox','roiPrf','nsplits');
    
    if length(rois)>1 %combine across ventral and dorsal ROIs
        oldNsd = nsd;
        nsd.voxResidual{1} = [];
        nsd.voxLenResidual{1} = [];
        nsd.voxResidualSplit{1} = [];
        nsd.voxLenResidualSplit{1} = [];
        nsd.r2{1} = [];
        nsd.r2len{1} = [];
        nsd.r2split{1} = [];
        nsd.r2lenSplit{1} = [];
        nsd.voxCoef{1} = [];
        nsd.voxLenCoef{1} = [];
        nsd.voxPredLenCoef{1} = [];
        nsd.voxLenPredLenCoef{1} = [];
        nsd.voxResidLenCoef{1} = [];
        nsd.voxLenResidLenCoef{1} = [];
        
        nsd.voxPredLenR2{1} = [];
        nsd.voxLenPredLenR2{1} = [];
        nsd.voxResidLenR2{1} = [];
        nsd.voxLenResidLenR2{1} = [];
        
        nsd.pearsonRlen{1} = [];
        nsd.pearsonR{1} = [];
        
        nsd.roiInd{1} = [];
        
%         oldSynth = synth;
%         synth.pearsonR{1} = [];
%         synth.pearsonRori{1} = [];
%         synth.r2ori{1} = [];
%         synth.voxResidual{1} = [];
%         synth.voxOriResidual{1} = [];
%         synth.voxCoef{1} = [];
%         synth.voxOriCoef{1} = [];
        
        for iroi=1:length(rois)
            nsd.voxResidual{1} = cat(2,nsd.voxResidual{1},oldNsd.voxResidual{iroi});
            nsd.voxLenResidual{1} = cat(2,nsd.voxLenResidual{1},oldNsd.voxLenResidual{iroi});
            nsd.voxResidualSplit{1} = cat(2,nsd.voxResidualSplit{1},oldNsd.voxLenResidual{iroi});
            nsd.voxLenResidualSplit{1} = cat(2,nsd.voxLenResidualSplit{1},oldNsd.voxLenResidual{iroi});
            
            nsd.pearsonRlen{1} = cat(2,nsd.pearsonRlen{1},oldNsd.pearsonRlen{iroi});
            nsd.pearsonR{1} = cat(2,nsd.pearsonR{1},oldNsd.pearsonR{iroi});
            nsd.r2{1} = cat(2,nsd.r2{1},oldNsd.r2{iroi});
            nsd.r2len{1} = cat(2,nsd.r2len{1},oldNsd.r2len{iroi});
            nsd.r2split{1} = cat(2,nsd.r2split{1},oldNsd.r2split{iroi});
            nsd.r2lenSplit{1} = cat(2,nsd.r2lenSplit{1},oldNsd.r2lenSplit{iroi});
            nsd.voxCoef{1} = cat(2,nsd.voxCoef{1},oldNsd.voxCoef{iroi});
            nsd.voxLenCoef{1} = cat(2,nsd.voxLenCoef{1},oldNsd.voxLenCoef{iroi});
            nsd.voxPredLenCoef{1} = cat(2,nsd.voxPredLenCoef{1},oldNsd.voxPredLenCoef{iroi});
            nsd.voxLenPredLenCoef{1} = cat(2,nsd.voxLenPredLenCoef{1},oldNsd.voxLenPredLenCoef{iroi});
            nsd.voxResidLenCoef{1} = cat(2,nsd.voxResidLenCoef{1},oldNsd.voxResidLenCoef{iroi});
            nsd.voxLenResidLenCoef{1} = cat(2,nsd.voxLenResidLenCoef{1},oldNsd.voxLenResidLenCoef{iroi});
            
            nsd.voxPredLenR2{1} = cat(2,nsd.voxPredLenR2{1},oldNsd.voxPredLenR2{iroi});
            nsd.voxLenPredLenR2{1} = cat(2,nsd.voxLenPredLenR2{1},oldNsd.voxLenPredLenR2{iroi});
            nsd.voxResidLenR2{1} = cat(2,nsd.voxResidLenR2{1},oldNsd.voxResidLenR2{iroi});
            nsd.voxLenResidLenR2{1} = cat(2,nsd.voxLenResidLenR2{1},oldNsd.voxLenResidLenR2{iroi});
            
            nsd.roiInd{1} = cat(1,nsd.roiInd{1}, oldNsd.roiInd{iroi});
            
%             synth.voxResidual{1} = cat(2,synth.voxResidual{1},oldSynth.voxResidual{iroi});
%             synth.voxOriResidual{1} = cat(2,synth.voxOriResidual{1},oldSynth.voxOriResidual{iroi});
%             synth.pearsonRori{1} = cat(2,synth.pearsonRori{1},oldSynth.pearsonRori{iroi});
%             synth.pearsonR{1} = cat(2,synth.pearsonR{1},oldSynth.pearsonR{iroi});
%             synth.voxCoef{1} = cat(1,synth.voxCoef{1},oldSynth.voxCoef{iroi});
%             synth.voxOriCoef{1} = cat(1,synth.voxOriCoef{1},oldSynth.voxOriCoef{iroi});
%             
        end
        oldPrf = roiPrf; clear roiPrf;
        roiPrf{1}.ecc=[];
        roiPrf{1}.ang=[];
        roiPrf{1}.sz=[];
        %         roiPrf{1}.exponent=[];
        %         roiPrf{1}.gain=[];
        roiPrf{1}.r2=[];
        roiPrf{1}.x=[];
        roiPrf{1}.y=[];
        for iroi=1:length(rois)
            roiPrf{1}.ecc = cat(1,roiPrf{1}.ecc,oldPrf{iroi}.ecc);
            roiPrf{1}.ang = cat(1,roiPrf{1}.ang,oldPrf{iroi}.ang);
            roiPrf{1}.sz = cat(1,roiPrf{1}.sz,oldPrf{iroi}.sz);
            %             roiPrf{1}.exponent = cat(1,roiPrf{1}.exponent,oldPrf{iroi}.exponent);
            %             roiPrf{1}.gain = cat(1,roiPrf{1}.gain,oldPrf{iroi}.gain);
            roiPrf{1}.r2 = cat(1,roiPrf{1}.r2,oldPrf{iroi}.r2);
            roiPrf{1}.x = cat(1,roiPrf{1}.x,oldPrf{iroi}.x);
            roiPrf{1}.y = cat(1,roiPrf{1}.y,oldPrf{iroi}.y);
        end
        rois = 1;
    end
    
    %% AVERAGE SPLITS
    nsd.voxCoef{1}(nsplits+1,:,:) = mean(nsd.voxCoef{1},1);
    nsd.voxLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxLenCoef{1},1);
    nsd.voxPredLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxPredLenCoef{1},1);
    nsd.voxLenPredLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxLenPredLenCoef{1},1);
    nsd.voxResidLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxResidLenCoef{1},1);
    nsd.voxLenResidLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxLenResidLenCoef{1},1);
    
    nsd.pearsonRlen{1}(nsplits+1,:) = mean(nsd.pearsonRlen{1},1);
    nsd.pearsonR{1}(nsplits+1,:) = mean(nsd.pearsonR{1},1);
    nsd.r2{1}(nsplits+1,:) = mean(nsd.r2{1},1);
    nsd.r2len{1}(nsplits+1,:) = mean(nsd.r2len{1},1);
    nsd.r2split{1}(nsplits+1,:) = mean(nsd.r2split{1},1);
    nsd.r2lenSplit{1}(nsplits+1,:) = mean(nsd.r2lenSplit{1},1);
    
%     synth.voxResidual{1}(nsplits+1,:,:) = mean(synth.voxResidual{1},1);
%     synth.voxOriResidual{1}(nsplits+1,:,:) = mean(synth.voxOriResidual{1},1);
%     synth.pearsonRori{1}(nsplits+1,:) = mean(synth.pearsonRori{1},1);
%     synth.pearsonR{1}(nsplits+1,:) = mean(synth.pearsonR{1},1);
    
    nsplits = nsplits+1;
    
    %% COMPUTE PREFERRED LEVEL - CENTER OF MASS
    clear vigPrefLevel fullPrefLevel residPrefLevel residLenPrefLevel predCoefLevel predLenCoefLevel
    for  iroi=1:length(rois)
        for isplit=1:nsplits
            %constrained model
            vigCoef = squeeze(nsd.voxCoef{iroi}(isplit,:,1:end-1));
            vigCoef = vigCoef';
            vigPrefLevel{iroi}(isplit,:) = gratingPrefFreq(vigCoef,gratings.sumLenEnergy);
            
            %full orientation model
            fullCoef = squeeze(nsd.voxLenCoef{iroi}(isplit,:,1:end-1));
            fullPrefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelLenEnergy);
            
            %full model on residuals of constrained model
            fullCoef = squeeze(nsd.voxResidLenCoef{iroi}(isplit,:,1:end-1));
            residPrefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelLenEnergy);
            
            %full model on residual of full model
            fullCoef = squeeze(nsd.voxLenResidLenCoef{iroi}(isplit,:,1:end-1));
            residLenPrefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelLenEnergy);
            
            %full model on constrained prediction
            fullCoef = squeeze(nsd.voxPredLenCoef{iroi}(isplit,:,1:end-1));
            predCoefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelLenEnergy);
            
            %full model on full model prediction
            fullCoef = squeeze(nsd.voxLenPredLenCoef{iroi}(isplit,:,1:end-1));
            predLenCoefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelLenEnergy);
            
        end
        %%%%%%%%%%%%%% SYNTH
        %constrained model
%         synthVigCoef = squeeze(synth.voxCoef{iroi}(:,1:end-1));
%         synthVigPrefLevel{iroi} = gratingPrefFreq(synthVigCoef,gratings.sumOriEnergy);
%         
%         %full orientation model
%         synthFullCoef = squeeze(synth.voxOriCoef{iroi}(:,1:end-1));
%         synthFullPrefLevel{iroi} = gratingPrefFreq(synthFullCoef,gratings.modelOriEnergy);
    end
    
    
    %% COMPUTE PREFERRED ORIENTATION - CIRCULAR CENTER OF MASS
    
    clear fullPrefLen residPrefLen residLenPrefLen predCoefLen predLenCoefLen lenDeviation vertDeviation cardDeviation
    clear fullLenModul fullPrefAmp fullAntiAmp
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            %full orientation model
            fullCoef = squeeze(nsd.voxLenCoef{iroi}(isplit,:,1:end-1));
            fullPrefLen{iroi}(isplit,:) = gratingPrefLen(fullCoef,gratings.modelLenEnergy);
            
            %full model on residuals of constrained model
            fullCoef = squeeze(nsd.voxResidLenCoef{iroi}(isplit,:,1:end-1));
            residPrefLen{iroi}(isplit,:) = gratingPrefLen(fullCoef,gratings.modelLenEnergy);
            
            %full model on residual of full model
            fullCoef = squeeze(nsd.voxLenResidLenCoef{iroi}(isplit,:,1:end-1));
            residLenPrefLen{iroi}(isplit,:) = gratingPrefLen(fullCoef,gratings.modelLenEnergy);
            
            %full model on constrained prediction
            fullCoef = squeeze(nsd.voxPredLenCoef{iroi}(isplit,:,1:end-1));
            predCoefLen{iroi}(isplit,:) = gratingPrefLen(fullCoef,gratings.modelLenEnergy);
            
            %full model on full model prediction
            fullCoef = squeeze(nsd.voxLenPredLenCoef{iroi}(isplit,:,1:end-1));
            predLenCoefLen{iroi}(isplit,:) = gratingPrefLen(fullCoef,gratings.modelLenEnergy);
            
        end
        %cross-validated measure of orientation modulation
        for isplit=1:2
            fullCoef = squeeze(nsd.voxLenCoef{iroi}(isplit,:,1:end-1));
            [fullLenModul{iroi}(isplit,:), fullPrefAmp{iroi}(isplit,:), fullAntiAmp{iroi}(isplit,:)] = gratingLenModulation(fullPrefLen{iroi}(3-isplit,:),fullCoef,gratings.modelLenEnergy);
        end
        fullLenModul{iroi}(nsplits,:) = mean(fullLenModul{iroi}(1:2,:),1);
        fullPrefAmp{iroi}(nsplits,:) = mean(fullPrefAmp{iroi}(1:2,:),1);
        fullAntiAmp{iroi}(nsplits,:) = mean(fullAntiAmp{iroi}(1:2,:),1);
        
        %%%%%%%%%%%%% SYNTH
        %full orientation model
%         synthFullCoef = squeeze(synth.voxLenCoef{iroi}(:,1:end-1));
%         synthFullPrefLen{iroi} = gratingPrefLen(synthFullCoef,gratings.modelLenEnergy);
    end
    
    %% SCATTER PLOT ORIENTATION DEVIATION FROM RADIAL VS ECCENTRICITY
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            %multiply by 2 for circular dist, then divide by 2
            lenDeviation{iroi}(isplit,:) = 0.5*circ_dist(2*(pi/2-fullPrefLen{iroi}(isplit,:)),2*(mod(roiPrf{iroi}.ang*pi/180,pi)'));
        end
    end
    %% SCATTER PLOT ORIENTATION DEVIATION FROM VERTICAL VS ECCENTRICITY
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            vertDeviation{iroi}(isplit,:) = 0.5*circ_dist(2*fullPrefLen{iroi}(isplit,:),0);
        end
    end
    
    %% SCATTER PLOT ORIENTATION DEVIATION FROM CARDINAL VS ECCENTRICITY
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            vertAngle1 = roiPrf{iroi}.ang*pi/180>=pi/4 & roiPrf{iroi}.ang*pi/180<3*pi/4;
            vertAngle2 = roiPrf{iroi}.ang*pi/180>=5*pi/4 & roiPrf{iroi}.ang*pi/180<7*pi/4;
            vertAngle = vertAngle1 | vertAngle2;
            horizAngle1 = roiPrf{iroi}.ang*pi/180<pi/4;
            horizAngle2 = roiPrf{iroi}.ang*pi/180>=3*pi/4 & roiPrf{iroi}.ang*pi/180<5*pi/4;
            horizAngle3 = roiPrf{iroi}.ang*pi/180>=7*pi/4;
            horizAngle = horizAngle1 | horizAngle2 | horizAngle3;
            %multiply by 2 for circular dist to work on pi cycle instead of 2pi.
            cardDeviation{iroi}(isplit,vertAngle) = 0.5*circ_dist(2*fullPrefLen{iroi}(isplit,vertAngle),0);
            cardDeviation{iroi}(isplit,horizAngle) = 0.5*circ_dist(2*fullPrefLen{iroi}(isplit,horizAngle),pi);%using pi instead of pi/2
        end
    end
    
    %% Distribution of cross-validated SYNTHETIC  R^2 across voxels
%     for  iroi=1:length(rois)%rois=1
%         for isplit=1:nsplits
%             [nsdSynthImprov_corr(iregion,isplit), nsdSynthImprov_pval(iregion,isplit)] = corr((nsd.pearsonRori{iroi}(isplit,:)-nsd.pearsonR{iroi}(isplit,:))', (synth.pearsonRori{iroi}(isplit,:) - synth.pearsonR{iroi}(isplit,:))', 'rows','complete');
%         end
%     end
    
    %%
    %save preferred orientation and level for this ROI
    allRoiPrf{iregion} = roiPrf{iroi};%iroi=1
    roiLevVig{iregion} = vigPrefLevel{iroi};
    roiLevFull{iregion} = fullPrefLevel{iroi};
    roiLen{iregion} = fullPrefLen{iroi};
    residLen{iregion} = residPrefLen{iroi};
    residLenLen{iregion} = residLenPrefLen{iroi};
    predLen{iregion} = predCoefLen{iroi};
    predLenLen{iregion} = predLenCoefLen{iroi};
    lenModulation{iregion} = fullLenModul{iroi};
    lenPrefAmp{iregion} = fullPrefAmp{iroi};
    lenAntiAmp{iregion} = fullAntiAmp{iroi};
    
    
    roiInd{iregion} = nsd.roiInd{iroi};
    roiNsdCorr{iregion} = nsd.pearsonR{iroi};
    roiNsdLenCorr{iregion} = nsd.pearsonRlen{iroi};
    roiNsdLenR2{iregion} = nsd.r2lenSplit{iroi};
    roiNsdR2{iregion} = nsd.r2split{iroi};
    
    roiNsdLenPredLenR2{iregion} = nsd.voxLenPredLenR2{iroi};
    roiNsdLenResidLenR2{iregion} = nsd.voxLenResidLenR2{iroi};
    roiNsdPredLenR2{iregion} = nsd.voxPredLenR2{iroi};
    roiNsdResidLenR2{iregion} = nsd.voxResidLenR2{iroi};
    
%     roiSynthCorr{iregion} = synth.pearsonR{iroi};
%     roiSynthOriCorr{iregion} = synth.pearsonRori{iroi};
%     roiSynthOri{iregion} = synthFullPrefOri{iroi};
%     roiSynthLevVig{iregion} = synthVigPrefLevel{iroi};
%     roiSynthLevFull{iregion} = synthFullPrefLevel{iroi};
%     roiOriDeviation{iregion} = oriDeviation{iroi};
%     roiVertDeviation{iregion} = vertDeviation{iroi};
%     roiCardDeviation{iregion} = cardDeviation{iroi};
end


%save all ROIs to create overlay
roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
combinedRoiNames = {'V1','V2','V3','hV4'};

save([prffolder 'voxModelPref_sub' num2str(isub) '.mat'],'allRoiPrf','roiLevVig','roiLevFull',...
    'roiLen','roiNsdCorr','roiNsdLenCorr','roiNsdLenR2','roiNsdR2',...
    'roiNsdResidLenR2','roiNsdLenResidLenR2','roiNsdPredLenR2','roiNsdLenPredLenR2',...
    'residLen','residLenLen','predLen','predLenLen',...
    'lenModulation','lenPrefAmp','lenAntiAmp',...
    'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis','nsplits');

%%
    function prefAngle = gratingPrefLen(fullCoef,modelLenEnergy);
        [numFreqs, numAngles, numLevels, numCurvs] = size(modelLenEnergy);
        numvox = size(fullCoef,1);
        coefMat = reshape(fullCoef,numvox,numLevels,numCurvs);%vox x levels x orientations
        coefVec = reshape(fullCoef,numvox,numLevels*numCurvs);%vox x coefficients
        modelEnergy = reshape(modelLenEnergy,numFreqs*numAngles,numLevels*numCurvs);
        modelEnergy = modelEnergy';%(ilev*orientation, isf*iangle)
        
        voxGratingResp = coefVec*modelEnergy;%vox,isf*iangle
        angles = linspace(0,180,numAngles+1);
        angles = angles(1:numAngles);
        switch prefAnalysis
            case 1
                % 1 - simple max
                [m maxInd] = max(voxGratingResp,[],2);%vox
                [prefFreqNum, prefAngleNum] = ind2sub([numFreqs, numAngles],maxInd);%vox
                prefAngle = prefAngleNum;
            case 2
                % 2 - average across gratings spatial frequency, and then max
                voxAngleResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),2));
                [m prefAngleNum] = max(voxAngleResp,[],2);%ivox
                keyboard;
                prefAngle = prefAngleNum;
            case 3
                % 3 - average across gratings spatial frequency, and then circular weighted mean
                voxAngleResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),2));
                %subtract minimum response (might be negative)
                voxAngleResp = voxAngleResp - min(voxAngleResp,[],2);
                theta = linspace(0,2*pi,numAngles+1);%for circular calculation
                theta = theta(1:end-1);
                keyboard;
                for ivox=1:numvox
                    prefAngle(ivox) = circ_mean(theta',voxAngleResp(ivox,:)');
                end
                
                % prefAngle = mod(prefAngle,2*pi);%from [-pi, pi] to [0 2pi]
                % prefAngle = prefAngle./2;%range 0 to pi.
        end
    end

%%
    function [lenModul, prefLenAmp, antiPrefLenAmp, prefLenInd] = gratingLenModulation(prefLen,fullCoef,modelLenEnergy);
        [numFreqs, numAngles, numLevels, numCurvs] = size(modelLenEnergy);
        numvox = size(fullCoef,1);
        coefMat = reshape(fullCoef,numvox,numLevels,numCurvs);%vox x levels x orientations
        coefVec = reshape(fullCoef,numvox,numLevels*numCurvs);%vox x coefficients
        modelEnergy = reshape(modelLenEnergy,numFreqs*numAngles,numLevels*numCurvs);
        modelEnergy = modelEnergy';%(ilev*orientation, isf*iangle)
        
        voxGratingResp = coefVec*modelEnergy;%vox,isf*iangle
        angles = linspace(0,180,numAngles+1);
        angles = angles(1:numAngles);
        prefLenInd = zeros(numvox,1);
        minDist = zeros(numvox,1);
        prefLenAmp = zeros(numvox,1);
        antiPrefLenAmp = zeros(numvox,1);
        
        randorder = randperm(numvox);
        for ivox=1:numvox
            [minDist(ivox), prefLenInd(ivox)] = min(abs(circ_dist(angles*pi/180,prefLen(ivox))));%multiplied by pi/180 to be circular around 2pi
        end
        antiPrefLenInd = 1+mod(prefLenInd - 1 - length(angles)/2,length(angles));
        
        %how to deal with frequency: 1 - average across frequencies. 2 -
        %use preferred frequency
        
        voxGratingResp = reshape(voxGratingResp, numvox,numFreqs,numAngles);
        for ivox=1:numvox
            prefLenAmp(ivox) = mean(voxGratingResp(ivox,:,prefLenInd(ivox)),2);%mean across frequencies
            antiPrefLenAmp(ivox) = mean(voxGratingResp(ivox,:,antiPrefLenInd(ivox)),2);
        end
        lenModul = (prefLenAmp - antiPrefLenAmp)./(prefLenAmp + antiPrefLenAmp);
        
    end

%%
    function prefFreqNum = gratingPrefFreq(fullCoef,modelLenEnergy);
        %         global prefAnalysis
        [numFreqs, numAngles, numLevels, numCurvs] = size(modelLenEnergy);
        numvox = size(fullCoef,1);
        coefMat = reshape(fullCoef,numvox,numLevels,numCurvs);%vox x levels x orientations
        coefVec = reshape(fullCoef,numvox,numLevels*numCurvs);%vox x coefficients
        modelEnergy = reshape(modelLenEnergy,numFreqs*numAngles,numLevels*numCurvs);
        modelEnergy = modelEnergy';%(ilev*orientation, isf*iangle)
        
        voxGratingResp = coefVec*modelEnergy;%vox,isf*iangle
        switch prefAnalysis
            case 1
                % 1 - simple max
                [m maxInd] = max(voxGratingResp,[],2);%vox
                [prefFreqNum, prefAngleNum] = ind2sub([numFreqs, numAngles],maxInd);%vox
            case 2
                % 2 - average across gratings orientation, and then max
                voxFreqResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),3));
                [m prefFreqNum] = max(voxFreqResp,[],2);%ivox
            case 3
                % 3 - average across gratings orientation, and then weighted mean
                voxFreqResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),3));
                %subtract minimum response
                voxFreqResp = voxFreqResp - min(voxFreqResp,[],2);
                prefFreqNum = (voxFreqResp * [1:numFreqs]')./sum(voxFreqResp,2);
        end
    end
end