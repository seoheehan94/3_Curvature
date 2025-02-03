clear all
close all

SFlist=[1 2 4 8];

% thetaList=[0 pi/6 pi/3 pi/2 2*pi/3 5*pi/6];
thetaList = linspace(0,157.5,8);
thetaList = thetaList/1200;
thetaList(3)=thetaList(3)*0.9;
thetaList(4)=thetaList(4)*0.9;
thetaList(6)=thetaList(6)*1.2;
thetaList(7)=thetaList(7)*2.4;
thetaList(8)=thetaList(8)*5;

picFolder = '/bwlab/Users/SeoheeHan/curvatureProject/IndoorArchitecture/images/';
picList=dir([picFolder, '/*']); 
picList = picList(~ismember({picList(:).name},{'.','..'}));
nPics=length(picList);

outputMatrix=zeros(nPics,length(SFlist),length(thetaList));

% vecLD = importSVG(['images/', picList(1).name]);
% imgLD = renderLinedrawing(newVecLD);
% imgLD = squeeze(imgLD(:,:,1)); 
% Img = 1-imgLD;
for pic=1:nPics
    pic
    ImgName=[picFolder, picList(pic).name];
    Img = im2double(imread(ImgName));
    % Img = double(imread(ImgName));
    % Img = im2double(imread(ImgName));
    % Img = 1-im2double(imread(ImgName));

    imsize = size(Img);
    recMap=NaN(imsize(1),imsize(2), length(SFlist), length(thetaList));

    for s=1:length(SFlist)
        SF=SFlist(s);
        for t=1:length(thetaList)
            theta=thetaList(t);
            [outputMatrix(pic,s,t), recMap(:,:,s,t)] =applyAngleFilters_AllRotations(Img,SF,theta,0);
        end
    end
    recMap_mean=squeeze(mean(recMap,3));
    recMap_max = squeeze(max(recMap,[],3));
end

%normalize
normData=zeros(size(outputMatrix));
for t=1:length(thetaList)
    for p=1:nPics
        for sf=1:length(SFlist)
            temp=outputMatrix(:,sf,t);
            normData(:,sf,t)=(temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
        end
    end
end

finalData=squeeze(mean(normData,2));
finalData=mean(finalData,2);

picNames = {picList(:).name}';
totalData = table(picNames, finalData, 'VariableNames', {'imgName', 'meanCurv'});

save('/bwlab/Users/SeoheeHan/curvatureProject/IndoorArchitecture/curvData.mat', 'totalData');
% load('/bwlab/Users/SeoheeHan/curvatureProject/IndoorArchitecture/curvData.mat')
%rank images and output
% [val,ind]=sort(finalData(:,4),'ascend');


