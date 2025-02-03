function [outputMatrix, recMap_mean,recMap_max]=applyFilters(ImgName, type)

SFlist=[1 2 4 8];

thetaList = linspace(0,157.5,8);
thetaList = thetaList/1200;
thetaList(3)=thetaList(3)*0.9;
thetaList(4)=thetaList(4)*0.9;
thetaList(6)=thetaList(6)*1.2;
thetaList(7)=thetaList(7)*2.4;
thetaList(8)=thetaList(8)*5;

outputMatrix=zeros(length(SFlist),length(thetaList));

% vecLD = importSVG(['images/', picList(1).name]);
% imgLD = renderLinedrawing(newVecLD);
% imgLD = squeeze(imgLD(:,:,1));
% Img = 1-imgLD;

if ~islogical(ImgName)
    Img = im2double(ImgName);
    % Img = double(imread(ImgName));
    % Img = 1-im2double(imread(ImgName));
else
    Img = ImgName;
end

imsize = size(Img);
recMap=NaN(imsize(1),imsize(2), length(SFlist), length(thetaList));

for s=1:length(SFlist)
    SF=SFlist(s);
    for t=1:length(thetaList)
        theta=thetaList(t);
        [outputMatrix(s,t), recMap(:,:,s,t)] =applyFilters_AllRotations(Img,type, SF,theta,0);
    end
end
recMap_mean=squeeze(mean(recMap,3));
recMap_max = squeeze(max(recMap,[],3));


