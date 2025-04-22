function curvMap = generateCurvatureMap(Img, prefAnalysis)

imgSize = size(Img);
thisMap = zeros(imgSize(1), imgSize(2));

SFlist=[1 2 4 8];

thetaList = linspace(0,157.5,8);
thetaList = thetaList/1200;
thetaList(3)=thetaList(3)*0.9;
thetaList(4)=thetaList(4)*0.9;
thetaList(6)=thetaList(6)*1.2;
thetaList(7)=thetaList(7)*2.4;
thetaList(8)=thetaList(8)*5;

imsize = size(Img);
thisMap=NaN(imsize(1),imsize(2), length(SFlist), length(thetaList));

for s=1:length(SFlist)
    SF=SFlist(s);
    for t=1:length(thetaList)
        theta=thetaList(t);
        [~, thisMap(:,:,s,t)] =applyAngleFilters_AllRotations(Img,'curve',SF,theta,0);
    end
end

switch (lower(prefAnalysis))
    case 'mean'
        curvMap=squeeze(mean(thisMap,3));

    case 'max'
        curvMap = squeeze(max(thisMap,[],3));
       
        %softmax
end

curvMap = curvMap.^(1/8);


