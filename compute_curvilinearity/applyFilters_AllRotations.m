function [output, recMap]=applyFilters_AllRotations(Img,type,SF,theta,resize)

%constant filter parameters
% mA=4;%control sigma
% xLength=300;
% yLength=300;
mA=4;%control sigma
xLength=30;
yLength=30;

type = lower(type);

if size(Img,3)>1
    Img = mean(Img,3);
end


%resize if necessary
 if resize==1
    Img=imresize(Img,.5,'bilinear');
 end

%find edges
if ~islogical(Img)
    Img=edge(Img,'canny',.2);
end


%parameters for invariace
rotationList=0:pi/8:(2*pi)-pi/8;
accumulator=zeros([size(Img) length(rotationList)]);
waveLength=1/SF;

for rot=1:length(rotationList)
    alpha=rotationList(rot);

    if strcmp(type,'angle')

        [SpaceKernel, FreKernel, SpaceKernel_realsize] = AngleFilter(waveLength, alpha, theta, mA, xLength, yLength);
    elseif strcmp(type,'curve')
        [SpaceKernel, FreKernel, SpaceKernel_realsize] = CurveFilter(waveLength, alpha, theta, mA, xLength, yLength);
    end
    filter=real(SpaceKernel);

    %apply filter
    filteredImg=conv2(double(Img),filter,'same');
    filteredImg=abs(filteredImg./sum(sum(abs(filter))));
    
    accumulator(:,:,rot)=(abs(filteredImg).^8);
end

accumulator=squeeze(mean(accumulator,3));
recMap = NaN(size(Img));
recMap(Img == 1) = accumulator(Img == 1);
output=mean(accumulator(Img==1));
