
function [SpaceKernel, FreKernel, SpaceKernel_realsize] = AngleFilter(kA, alphaA, theta, mA, xA, yA)

%%% the function is used to generate banana wavelet kernels.  The kernels
%%% can be used to filter a image to quantify curvatures.
%%% 
%%% usage:  [SpaceKernel, FreKernel, SpaceKernel_realsize] = bananakernel(kA, alphaA, bA, mA, xA, yA)
%%%
%%% kA:          length of the wave vector K
%%% alphaA:      direction of the wave vector
%%% bA:          bending value b
%%% mA:          magnitude value m
%%% xA:          x-coordinate x
%%% yA:          y-coordinate y
%%%
%%% for references:
%%% preFactorA:  pre-factor p
%%% DCPartRealA: real dc-part
%%% DCPartImagA: imaginary dc-part
%%% gaussPartA:  Gaussian part

%%% the function is written by Xiaomin Yue at 9/10/2012, based on Norbert
%%% Kruger's method
%%%
%%%

bA=1/(tan((pi/2)-(theta/2)));
sigmaXbend = 2;
sigmaYbend = 2;

[xA, yA] = meshgrid(-xA:xA-1,-yA:yA-1);

xRotL = cos(alphaA)*xA + sin(alphaA)*yA; 
yRotL = cos(alphaA)*yA - sin(alphaA)*xA;

%xRotBendL = xRotL + bA * (yRotL).^2;
xRotBendL = abs(xRotL) + bA * (yRotL);
yRotBendL = yRotL;

%correct for distortion
temp=(cos(kA*xRotBendL));


%1)find target wavelength

xRotLOrig = cos(0)*xA + sin(0)*yA; 
pre=(cos(kA*xRotLOrig));

throughCount=0;
i=1;
while throughCount<2
    i=i+1;
    if pre(1,i-1)>pre(1,i) && pre(1,i+1)> pre(1,i)
        throughCount=throughCount+1;
        if throughCount==1
            A=i;
        else
            B=i;
        end
    end
end
targetLambda=B-A;



%find effective wavelength by rotating to perpendicular direction and
%figuring out line thickness
alphaRot=(pi/2)-(theta/2);
xRotPerp = cos(alphaRot)*xA + sin(alphaRot)*yA; 
yRotPerp = cos(alphaRot)*yA - sin(alphaRot)*xA;
xRotBendPerp = abs(xRotPerp) + bA * (yRotPerp);
perp=(cos(kA*xRotBendPerp));


throughCount=0;
i=length(temp(:,end))-1;
%i=1;
while throughCount<2
    i=i-1;
  %  i=i+1;
    if perp(i-1,end)>perp(i,end) && perp(i+1,end)> perp(i,end)
    %if temp(i-1,1)>temp(i,1) && temp(1,i+1)> temp(1,i)
        throughCount=throughCount+1;
        if throughCount==1
            A=i;
        else
            B=i;
        end
    end
end

effectiveLambda=abs(A-B);

kA=kA*(effectiveLambda/targetLambda);


%gaussian mask

gaussMask=zeros(size(xRotL));
gsize = size(gaussMask);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
sigma=gsize(1)*10;
center=[gsize(1)/2 gsize(2)/2 ];
gaussMask = gaussC(R,C, sigma, center);


sigmaYbend=10;
%%% make the DC free
tmpgaussPartA = exp(-0.5*(kA)^2*((xRotBendL/sigmaXbend).^2 + (yRotBendL/(mA*sigmaYbend)).^2));
tmprealteilL  = 1*tmpgaussPartA.*(cos(kA*xRotBendL) - 0);
tmprealteilL = tmprealteilL.*gaussMask;
tmpimagteilL  = 1*tmpgaussPartA.*(sin(kA*xRotBendL) - 0);
tmpimagteilL = tmpimagteilL.*gaussMask;

numeratorRealL = sum(tmprealteilL(:));
numeratorImagL = sum(tmpimagteilL(:));
denominatorL   = sum(tmpgaussPartA(:));

DCPartRealA = numeratorRealL/denominatorL;
DCPartImagA = numeratorImagL/denominatorL;

DCValueAnalysis = exp(-0.5 * sigmaXbend * sigmaXbend);
if DCPartRealA < DCValueAnalysis
    DCPartRealA = DCValueAnalysis;
    DCPartImagA = 0;
end    

%%% generate a space kernel
preFactorA = kA^2;
gaussPartA = exp(-0.5*(kA)^2*((xRotBendL/sigmaXbend).^2 + (yRotBendL/(mA*sigmaYbend)).^2));
realteilL  = preFactorA*gaussPartA.*(cos(kA*xRotBendL) - DCPartRealA);
realteilL = realteilL.*gaussMask;
imagteilL  = preFactorA*gaussPartA.*(sin(kA*xRotBendL) - DCPartImagA);
imagteilL = imagteilL.*gaussMask;

% %%% set values to zeros outside of kernel.
% ZeroOutSideRadius2L = imcircle(ceil(2*4*sigmaYbend*mA*(1/kA)));
% KernelSize          = size(ZeroOutSideRadius2L);
% if (mod(KernelSize(1),2) ~= 0)
%     tmpZeros = zeros(KernelSize+1);
%     tmpZeros(1:end-1, 1:end-1) = ZeroOutSideRadius2L;
%     ZeroOutSideRadius2L = tmpZeros;
%     KernelSize = size(ZeroOutSideRadius2L);
% end    
% RadiusMask          = zeros(size(xA));
% RadiusMaskSize      = size(RadiusMask);
% xx = RadiusMaskSize(1)/2-(KernelSize(1)/2):RadiusMaskSize(1)/2+(KernelSize(1)/2)-1;
% yy = RadiusMaskSize(2)/2-(KernelSize(2))/2:RadiusMaskSize(2)/2+(KernelSize(2)/2)-1;
% RadiusMask(xx,yy) = ZeroOutSideRadius2L;

realteilL_masked  = realteilL;%.*RadiusMask;
imagteilL_masked  = imagteilL;%.*RadiusMask;

%%% normalize the kernel 
normRealL   = sqrt(sum(realteilL_masked(:).^2));
normImagL   = sqrt(sum(imagteilL_masked(:).^2));
normFactorL = kA^2;

norm_realteilL = realteilL_masked*normFactorL/(0.5*(normRealL + normImagL));
norm_imagteilL = imagteilL_masked*normFactorL/(0.5*(normRealL + normImagL));

SpaceKernel = complex(norm_realteilL, norm_imagteilL);
FreKernel   = fft2(SpaceKernel);

SpaceKernel_realsize = complex(norm_realteilL, norm_imagteilL);





