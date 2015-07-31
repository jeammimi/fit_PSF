% gaussianfft(bits, width(um), phase koeff)
%clear;
load('SP_phaseMask.mat');
% RGB=imread('DH_phase.png');
% DH_phase = double(rgb2gray(RGB))/226*2*pi;
n=9;
x = (-2^(n-1):1:2^(n-1)-1)/2^n; % x vector
y= x';
truncatecirle=zeros(2^n);
for i=1:length(x)
    for j=1:length(y)
        if (i-2^(n-1))^2+(j-2^(n-1))^2 < (2^(n-2))^2
            truncatecirle(i,j) = 1;
        end
    end
end

sphere=zeros(2^n);
for i=1:length(x)
    for j=1:length(y)
        sphere(i,j)=((x(i))^2+(y(j))^2); % spherical lens
    end
end

n=9;
Icam=zeros(2^n,2^n);
tx = (-2^(n-1):1:2^(n-1)-1)/2^n;
ty=tx';
PSFzrange=5000;
PSFzframes=51;
cropsize=64;
PSFzscale=-PSFzrange/2:PSFzrange/(PSFzframes-1):PSFzrange/2;
inputph=SP_phase;
%hold on;
% for z1=0:zstep:zmax
PSFarray=zeros(2^n, 2^n, PSFzframes);
PSFarrayfft=zeros(2^n, 2^n, PSFzframes);
PSFarraysm=zeros(cropsize+1, cropsize+1, PSFzframes);
PSFarraysmnorm=zeros(cropsize+1, cropsize+1, PSFzframes);
zp=-PSFzrange/2:PSFzrange/(PSFzframes-1):PSFzrange/2;
#for ind=1:PSFzframes
I2=gaussianfft2(n,0.1,0,200,zp(25),inputph,sphere,truncatecirle);
#    PSFarray(:,:,ind)=I2;
    %    PSFarraysm(:,:,ind)=imcrop(I2,[230,230,64,64]);
    % PSFarraysm(:,:,ind)=imcrop(I2,[2^(n-1)-cropsize/2,2^(n-1)-cropsize/2,cropsize,cropsize]);
    %    PSFarrayfft(:,:,ind)=abs(fftshift(ifft2(I2)));
    % tmp=PSFarraysm(:,:,ind);
    % PSFarraysmnorm(:,:,ind)=tmp/max(tmp(:));
#end
