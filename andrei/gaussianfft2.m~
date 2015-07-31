function I2 = gaussianfft2( n,d,xc,yc,z,SLM,sphere,truncatecirle)
%gaussianfft Summary of this function goes here
%I use the gaussian source with n bits and w width
%after fourier conversion I add a quadratic phase to emulate out of focus
%in microsope and make back fourier transform to see the shape broadening.
% The output will be the fwhm of the resulting gaussian profile.
%   Detailed explanation goes here
% n - number of bits
% d - source diameter in um
% z - z-distance out of focus
% xc - shift in um
% yc ===
%load('/Users/andrey/Downloads/SP_phaseMask.mat');

%n=9;
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

load('SP_phaseMask.mat');
SLM=SP_phase;

FOV=50; % Field of view of the objective (um)
x = (-2^(n-1):1:2^(n-1)-1)/2^n; % x vector
y= x';
we=d/FOV;
xc=xc/FOV;
yc=yc/FOV;
% E0 = (exp(-y.^2)*exp(-x.^2))/(d/FOV)^2;   % Mag
E0 = exp(-(y-yc).^2/(we^2/2))*(exp(-(x-xc).^2/(we^2/2)));
I0=E0.^2; %intensity
% k = z/5; %derived from Rayleigh length for wl=600nm, w0=170nm.
k = z/10; %adapted using Grover DH

nlph=SLM+sphere*k;
%E0=gpuArray(E0);
FFE0=fft2(E0);
FPlane=abs(fftshift(FFE0));
% truncate by circle
FPlane=FPlane.*truncatecirle;
% Truncate by rectangle
% for i=1:length(x)
%     for j=1:length(y)
%         if abs((j-2^(n-1))) > 2^(n-3);
%             ph(i,j) = 0;
%         end
%     end
% end

ph=fftshift(nlph);
%Back Fourier including nonlinear phase
E2=ifft2(fftshift(FPlane).*exp(1i*(angle(FFE0)+ph))); 
amp=abs(E2);
I2=amp.^2;
%I2=gather(I2);
tum=x.*512*.1; % t in um
gss=0;%fit(tum.',I2.','gauss1');
gssc=0;%coeffvalues(gss);
width=0;%gssc(3)*2*sqrt(2*log(2)); %width in nm
I2=amp.^2;
end

