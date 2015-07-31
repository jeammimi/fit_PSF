% initialization routines
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

% load('zcorrtable.mat');
% load('z3corrtable.mat');
% load('z2corrtable.mat');
% 
% zcorrtable3(:,2)=zcorrtable(:,2) + z2corrtable(:,2) + z3corrtable(:,2);
% zcorrtable3(:,1)=zcorrtable(:,1);

% load('zcortable7may.mat');

SP_PSF_generation

% Correction table for xcorr2
for i=1:PSFzframes
tmp=(xcorr2(PSFarraysm(:,:,i),PSFarraysm(:,:,i)));
CC(i)=max(tmp(:));
end
CC1=CC/max(CC(:));
CC2=CC1*0.5+ones(1,PSFzframes)*.5;
xcorrnorm1=CC1.^(-1);
xcorrnorm2=CC2.^(-1);
xvector=ones(1,576);
xcorr2normmartix=xcorrnorm1'*xvector;
xcorr2normmartix2=xcorrnorm2'*xvector;
