% This piece of code can take a phase profile and give you the resulting
% electric field in the fourier plane. 

% close all;
clear all;
clc;
%% Setting up system parameters

s_factor = 1;       % Set this to 1 for default (actual) values. Set this to larger than 1 for scaling the fft2. This ...
                    % ...will introduce increasing errors as you increase its value largr.

D = pwd;            % Working Directory of this file
Lambda = 532e-9;    % Wavelength
k = 2*pi/Lambda;    % Wave number
f_lens = 50e-3;     % Focal length of the Fourier Transforming Lens
pp_slm = s_factor*8e-6; %spatial period of sampling (pixels)
f_max = 1/pp_slm;      % spatial freq of sampling (pixels)

%% Reading image file and converting to proper phase scaling
S = fullfile(D,'./data/1920x1080/DOEs/grid.png');   %input the phase profile
phi_slm = double(imread(S))/255*2*pi;


%% Generating 2D Eletric Field
Nx = 1920;  % X axis
Ny = 1080;  % Y axis

x_slm = pp_slm*((1:Nx)- Nx/2);
y_slm = pp_slm*((1:Ny)- Ny/2);

a_slm = 1;
field_slm = a_slm.*exp(1i*phi_slm);


%% The Fourier transform (2D) 
% The electric field in the Fourier plane
field_fourier = fftshift(fft2(field_slm)/length(field_slm));


% spatial coordinates in the reconstruction plane
fx = ((1:Nx)- Nx/2)*(1/pp_slm/Nx);
x = s_factor*fx*Lambda*f_lens;       % converting from spatial frequency to coordinates

fy = ((1:Ny)- Ny/2)*(1/pp_slm/Ny);
y = s_factor*fy*Lambda*f_lens;       % converting from spatial frequency to coordinates

%% Display image
figure('NumberTitle','off','Name','2D Computation: Fourier Transform');

s(1) = subplot(2,2,1);
imagesc(x_slm, y_slm, abs(field_slm));colorbar();
xlabel('x_{SLM}');
ylabel('y_{SLM}');
title(s(1),'Amplitude profile SLM plane (normalized)');

s(2) = subplot(2,2,2);
imshow(phi_slm,[0 2*pi]); colorbar(); title(s(2),'Input Phase profile');

s(3) = subplot(2,2,3);
imagesc(x/(1e-3), y/(1e-3), (abs(field_fourier).^2)); colorbar();
xlim([-0.8 0.8]);
ylim([-1 1]);
title(s(3),'Amplitude profile Image plane');
xlabel('x_{image} / (mm)');
ylabel('y_{image} / (mm)');



