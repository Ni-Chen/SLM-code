close all;
clear all; 
clc;

%% Setting up system parameters
s_factor = 1;       % Set this to 1 for default (actual) values. Set this to larger than 1 for scaling the fft2. This ...
                    % ...will introduce increasing errors as you increase its value largr.
Lambda = 532e-9;    %Wavelength
k = 2*pi/Lambda;    %Wave number
f_lens = (400)*1e-3; %Focal length of the Fourier Transforming Lens
pp_slm = s_factor*12.5e-6; %spatial period of sampling (pixels)
Nx = 1024;
x_slm = pp_slm*((1:Nx)- Nx/2);

%% Input phase profile here: (example of sinusoidal grating shown)
% other 1D phase profiles of length 1024 can be included
c = 2*pi;           % contrast
period_px = 50;      % period in terms of pixels
f_0 = 1/(pp_slm*period_px);   % spatial frequency of the grating 
phi_sin = c/2 * sin(2*pi*f_0*x_slm)+c/2;
figure('Name','Sinusoidal grating Phase Profile (ColorMap)','NumberTitle','off');
imagesc(unwrap(phi_sin));
colormap gray;
h = colorbar;
ylabel(h,'\phi (units of \pi)');
axis off;

field_slm = exp(1i*phi_sin);		% The electric field

%% Computing the Fresnel Diffraction Integral
y = fftshift(fft(field_slm,Nx));

% Setting the axis
fx = ((1:Nx)- Nx/2)*(1/pp_slm/Nx);
x = s_factor*fx*Lambda*f_lens;       % converting from spatial frequency to coordinates

% Plotting the Phase profile in 1D
figure('NumberTitle','off','Name','1D Computation: Fourier Transform');
s(1) = subplot(2,2,1);
plot(x_slm, unwrap(angle(field_slm)));
% plot(x_slm, angle(field_slm)/pi);
ylabel('Phase \Phi (units of \pi)');
xlabel('X axis');

s(2) = subplot(2,2,2);
plot(x_slm, abs(field_slm));
ylabel('Amplitude');
xlabel('X axis');
title(s(1),'Phase profile at SLM');
title(s(2),'Amplitude profile at SLM (normalized)');

% Plotting the Image plane Field
s(4) = subplot(2,2,4);
plot(x/(1e-3),abs(y));
ylabel('Amplitude');
xlabel('X axis (mm)');

s(3) = subplot(2,2,3);
plot(x/(1e-3), unwrap(angle(y)));
ylabel('Phase \Phi (units of \pi)');
xlabel('X axis (mm)');
title(s(3),'Phase profile Image plane');
title(s(4),'Amplitude profile Image plane');

