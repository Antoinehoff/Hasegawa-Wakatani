function [ conv ] = conv_thm_2D_AA_sym( F, G, AA )
%% 
[Nx, Ny] = size(F);

%% FFT of the functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = ifft2(fftshift(F),Nx,Nx,'symmetric')); % Filter to ensure that the functions are real
g = ifft2(fftshift(G),Nx,Nx,'symmetric'));

%% Convolution theorem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv_pad = ifftshift(fft2(f.*g));

%% Cut the domain back to input size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv_pad = conv_pad(1:Nx,1:Ny);

%% Cropping for frequencies > 2/3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv = AA.* conv_pad;

