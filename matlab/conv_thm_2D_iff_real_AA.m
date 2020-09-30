function [ conv ] = conv_thm_2D_iff_real_AA( F, G, AA )
%% FFT of the functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = real(ifft2(fftshift(F),'symmetric')); % Filter to ensure that the functions are real
g = real(ifft2(fftshift(G),'symmetric'));

%% Convolution theorem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv_pad = ifftshift(fft2(f.*g));

%% Cropping for frequencies > 2/3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv = AA.* conv_pad;

