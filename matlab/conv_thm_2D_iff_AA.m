function [ conv ] = conv_thm_2D_iff_AA( ff, gg, AA )
%% FFT of the functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_pad = (ifft2(fftshift(ff)));% Normalized to match the fourier transform
G_pad = (ifft2(fftshift(gg)));

%% Convolution theorem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv_pad = ifftshift(fft2(F_pad.*G_pad));

%% Cropping for frequencies > 2/3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv = AA.* conv_pad;

