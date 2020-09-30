function [ conv ] = conv_thm_2D_iff_real( ff, gg )
%% FFT of the functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_pad = real(ifft2(fftshift(ff)));% Normalized to match the fourier transform
G_pad = real(ifft2(fftshift(gg)));

%% Convolution theorem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv = ifftshift(fft2(F_pad.*G_pad));
end

