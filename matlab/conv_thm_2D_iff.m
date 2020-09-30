function [ conv ] = conv_thm_2D_iff( ff, gg, Pad )
%% Dealiasing through zero-padding FF and GG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Nx,Ny] = size(ff);
ixmin = Nx*(Pad-1)/2; iymin = Ny*(Pad-1)/2;

f_pad = zeros(Pad*Nx, Pad*Ny); 
f_pad(ixmin+(1:Nx),iymin+(1:Ny)) = ff;   

g_pad = zeros(Pad*Nx, Pad*Ny);
g_pad(ixmin+(1:Nx),iymin+(1:Ny)) = gg;

%% FFT of the functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_pad = ifft2(fftshift(f_pad));% Normalized to match the fourier transform
G_pad = ifft2(fftshift(g_pad));

%% Convolution theorem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv_pad = ifftshift(fft2(F_pad.*G_pad));

%% Cropping the padding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv = conv_pad(ixmin+(1:Nx),iymin+(1:Ny));
end

