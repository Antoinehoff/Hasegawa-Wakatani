function [ conv ] = conv_thm_2D_AA_sym( F, G, AA )
%% 
[Nx, Ny] = size(F);

if Ny == Nx/2
    shift = @(x) fftshift(x,1);
   ishift = @(x) ifftshift(x,1);
else
    shift = @(x) fftshift(x);
   ishift = @(x) ifftshift(x);
end
%% FFT of the functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = ifft2(shift(F),Nx,Nx,'symmetric'); % Symmetric ensures Real
g = ifft2(shift(G),Nx,Nx,'symmetric');

%% Convolution theorem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv_pad = ishift(fft2(f.*g));

%% Cut the domain back to input size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv_pad = conv_pad(1:Nx,1:Ny);

%% Cropping for frequencies > 2/3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv = AA.* conv_pad;

