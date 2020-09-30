function [ conv ] = conv_thm_1D( F, G, L, Pad )
%conv_thm_2D computes the convolution between F and G with a zero pading Pad
% kspace -> real -> product -> kspace

Nr = numel(F);
Mr = Nr * Pad;

f  = ifft((F),Mr);
g  = ifft((G),Mr);
conv_pad = fft(f.*g); % convolution becomes product

if Pad == 2
    cmin = (Pad-1)/2; cmax = (Pad+1)/2;
    conv = L*Pad*(conv_pad(cmin*Nr+1:cmax*Nr)); % remove padding
elseif Pad == 1
    conv = L*fftshift(conv_pad);
end

end

