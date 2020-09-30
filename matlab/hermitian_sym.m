function [ FF_full ] = hermitian_sym( FF_half )
%HERMITIAN_SYM Reconstruct the full field from a half domain including 0
[Nx,Ny] = size(FF_half);
Ny      = 2*(Ny-1);

FF_full = zeros(Nx,Ny);

FF_full(:,1:Ny/2+1) = FF_half;

V1 = FF_half(1,2:Ny/2);
V2 = FF_half(Nx/2+1,2:Ny/2);

S1 = FF_half(2:Nx/2,2:Ny/2);
S2 = FF_half(2+Nx/2:end,2:Ny/2);

FF_full(1,Ny/2+2:end) = flip((V1));
FF_full(Nx/2+1,Ny/2+2:end) = flip(conj(V2));
FF_full(2:Nx/2,Ny/2+2:end) = flip(flip(conj(S2),2),1);
FF_full(Nx/2+2:end,Ny/2+2:end) = flip(flip(conj(S1),2),1);

end

