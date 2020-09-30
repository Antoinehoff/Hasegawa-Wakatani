%% Hasegawa-Wakatani linear dispersion relation
function [maxG, GG, KX_, KY_] = lin_disp_rel(alpha,kappa)
    gamma = @(kx,ky) max(imag(roots([(kx^2+ky^2), 1i*alpha*((kx^2+ky^2)+1), -alpha^2 + alpha - 1i*ky*kappa])));

    %% Frequency space %%%%%%%%%%%%%%%%%%%%%%%%%%
    Nx_ = 100; Ny_ = 100;
    dkx_= 0.1; dky_ = 0.1;
    kx_ = dkx_*(0:Nx_-1);
    ky_ = dky_*(0:Ny_-1);
    [KY_,KX_] = meshgrid(ky_,kx_);

    %% Evaluation of gamma on the grid %%%%%%%%%%
    GG = zeros(Nx_,Ny_);

    if alpha == 0 || kappa == 0
        maxG = 0;
    else
        for ikx = 1:Nx_
            for iky = 1:Ny_
                GG(ikx,iky) = gamma(kx_(ikx),ky_(iky));
            end
        end
        maxG = (max(max(GG)));
    end
end