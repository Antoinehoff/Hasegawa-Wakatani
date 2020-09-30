%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IFFT
zeta   = zeros(Nx,Nx,Ns);
ni     = zeros(Nx,Nx,Ns);
phi    = zeros(Nx,Nx,Ns);
if HALF
    shift =  @(x) ifftshift(x,1);
else
    shift = @(x) ifftshift(x);
end

for it = 1:numel(PP(1,1,:))
    ZZ_ = ZZ(:,:,it); NN_ = NN(:,:,it); PP_ = PP(:,:,it);
    F_          = (ifft2(shift(ZZ_),Nx,Nx,'symmetric'));
    zeta(:,:,it)= ifftshift(F_);
    F_          = (ifft2(shift(NN_),Nx,Nx,'symmetric'));
    ni(:,:,it)   = ifftshift(F_);
    F_          = (ifft2(shift(PP_),Nx,Nx,'symmetric'));
    phi(:,:,it) = ifftshift(F_);
end

%% Post processing
phi_ST = zeros(Ny,Ns);   % Space-Time diagram of ES potential
ni_ST  = zeros(Ny,Ns);   % Space-Time diagram of densitz
phi_00 = zeros(1,Ns);    % Time evolution of ES potential at origin
ni_00  = zeros(1,Ns);    % Time evolution of density at origin
Sn_norm= zeros(1,Ns);    % Time evolution of the amp of density nonlin term
Sz_norm= zeros(1,Ns);    % Time evolution of the amp of vorti. nonlin term
E_pot  = zeros(1,Ns);    % Potential energy n^2
E_kin  = zeros(1,Ns);    % Kinetic energy grad(phi)^2
W      = zeros(1,Ns);    % Enstrophy
G_n    = zeros(1,Ns);    % Background density energy source
G_a    = zeros(1,Ns);    % Adiabadicity energy sink
D_E    = zeros(1,Ns);    % Dissipative energy term
D_W    = zeros(1,Ns);    % Dissipative vorticity term
ExB    = zeros(1,Ns);    % ExB drift intensity \propto |\grad \phi|
CFL    = zeros(1,Ns);    % CFL time step
Ddx = 1i*KX; Ddy = 1i*KY; lapl   = Ddx.^2 + Ddy.^2; 
for it = 1:numel(PP(1,1,:))
    ZZ_ = ZZ(:,:,it); NN_ = NN(:,:,it); PP_ = PP(:,:,it);
    phi_ST(:,it)= phi(:,y==0,it);
    ni_ST(:,it) = ni(:,y==0,it);
    phi_00(it)  = phi(x==0,y==0,it);
    ni_00(it)   = ni(x==0,y==0,it);
    Sn_norm(it) = sum(sum(abs(SN(:,:,it))));
    Sz_norm(it) = sum(sum(abs(SZ(:,:,it))));
    E_pot(it)   = pi/Lx/Ly*sum(sum(abs(NN_).^2))/Nkx/Nky; % integrate through Parseval id
    E_kin(it)   = pi/Lx/Ly*sum(sum(abs(Ddx.*PP_).^2+abs(Ddy.*PP_).^2))/Nkx/Nky;
    W(it)       = pi/Lx/Ly*sum(sum(abs(NN_ - lapl.*PP_).^2))/Nkx/Nky;
    G_n(it)     =real(-2*pi*kappa/Lx/Ly*sum(sum((NN_.*conj(Ddy.*PP_))))/Nkx/Nky);
    G_a(it)     = 2*pi*alpha/Lx/Ly*sum(sum(abs(NN_-PP_).^2))/Nkx/Nky;
    D_E(it)     = 2*pi*mu   /Lx/Ly*sum(sum(abs(lapl.*NN_).^2 + abs(Ddx.*(lapl.*PP_)).^2 + abs(Ddy.*(lapl.*PP_)).^2))/Nkx/Nky;
    D_W(it)     = real(2*pi*mu   /Lx/Ly*sum(sum(abs(lapl.*NN_).^2 + abs(lapl.^2.*PP_).^2 - 2*(lapl.*NN_).*conj(lapl.^2.*PP_)))/Nkx/Nky);
    ExB(it)     = max(max(max(abs(phi(3:end,:,it)-phi(1:end-2,:,it))/(2*dx))),max(max(abs(phi(:,3:end,it)-phi(:,1:end-2,it))'/(2*dy))));
    CFL(it)     = 4*min([dx^2/mu,dx/kappa,ExB(it)]);
end
E_kin_ky = mean(mean(abs(Ddx.*PP(:,:,it)).^2+abs(Ddy.*PP(:,:,it)).^2,3),1);
E_kin_kx = mean(mean(abs(Ddx.*PP(:,:,it)).^2+abs(Ddy.*PP(:,:,it)).^2,3),2);
dEdt     = diff(E_pot+E_kin)./diff(Ts);
dWdt     = diff(W)./diff(Ts);
%% PLOTS
FTYPE = '.fig';
%% Time evolutions
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',SID)];
subplot(221)
    semilogy(Ts,abs(ni_00),'-','DisplayName','$n$')
    grid on; xlabel('$t$'); ylabel('$|n(x=0,y=0)|$');
subplot(222)
    semilogy(Ts,abs(phi_00),'-','DisplayName','$\phi$')
    grid on; xlabel('$t$'); ylabel('$|\phi(x=0,y=0)|$');
subplot(223)
    semilogy(Ts,E_kin+E_pot,'-','DisplayName','$\sum|ik\tilde\phi_i|^2+\sum|\tilde n_i|^2$')
    hold on;
    if LINEAR % Plot linear growth rate
        [Gmax, GG] = HW_lin_disp_rel(alpha,kappa,mu,KX,KY);
        semilogy(Ts,(E_pot(end)+E_kin(end)).*exp(2*Gmax.*(Ts-Ts(end))),'--k','DisplayName','$\exp(2\gamma_{\max}t)$')
    end
    grid on; xlabel('$t$'); ylabel('$E$'); legend('show');
subplot(224)
    semilogy(Ts,Sn_norm,'-','DisplayName','$\sum|S_n|$'); 
    hold on;
    semilogy(Ts,Sz_norm,'-','DisplayName','$\sum|S_z|$');
    grid on; xlabel('$t$'); ylabel('$S$'); legend('show');
save_figure

%% Energy balance
Gn_mid = (G_n(1:end-1)+G_n(2:end))/2.0;
Ga_mid = (G_a(1:end-1)+G_a(2:end))/2.0;
DE_mid = (D_E(1:end-1)+D_E(2:end))/2.0;
DW_mid = (D_W(1:end-1)+D_W(2:end))/2.0;
fig = figure; FIGNAME = ['Energy_balance',sprintf('_%.2d',SID)];
    subplot(221); title('Energy balance terms')
        plot(Ts(2:end),dEdt,'DisplayName','$\partial_t E$'); hold on;
        plot(Ts, G_n, 'DisplayName', '$\Gamma_n$');
        plot(Ts, G_a, 'DisplayName', '$\Gamma_a$');
        plot(Ts, D_E, 'DisplayName', '$D_E$');
        grid on; xlabel('$t$');  legend('show');
    subplot(223); title('Enstrophy balance terms')
        plot(Ts(2:end),dWdt,'DisplayName','$\partial_t W$'); hold on;
        plot(Ts, G_n, 'DisplayName', '$\Gamma_n$');
        plot(Ts, D_W, 'DisplayName', '$D_W$');
        grid on; xlabel('$t$');  legend('show');
    subplot(122); title('Conservation rel. error')
        semilogy(Ts(2:end),100*abs(dEdt - (Gn_mid-Ga_mid-DE_mid))./abs(dEdt),'DisplayName', 'Energy'); hold on;
        plot(Ts(2:end),100*abs(dWdt - (Gn_mid-DW_mid))./abs(dWdt),'DisplayName', 'Enstrophy')
        grid on; xlabel('$t$'); ylabel('$\epsilon[\%]$'); legend('show');
save_figure

%% Spectra energy
fig = figure; FIGNAME = ['Energy_kin_ky',sprintf('_%.2d',SID)];
semilogy(kx(end/2+1:end),E_kin_kx(end/2+1:end),'o','DisplayName','$\sum_y\langle|ik\tilde\phi_i|^2\rangle_t$')
hold on;
loglog(ky(end/2+1:end),E_kin_ky(end/2+1:end),'o','DisplayName','$\sum_x\langle|ik\tilde\phi_i|^2\rangle_t$')
grid on; xlabel('$k$');  legend('show');
save_figure

%% CFL condition
fig = figure; FIGNAME = ['CFL',sprintf('_%.2d',SID)];
semilogy(Ts,dy./ExB,'-','DisplayName','$|\nabla \phi|\Delta y$');
hold on;
plot(Ts,dx*dy/mu*ones(1,numel(Ts)),'-','DisplayName','$\Delta x \Delta y / \mu$');
plot(Ts,dy/kappa*ones(1,numel(Ts)),'-','DisplayName','$\Delta y/ \kappa$');
plot(Ts,dt*ones(1,numel(Ts)),'--k','DisplayName','$\Delta t$');
grid on; xlabel('$t$'); ylabel('$\Delta t$'); legend('show');
save_figure

%% Space-Time diagram at ky = 0
plt = @(x) real(x);
%% phi
fig = figure; FIGNAME = ['phi_ST',sprintf('_%.2d',SID)];
[TY,TX] = meshgrid(Ts,y);
pclr = pcolor(TX,TY,(plt(phi_ST))); set(pclr, 'edgecolor','none'); colorbar;
xlabel('$x\,(y=0)$'); ylabel('$t$'); title('$\phi$');
save_figure

if 0
%% Show frame
it = 100;
pclr = pcolor(KX,KY,(real(PP(:,:,it)))); set(pclr, 'edgecolor','none'); colorbar;
xlabel('$kx$'); ylabel('$ky$'); title(sprintf('t=%.3d',Ts(it)));
end
%%
if 0
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DELAY = 0.1; skip_ = 1;
%% Vorticity
GIFNAME = ['zeta',sprintf('_%.2d',SID)]; FIELDNAME = '$\zeta$';
FIELD = real(zeta(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
%% Density
GIFNAME = ['n',sprintf('_%.2d',SID)]; FIELDNAME = '$n$';
FIELD = real(ni(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',SID)]; FIELDNAME = '$\phi$';
FIELD = real(phi(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
%% FPhi
GIFNAME = ['Fphi',sprintf('_%.2d',SID)]; FIELDNAME = '$\phi$';
FIELD = real(PP(:,:,1:skip_:end)); X = KX; Y = KY; T = Ts(1:skip_:end);
create_gif
end