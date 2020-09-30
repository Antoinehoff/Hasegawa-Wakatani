close all; clear all;
addpath(genpath('../matlab')) % ... add
default_plots_options
SIMID = 'test_init_';
%% Basic functions in real space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gain   = 0.5; sigma_exp = 2.0;
f      = @( x, y) Gain*exp(-(x.^2+y.^2)/sigma_exp^2);
FT_f   = @(kx,ky) Gain*sigma_exp*exp(-sigma_exp^2*(kx.^2+ky.^2)/4)/sqrt(2);
ggf    = @( x, y) Gain*4.0/sigma_exp^2*(x.^2/sigma_exp^2 + y.^2/sigma_exp^2 - 1).*f(x,y);
FT_ggf = @(kx,ky) Gain*-(kx.^2+ky.^2).*FT_f(kx,ky);

%% Grid in real space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HALF = 1;
Nx   = 256;   Ny   = Nx;
Lx   = 40;    Ly   = Lx;
dx   = Lx/(Nx); dy   = Ly/(Ny);
x = dx*(-Nx/2:(Nx/2-1));    % x = linspace(-Lx/2,Lx/2,Nx);
y = dy*(-Ny/2:(Ny/2-1));% y = linspace(-Ly/2,Ly/2,Ny);
[YY,XX] = meshgrid(y,x);

%% Frequency grid linked to the real space grid %%%%%%%%%%%%%%%%%%%%%%%%%%
Nkx= Nx;     Nky  = Ny*(1-0.5*HALF);
dkx= (2*pi/Nx/dx); dky = (2*pi/Ny/dy);
kx = dkx*(-(Nx/2):(Nx/2-1));
ky = dky*(-(Ny/2*(1-HALF)):(Ny/2-1));
[KY,KX] = meshgrid(ky,kx);
RES = [num2str(Nkx),'x',num2str(Nky),'_'];
GRID= ['L_',num2str(Lx),'_'];

%% Model and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LINEAR   = 0; % Cancel non-linearity, analysis computes analytical growth rate
MHW      = 0; % modified Hasegawa-Wakatani system
alpha    = 1.0;
kappa    = 1.0;
mu       = 1e-4;
PARAMS   = ['a_',num2str(alpha),'_mu_',num2str(mu),'_'];
if LINEAR; PARAMS = ['lin_',PARAMS]; end
if MHW;    MODEL  = 'MHW_'; else; MODEL = 'HW_'; end

% Hasegawa-Wakatani system
HWz = @(Z,N,P,S) -S + alpha*(P-N) - mu*((KX.^2+KY.^2).^2).*Z;
HWn = @(Z,N,P,S) -S + alpha*(P-N) - mu*((KX.^2+KY.^2).^2).*N - kappa*1i*KY.*P;

% Poisson bracket computation
kxmax = max((kx)); kymax = max((ky));
AA = real( kx>-2/3*kxmax & kx<2/3*kxmax )' ... % Anti aliasing filter
    *real( ky>-2/3*kymax & ky<2/3*kymax );
% Pb = @(f,g) conv_thm_2D_iff_real_AA(1i*KX.*f, 1i*KY.*g, AA)*Nx*dx*Ny*dy ...
%            -conv_thm_2D_iff_real_AA(1i*KY.*f, 1i*KX.*g, AA)*Nx*dx*Ny*dy;
Poisson_bracket = @(f,g) ...
    conv_thm_2D_AA_sym(1i*KX.*f, 1i*KY.*g, AA)*Nx*dx*Ny*dy ...
    -conv_thm_2D_AA_sym(1i*KY.*f, 1i*KX.*g, AA)*Nx*dx*Ny*dy;
       
%% Time parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmax  = 50;
dt    = 1e-3;
SPS   = 2; % Sampling per time unit
START = 1; % To start the simulation from t=0
SID   = 00;% To load data and run from time trst to Tmax
t_rst = 29; 
%% Recap of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BASIC.SIMID = [SIMID,MODEL,RES,GRID,PARAMS(1:end-1)]; disp(BASIC.SIMID)
disp(['START = ',num2str(START),', SID = ',num2str(SID)])
