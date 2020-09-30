%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if START % Start from 0
    NN_0 =  FT_f(KX,KY);
    ZZ_0 = -(KX.^2+KY.^2).*NN_0; 
    PP_0 = NN_0.*((KX~=0).*(KY~=0));
    NN_0 = NN_0(1:Nkx,(Ny/2*(HALF)+1):Ny); 
    PP_0 = PP_0(1:Nkx,(Ny/2*(HALF)+1):Ny); 
    ZZ_0 = ZZ_0(1:Nkx,(Ny/2*(HALF)+1):Ny);
    % init counters
    it0  = 2; t0 = 0; SID = 00;
else % restart from a save
    load_data
    SID = SID + 1;
end

%% Time domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = t0:dt:t0+Tmax;
Nt   = numel(time);
Skip = ceil(1/SPS/dt);
Ns   = ceil(Nt/Skip);
Ts   = time(1:Skip:end);

%% Initialization and time integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZZ = zeros(Nkx,Nky,Ns); ZZ(:,:,1) = ZZ_0; % records
NN = zeros(Nkx,Nky,Ns); NN(:,:,1) = NN_0; % records
PP = zeros(Nkx,Nky,Ns); PP(:,:,1) = PP_0; % records
SN = zeros(Nkx,Nky,Ns); SN(:,:,1) = 0; % records
SZ = zeros(Nkx,Nky,Ns); SZ(:,:,1) = 0; % records
ZZ_1 = ZZ_0; NN_1 = NN_0; PP_1  = PP_0; % Next time step variables

% Substract zonal component of adiabatic term if Modified-Hasegawa-Wakatani
if MHW
   KY_eq_0 = (KY==0);
   RHSz = @(Z,N,P,S) HWz(Z,N,P,S) - alpha * (P-N) .* KY_eq_0;
   RHSn = @(Z,N,P,S) HWn(Z,N,P,S) - alpha * (P-N) .* KY_eq_0;
else
   RHSz = @(Z,N,P,S) HWz(Z,N,P,S);
   RHSn = @(Z,N,P,S) HWn(Z,N,P,S);
end

% Poisson bracket function (cancel it if linear sim)
if LINEAR
    compute_S = @(f,g) 0;
else
    compute_S = @(f,g) Poisson_bracket(f,g);
end

%% RK4 explicit solver
disp('Run RK4 solver...')
nbytes = fprintf(2,'step %d/%d',0,Nt-1);
    for it = it0:numel(time)
        
        %Computing non linear terms
        Spz = compute_S(PP_0,ZZ_0);
        Spn = compute_S(PP_0,NN_0);

        %Updating vorticity ZZ_1
        k_1 = RHSz(ZZ_0             , NN_0, PP_0, Spz);
        k_2 = RHSz(ZZ_0 + 0.5*dt*k_1, NN_0, PP_0, Spz);
        k_3 = RHSz(ZZ_0 + 0.5*dt*k_2, NN_0, PP_0, Spz);
        k_4 = RHSz(ZZ_0 +     dt*k_3, NN_0, PP_0, Spz);    

        ZZ_1 = ZZ_0 + dt*(1./6.)*(k_1 + 2.*k_2 + 2.*k_3 + k_4);

        % Updating density NN_1
        k_1 = RHSn(ZZ_0, NN_0             , PP_0, Spn);
        k_2 = RHSn(ZZ_0, NN_0 + 0.5*dt*k_1, PP_0, Spn);
        k_3 = RHSn(ZZ_0, NN_0 + 0.5*dt*k_2, PP_0, Spn);
        k_4 = RHSn(ZZ_0, NN_0 +     dt*k_3, PP_0, Spn);    

        NN_1 = NN_0 + dt*(1./6.)*(k_1 + 2.*k_2 + 2.*k_3 + k_4);
        
        % Solving Poisson equation
        PP_1 = -ZZ_1./(KX.^2+KY.^2);
        PP_1(abs(PP_1)==inf) = 0.0; %cancel singularity

        % Old -> new
        ZZ_0 = ZZ_1; NN_0 = NN_1; PP_0 = PP_1;

        % Recordings and check
        if mod(it,Skip) == 0
            PP(:,:,floor(it/Skip)+1) = PP_0;
            NN(:,:,floor(it/Skip)+1) = NN_0;
            ZZ(:,:,floor(it/Skip)+1) = ZZ_0;
            SN(:,:,floor(it/Skip)+1) = Spn;
            SZ(:,:,floor(it/Skip)+1) = Spz;
            % terminal info
            while nbytes > 0
              fprintf('\b')
              nbytes = nbytes - 1;
            end
            nbytes = fprintf(2,'step %d/%d',it,Nt-1);
            if sum(sum(isnan(ZZ_0)+isnan(NN_0)+isnan(PP_0))) > 0
                disp('')
                disp('Nan detected')
                FILTER = 1;
                break
            end
        end
        FILTER = 0;
    end
disp(' ')
disp('Done')
%% Filter the nan
if FILTER
    for it = 1:Ns
        if sum(sum(isnan(ZZ(:,:,it)))) > 0 ||...
           sum(sum(isnan(NN(:,:,it)))) > 0 ||...
           sum(sum(isnan(PP(:,:,it)))) > 0
           itcut = it-1;
           break
        end
    end
    ZZ = ZZ(:,:,1:itcut); NN = NN(:,:,1:itcut); PP = PP(:,:,1:itcut);
    Ns = numel(NN(1,1,:)); Ts = Ts(1:Ns);
end