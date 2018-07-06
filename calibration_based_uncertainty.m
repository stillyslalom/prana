function [Ux,Uy,UxLB,UxUB,UyLB,UyUB]=calibration_based_uncertainty(fitmodelname,metric,correlation_method)
% Compare input string to evaluate which calibration based method
% coeeficients should be used to evaluate uncertainty
if strcmp(fitmodelname,'PPR_Charonkomodel')
    
    ppr=metric;
    if strcmp(correlation_method,'SCC')
        Ap=0.226;Bp=1;Cp=0.08; Mp=13.1;Np=1;Sp=0.317;
        Uppr=sqrt((Mp*exp(-0.5*((ppr-Np)/Sp).^2)).^2+(Ap*ppr.^(-Bp)).^2+Cp.^2);
    elseif strcmp(corrleation_method,'RPC')
        Ap=1.41;Bp=1;Cp=1.7e-5; Mp=9.76;Np=1;Sp=1.139;
        Uppr=sqrt((Mp*exp(-0.5*((ppr-Np)/Sp).^2)).^2+(Ap*ppr.^(-Bp)).^2+Cp.^2);
    end
    
    Ux=Uppr/sqrt(2);
    Uy=Uppr/sqrt(2);
    UxLB=0;
    UyLB=0;
    UxUB=0;
    UyUB=0;

elseif strcmp(fitmodelname,'PPR_Xuemodel')
    
    ppr=metric;
    
    if strcmp(correlation_method,'SCC')
        
        % Uncertainty Lower Bound
        Apl=0.1043;Bpl=0.6786;Cpl=0;Mpl=0.278;Npl=1;Spl=0.1927;
        UpprLB=sqrt((Mpl*exp(-0.5*((ppr-Npl)/Spl).^2)).^2+(Apl*ppr.^(-Bpl)).^2+Cpl.^2);
        % Uncertainty Upper Bound
        Apu=0.6888;Bpu=0.846;Cpu=0;Mpu=10.59;Npu=1;Spu=0.1925;
        UpprUB=sqrt((Mpu*exp(-0.5*((ppr-Npu)/Spu).^2)).^2+(Apu*ppr.^(-Bpu)).^2+Cpu.^2);
    elseif strcmp(correlation_method,'RPC')
        
        % Uncertainty Lower Bound
        Apl=0.09359;Bpl=0.5597;Cpl=0;Mpl=0.09828;Npl=1;Spl=0.2258;
        UpprLB=sqrt((Mpl*exp(-0.5*((ppr-Npl)/Spl).^2)).^2+(Apl*ppr.^(-Bpl)).^2+Cpl.^2);
        % Uncertainty Upper Bound
        Apu=0.4583;Bpu=0.5696;Cpu=0;Mpu=25.11;Npu=1;Spu=0.2874;
        UpprUB=sqrt((Mpu*exp(-0.5*((ppr-Npu)/Spu).^2)).^2+(Apu*ppr.^(-Bpu)).^2+Cpu.^2);
        
    end
    
    
    Ux=0;
    Uy=0;
    UxLB=UpprLB/sqrt(2);
    UyLB=UpprLB/sqrt(2);
    UxUB=UpprUB/sqrt(2);
    UyUB=UpprUB/sqrt(2);

elseif strcmp(fitmodelname,'MI_Xuemodel')
    
    mi=metric;
    
    if strcmp(correlation_method,'SCC')
        % Uncertainty Lower Bound
        Al=0.1525;Bl=0.7294;Cl=0.0164;Ml=0.3874;Nl=0;Sl=1.048;
        UmiLB=sqrt((Ml*exp(-0.5*((mi-Nl)/Sl).^2)).^2+(Al*mi.^(-Bl)).^2+Cl.^2);
        % Uncertainty Upper Bound
        Au=1.101;Bu=1.073;Cu=0.0980;Mu=65.59;Nu=0;Su=0.9197;
        UmiUB=sqrt((Mu*exp(-0.5*((mi-Nu)/Su).^2)).^2+(Au*mi.^(-Bu)).^2+Cu.^2);
        
    elseif strcmp(correlation_method,'RPC')
        
        % Uncertainty Lower Bound
        Al=0.0479;Bl=0.3252;Cl=0;Ml=0.0012;Nl=0;Sl=0.0067;
        UmiLB=sqrt((Ml*exp(-0.5*((mi-Nl)/Sl).^2)).^2+(Al*mi.^(-Bl)).^2+Cl.^2);
        % Uncertainty Upper Bound
        Au=0.1289;Bu=0.4222;Cu=0.0662;Mu=0.3367;Nu=0;Su=3.575;
        UmiUB=sqrt((Mu*exp(-0.5*((mi-Nu)/Su).^2)).^2+(Au*mi.^(-Bu)).^2+Cu.^2);
        
    end
    
    Ux=0;
    Uy=0;
    UxLB=UmiLB/sqrt(2);
    UyLB=UmiLB/sqrt(2);
    UxUB=UmiUB/sqrt(2);
    UyUB=UmiUB/sqrt(2);

end

% Further calibration based methods can be similarly added

end
