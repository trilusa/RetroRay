 function f3db = calcPDBandwidth(D_d)
     %D_d is vector of diode dimaters
    A = pi * (D_d/2).^2;  %m^2
    
    ep0 = 8.85e-12;  %pemitivity of free space (F/m)
    epSi=11.9;      %dielectric constant of Si (ratio)
    mu = .14;       % mobiltiy of electrons at 300K m^2/V
    rho = 2.3;       % resistivity of Silicon Ohm*meter (2.3e3)
    Vbuiltin = 0;     %builtin Voltage of silicon, assume neglible (V)
    Vbias = 5;     %applied reverse bias. Will set at 10V for now

    Wd = sqrt(2*ep0*epSi*mu*rho*(Vbias+Vbuiltin));  %delpetion width
    Cjunc = (ep0 * epSi  * A) / Wd;  % junction capactinace (F)
    Cstray = 0;                        %neglible stray capaciteance
    Ctotal = Cjunc + Cstray;


    Rseries = 1000;         %should be 0 but 10-1000 in practice, meausred (ohms)
    Rload = 50;             %load resistance, recomended to be 50 ohms
    Rtotal = Rseries + Rload;  % treating shunt resistance as neglible

    tRC = 2.2*Rtotal*Cjunc; %rise time with just RC (seconds)
    tdrift = 0;             %assuming at large Area tRC will dominate
    tdiffuse = 0;            %assuming at large Area RC will dominate
    tRise = sqrt( tdrift.^2 + tdiffuse.^2 + tRC.^2 );
    f3db = .35 ./ tRise;
    % f3db = 1/(2*pi*Rtotal*Ctotal) also works
 end
%Specs for PD we ordered
 % 380 pF @ 5 V
%3.8000e-10
% 150 ns / 150 nsd @ 5 V
%1.5000e-07