% ==============================================================================================
% This code generates the small-signal model of a grid-forming PEM hydrogen
% electrolyzer system. It utilizes the parameters of the PEM electrolyzer,
% the DC-DC buck converter, the DC-AC inverter, and the operating point 
% around which the small-signal model is linearized.
%
% The output of this script is a linearized state-space representation 
% suitable for analyzing the dynamic behavior of the PEMEL unit.
% This model is particularly useful for studying integration scenarios 
% involving other loads and energy resources within power systems.
% ==============================================================================================

%% PEMEL parameters initialization
format long
% ==============================================================================================
% Environment conditions
env_p = 0.101325; % [MPa] Environment Pressure
env_T = 20; % [degC] Environment Temperature
% ==============================================================================================
% PEMEL parameters
cells_per_stack=48; % Number of cells in a stack
Stacks_No=5; % Number of stacks 
stack_num_cells = cells_per_stack*Stacks_No; % [-] Entire number of cells
A_PEM= 180; % [cm^2] Membrane sectional area in cm^2
L_PEM = 0.125; % [cm] Membrane thickness
stack_t_gdl_A = 25; % [um] Anode gas diffusion layer thickness
stack_t_gdl_C = 250; % [um] Cathode gas diffusion layer thickness
stack_mea_rho = 1800; % [kg/m^3] Overall density of membrane electrode assembly
stack_mea_cp = 870; % [J/(kg*K)] Overall specific heat of membrane electrode assembly
J_cell_0=8e-05;%[A/cm^2] % Exchange current density
T_ref_st=273.15+60; %[k] Electrolyzer operational temperature reference 
h_t=(9.4+4.7)*(L_PEM*1e4 + stack_t_gdl_A + stack_t_gdl_C)*1e-6*sqrt(A_PEM*1e-4)*stack_num_cells+4*5.9*A_PEM*1e-4; % Entire stacks assembly heat convection factor
R_th_cell=1/h_t; % [K/W]% thermal resistance of the electrolysis stacks assembly 
C_P_cell=stack_mea_rho * A_PEM*1e-4 * (L_PEM*1e4 + stack_t_gdl_A + stack_t_gdl_C)*1e-6 * stack_num_cells*stack_mea_cp;% [J/K] thermal capacitance of the Entire stacks assembly
M_membrane=1.1; % kg/mol Equivalent molar weight of dry membrane
rho_membrane =2000;% kg/m^3 % Density of dry membrane
tau_e=0.1; % molar flow response delay with current variations
C_dl=0.051000; %[F]  double layer charging effect capacitance
% ==============================================================================================
% Phyiscal Properties of gas species and liquid water with temperature
% variation
T_TLU = [-56.55, -50:10:-10, -5:1:5, 10:10:350]';%, 'degC'}; % Temperature varition vector for lookup table of gases properties

% Hydrogen property tables
H2_R = 4124.48151675695; % [J/(kg*K)] Specific gas constant
H2_h = [2783.66044045879, 2873.94181645932, 3012.62601981068, 3152.22537301871, 3292.62486344326, 3433.72138359680, 3504.50184295996, 3518.67512520997, 3532.85394543712, 3547.03822323045, 3561.22787913325, 3575.42283463482, 3589.62301216223, 3603.82833507219, 3618.03872764290, 3632.25411506595, 3646.47442343830, 3717.64725530348, 3860.32198734890, 4003.38288163040, 4146.77354613486, 4290.44463593559, 4434.35318483605, 4578.46197843019, 4722.73896833090, 4867.15672726297, 5011.69194456718, 5156.32496143709, 5301.03934494342, 5445.82149962224, 5590.66031514150, 5735.54684833147, 5880.47403768083, 6025.43644826542, 6170.43004499140, 6315.45199199528, 6460.50047604533, 6605.57455182556, 6750.67400704912, 6895.79924543541, 7040.95118568927, 7186.13117473630, 7331.34091359044, 7476.58239435534, 7621.85784698693, 7767.16969456766, 7912.52051596217, 8057.91301483815, 8203.34999414341, 8348.83433523064, 8494.36898091475, 8639.95692183309]; % [kJ/kg] Specific enthalpy vector with temperature variation
% H2_D = 74; % [mm^2/s] Diffusivity
C_Pm_H2=14.30;%KJ/(kg.K) heat capacity of H2 gas

% Oxygen property tables
O2_R = 259.836612622973; % [J/(kg*K)] Specific gas constant
O2_h = [196.314045635230, 202.302438242727, 211.445560194021, 220.590899228792, 229.740231593225, 238.895348952500, 243.475641264001, 244.391947430149, 245.308340402054, 246.224821996439, 247.141394030059, 248.058058319568, 248.974816681386, 249.891670931563, 250.808622885637, 251.725674358497, 252.642827164235, 257.230174605453, 266.413508632715, 275.609852765421, 284.820965750972, 294.048557917965, 303.294277512928, 312.559698673577, 321.846311315623, 331.155513043735, 340.488603075706, 349.846778083808, 359.231129801281, 368.642644208585, 378.082202097877, 387.550580810768, 397.048456949970, 406.576409877154, 416.134925824813, 425.724402467528, 435.345153816409, 444.997415318718, 454.681349062211, 464.397049000008, 474.144546126762, 483.923813550211, 493.734771414086, 503.577291638585, 513.451202453549, 523.356292706982, 533.292315937896, 543.258994207727, 553.256021688858, 563.283068012223, 573.339781378746, 583.425791441447]; % [kJ/kg] Specific enthalpy vectdor with temperature variation
% O2_D = 18; % [mm^2/s] Diffusivity
C_Pm_O2=0.918;%KJ/(kg.K) heat capacity of O2 gas

% Water vapor saturation pressure table with temperature varition
p_ws_TLU = [
        1.71168953425982e-06
        3.93770601979424e-06
        1.28411717710344e-05
        3.80051394873809e-05
        0.000103239029002090
        0.000259873810798063
        0.000401741022116384
        0.000437454836117844
        0.000476040820618202
        0.000517704666185737
        0.000562664883741207
        0.000611212677444345
        0.000657088049117899
        0.000705987905898279
        0.000758082381149132
        0.000813549384183233
        0.000872574861129522
        0.00122818386934022
        0.00233921476677690
        0.00424668834054807
        0.00738442748706953
        0.0123512704340234
        0.0199458019246787
        0.0312006356960619
        0.0474147199263783
        0.0701823607447713
        0.101417977921310
        0.143375967241115
        0.198665399739302
        0.270259606559999
        0.361500961984849
        0.476101381081492
        0.618139196722055
        0.792053183687694
        1.00263456881210
        1.25501792086105
        1.55467186826983
        1.90739066433297
        2.31928772772542
        2.79679245576864
        3.34665187151016
        3.97593907083533
        4.69207105435535
        5.50283947408830
        6.41645928168037
        7.44164254361690
        8.58770832955728
        9.86474556030261
        11.2838558865497
        12.8575218898034
        14.6001810568052
        16.5291642526045];% [MPa]  Water vapor saturation pressure vector


% Water vapor enthalpy:
h_w_vap_TLU=[2836.88241275372; 2837.81392500514; 2838.63937175807; 2838.7309929628; 2838.06905313927; 2836.62597341095; 2835.60023952573; 2835.37006381376; 2835.13143158077; 2834.88429145066; 2834.62859062897; 2500.93420564316; 2498.55329907119; 2496.17495082036; 2493.79885912205; 2491.42474723191; 2489.05236104642; 2477.20875029194; 2453.54955988604; 2429.83856603313; 2406.00136954922; 2381.97406342174; 2357.69101156389; 2333.08088387814; 2308.06565480412; 2282.56034405021; 2256.47287422313; 2229.70428017508; 2202.14968030993; 2173.69998896199; 2144.24368406447; 2113.66758247452; 2081.85585110571; 2048.68725710494; 2014.03141249899; 1977.7449923284; 1939.66849592886; 1899.62342931686; 1857.40930169838; 1812.79975416601; 1765.53721731251; 1715.32525772152; 1661.81700477787; 1604.59703745016; 1543.1534616334; 1476.83692476464; 1404.80240352424; 1325.92091736114; 1238.61667822208; 1140.5102987018; 1027.62017777647; 892.733785613825];%kJ/kg evaporation enthalpy
h_w_TLU=[2396.55944251649; 2408.68643343608; 2427.1988031141; 2445.702165897; 2464.18429108356; 2482.62529466839; 2491.82135326629; 2493.6580151792; 2495.49372578218; 2497.32843835745; 2499.16210446923; 2500.99462758899; 2502.83092214066; 2504.6665223621; 2506.50140563013; 2508.335548891; 2510.16892865793; 2519.32352241995; 2537.56068088674; 2555.67742615292; 2573.6403223998; 2591.4109179101; 2608.9455037701; 2626.19492016262; 2643.10444358428; 2659.61377634497; 2675.58278696023; 2696.1846256545; 2716.49989553741; 2736.6235210957; 2756.61216047011; 2776.50561773853; 2796.33378922508; 2816.11979583049; 2835.88181423877; 2855.63430574511; 2875.38889327595; 2895.15500577172; 2914.94035964689; 2934.75132290228; 2954.59319330244; 2974.47041285505; 2994.38673458964; 3014.34535328009; 3034.34900867103; 3054.40006755764; 3074.50058946898; 3094.65237953784; 3114.85703128106; 3135.11596137801; 3155.43043805905; 3175.80160435813];%kJ/kg water vapor enthalpy
C_Pm_H2O_lid=4.1816;%[KJ/(kg.K)] heat capacity of H2O liquid
C_Pm_H2O_gas=1.87;%[KJ/(kg.K)] heat capacity of H2O gas
rho_H2O=997;%[Kg/m^3] liquid water density

% ==============================================================================================
% Universal constants
F     = 96485.33212; %     [C/mol]  ; % Faraday constant
R     = 8.31446261815324;%, [J/K/mol] ; % Universal gas constant
v_H2O=1; %[-] stoichiometric coefficient of water vapor in the electrolysis
v_H2=1; %[-] stoichiometric coefficient of Hydrogen gas in the electrolysis
v_O2=0.5; %[-] stoichiometric coefficient of Oxygen gas in the electrolysis
R_w=461.523;%[J/Kg/K] water vapor gas constant
MW_H2O = R/R_w; %[kg/mols] Molar mass of water vapor
MW_O2  = R/O2_R;%[kg/mols] Molar mass of Oxygen gas
MW_H2  = R/H2_R;%[kg/mols] Molar mass of Hydrogen gas
g=9.8; %[N/kg] free-fall acceleration


% ==============================================================================================
% Oxygen Seperator vessel (OSV) parameters
A_SEPO=0.25;%[m^2] %  vessel base area
L_SEPO_h=0.5;%[m] % Height of vessel
L_SEPO_min=0.25;%[m] % minimum liquid water level inside the vessel
V_SEPO=A_SEPO*L_SEPO_h; % Vessel volume 
max_flow_tank=0.05;% [kg/s] % max inflow in water pumps (OSV/HSV/ reciruclation circuits)
% the PI is set to saturate at low value of 0.25 

% valve flow resistance between OSV and atmosphere
R_SEP_AND=50000;
R_cond_SEPO=0.167470259285884;%[K/W] %thermal resistance of the OSV with the enclosure
C_th_SEPO_v=45.902149439227564;%[kJ/k] %Thermal capacitance of void OSV
w_eva=100; % water evaporation filter cutt-off frequency
wc_RES=1/1;
KpL=5000*4;
KiL=1000*4;
% ==============================================================================================
% Hydrogen seperator vessel (HSV) parameters
A_SEPH=0.25;%[m^2] % vessel base area
L_SEPH_h=0.5;%[m] % Height of vessel
L_SEPH_min=0.25;%[m] minimum liquid water level inside the vessel
V_SEPH=A_SEPH*L_SEPH_h; % Vessel volume
% valve flow resistance between HSV and atmosphere
R_SEP_CAT=50000;
C_th_SEPH_v=C_th_SEPO_v; %Thermal capacitance of void HSV
R_cond_SEPH=R_cond_SEPO;
Kp_rsh=KpL;
Ki_rsh=KiL;
wc_RSH=wc_RES;
% ==============================================================================================
% Heat exchanger and anode reciruclation circuit
exchanger_L = 1; % [m] Overall radiator length
exchanger_W = 0.025; % [m] Overall radiator width
exchanger_N_tubes = 25; % [-] Number of coolant tubes
exchanger_tube_H = 0.0015; % [m] Height of each coolant tube
R_ech=0.357772792810198;%[K.s/kJ]% Heat Exchanger thermal resistance
C_th_RECO=0.51*0.001*(3.14*0.025)*0.5*7900+3.14*(0.025/2)^2*0.5*rho_H2O*C_Pm_H2O_lid+...
    C_Pm_H2O_lid*rho_H2O*exchanger_L*exchanger_W*exchanger_tube_H*exchanger_N_tubes;%[kJ/K]; %  thermal capacitance of anode recirculation
R_cond_RECO=3.397027600849257;% thermal resistance between recircuilation circuit and enclosure
C_th_COOL=0.51*0.001*(3.14*0.025)*0.35*7900+3.14*(0.025/2)^2*0.35*rho_H2O*C_Pm_H2O_lid+...
    C_Pm_H2O_lid*rho_H2O*exchanger_L*exchanger_W*exchanger_tube_H*exchanger_N_tubes;%[kJ/K]; %  thermal capacitance of coolant circuit
C_th_COLD=0.51*0.003*6*(0.3*0.3)*7900+0.3*0.3*0.3*rho_H2O*C_Pm_H2O_lid;
max_flow_and=0.05;%[kg/s] % max flow of the reciruclation circuit
cooling_flow=0.2;%[kg/s] 
wc_RECO=1/1;
KpRECO=0.4*4; % proportional gain of the PI controller governing the recirculation flow
KiRECO=0.05*4; % interal gain of the PI controller governing the recirculation flow
% ==============================================================================================
% Cathode reciruclation circuit
C_th_RECH=0.51*0.001*(3.14*0.025)*0.35*7900+3.14*(0.025/2)^2*0.35*rho_H2O*C_Pm_H2O_lid;%kJ/k;
R_cond_RECH=4.852896572641796; % thermal resistance of the recirculation circuit at the cathode side
% ==============================================================================================
% Purification unit
V_DRY=1e-2; % volume of purification unit
R_cond_DRY=0.441482518449157;%K/W thermal resistance between purification unit and unit enclosure
C_th_DRY=3.366173053912660;%kJ/K thermal capacitance of the purification unit
ValvArea=7.0686e-05;% This value represent the maximum area of the valve seperating the Hydrogen storage from the purifier
R_Dry_Fix=1.399428921390710;
w_ads=1000;% Fast condesation rate of time constant of 1e-3 
% ==============================================================================================
% Enclosure of the unit
C_p_air=1.012;%kJ/(kg.K) % air density
C_th_Enc=1.007250000000000e+02; % Enclosure thermal capacitance
R_cond_enc=0.037313432835821; % Enclosure thermal resistance with the atmosphere
m_dot_fan=0.769300000000000; % air flow of the fan
% ==============================================================================================
% Storage tank control
wc_STRG=1/0.1; % valve dynamics filter
KpSTRG=0.1e-6; % proportional gain of the storage tank valve control
KiSTRG=0.1e-6; % Integral gain of the storage tank valve control


%% PEIs parameters:
wc=25; %power measurement filter
vnom=380;% nominal voltage
Vdc=622.26;% dc link reference voltage
Rf=0.1; % LC filter resistance 
Lf=1.35e-3;% LC filter inductance 
Cf=50e-6; % LC filter capacitance
Kpi=10.5; % ac inverter current controller proportional gain 
Kii=16000;% ac inverter current controller integral gain 
Kiv=390;  % ac inverter voltage controller integral gain
Kpv=0.05; % ac inverter voltage controller proportional gain
wnom=2*pi*50; % ac nominal rotational frequency of the grid forming operation
H=0.75; % feed forward gain for the voltage controller
Kdcvp=0.85*50; % DC link voltage controller proportional gain
Kdcvi=50; % DC link voltage controller intgeral gain
Kdcip=0.5; % DC link current controller proportional gain
Kdcii=0.2; % DC link current controller intgeral gain
Ls=0.5e-3; % DC buck converter inductance
fsw=8000; % switching frequency
Rc=0.03; % coubling line resistance
Lc=0.35e-3; % coubling line inductance
Tsim=1e-6; % simulation time step
m1=1.33e-4; % active power droop gain
n1=1.33e-3; % reactive power droop gain
nd1=n1*0.1; % reactive power derivative droop gain
md1=m1*0.1; % reactive power droop gain
rn=10000; % virtual resistance
C_el=1e-3; % High voltage DC link capacitance
Cs=1e-2; % Low voltage DC capacitance
R_load=10.62; % grid load resistance
L_load=0.02032*5; % grid load inductance


%% Operating Point Initialization
% ==============================================================================================
% This part of the script is developed to enable the script to obtain the
% operating point details from a running nonlinear time domain simulation.
% The nonlinear time domain model originally used in the paper is preserved
% due to funding restrictions. A sample operating point is provided in the
% github link eo enable running the manuscript. The sample operating point
% is found in the Github link (https://github.com/bhamad97/PEMELSSM.git)
% as matlab file OperatingPointGithub.mat
% ==============================================================================================

load("OperatingPointGithub.mat"); % the file is supposed to be in the current directory

ssp_I_cell=out.I_cell.Data(end);% cell currebt
ssp_T_stack=out.T_stack.Data(end); % stack temperature
J_0_ref=J_cell_0; % exchange current density
ssp_I_cell_0=J_0_ref*A_PEM; 
ssp_R_ohm=out.R_mem.Data(end); % ohmic resistance
ssp_E_stk=out.E_stack.Data(end); % Stack voltage
ssp_E_cell=ssp_E_stk/cells_per_stack; % cell voltage
ssp_E_rev=out.E_rev.Data(end); % reverisble voltage of the cell
ssp_E_act=out.E_act.Data(end); % activation overpotential
ssp_E_ss=R/F*ssp_T_stack*log(ssp_I_cell/ssp_I_cell_0); % activation overpotential steady state
ssp_P_SEPH_H2=out.P_SEPH_H2.Data(end); % hydrogen partial pressure in HSV
ssp_P_SEPH_H2O=out.P_SEPH_H2O.Data(end); % water gas partial pressure in HSV
ssp_P_SEPH=ssp_P_SEPH_H2+ssp_P_SEPH_H2O; % HSV pressure 
ssp_m_H2O_SEPH=out.m_SEPH_H2O.Data(end); % water vapor mass in HSV
ssp_m_H2_SEPH=out.m_SEPH_H2.Data(end); % Hydrogen mass in HSV
ssp_sigma_PEM=out.sigma_PEM.Data(end); % conductivity of the membrane
ssp_P_SEPO_H2O=out.p_SEPO_H2O.Data(end); % water vapor partial pressure in OSV
ssp_m_H2O_SEPO=out.m_SEPO_H2O.Data(end); % water vapor mass in OSV
ssp_m_O2_SEPO=out.m_SEPO_O2.Data(end); % Oxygen gas mass in OSV
ssp_P_SEPO_O2=out.p_SEPO_O2.Data(end); % Oxygen partial pressure in OSV
ssp_P_SEPO=ssp_P_SEPO_O2+ssp_P_SEPO_H2O; % OSV pressure
ssp_a_SEPO_H2O=out.R_HO.Data(end); % Humidity of OSV
ssp_P_OXYG=env_p*1e6; % Atmosphere pressure
ssp_H_H2O_0=h_w_TLU(12);% water vapor enthalpy at zero degree
ssp_H_H2_0=H2_h(12);% Hydrogen gas enthalpy at zero degree
ssp_H_O2_0=O2_h(12);% Oxygen gas enthalpy at zero degree
ssp_H_H_2_g=C_Pm_H2*(ssp_T_stack-273.15)+ssp_H_H2_0;% Hydrogen gas enthalpy
ssp_H_O_2_g=C_Pm_O2*(ssp_T_stack-273.15)+ssp_H_O2_0;% Oxygen gas enthalpy
ssp_H_H_2O_g=C_Pm_H2O_gas*(ssp_T_stack-273.15)+ssp_H_H2O_0;% water vapor enthalpy
ssp_m_d_RECH=out.m_dot_RECH.Data(end); % cathode recirculation water flow rate
ssp_m_d_RECO=out.m_dot_RECO.Data(end); % anode recirculation water flow rate
ssp_m_d_OXYG=out.m_dot_OXYG.Data(end); % gas mass flow rate leaving the OSV to the atmosphere
ssp_m_d_SPVA=out.m_dot_SPVA.Data(end); % gas mass flow rate leaving the OSV to the purification unit
ssp_m_d_eos=out.m_dot_H2O_eos.Data(end); % electro-osmotic drag flow rate
ssp_m_d_diff=out.m_dot_H2O_diff.Data(end); % diffusion flow rate
ssp_T_RECH=out.T_RECH.Data(end); % Cathode recirculation temperature
ssp_T_RECO=out.T_RECO_td.Data(end); % anode recirculation temperature
ssp_T_ENC=out.T_Encl.Data(end); % Enclosure temperature
ssp_T_RES=env_T+273.15; % water reservoir temperature at the OSV side
ssp_T_RSH=env_T+273.15; % water reservoir temperature at the HSV side
ssp_xi=out.xi.Data(end); % cell molar flow rate
ssp_T_SEPO=out.T_SEPO.Data(end); % OSV temperature
ssp_L_SEPO=out.L_SEPO.Data(end); % OSV water level
ssp_T_SEPH=out.T_SEPH.Data(end); % HSV temperature
ssp_L_SEPH=out.L_SEPH.Data(end); % HSV water level
ssp_Htot_SEPO=out.Htot_SEPO_td.Data(end); % total enthalpy flow into OSV
ssp_C_th_SEPO=out.C_th_SEPO_td.Data(end); % thermal capacitance of the OSV
ssp_Htot_SEPH=out.Htot_SEPH_td.Data(end); % total enthalpy flow into HSV
ssp_C_th_SEPH=out.C_th_SEPH_td.Data(end); % thermal capacitance of the HSV
ssp_P_SPVA_H2O=out.P_SPVA_H2O.Data(end); % water vapor partial pressure in purification unit
ssp_P_SPVA_H2=out.P_SPVA_H2.Data(end); % Hydrogen partial pressure in purification unit
ssp_P_SPVA=ssp_P_SPVA_H2+ssp_P_SPVA_H2O; % pressure of the purification unit
P_PURI=env_p*1e6; % pressure of the storage unit
ssp_m_d_PURI=out.m_dot_PURI.Data(end); % mass flow to the storage unit
ssp_m_H2_DRY=out.m_DRY_H2.Data(end); % Hydrogen gas mass in the purification unit
ssp_m_H2O_DRY=out.m_DRY_H2O.Data(end); % Water vapor mass in the purification unit
ssp_T_DRY=out.T_DRY.Data(end); % Purification unit temperature
ssp_m_d_ads=out.ssp_m_d_ads.Data(end); %  adsorpation rate
ssp_A_rest=out.A_rest.Data(end); % storage valve opeing area
ssp_Psat=out.P_sat_SEPO.Data(end); % Saturation pressure
ssp_m_d_evo=out.m_dot_evaO.Data(end); % evaporation pressure in the OSV
ssp_m_d_evh=out.m_dot_evaH.Data(end); % evaporation pressure in the HSV
ssp_Dw=out.D_H2O_m.Data(end); % water diffusion factor
ssp_CH2O_chh=out.Conc_H2O_ccl.Data(end); % water vapor concentration at the anode electrolyte
ssp_CH2O_cho=out.Conc_H2O_acl.Data(end); % water vapor concentration at the cathode electrolyte
lambda=out.Lambda.Data(end); % water content of the membrane
ssp_PI_L_SEPO=out.PI_L_SEPO.Data(end); % OSV water level PI controller output
ssp_PI_L_SEPH=out.PI_L_SEPH.Data(end); % HSV water level PI controller output
ssp_PI_T_SEPO=out.PI_T_SEPO.Data(end); % Temperature level PI controller output
ssp_PI_A_rest=out.PI_A_rest.Data(end); % Valve opening PI controller output
% the following parameters describe if the corresponding PI controller is
% active or reached its saturation
ac_m_rso=double(ssp_PI_L_SEPO>0 && ssp_PI_L_SEPO<1); % activity of OSV water level PI controller
ac_m_rsh=double(ssp_PI_L_SEPH>-1 && ssp_PI_L_SEPH<1); % activity of HSV water level PI controller
ac_T_spo=double(ssp_PI_T_SEPO>0 && ssp_PI_T_SEPO<1); % activity of temperature PI controller
ac_A_rst=double(ssp_PI_A_rest>0 && ssp_PI_A_rest<0.9); % activity of Valve opening PI controller
ssp_T_COOL=out.T_COOL_td.Data(end); % Coolant temperature
ssp_T_COLD=out.T_COLD_td.Data(end); % cooling unit temperature
Dbar=ssp_E_stk/Vdc; % DC buck duty cycle 
Iod1=out.I_od.Data(end); % direct component of the VSC output current 
Ioq1=out.I_oq.Data(end); % quadrature component of the VSC output current 
Vod1=out.V_od.Data(end); % direct component of the VSC output voltage
Voq1=0; % quadrature component of the VSC output voltage
Pssp=Vod1*Iod1; % VSC output active power
Qssp=-Vod1*Ioq1; % VSC output reactive power
Ild1=out.I_ld.Data(end); % LC filter d-axis current
Ilq1=out.I_lq.Data(end); % LC filter q-axis current
w0=wnom-m1*Pssp; % VSC angular frequency
Iloadq=out.I_loadq.Data(end); % load quadrature current
Iloadd=out.I_loadd.Data(end); % load direct current
Vid1=Vod1+Rf*Ild1-w0*Ilq1*Lf; % Inverter terminal d-axis voltage
Viq1=Voq1+Rf*Ilq1+w0*Ild1*Lf; % Inverter terminal q-axis voltage
Vbd1=Vod1+w0*Ioq1*Lc-Rc*Iod1; % Point of connection d-axis voltage
Vbq1=Voq1-w0*Iod1*Lc-Rc*Ioq1; % Point of connection q-axis voltage
Pinvbar=-ssp_E_stk*Stacks_No*ssp_I_cell; % power fed to the inverter by the DC link
Ib=out.i_b.Data(end); % DC buck low voltage side current feeding the electrolyzer stack assembly  

% the operating point of the system is given as follows
Lin_op=[ssp_E_act;ssp_xi;ssp_T_stack;ssp_T_SEPO;ssp_T_SEPH;ssp_T_RECO;ssp_T_COOL;ssp_T_COLD;ssp_T_RECH;ssp_T_DRY;ssp_T_ENC;ssp_m_H2O_SEPO;ssp_m_O2_SEPO;ssp_m_H2_SEPH;ssp_m_H2O_SEPH;ssp_m_H2_DRY;ssp_m_H2O_DRY;ssp_L_SEPO*rho_H2O*A_SEPO;ssp_L_SEPH*rho_H2O*A_SEPH;ssp_PI_L_SEPO*max_flow_tank;ssp_PI_L_SEPH*max_flow_tank;ssp_m_d_RECO;ssp_A_rest;0;0;0;ssp_E_stk;Ib;Vdc;0;0;Pssp;Qssp;0;0;0;0;Ild1;Ilq1;Vod1;Voq1;Iod1;Ioq1;Iloadd;Iloadq];
%==========================================================================
%% Linearization elements of the PEMEL unit
%==========================================================================
% The elements of the linearization are provided according to details in
% the main manuscript and the supplementary document. 
%==========================================================================
% Electrolysis Stack Subsystem
%==========================================================================
% PEM Electrochemical subsystem linearization elements
z1=1/C_dl*(1-ssp_E_act/ssp_E_ss);
z2=-ssp_I_cell/C_dl/ssp_E_ss;
z3=ssp_I_cell/C_dl*ssp_E_act/ssp_E_ss^2;
z4=R/F*log(ssp_I_cell/ssp_I_cell_0); 
z5=R/F*ssp_T_stack/ssp_I_cell; 
z6=1/ssp_R_ohm;
z7=-1/ssp_R_ohm;
z8=-1/ssp_R_ohm;
z9=-(ssp_E_cell-ssp_E_rev-ssp_E_act)/(ssp_R_ohm^2);
z10=(1.5421*10^(-3)+9.523*10^(-5)*(log(ssp_T_stack)+1)+2*9.84e-8*ssp_T_stack);
z11=R/2/F*log(ssp_P_SEPH_H2*sqrt(ssp_P_SEPO_O2)/ssp_a_SEPO_H2O/(env_p*1e6)^1.5); 
z12=R*ssp_T_stack/2/F/ssp_P_SEPH_H2; 
z13=R*ssp_T_stack/4/F/ssp_P_SEPO_O2; 
z14=-R*ssp_T_stack/2/F/ssp_P_SEPO_H2O; 
z15=0.005139*exp(1268*(1/303-1/ssp_T_stack)); 
z16=((0.005139*lambda-0.00326)*exp(1268*(1/303-1/ssp_T_stack))*1268/ssp_T_stack^2); 
z17=0.5*(17.81/(ssp_Psat)-2*39.85*ssp_P_SEPO_H2O/(ssp_Psat)^2+3*36*ssp_P_SEPO_H2O^2/(ssp_Psat)^3);
z18=0.5*(17.81/(ssp_Psat)-2*39.85*ssp_P_SEPH_H2O/(ssp_Psat)^2+3*36*ssp_P_SEPH_H2O^2/(ssp_Psat)^3);
z19=-L_PEM/(A_PEM*ssp_sigma_PEM^2);

a_pem_e_21=z8/(2*F*tau_e);
a_pem_e_22=-1/tau_e;
b_pem_e_21=1/(2*F*tau_e)*(z7*(z10+z11)+z9*(z19*z16));
b_pem_e_22=z7/(2*F*tau_e)*z12;
b_pem_e_23=z7/(2*F*tau_e)*z13;
b_pem_e_24=1/(2*F*tau_e)*(z7*z14+z9*(z19*(z15*z17)));
b_pem_e_25=1/(2*F*tau_e)*z9*(z19*(z15*z18));
b_pem_e_26=z6/(2*F*tau_e);
c_pem_e_21=z8;
d_pem_e_21=b_pem_e_21*(2*F*tau_e);
d_pem_e_22=b_pem_e_22*(2*F*tau_e);
d_pem_e_23=b_pem_e_23*(2*F*tau_e);
d_pem_e_24=b_pem_e_24*(2*F*tau_e);
d_pem_e_25=b_pem_e_25*(2*F*tau_e);
d_pem_e_26=b_pem_e_26*(2*F*tau_e);
a_pem_e_11=z2+z1*z8+z3*(z5*z8);
b_pem_e_11=z3*(z4+z5*(z7*(z10+z11)+z9*(z19*z16)))+z1*(z7*(z10+z11)+z9*(z19*z16));
b_pem_e_12=z3*z5*z7*z12+z1*z7*z12;
b_pem_e_13=z3*z5*z7*z13+z1*z7*z13;
b_pem_e_14=z3*z5*(z7*z14+z9*z19*(z15*z17))+z1*(z7*z14+z9*z19*(z15*z17));
b_pem_e_15=z3*z5*z9*(z19*(z15*z18))+z1*z9*(z19*(z15*z18));
b_pem_e_16=z3*z5*z6+z1*z6;
%==========================================================================
% Stack thermal linearization elements:
e1=1/1000/R_th_cell; 
e2=0.0134*stack_num_cells*MW_H2O*ssp_xi;
e3=stack_num_cells*MW_H2O*(0.0134*ssp_T_stack+0.03);
e4=-(ssp_T_stack-273.15)*C_Pm_H2O_lid;
e5= C_Pm_H2O_gas*(ssp_T_stack-273.15)+ssp_H_H2O_0;
e6=C_Pm_H2O_lid*(ssp_m_d_RECO-ssp_m_d_eos-ssp_xi*MW_H2O*v_H2O*stack_num_cells);
e7=C_Pm_H2O_gas*ssp_m_d_diff;
e8=C_Pm_O2*ssp_xi*MW_O2*v_O2*stack_num_cells;
e9=C_Pm_H2O_lid*(ssp_m_d_RECH-ssp_m_d_diff);
e10=C_Pm_H2O_gas*ssp_m_d_eos;
e11=C_Pm_H2*ssp_xi*MW_H2*v_H2*stack_num_cells;
a_pem_t_11=1000/C_P_cell*(-e9-e10-e5*e2-e6-e7-e4*e2-e1-e11-e8);
b_pem_t_18=1000/C_P_cell*C_Pm_H2O_lid*ssp_m_d_RECH;
b_pem_t_19=1000/C_P_cell*C_Pm_H2O_lid*ssp_m_d_RECO;
b_pem_t_110=1/C_P_cell/R_th_cell;
b_pem_t_11=1000/C_P_cell*C_Pm_H2O_lid*(ssp_T_RECH-ssp_T_stack);
b_pem_t_12=1000/C_P_cell*C_Pm_H2O_lid*(ssp_T_RECO-ssp_T_stack);
b_pem_t_13=-1000/C_P_cell*e5-1000/C_P_cell*e4;
b_pem_t_15=1000*stack_num_cells/C_P_cell*(-e4*MW_H2O*v_H2O-v_H2*MW_H2*ssp_H_H_2_g-v_O2*MW_O2*ssp_H_O_2_g)-1000/C_P_cell*e3*(e4+ssp_H_H_2O_g);
b_pem_t_16=stack_num_cells/C_P_cell*ssp_E_cell;
b_pem_t_17=stack_num_cells/C_P_cell*ssp_I_cell;
%==========================================================================
% Stack gaseous diffusion linearization elements
d_pem_g_1=-stack_num_cells*ssp_Dw*A_PEM*0.0001*MW_H2O*rho_membrane/(0.01*L_PEM)/M_membrane*z17;
d_pem_g_2=stack_num_cells*ssp_Dw*A_PEM*0.0001*MW_H2O*rho_membrane/(0.01*L_PEM)/M_membrane*z18;
d_pem_g_3=stack_num_cells*A_PEM*0.0001*MW_H2O/(0.01*L_PEM)*(ssp_CH2O_chh-ssp_CH2O_cho)*1.25e-10*exp(2416/303.15-2416/ssp_T_stack)/ssp_T_stack^2*2416;
%==========================================================================
% small signal subsystem of the Oxygen separator vessel
%==========================================================================
% gaseous domain linearization elements
a_SEPO_g_11=-ssp_m_O2_SEPO*ssp_m_d_OXYG/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2;
a_SEPO_g_12=ssp_m_H2O_SEPO*ssp_m_d_OXYG/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2;
a_SEPO_g_21= ssp_m_O2_SEPO*ssp_m_d_OXYG/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2;
a_SEPO_g_22=-(ssp_m_H2O_SEPO)*ssp_m_d_OXYG/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2;
b_SEPO_g_13=-ssp_m_H2O_SEPO/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO);
b_SEPO_g_23=-ssp_m_O2_SEPO/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO);
b_SEPO_g_24=stack_num_cells*MW_O2*v_O2;
c_SEPO_g_11=R*ssp_T_SEPO/(MW_H2O*(V_SEPO-A_SEPO*ssp_L_SEPO));
c_SEPO_g_22=R*ssp_T_SEPO/(MW_O2*(V_SEPO-A_SEPO*ssp_L_SEPO));
d_SEPO_g_11=R*ssp_T_SEPO*ssp_m_H2O_SEPO/(MW_H2O*(V_SEPO-A_SEPO*ssp_L_SEPO)^2*rho_H2O);
d_SEPO_g_21=R*ssp_T_SEPO*ssp_m_O2_SEPO/(MW_O2*(V_SEPO-A_SEPO*ssp_L_SEPO)^2*rho_H2O);
d_SEPO_g_12=R*ssp_m_H2O_SEPO/(MW_H2O*(V_SEPO-A_SEPO*ssp_L_SEPO));
d_SEPO_g_22=R*ssp_m_O2_SEPO/(MW_O2*(V_SEPO-A_SEPO*ssp_L_SEPO));
c_SEPO_g_33=-w_eva;
d_SEPO_g_31=-w_eva*ssp_Psat*MW_H2O/R/ssp_T_SEPO/rho_H2O;
d_SEPO_g_32=-w_eva*ssp_Psat*MW_H2O/R/(ssp_T_SEPO^2)*(V_SEPO-A_SEPO*ssp_L_SEPO);
%==========================================================================
% thermal domain linearization elements:
b_SEPO_t_11=1/ssp_C_th_SEPO*(e6+e4*e2+e7+e8);
a_SEPO_t_11=-1/ssp_C_th_SEPO*(ssp_m_d_RECO*C_Pm_H2O_lid+0.001/R_cond_SEPO+ssp_m_d_OXYG*(ssp_m_O2_SEPO*C_Pm_O2+ssp_m_H2O_SEPO*C_Pm_H2O_gas)/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)+ssp_m_d_evo*(C_Pm_H2O_gas-C_Pm_H2O_lid));
b_SEPO_t_12=0.001/ssp_C_th_SEPO/R_cond_SEPO;
b_SEPO_t_13=1/ssp_C_th_SEPO*ssp_T_stack*C_Pm_H2O_lid-1/ssp_C_th_SEPO*ssp_T_SEPO*C_Pm_H2O_lid;
b_SEPO_t_14=1/ssp_C_th_SEPO*(ssp_T_RES-273.15)*C_Pm_H2O_lid;
b_SEPO_t_15=-1/ssp_C_th_SEPO*((ssp_m_O2_SEPO*C_Pm_O2+ssp_m_H2O_SEPO*C_Pm_H2O_gas)*(ssp_T_SEPO-273.15)+ssp_m_O2_SEPO*ssp_H_O2_0+ssp_m_H2O_SEPO*ssp_H_H2O_0)/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO);
b_SEPO_t_18=1/ssp_C_th_SEPO*e4*(e3+stack_num_cells*MW_H2O*v_H2O)+1/ssp_C_th_SEPO*ssp_H_O_2_g*stack_num_cells*MW_O2*v_O2;
b_SEPO_t_16=1/ssp_C_th_SEPO*e5;
b_SEPO_t_19=1/ssp_C_th_SEPO*ssp_m_d_OXYG*ssp_m_O2_SEPO/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2*(C_Pm_O2*(ssp_T_SEPO-273.15)+ssp_H_O2_0-C_Pm_H2O_gas*(ssp_T_SEPO-273.15)-ssp_H_H2O_0);
b_SEPO_t_110=+1/ssp_C_th_SEPO*ssp_m_d_OXYG*ssp_m_H2O_SEPO/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2*(C_Pm_H2O_gas*(ssp_T_SEPO-273.15)+ssp_H_H2O_0-C_Pm_O2*(ssp_T_SEPO-273.15)-ssp_H_O2_0);
b_SEPO_t_111=0;
b_SEPO_t_17=-1/ssp_C_th_SEPO*(-C_Pm_H2O_lid*(ssp_T_SEPO-273.15)+C_Pm_H2O_gas*(ssp_T_SEPO-273.15)+ssp_H_H2O_0);
%==========================================================================
% liquid fluidic domain linearization elements:
a_SEPO_f_21=-wc_RES*max_flow_tank*KpL/A_SEPO/rho_H2O*ac_m_rso;  
a_SEPO_f_22=-wc_RES;
a_SEPO_f_23=wc_RES*max_flow_tank*KiL*ac_m_rso;
a_SEPO_f_31=-1/A_SEPO/rho_H2O*ac_m_rso;  
b_SEPO_f_11=-MW_H2O*ssp_xi*0.0134*stack_num_cells;
b_SEPO_f_12=-stack_num_cells*MW_H2O*(0.03+0.0134*ssp_T_stack+v_H2O);
d_SEPO_f_11=1/(2*R_SEP_AND*sqrt(ssp_P_SEPO-ssp_P_OXYG))*sign(ssp_m_d_OXYG);
d_SEPO_f_12=d_SEPO_f_11;
%==========================================================================
% small signal subsystem of the Hydrogen separator vessel
%==========================================================================
% gaseous domain linearization elements
a_SEPH_g_11=-ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;
a_SEPH_g_12=ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;
a_SEPH_g_21=ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;
a_SEPH_g_22=-ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;
b_SEPH_g_11=-(ssp_m_H2_SEPH)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);
b_SEPH_g_21=-(ssp_m_H2O_SEPH)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);
b_SEPH_g_14=stack_num_cells*MW_H2*v_H2;
b_SEPH_g_24=stack_num_cells*MW_H2O*(0.03+0.0134*ssp_T_stack);
b_SEPH_g_23=stack_num_cells*MW_H2O*ssp_xi*0.0134;
c_SEPH_g_11=R*ssp_T_SEPH/(MW_H2*(V_SEPH-A_SEPH*ssp_L_SEPH));
c_SEPH_g_22=R*ssp_T_SEPH/(MW_H2O*(V_SEPH-A_SEPH*ssp_L_SEPH));
d_SEPH_g_11=R*ssp_T_SEPH*ssp_m_H2_SEPH/(MW_H2*(V_SEPH-A_SEPH*ssp_L_SEPH)^2*rho_H2O);
d_SEPH_g_21=R*ssp_T_SEPH*ssp_m_H2O_SEPH/(MW_H2O*(V_SEPH-A_SEPH*ssp_L_SEPH)^2*rho_H2O);
d_SEPH_g_12=R*ssp_m_H2_SEPH/(MW_H2*(V_SEPH-A_SEPH*ssp_L_SEPH));
d_SEPH_g_22=R*ssp_m_H2O_SEPH/(MW_H2O*(V_SEPH-A_SEPH*ssp_L_SEPH));
c_SEPH_g_33=-w_eva;
d_SEPH_g_31=-w_eva*ssp_Psat*MW_H2O/R/ssp_T_SEPH/rho_H2O;
d_SEPH_g_32=-w_eva*ssp_Psat*MW_H2O/R/ssp_T_SEPH^2*(V_SEPH-A_SEPH*ssp_L_SEPH);
%==========================================================================
% thermal domain linearization elements:
a_SEPH_t_11=-1/ssp_C_th_SEPH*(ssp_m_d_RECH*C_Pm_H2O_lid+0.001/R_cond_SEPH+ssp_m_d_SPVA*(ssp_m_H2O_SEPH*C_Pm_H2O_gas+ssp_m_H2_SEPH*C_Pm_H2)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH))-1/ssp_C_th_SEPH*(C_Pm_H2O_gas-C_Pm_H2O_lid)*ssp_m_d_evh;
b_SEPH_t_11=1/ssp_C_th_SEPH*(e9+e10+e11+e5*e2);
b_SEPH_t_12=0.001/ssp_C_th_SEPH/R_cond_SEPH;
b_SEPH_t_13=1/ssp_C_th_SEPH*ssp_T_stack*C_Pm_H2O_lid-1/ssp_C_th_SEPH*ssp_T_SEPH*C_Pm_H2O_lid;
b_SEPH_t_14=1/ssp_C_th_SEPH*(ssp_T_RSH-273.15)*C_Pm_H2O_lid;
b_SEPH_t_15=-1/ssp_C_th_SEPH*((ssp_m_H2O_SEPH*C_Pm_H2O_gas+ssp_m_H2_SEPH*C_Pm_H2)*(ssp_T_SEPH-273.15)+ssp_m_H2_SEPH*ssp_H_H2_0+ssp_m_H2O_SEPH*ssp_H_H2O_0)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);
b_SEPH_t_16=-1/ssp_C_th_SEPH*(ssp_T_stack-273.15)*C_Pm_H2O_lid;
b_SEPH_t_17=-1/ssp_C_th_SEPH*(-C_Pm_H2O_lid*(ssp_T_SEPH-273.15)+C_Pm_H2O_gas*(ssp_T_SEPH-273.15)+ssp_H_H2O_0);
b_SEPH_t_18=1/ssp_C_th_SEPH*e5*e3+1/ssp_C_th_SEPH*(ssp_H_H_2_g)*stack_num_cells*MW_H2*v_H2;
b_SEPH_t_19=1/ssp_C_th_SEPH*ssp_m_d_SPVA*ssp_m_H2O_SEPH/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2*(C_Pm_H2O_gas*(ssp_T_SEPH-273.15)+ssp_H_H2O_0-C_Pm_H2*(ssp_T_SEPH-273.15)-ssp_H_H2_0);
b_SEPH_t_110=1/ssp_C_th_SEPH*ssp_m_d_SPVA*ssp_m_H2_SEPH/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2*(C_Pm_H2*(ssp_T_SEPH-273.15)+ssp_H_H2_0-C_Pm_H2O_gas*(ssp_T_SEPH-273.15)-ssp_H_H2O_0);
b_SEPH_t_111=0;
%==========================================================================
% liquid fluidic domain linearization elements:
a_SEPH_f_21=-wc_RSH*max_flow_tank*Kp_rsh/A_SEPH/rho_H2O*ac_m_rsh;  
a_SEPH_f_22=-wc_RSH;
a_SEPH_f_23=wc_RSH*max_flow_tank*Ki_rsh*ac_m_rsh;
a_SEPH_f_31=-1/A_SEPH/rho_H2O*ac_m_rsh;  
d_SEPH_f_11=1/(2*sqrt(abs(ssp_P_SEPH_H2+ssp_P_SEPH_H2O-ssp_P_SPVA_H2O-ssp_P_SPVA_H2)*R_SEP_CAT))*sign(ssp_m_d_SPVA);
d_SEPH_f_12=d_SEPH_f_11;
d_SEPH_f_13=-d_SEPH_f_11;
d_SEPH_f_14=-d_SEPH_f_11;
%==========================================================================
% small-signal subsysem of the RECO
%==========================================================================
% thermal domain linearization elements:
a_RECO_t_11=-1/C_th_RECO*(C_Pm_H2O_lid*ssp_m_d_RECO+0.001/R_cond_RECO+1/R_ech);
a_RECO_t_12=1/C_th_RECO*1/R_ech;
a_RECO_t_21=1/C_th_COOL/R_ech;
a_RECO_t_22=-1/C_th_COOL*(C_Pm_H2O_lid*cooling_flow+1/R_ech);
a_RECO_t_23=1/C_th_COOL*C_Pm_H2O_lid*cooling_flow;
a_RECO_t_32=1/C_th_COLD*C_Pm_H2O_lid*cooling_flow;
a_RECO_t_33=-1/C_th_COLD*C_Pm_H2O_lid*cooling_flow;
b_RECO_t_11=1/C_th_RECO*C_Pm_H2O_lid*ssp_m_d_RECO;
b_RECO_t_12=0.001/C_th_RECO/R_cond_RECO;
b_RECO_t_13=1/C_th_RECO*C_Pm_H2O_lid*(ssp_T_SEPO-ssp_T_RECO);
%==========================================================================
% liquid fluidic domain linearization elements:
a_RECO_f_11=-wc_RECO;
a_RECO_f_12=wc_RECO*max_flow_and*KiRECO*ac_T_spo;
b_RECO_f_11=wc_RECO*max_flow_and*KpRECO*ac_T_spo;
%==========================================================================
% small-signal subsysem of the RECH
%==========================================================================
% thermal domain linearization elements:
a_RECH_t_11=-1/C_th_RECH*(C_Pm_H2O_lid*ssp_m_d_RECH+0.001/R_cond_RECH);
b_RECH_t_11=1/C_th_RECH*C_Pm_H2O_lid*ssp_m_d_RECH;
b_RECH_t_12=0.001/C_th_RECH/R_cond_RECH;
b_RECH_t_13=1/C_th_RECH*C_Pm_H2O_lid*(ssp_T_SEPH-ssp_T_RECH);
%==========================================================================
% small-signal subsystem of the purification unit
%==========================================================================
% thermal domain linearization elements:
a_DRY_t_11=1/C_th_DRY*((C_Pm_H2O_gas-2*C_Pm_H2O_lid)*ssp_m_d_ads+...
                        -ssp_m_d_PURI*(ssp_m_H2O_DRY*C_Pm_H2O_gas+ssp_m_H2_DRY*C_Pm_H2)/(ssp_m_H2O_DRY+ssp_m_H2_DRY)+...
                       -0.001/R_cond_DRY...
                       );
b_DRY_t_11=1/C_th_DRY*ssp_m_d_SPVA*(ssp_m_H2O_SEPH*C_Pm_H2O_gas+ssp_m_H2_SEPH*C_Pm_H2)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);
b_DRY_t_12=0.001/C_th_DRY/R_cond_DRY;
b_DRY_t_13=1/C_th_DRY*ssp_m_H2O_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2*(C_Pm_H2O_gas*(ssp_T_DRY-273.15)+ssp_H_H2O_0-C_Pm_H2*(ssp_T_DRY-273.15)-ssp_H_H2_0);
b_DRY_t_14=1/C_th_DRY*ssp_m_H2_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2*(C_Pm_H2*(ssp_T_DRY-273.15)+ssp_H_H2_0-C_Pm_H2O_gas*(ssp_T_DRY-273.15)-ssp_H_H2O_0);
b_DRY_t_15=-1/C_th_DRY*ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2*(C_Pm_H2O_gas*(ssp_T_SEPH-273.15)+ssp_H_H2O_0-C_Pm_H2*(ssp_T_SEPH-273.15)-ssp_H_H2_0);
b_DRY_t_16=-1/C_th_DRY*ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2*(C_Pm_H2*(ssp_T_SEPH-273.15)+ssp_H_H2_0-C_Pm_H2O_gas*(ssp_T_SEPH-273.15)-ssp_H_H2O_0);
b_DRY_t_17=-1/C_th_DRY*((ssp_m_H2O_DRY*C_Pm_H2O_gas+ssp_m_H2_DRY*C_Pm_H2)*(ssp_T_DRY-273.15)+ssp_m_H2_DRY*ssp_H_H2_0+ssp_m_H2O_DRY*ssp_H_H2O_0)/(ssp_m_H2O_DRY+ssp_m_H2_DRY);
b_DRY_t_18=1/C_th_DRY*b_SEPH_t_15*-ssp_C_th_SEPH;
b_DRY_t_19=1/C_th_DRY*((C_Pm_H2O_gas-2*C_Pm_H2O_lid)*(ssp_T_DRY-273.15)+ssp_H_H2O_0);
%==========================================================================
% gaseous domain linearization elements
a_DRY_g_11=-ssp_m_H2O_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2;
a_DRY_g_12=ssp_m_H2_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2;
a_DRY_g_21=ssp_m_H2O_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2;
a_DRY_g_22=-ssp_m_H2_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2;
b_DRY_g_11=ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;
b_DRY_g_12=-ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;
b_DRY_g_21=-ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;
b_DRY_g_22=ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;
b_DRY_g_13=ssp_m_H2_SEPH/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);
b_DRY_g_23=ssp_m_H2O_SEPH/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);
b_DRY_g_14=-ssp_m_H2_DRY/(ssp_m_H2O_DRY+ssp_m_H2_DRY);
b_DRY_g_24=-ssp_m_H2O_DRY/(ssp_m_H2O_DRY+ssp_m_H2_DRY);
c_DRY_g_11=R*ssp_T_DRY/(MW_H2*V_DRY);
c_DRY_g_22=R*ssp_T_DRY/(MW_H2O*V_DRY);
d_DRY_g_11=R*ssp_m_H2_DRY/(MW_H2*V_DRY);
d_DRY_g_21=R*ssp_m_H2O_DRY/(MW_H2O*V_DRY);
%==========================================================================
% small-signal subsystem of the unit enclosure
%==========================================================================
% thermal domain linearization elements:
a_Enc_t_11=-0.001/C_th_Enc*(1/R_th_cell+1/R_cond_RECH+1/R_cond_DRY+1/R_cond_RECO+1/R_cond_SEPH+1/R_cond_SEPO+1/R_cond_enc+1000*C_p_air*m_dot_fan);
b_Enc_t_11=0.001/C_th_Enc*1/R_th_cell;
b_Enc_t_12=0.001/C_th_Enc*1/R_cond_SEPO;
b_Enc_t_13=0.001/C_th_Enc*1/R_cond_SEPH;
b_Enc_t_14=0.001/C_th_Enc*1/R_cond_RECO;
b_Enc_t_15=0.001/C_th_Enc*1/R_cond_RECH;
b_Enc_t_16=0.001/C_th_Enc*1/R_cond_DRY;
b_Enc_t_17=1/C_th_Enc*(C_p_air*m_dot_fan+1/R_cond_enc);
%==========================================================================
% small-signal subsystem of the Hydrogen storage tank##
%==========================================================================
a_STRG_g_11=-wc_STRG;
a_STRG_g_12=wc_STRG*KiSTRG*ValvArea*ac_A_rst;
b_STRG_g_11=wc_STRG*KpSTRG*ValvArea*ac_A_rst;
b_STRG_g_12=wc_STRG*KpSTRG*ValvArea*ac_A_rst;
c_STRG_g_11=sqrt(ssp_P_SPVA-P_PURI)/R_Dry_Fix*sign(ssp_m_d_PURI);
d_STRG_g_11=ssp_A_rest/R_Dry_Fix/2/sqrt(ssp_P_SPVA-P_PURI)*sign(ssp_m_d_PURI);
d_STRG_g_12=d_STRG_g_11;
%==========================================================================
%% PEMEL small-signal model construction
%==========================================================================
% the following represent the rows in the state matrix A_el of the PEMEL 
a01=[a_pem_e_11	0	        b_pem_e_11	b_pem_e_13*d_SEPO_g_22+b_pem_e_14*d_SEPO_g_12	b_pem_e_12*d_SEPH_g_12+b_pem_e_15*d_SEPH_g_22	0	0	0	0	0	0	b_pem_e_14*c_SEPO_g_11	b_pem_e_13*c_SEPO_g_22	b_pem_e_12*c_SEPH_g_11	b_pem_e_15*c_SEPH_g_22	0	0	b_pem_e_13*d_SEPO_g_21+b_pem_e_14*d_SEPO_g_11	b_pem_e_12*d_SEPH_g_11+b_pem_e_15*d_SEPH_g_21	0	0	0	0	0	0	0	0];
a02=[a_pem_e_21	a_pem_e_22	b_pem_e_21	b_pem_e_23*d_SEPO_g_22+b_pem_e_24*d_SEPO_g_12	b_pem_e_22*d_SEPH_g_12+b_pem_e_25*d_SEPH_g_22	0	0	0	0	0	0	b_pem_e_24*c_SEPO_g_11	b_pem_e_23*c_SEPO_g_22	b_pem_e_22*c_SEPH_g_11	b_pem_e_25*c_SEPH_g_22	0	0	b_pem_e_23*d_SEPO_g_21+b_pem_e_24*d_SEPO_g_11	b_pem_e_22*d_SEPH_g_11+b_pem_e_25*d_SEPH_g_21	0	0	0	0	0	0	0	0];
a03=[b_pem_t_16*c_pem_e_21	b_pem_t_15	a_pem_t_11+b_pem_t_13*d_pem_g_3+b_pem_t_16*d_pem_e_21	b_pem_t_13*d_pem_g_1*d_SEPO_g_12+b_pem_t_16*d_pem_e_23*d_SEPO_g_22+b_pem_t_16*d_pem_e_24*d_SEPO_g_12	b_pem_t_13*d_pem_g_2*d_SEPH_g_22+b_pem_t_16*d_pem_e_22*d_SEPH_g_12+b_pem_t_16*d_pem_e_25*d_SEPH_g_22	b_pem_t_19	0	0	b_pem_t_18	0	b_pem_t_110	b_pem_t_13*d_pem_g_1*c_SEPO_g_11+b_pem_t_16*d_pem_e_24*c_SEPO_g_11	b_pem_t_16*d_pem_e_23*c_SEPO_g_22	b_pem_t_16*d_pem_e_22*c_SEPH_g_11	b_pem_t_13*d_pem_g_2*c_SEPH_g_22+b_pem_t_16*d_pem_e_25*c_SEPH_g_22	0	0	b_pem_t_13*d_pem_g_1*d_SEPO_g_11+b_pem_t_16*d_pem_e_23*d_SEPO_g_21+b_pem_t_16*d_pem_e_24*d_SEPO_g_11	b_pem_t_13*d_pem_g_2*d_SEPH_g_21+b_pem_t_16*d_pem_e_22*d_SEPH_g_11+b_pem_t_16*d_pem_e_25*d_SEPH_g_21	0	0	b_pem_t_12+b_pem_t_11   0	0	0	0	0];
a04=[0	b_SEPO_t_18	b_SEPO_t_11+b_SEPO_t_16*d_pem_g_3	a_SEPO_t_11+b_SEPO_t_15*d_SEPO_f_11*d_SEPO_g_12+b_SEPO_t_15*d_SEPO_f_12*d_SEPO_g_22+b_SEPO_t_16*d_pem_g_1*d_SEPO_g_12+b_SEPO_t_17*d_SEPO_g_32	b_SEPO_t_16*d_pem_g_2*d_SEPH_g_22	0	0	0	0	0	b_SEPO_t_12	b_SEPO_t_19+b_SEPO_t_15*d_SEPO_f_11*c_SEPO_g_11+b_SEPO_t_16*d_pem_g_1*c_SEPO_g_11+b_SEPO_t_17*c_SEPO_g_33	b_SEPO_t_110+b_SEPO_t_15*d_SEPO_f_12*c_SEPO_g_22	0	b_SEPO_t_16*d_pem_g_2*c_SEPH_g_22	0	0	b_SEPO_t_111+b_SEPO_t_15*d_SEPO_f_11*d_SEPO_g_11+b_SEPO_t_15*d_SEPO_f_12*d_SEPO_g_21+b_SEPO_t_16*d_pem_g_1*d_SEPO_g_11+b_SEPO_t_17*d_SEPO_g_31	b_SEPO_t_16*d_pem_g_2*d_SEPH_g_21	b_SEPO_t_14	0	b_SEPO_t_13	0	0	0	0	0];
a05=[0	b_SEPH_t_18	b_SEPH_t_11+b_SEPH_t_16*d_pem_g_3	b_SEPH_t_16*d_pem_g_1*d_SEPO_g_12	a_SEPH_t_11+b_SEPH_t_15*d_SEPH_f_11*(d_SEPH_g_12+d_SEPH_g_22)+b_SEPH_t_16*d_pem_g_2*d_SEPH_g_22+b_SEPH_t_17*d_SEPH_g_32	0	0	0	0	b_SEPH_t_15*d_SEPH_f_14*(d_DRY_g_11+d_DRY_g_21)	b_SEPH_t_12	b_SEPH_t_16*d_pem_g_1*c_SEPO_g_11	0	b_SEPH_t_19+b_SEPH_t_15*d_SEPH_f_11*c_SEPH_g_11	b_SEPH_t_110+b_SEPH_t_15*d_SEPH_f_12*c_SEPH_g_22+b_SEPH_t_16*d_pem_g_2*c_SEPH_g_22+b_SEPH_t_17*c_SEPH_g_33	b_SEPH_t_15*d_SEPH_f_14*c_DRY_g_11	b_SEPH_t_15*d_SEPH_f_14*c_DRY_g_22	b_SEPH_t_16*d_pem_g_1*d_SEPO_g_11	b_SEPH_t_111+b_SEPH_t_15*d_SEPH_f_11*(d_SEPH_g_11+d_SEPH_g_21)+b_SEPH_t_16*d_pem_g_2*d_SEPH_g_21+b_SEPH_t_17*d_SEPH_g_31	0	b_SEPH_t_14	b_SEPH_t_13	0	0	0	0   0];
a06=[0	0	0	b_RECO_t_11	0	a_RECO_t_11	a_RECO_t_12	0	0	0	b_RECO_t_12	0	0	0	0	0	0	0	0	0	0	b_RECO_t_13	0	0	0	0	0];
a07=[0	0	0	0	0	a_RECO_t_21	a_RECO_t_22	a_RECO_t_23	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
a08=[0	0	0	0	0	0	a_RECO_t_32	a_RECO_t_33	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
a09=[0	0	0	0	b_RECH_t_11	0	0	0	a_RECH_t_11	0	b_RECH_t_12	0	0	0	0	0	0	0	0	0	0	b_RECH_t_13	0	0	0	0   0];
a10=[0	0	0	0	b_DRY_t_11+b_DRY_t_18*d_SEPH_f_11*(d_SEPH_g_12+d_SEPH_g_22)	0	0	0	0	a_DRY_t_11+b_DRY_t_17*d_STRG_g_11*d_DRY_g_11+b_DRY_t_17*d_STRG_g_12*d_DRY_g_21+b_DRY_t_18*d_SEPH_f_14*(d_DRY_g_11+d_DRY_g_21)	b_DRY_t_12	0	0	b_DRY_t_15+b_DRY_t_18*d_SEPH_f_11*c_SEPH_g_11	b_DRY_t_16+b_DRY_t_18*d_SEPH_f_12*c_SEPH_g_22	b_DRY_t_13+b_DRY_t_17*d_STRG_g_11*c_DRY_g_11+b_DRY_t_18*d_SEPH_f_14*c_DRY_g_11	b_DRY_t_14+b_DRY_t_17*d_STRG_g_12*c_DRY_g_22+b_DRY_t_18*d_SEPH_f_14*c_DRY_g_22+b_DRY_t_19*w_ads	0	b_DRY_t_18*d_SEPH_f_11*(d_SEPH_g_11+d_SEPH_g_21)	0	0	0	b_DRY_t_17*c_STRG_g_11	0	0	0	0];
a11=[0	0	b_Enc_t_11	b_Enc_t_12	b_Enc_t_13	b_Enc_t_14	0	0	b_Enc_t_15	b_Enc_t_16	a_Enc_t_11	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
a12=[0	0	1*d_pem_g_3	1*d_pem_g_1*d_SEPO_g_12+1*d_SEPO_g_32+b_SEPO_g_13*d_SEPO_f_11*d_SEPO_g_12+b_SEPO_g_13*d_SEPO_f_12*d_SEPO_g_22	1*d_pem_g_2*d_SEPH_g_22	0	0	0	0	0	0	a_SEPO_g_11+1*d_pem_g_1*c_SEPO_g_11+1*c_SEPO_g_33+b_SEPO_g_13*d_SEPO_f_11*c_SEPO_g_11	a_SEPO_g_12+b_SEPO_g_13*d_SEPO_f_12*c_SEPO_g_22	0	1*d_pem_g_2*c_SEPH_g_22	0	0	1*d_pem_g_1*d_SEPO_g_11+1*d_SEPO_g_31+b_SEPO_g_13*d_SEPO_f_11*d_SEPO_g_11+b_SEPO_g_13*d_SEPO_f_12*d_SEPO_g_21	1*d_pem_g_2*d_SEPH_g_21	0	0	0	0	0	0	0	0];
a13=[0	b_SEPO_g_24	0	b_SEPO_g_23*d_SEPO_f_11*d_SEPO_g_12+b_SEPO_g_23*d_SEPO_f_12*d_SEPO_g_22	0	0	0	0	0	0	0	a_SEPO_g_21+b_SEPO_g_23*d_SEPO_f_11*c_SEPO_g_11	a_SEPO_g_22+b_SEPO_g_23*d_SEPO_f_12*c_SEPO_g_22	0	0	0	0	b_SEPO_g_23*d_SEPO_f_11*d_SEPO_g_11+b_SEPO_g_23*d_SEPO_f_12*d_SEPO_g_21	0	0	0	0	0	0	0	0	0];
a14=[0	b_SEPH_g_14	0	0	b_SEPH_g_11*d_SEPH_f_11*(d_SEPH_g_12+d_SEPH_g_22)	0	0	0	0	b_SEPH_g_11*d_SEPH_f_14*(d_DRY_g_11+d_DRY_g_21)	0	0	0	a_SEPH_g_11+b_SEPH_g_11*d_SEPH_f_11*c_SEPH_g_11	a_SEPH_g_12+b_SEPH_g_11*d_SEPH_f_12*c_SEPH_g_22	b_SEPH_g_11*d_SEPH_f_14*c_DRY_g_11	b_SEPH_g_11*d_SEPH_f_14*c_DRY_g_22	0	b_SEPH_g_11*d_SEPH_f_11*(d_SEPH_g_11+d_SEPH_g_21)	0	0	0	0	0	0	0	0];
a15=[0	b_SEPH_g_24	b_SEPH_g_23	0	b_SEPH_g_21*d_SEPH_f_11*(d_SEPH_g_12+d_SEPH_g_22)+1*d_SEPH_g_32	0	0	0	0	b_SEPH_g_21*d_SEPH_f_14*(d_DRY_g_11+d_DRY_g_21)	0	0	0	a_SEPH_g_21+b_SEPH_g_21*d_SEPH_f_11*c_SEPH_g_11	a_SEPH_g_22+b_SEPH_g_21*d_SEPH_f_12*c_SEPH_g_22+1*c_SEPH_g_33	b_SEPH_g_21*d_SEPH_f_14*c_DRY_g_11	b_SEPH_g_21*d_SEPH_f_14*c_DRY_g_22	0	b_SEPH_g_21*d_SEPH_f_11*(d_SEPH_g_11+d_SEPH_g_21)+1*d_SEPH_g_31	0	0	0	0	0	0	0	0];
a16=[0	0	0	0	b_DRY_g_13*d_SEPH_f_11*(d_SEPH_g_12+d_SEPH_g_22)	0	0	0	0	b_DRY_g_13*d_SEPH_f_14*(d_DRY_g_11+d_DRY_g_21)+b_DRY_g_14*d_STRG_g_11*d_DRY_g_11+b_DRY_g_14*d_STRG_g_12*d_DRY_g_21	0	0	0	b_DRY_g_11+b_DRY_g_13*d_SEPH_f_11*c_SEPH_g_11	b_DRY_g_12+b_DRY_g_13*d_SEPH_f_12*c_SEPH_g_22	a_DRY_g_11+b_DRY_g_13*d_SEPH_f_14*c_DRY_g_11+b_DRY_g_14*d_STRG_g_11*c_DRY_g_11	a_DRY_g_12+b_DRY_g_13*d_SEPH_f_14*c_DRY_g_22+b_DRY_g_14*d_STRG_g_12*c_DRY_g_22	0	b_DRY_g_13*d_SEPH_f_11*(d_SEPH_g_11+d_SEPH_g_21)	0	0	0	b_DRY_g_14*c_STRG_g_11	0	0	0	0];
a17=[0	0	0	0	b_DRY_g_23*d_SEPH_f_11*(d_SEPH_g_12+d_SEPH_g_22)	0	0	0	0	b_DRY_g_23*d_SEPH_f_14*(d_DRY_g_11+d_DRY_g_21)+b_DRY_g_24*d_STRG_g_11*d_DRY_g_11+b_DRY_g_24*d_STRG_g_12*d_DRY_g_21	0	0	0	b_DRY_g_21+b_DRY_g_23*d_SEPH_f_11*c_SEPH_g_11	b_DRY_g_22+b_DRY_g_23*d_SEPH_f_12*c_SEPH_g_22	a_DRY_g_21+b_DRY_g_23*d_SEPH_f_14*c_DRY_g_11+b_DRY_g_24*d_STRG_g_11*c_DRY_g_11	a_DRY_g_22+b_DRY_g_23*d_SEPH_f_14*c_DRY_g_22+b_DRY_g_24*d_STRG_g_12*c_DRY_g_22-w_ads	0	b_DRY_g_23*d_SEPH_f_11*(d_SEPH_g_11+d_SEPH_g_21)	0	0	0	b_DRY_g_24*c_STRG_g_11	0	0	0	0];
a18=[0	b_SEPO_f_12	b_SEPO_f_11	-1*d_SEPO_g_32	0	0	0	0	0	0	0	-1*c_SEPO_g_33	0	0	0	0	0	-1*d_SEPO_g_31	0	1	0	0	0	0	0	0	0];
a19=[0	0	-1*d_pem_g_3	-1*d_pem_g_1*d_SEPO_g_12	-1*d_pem_g_2*d_SEPH_g_22-1*d_SEPH_g_32	0	0	0	0	0	0	-1*d_pem_g_1*c_SEPO_g_11	0	0	-1*d_pem_g_2*c_SEPH_g_22+-1*c_SEPH_g_33	0	0	-1*d_pem_g_1*d_SEPO_g_11	-1*d_pem_g_2*d_SEPH_g_21+-1*d_SEPH_g_31	0	1	0	0	0	0	0	0];
a20=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	a_SEPO_f_21	0	a_SEPO_f_22	0	0	0	a_SEPO_f_23	0	0	0];
a21=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	a_SEPH_f_21	0	a_SEPH_f_22	0	0	0	a_SEPH_f_23	0	0];
a22=[0	0	b_RECO_f_11	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	a_RECO_f_11	0	0	0	a_RECO_f_12	0];
a24=[0	0	0	0	0	0	0	0	0	b_STRG_g_11*(d_DRY_g_11+d_DRY_g_21)	0	0	0	0	0	b_STRG_g_11*c_DRY_g_11	b_STRG_g_12*c_DRY_g_22	0	0	0	0	0	a_STRG_g_11	0	0	0	a_STRG_g_12];
a25=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	a_SEPO_f_31	0	0	0	0	0	0	0	0	0];
a26=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	a_SEPH_f_31	0	0	0	0	0	0	0	0];
a27=ac_T_spo*[0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
a29=ac_A_rst*[0	0	0	0	0	0	0	0	0	(d_DRY_g_11+d_DRY_g_21)	0	0	0	0	0	c_DRY_g_11	c_DRY_g_22	0	0	0	0	0	0	0	0	0	0];

% the state-space representation of the PEMEL unit without its PEI interafaces is given as follows:
A_ele=[a01;a02;a03;a04;a05;a06;a07;a08;a09;a10;a11;a12;a13;a14;a15;a16;a17;a18;a19;a20;a21;a22;a24;a25;a26;a27;a29];
B_ele=[b_pem_e_16;b_pem_e_26;b_pem_t_17+b_pem_t_16*d_pem_e_26;zeros(size(A_ele,2)-3,1)];
C_icl=[c_pem_e_21	0	d_pem_e_21	d_pem_e_23*d_SEPO_g_22+d_pem_e_24*d_SEPO_g_12	d_pem_e_22*d_SEPH_g_12+d_pem_e_25*d_SEPH_g_22	0	0	0	0	0	0	d_pem_e_24*c_SEPO_g_11	d_pem_e_23*c_SEPO_g_22	d_pem_e_22*c_SEPH_g_11	d_pem_e_25*c_SEPH_g_22	0	0	d_pem_e_23*d_SEPO_g_21+d_pem_e_24*d_SEPO_g_11	d_pem_e_22*d_SEPH_g_11+d_pem_e_25*d_SEPH_g_21	0	0	0	0	0	0	0	0];
D_icl=d_pem_e_26;
%==========================================================================
%% Inverter and Buck converters linearization:
%==========================================================================
% block submatrcies of the dc converter are given as follows:
A_dc_lc=[0            1/C_el                0
        -1/Ls         0                     Dbar/Ls
        0            -Dbar/Cs               Pinvbar/(Vdc^2*Cs)];
B_dc_lc_1=[-Stacks_No/C_el
           0
           0];
B_dc_lc_2=[0;Vdc/Ls;-Ib/Cs];
B_dc_lc_3=[0;0;-1/Cs/Vdc];
A_dc_ct=[0       0
         Kdcvi   0];
B_dc_ct_1=[0 0 1
           0 0 Kdcvp];
B_dc_ct_2=[0
           -Stacks_No];
C_dc_ct=[Kdcip*Kdcvi Kdcii];
D_dc_ct_1=[0 0 Kdcip*Kdcvp];
D_dc_ct_2=-Stacks_No*Kdcip;
%==========================================================================
% The state space representation of the dc-dc buck converter subsystem
A_DC=[A_dc_lc+B_dc_lc_2*D_dc_ct_1          B_dc_lc_2*C_dc_ct
     B_dc_ct_1                          A_dc_ct];
B_DC_cl=[B_dc_lc_1+B_dc_lc_2*D_dc_ct_2
      B_dc_ct_2];
B_DC_p=[B_dc_lc_3
        zeros(2,1)];
C_DC_ec=[1/cells_per_stack 0 0 0 0];
C_DC_ed=[0 0 1 0 0];
C_DC=[C_DC_ec
    C_DC_ed];
%==========================================================================
% block submatrcies of the ac inverter are given as follows:
Ap=[-wc 0
    0 -wc];
Bp_1=[0     0   wc*Iod1           wc*Ioq1         wc*Vod1           wc*Voq1
      0     0   -wc*Ioq1          wc*Iod1         wc*Voq1           -wc*Vod1];
Cpw=[-m1+md1*wc 0];
Cpv=[0 -n1+nd1*wc
       0 0];
Dpw_1=[0     0   -md1*wc*Iod1     -md1*wc*Ioq1     -md1*wc*Vod1  -md1*wc*Voq1];

Dpv_1=[ 0     0   nd1*wc*Ioq1     -nd1*wc*Iod1    -nd1*wc*Voq1    nd1*wc*Vod1
        0     0   0                 0               0                0       ];
Bv1=[1 0
     0 1];       
Bv2=[0 0 -1 0 0 0
     0 0 0 -1 0 0];
Cv=[Kiv 0
    0 Kiv];      
Dv1=[Kpv 0
     0 Kpv];    
Dv2=[0 0 -Kpv       -wnom*Cf    H 0
     0 0 wnom*Cf    -Kpv        0 H];
Bc1=[1 0
     0 1];       
Bc2=[-1 0 0 0 0 0
     0 -1 0 0 0 0];
Cc=[Kii 0
    0   Kii];     
Dc1=[Kpi 0
     0   Kpi];    
Dc2=[-Kpi       -wnom*Lf 0 0 0 0
     wnom*Lf    -Kpi     0 0 0 0]; 
ALCL=[-Rf/Lf  w0     -1/Lf  0       0               0               
      -w0     -Rf/Lf 0      -1/Lf   0               0               
      1/Cf    0      0      w0      -1/Cf           0               
      0       1/Cf   -w0    0       0               -1/Cf           
      0       0      1/Lc   0       -Rc/Lc-rn/Lc    w0              
      0       0      0      1/Lc    -w0             -Rc/Lc-rn/Lc];   
BLCL1=[1/Lf 0
        0   1/Lf
        0   0
        0   0
        0   0
        0   0]; 
BLCL2=[0            0
       0            0
       0            0
       0            0
       -rn/Lc       0
       0            -rn/Lc];
BLCL3_1=[Ilq1
        -Ild1
        Voq1
        -Vod1
        Ioq1 
        -Iod1]; 
BLCL4_1=[0 1/Lf*Vid1/Vdc
         0 1/Lf*Viq1/Vdc
         0 0
         0 0
         0 0
         0 0]; 
BLCL5_1=[0                           0
         0                           0
         0                           0
         0                           0
         rn/Lc                       0
         0                           rn/Lc];
ALoad=[-rn/L_load-R_load/L_load    w0  
       -w0                         -rn/L_load-R_load/L_load];
Bload1=[0 0 0 0 rn/L_load    0
        0 0 0 0 0            rn/L_load];
Bload2=[Iloadq
        -Iloadd];
Bload3=[rn/L_load    0
        0            rn/L_load];
p_il=[Ild1 Ilq1];
p_ed=[Vid1/Vdc;Viq1/Vdc];
p_cc=[Vid1 Viq1 0 0 0 0];
C_AC_p=[p_il*Dc1*Dv1*Cpv         p_il*Dc1*Cv         p_il*Cc       p_il*(Dc1*Dv2+Dc2+Dc1*Dv1*Dpv_1)+p_cc     zeros(1,2)];
D_AC_p=p_il*p_ed;
%==========================================================================
%Inverter states [P Q phid phiq gamd gamq ild ilq vod voq iod ioq igd igq]
% Inverter subsystem state-space represention
Ainv=[Ap                                                                       zeros(2,2)           zeros(2,2)           Bp_1                                                     zeros(2,2)
        Bv1*Cpv                                                                zeros(2,2)           zeros(2,2)           Bv2+Bv1*Dpv_1                                            zeros(2,2)
        Bc1*Dv1*Cpv                                                            Bc1*Cv               zeros(2,2)           Bc1*Dv2+Bc2+Bc1*Dv1*Dpv_1                                zeros(2,2)
        BLCL1*Dc1*Dv1*Cpv+BLCL3_1*Cpw                                        BLCL1*Dc1*Cv         BLCL1*Cc             ALCL+BLCL1*(Dc1*Dv2+Dc2+Dc1*Dv1*Dpv_1)+BLCL3_1*Dpw_1     BLCL5_1
        Bload2*Cpw                                                             zeros(2,2)           zeros(2,2)           Bload1+Bload2*Dpw_1                                      ALoad];
Binv_is=[zeros(6,2)
        BLCL2
        Bload3];
Binv_ed=[zeros(2,2)
         zeros(2,2)
         zeros(2,2)
         BLCL4_1
         zeros(2,2)];
%==========================================================================
%% GFPEMEL entire small-signal model construction
%==========================================================================
% the following A, B pair represent the final accumulation of all previous
% subsystems of the Grid-Forming Proton Exchange Membrane Hydrogen
% Electrolyzer:
%==========================================================================
A_tot=[A_ele                                    B_ele*C_DC_ec                                           zeros(size(A_ele,2),size(Ainv,2))
       B_DC_cl*C_icl                            A_DC+B_DC_cl*D_icl*C_DC_ec+B_DC_p*D_AC_p*C_DC_ed        B_DC_p*C_AC_p
       zeros(size(Ainv,2),size(A_ele,2))      Binv_ed*C_DC                                            Ainv];
B_tot=[zeros(26,2)
       zeros(5,2)
       Binv_is];

% the eigenvalues of the system can be evaluated as follows:
Eigen=eig(A_tot);
