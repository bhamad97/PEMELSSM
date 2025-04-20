%% Code to define parameters for ssc_electrolyzer
format long

% Environment conditions
env_p = 0.101325; % [MPa] Pressure
env_T = 20; % [degC] Temperature
cells_per_stack=48;
Stacks_No=5;
stack_num_cells = cells_per_stack*Stacks_No; % [-] Number cells

%stack_area = 280; % [cm^2] Cell area
A_PEM= 180; % [cm^2] Cell area

% stack_t_membrane = 125; % [um] Membrane thickness
L_PEM = 0.125; % [cm]
% R_others= 0.005;



stack_t_gdl_A = 25; % [um] Anode gas diffusion layer thickness
stack_t_gdl_C = 250; % [um] Cathode gas diffusion layer thickness
% stack_w_channels = 1; % [cm] Gas channel width and height
% stack_num_channels = 8; % [-] Number of gas channels per cell
% stack_io = 1e-04; % [A/cm^2] Exchange current density
% stack_alpha = 0.7; % [-] Charge transfer coefficient
% stack_D_gdl_A = 0.07; % [cm^2/s] Water diffusivity in anode GDL
% stack_D_gdl_C = 0.07; % [cm^2/s] Water diffusivity in cathode GDL
% stack_membrane_rho = 2000; % [kg/m^3] Density of dry membrane
% stack_membrane_MW = 1.1; % [kg/mol] Equivalent weight of dry membrane
stack_mea_rho = 1800; % [kg/m^3] Overall density of membrane electrode assembly
stack_mea_cp = 870; % [J/(kg*K)] Overall specific heat of membrane electrode assembly


% water_pipe_D = 0.01; % m
% gas_pipe_D = 0.01; % m

% Heat exchanger dimensions
exchanger_L = 1; % [m] Overall radiator length
exchanger_W = 0.025; % [m] Overall radiator width
exchanger_H = 0.5; % [m] Overal radiator height
exchanger_N_tubes = 25; % [-] Number of coolant tubes
exchanger_tube_H = 0.0015; % [m] Height of each coolant tube
exchanger_fin_spacing = 0.002; % [-] Fin spacing
exchanger_eta_fin = 0.7; % [-] Fin efficiency

exchanger_gap_H = (exchanger_H - exchanger_N_tubes*exchanger_tube_H) / (exchanger_N_tubes - 1); % [m] Height between coolant tubes
exchanger_air_area_primary = 2 * (exchanger_N_tubes - 1) * exchanger_W * (exchanger_L + exchanger_gap_H); % [m^2] Primary air heat transfer surface area
exchanger_N_fins = (exchanger_N_tubes - 1) * exchanger_L / exchanger_fin_spacing; % [-] Total number of fins
exchanger_air_area_fins = 2 * exchanger_N_fins * exchanger_W * exchanger_gap_H; % [m^2] Total fin surface area
% exchanger_tube_Leq = 2*(exchanger_H + 20*exchanger_tube_H*exchanger_N_tubes); % [m] Additional equivalent tube length for losses due to manifold and splits

% Hydrogen property tables
H2_R = 4124.48151675695; % [J/(kg*K)] Specific gas constant
H2_T = [-56.55, -50:10:-10, -5:1:5, 10:10:350]; % [degC] Temperature vector
H2_h = [2783.66044045879, 2873.94181645932, 3012.62601981068, 3152.22537301871, 3292.62486344326, 3433.72138359680, 3504.50184295996, 3518.67512520997, 3532.85394543712, 3547.03822323045, 3561.22787913325, 3575.42283463482, 3589.62301216223, 3603.82833507219, 3618.03872764290, 3632.25411506595, 3646.47442343830, 3717.64725530348, 3860.32198734890, 4003.38288163040, 4146.77354613486, 4290.44463593559, 4434.35318483605, 4578.46197843019, 4722.73896833090, 4867.15672726297, 5011.69194456718, 5156.32496143709, 5301.03934494342, 5445.82149962224, 5590.66031514150, 5735.54684833147, 5880.47403768083, 6025.43644826542, 6170.43004499140, 6315.45199199528, 6460.50047604533, 6605.57455182556, 6750.67400704912, 6895.79924543541, 7040.95118568927, 7186.13117473630, 7331.34091359044, 7476.58239435534, 7621.85784698693, 7767.16969456766, 7912.52051596217, 8057.91301483815, 8203.34999414341, 8348.83433523064, 8494.36898091475, 8639.95692183309]; % [kJ/kg] Specific enthalpy vector
% H2_D = 74; % [mm^2/s] Diffusivity

% Oxygen property tables
O2_R = 259.836612622973; % [J/(kg*K)] Specific gas constant
O2_T = [-56.55, -50:10:-10, -5:1:5, 10:10:350]; % [degC] Temperature vector
O2_h = [196.314045635230, 202.302438242727, 211.445560194021, 220.590899228792, 229.740231593225, 238.895348952500, 243.475641264001, 244.391947430149, 245.308340402054, 246.224821996439, 247.141394030059, 248.058058319568, 248.974816681386, 249.891670931563, 250.808622885637, 251.725674358497, 252.642827164235, 257.230174605453, 266.413508632715, 275.609852765421, 284.820965750972, 294.048557917965, 303.294277512928, 312.559698673577, 321.846311315623, 331.155513043735, 340.488603075706, 349.846778083808, 359.231129801281, 368.642644208585, 378.082202097877, 387.550580810768, 397.048456949970, 406.576409877154, 416.134925824813, 425.724402467528, 435.345153816409, 444.997415318718, 454.681349062211, 464.397049000008, 474.144546126762, 483.923813550211, 493.734771414086, 503.577291638585, 513.451202453549, 523.356292706982, 533.292315937896, 543.258994207727, 553.256021688858, 563.283068012223, 573.339781378746, 583.425791441447]; % [kJ/kg] Specific enthalpy vectdor
% O2_D = 18; % [mm^2/s] Diffusivity

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

T_TLU = [-56.55, -50:10:-10, -5:1:5, 10:10:350]';%, 'degC'}; % Temperature vector

F     = 96485.33212; %     'C/mol'  }; % Faraday constant
R     = 8.31446261815324;%, 'J/K/mol'}; % Universal gas constant

J_cell_0=8e-05;%[A/cm^2]

T_ref=298; %[k]
T_ref_st=273.15+60; %[k]
T_ref_cool=env_T+273.15;
T_thres=10;

v_H2O=1; %[-]stoichiometric coefficient
v_H2=1; %[-]
v_O2=0.5; %[-]

R_w=461.523;%J/Kg/K

MW_H2O = R/R_w;%kg/mols
MW_O2  = R/O2_R;%kg/mols
MW_H2  = R/H2_R;%kg/mols

% exchanger_area=exchanger_air_area_primary + exchanger_eta_fin*exchanger_air_area_fins;
% Heat_tran_coeff_pem=300;
% R_th_cell=(exchanger_area)*Heat_tran_coeff_pem;
% h_stack=2*A_PEM*1e-4 +4*(stack_t_membrane + stack_t_gdl_A + stack_t_gdl_C)*1e-6
% the 4 is 4 W/m*K, the convection with air varies with the temperature
% and the direction of the heat transfer if it is upward downard or from
% the side. It depends as well on the shape of the material. Considering
% the stack to have 4 vertical surfaces and one surface upward and one
% downard, the convection happens at six surfaces. The air is considered
% the fluid surrounding the electrolyzer.
% https://quickfield.com/natural_convection.htm according to this website
% and for 80 degree stack temperature and 25 degree atm temperature, the
% follwing h are obtained:
% Vertical h considering length sqrt(A_PEM), h_v=5.9 W/m*K,
% Upward h (L_PEM*1e4 + stack_t_gdl_A + stack_t_gdl_C)*1e-6*sqrt(A_PEM*1e-4)*stack_num_cells
% (L_PEM*1e4 + stack_t_gdl_A + stack_t_gdl_C)*1e-6*2*stack_num_cells+2*sqrt(A_PEM*1e-4)
% Upward h=9.4 W/m*K
% downward h=4.7 W/m*K

h_t=(9.4+4.7)*(L_PEM*1e4 + stack_t_gdl_A + stack_t_gdl_C)*1e-6*sqrt(A_PEM*1e-4)*stack_num_cells+4*5.9*A_PEM*1e-4;
R_th_cell=1/h_t;
% for a 100 kW electrolyzer as based on the matlab ex, 22.5 kW as heat is
% generated. There is not specific paper that detailed the heat transfer
% coefficients or thermal mass of the electrolyzer with the components
% discussed in this modeling. However, the simplest and most general rules
% can apply. Assuming the air is filling the enclosure of the electrolysis
% system, the convective heat transfer 

%C_Pm_H2O=75.38;% J/(mol C) at 25
% T_vec=[273.1600 : 10 : 373.16]'; %#ok<NBRAK1> %K
% Pre_vec=[0.01, 0.1, 5 : 5 : 50];%MPa
% C_Pm_H2O=[4.2199, 4.2194, 4.1957, 4.1726, 4.1508, 4.13, 4.1103, 4.0916, 4.0738, 4.0569, 4.0408, 4.0254; 4.1955, 4.1951, 4.1769, 4.1592, 4.1422, 4.1259, 4.1103, 4.0954, 4.0812, 4.0676, 4.0545, 4.0421; 4.1843, 4.184, 4.169, 4.1543, 4.1401, 4.1264, 4.1133, 4.1007, 4.0885, 4.0768, 4.0656, 4.0548; 4.1801, 4.1798, 4.1667, 4.1538, 4.1413, 4.1293, 4.1177, 4.1064, 4.0956, 4.0851, 4.075, 4.0653; 4.1796, 4.1794, 4.1675, 4.1558, 4.1444, 4.1334, 4.1227, 4.1123, 4.1023, 4.0926, 4.0832, 4.0741; 4.1813, 4.1813, 4.1702, 4.1592, 4.1484, 4.138, 4.1279, 4.1181, 4.1086, 4.0993, 4.0903, 4.0816; 4.185, 4.185, 4.1742, 4.1636, 4.1532, 4.1432, 4.1334, 4.1239, 4.1146, 4.1056, 4.0969, 4.0884; 4.1901, 4.1901, 4.1795, 4.169, 4.1588, 4.1489, 4.1393, 4.1299, 4.1207, 4.1118, 4.1032, 4.0947; 4.1968, 4.1968, 4.1862, 4.1757, 4.1655, 4.1555, 4.1458, 4.1364, 4.1273, 4.1183, 4.1096, 4.1012; 4.2052, 4.2052, 4.1944, 4.1837, 4.1733, 4.1632, 4.1534, 4.1438, 4.1345, 4.1255, 4.1166, 4.108; 4.2136, 4.2136, 4.2045, 4.1935, 4.1828, 4.1724, 4.1623, 4.1524, 4.1429, 4.1335, 4.1245, 4.1157];%kJ/(K*kg)
% C_Pm_H2=(H2_h(2:end)-H2_h(1:end-1))./(H2_T(2:end)-H2_T(1:end-1));

C_Pm_H2=14.30;%KJ/(kg.K)
C_Pm_O2=0.918;%KJ/(kg.K)

% C_Pm_O2=(O2_h(2:end)-O2_h(1:end-1))./(O2_T(2:end)-O2_T(1:end-1));

C_Pm_H2O_avg=4.1816;%[KJ/(kg.K)]

%Water vapor enthalpy at 0 C 
% H_0_H2O=2.500994627588990e+03;%kJ/Kg

%Water vapor enthalpy:
h_w_vap_TLU=[2836.88241275372; 2837.81392500514; 2838.63937175807; 2838.7309929628; 2838.06905313927; 2836.62597341095; 2835.60023952573; 2835.37006381376; 2835.13143158077; 2834.88429145066; 2834.62859062897; 2500.93420564316; 2498.55329907119; 2496.17495082036; 2493.79885912205; 2491.42474723191; 2489.05236104642; 2477.20875029194; 2453.54955988604; 2429.83856603313; 2406.00136954922; 2381.97406342174; 2357.69101156389; 2333.08088387814; 2308.06565480412; 2282.56034405021; 2256.47287422313; 2229.70428017508; 2202.14968030993; 2173.69998896199; 2144.24368406447; 2113.66758247452; 2081.85585110571; 2048.68725710494; 2014.03141249899; 1977.7449923284; 1939.66849592886; 1899.62342931686; 1857.40930169838; 1812.79975416601; 1765.53721731251; 1715.32525772152; 1661.81700477787; 1604.59703745016; 1543.1534616334; 1476.83692476464; 1404.80240352424; 1325.92091736114; 1238.61667822208; 1140.5102987018; 1027.62017777647; 892.733785613825];%kJ/kg
h_w_TLU=[2396.55944251649; 2408.68643343608; 2427.1988031141; 2445.702165897; 2464.18429108356; 2482.62529466839; 2491.82135326629; 2493.6580151792; 2495.49372578218; 2497.32843835745; 2499.16210446923; 2500.99462758899; 2502.83092214066; 2504.6665223621; 2506.50140563013; 2508.335548891; 2510.16892865793; 2519.32352241995; 2537.56068088674; 2555.67742615292; 2573.6403223998; 2591.4109179101; 2608.9455037701; 2626.19492016262; 2643.10444358428; 2659.61377634497; 2675.58278696023; 2696.1846256545; 2716.49989553741; 2736.6235210957; 2756.61216047011; 2776.50561773853; 2796.33378922508; 2816.11979583049; 2835.88181423877; 2855.63430574511; 2875.38889327595; 2895.15500577172; 2914.94035964689; 2934.75132290228; 2954.59319330244; 2974.47041285505; 2994.38673458964; 3014.34535328009; 3034.34900867103; 3054.40006755764; 3074.50058946898; 3094.65237953784; 3114.85703128106; 3135.11596137801; 3155.43043805905; 3175.80160435813];%kJ/kg
% ads_C=(h_w_vap_TLU(12:25)-h_w_vap_TLU(13:26))./(T_TLU(12:25)-T_TLU(13:26));
% ads_C=2.428305185785030;%[kJ/(kg.K)]

% C_Pm_H2O_gas=(h_w_TLU(2:end)'-h_w_TLU(1:end-1)')./(O2_T(2:end)-O2_T(1:end-1));
C_Pm_H2O_gas=1.87;%KJ/(kg.K)
% the enthalpy of water starts from 0.06 at 0C and goes with a specific
% heat coefficient of around 4.1816. The water vapor 
% hfg_H2O_std=2396.55944251649;
% hfg_O2_std=2396.55944251649;
% hfg_H2_std=2396.55944251649;

epsi_0=8.8541878128e-12; % the permitivity in vaccuam 
epsilon_0=epsi_0*4;
% d_eps=0;
tau_e=0.1;

M_membrane=1.1; % kg/mol Equivalent weight of dry membrane
rho_membrane =2000;% kg/m^3 % Density of dry membrane

A_SEPO=0.25;%[m^2]
L_SEPO_max=0.5;%[m]
L_SEPO_min=0.25;%[m]
g=9.8;%N/kg

% water density Table:
% rhp_H2O=[999.7973, 999.8431, 1002.3216, 1004.8218, 1007.293, 1009.7358, 1012.1506, 1014.5379, 1016.8979, 1019.2311, 1021.538, 1023.8188; 999.6579, 999.7009, 1002.0303, 1004.382, 1006.7084, 1009.0101, 1011.2874, 1013.5407, 1015.7703, 1017.9766, 1020.1599, 1022.3206; 998.1632, 998.2045, 1000.4375, 1002.6924, 1004.9239, 1007.1324, 1009.3183, 1011.4819, 1013.6236, 1015.7438, 1017.8428, 1019.921; 995.6057, 995.6458, 997.81864, 1000.0128, 1002.1841, 1004.3332, 1006.4603, 1008.566, 1010.6505, 1012.7144, 1014.758, 1016.7815; 992.1724, 992.2119, 994.35133, 996.5113, 998.6484, 1000.7633, 1002.8563, 1004.928, 1006.9788, 1009.009, 1011.0192, 1013.0097; 988.02, 988.0299, 990.15684, 992.3035, 994.4268, 996.5273, 998.6056, 1000.6622, 1002.6976, 1004.7123, 1006.7067, 1008.6813; 983.18, 983.1901, 985.3217, 987.4721, 989.5983, 991.7008, 993.7803, 995.8374, 997.8726, 999.8864, 1001.8795, 1003.8522; 977.75, 977.7583, 979.9094, 982.0783, 984.2216, 986.34, 988.4343, 990.505, 992.5529, 994.5786, 996.5825, 998.5654; 971.77, 971.7836, 973.96734, 976.1678, 978.3409, 980.4876, 982.6087, 984.7049, 986.7768, 988.8254, 990.8511, 992.8545; 965.26, 965.3023, 967.53107, 969.7753, 971.9901, 974.1764, 976.3354, 978.4676, 980.5741, 982.6556, 984.7129, 986.7465; 958.62, 958.6315, 960.627, 962.9267, 965.1944, 967.4314, 969.6387, 971.8173, 973.9681, 976.0921, 978.19, 980.2628];%Kg/m^3
rho_H2O_avg=997;%Kg/m^3
V_SEPO=A_SEPO*L_SEPO_max;

A_SEPH=0.25;%[m^2]
L_SEPH_max=0.5;%[m]
L_SEPH_min=0.25;%[m]
V_SEPH=A_SEPH*L_SEPH_max;

max_flow_cat=0.2;% kg/s
max_flow_tank=0.2;% [kg/s]
L_SEPO_ref=(L_SEPO_max+L_SEPO_min)/2;% m
L_SEPH_ref=(L_SEPH_max+L_SEPH_min)/2;% m

R_SEP_AND=50000;
% for a nominal pressure drop of 0.001 MPa and a nominal mass flow of 0.002 kg/s
% This valve resistance for the gas flow from the vessel to atmoshpere
% given that the oxygen gas is evacuated.
% R_Valve_sep=1.5811e+04;
% % for a nominal pressure drop of 0.001 MPa and a nominal mass flow of 0.002 kg/s
% % This valve resistance for the water flow between the Hydrogen and
% % the Oxygen vessels.

R_SEP_CAT=50000;


% considering cylindrical vessel 
% For 80 degree stack temperature and 25 degree atm temperature, the
% follwing h are obtained:
% Vertical h considering length L_SEPO_max, h_v=4.2 W/m*K,
% Upward h A_SEPO
% 2*3.14*sqrt(A_SEPO/3.14)
% Upward h=6 W/m*K
% downward h=3 W/m*K

h_spo=(6+3)*A_SEPO+2*3.14*sqrt(A_SEPO/3.14)*L_SEPO_max*4.2;
R_cond_SEPO=1/h_spo;%K/W

% The required thickness of a pressure vessel made from 304L stainless steel
% for a pressure of 3 MPa and temperature of 80C can be calculated using the
% ASME Boiler and Pressure Vessel Code (BPVC) or similar standards.
% t= PR/((SE)0.6P), t is the thickness, P is the pressure to withstand
% S is allowable stress for 304L at 80C, E for fully radiographed welds or
% seamless construction E is 1, R is the radius and is here 0.2822 m. The
% thickness is then 8.22 mm with 2 mm is for  corrosion allowance. The
% specific heat capacity of 304L is 0.51 kJ/(Kg.k). The density of 304L is
% 7900 kg/m^3

C_th_SEPO_void=0.51*0.00822*(A_SEPO*2+2*3.14*sqrt(A_SEPO/3.14)*L_SEPO_max)*7900;%kJ/k



R_cond_SEPH=R_cond_SEPO;
% So R_cond_RECO represents the heat resistance of the recirculation
% circuit with the enclosure unit. It is different from the exchanger
% pipes, and it contains the channels surrounding the stack and pipes
% connecting osygen vessel with the stack. The channels length
% sqrt(A_PEM*1e-4). Neglecting the channel length and assuming 0.5 m
% connecting pipes, the resistance with the enclosure depends on the radius
% of the pipe. Assuming 1 cm radius of the pipe, the horizontal cylindar
% law of convective heat transfer can then be given as follows: 
% h=7.5 W/(m^2.K).
h_rco=7.5*3.14*0.025*0.5;
R_cond_RECO=1/h_rco;
h_rch=7.5*3.14*0.025*0.35;
R_cond_RECH=1/h_rch;

R_ech=1/(9.3169*0.3);%[K.s/kJ]% considering the values in the electrolyzer.ssc, 

C_th_RECO=0.51*0.001*(3.14*0.025)*0.5*7900+3.14*(0.025/2)^2*0.5*rho_H2O_avg*C_Pm_H2O_avg+...
    C_Pm_H2O_avg*rho_H2O_avg*exchanger_L*exchanger_W*exchanger_tube_H*exchanger_N_tubes;%[kJ/K];
% considering same material used in the manufacturing of the vessel to be
% used in the manufacturing of the pipes. 0.5 m is the length of the pipes,
% 25 mm is the diameter of the pipe and 1mm is the thickness of the pipe.
% For its counterpart for the RECH, the length of the pipe only changes to
% 0.35 m
C_th_RECH=0.51*0.001*(3.14*0.025)*0.35*7900+3.14*(0.025/2)^2*0.35*rho_H2O_avg*C_Pm_H2O_avg;%kJ/k;

C_th_COOL=0.51*0.001*(3.14*0.025)*0.35*7900+3.14*(0.025/2)^2*0.35*rho_H2O_avg*C_Pm_H2O_avg+...
    C_Pm_H2O_avg*rho_H2O_avg*exchanger_L*exchanger_W*exchanger_tube_H*exchanger_N_tubes;%[kJ/K];

% Assume the pipes of the cooling circuit inside the heat exchanger has the
% same size as the pipes of the recirculation circuit which is logical due
% to their similar functionality of exchanging heat. Due to having
% deionized water as coolant, then they would have the same heat capacity.
% Assuming also that the connecting pipes have a length of 0.35 m.

C_th_COLD=0.51*0.003*6*(0.3*0.3)*7900+0.3*0.3*0.3*rho_H2O_avg*C_Pm_H2O_avg;

% The cold unit is simply a tank of the coolant liquid with considerable
% heat capacity so it can maintain significant temperature difference
% facilitating heat transfer. Considering 304L to be the material of the
% tank and that the walling has a thickness of 3 mm, and the tank sizing is
% 0.3 x 0.3 x 0.3.

V_DRY=1e-2;
% A_DRY=V_DRY^(2/3);

% So R_cond_DRY represents the heat resistance of the purification vessel
% circuit with the enclosure unit. It is different from the exchanger
% pipes, and it contains the channels surrounding the stack and pipes
% connecting osygen vessel with the stack. The channels length
% sqrt(A_PEM*1e-4). Neglecting the channel length and assuming 0.5 m
% connecting pipes, the resistance with the enclosure depends on the radius
% of the pipe. Assuming 1 cm radius of the pipe, the horizontal cylindar
% law of convective heat transfer can then be given as follows: 
% h=7.5 W/(m^2.K).
% For 80 degree stack temperature and 25 degree atm temperature, the
% follwing h are obtained:
% Vertical h considering length V_DRY^(1/3), h_v=8 W/m*K,
% Upward h V_DRY^(2/3)=1e-4;
% 4*V_DRY^(1/3)=0.1857
% Upward h=11.2 W/m*K
% downward h=5.6 W/m*K
h_dry=8*V_DRY^(2/3)*4+(11.2+5.6)*V_DRY^(2/3);

R_cond_DRY=1/h_dry;%K/W
C_th_SEPH_void=C_th_SEPO_void;
C_th_DRY=0.51*0.003*V_DRY^(2/3)*6*7900;%kJ/K
% The thickness of the dryer vessel is 3mm.

max_flow_and=0.2;%[kg/s] 
cooling_flow=0.2;%[kg/s] 
cooling_heat=-5;%[kJ/s]





P_SPVA_ref_pres=1e6;% 3 MPa pressure of the dryer vessel

ValvArea=7.0686e-05;% This value represent the maximum area of the valve seperating the Hydrogen storage from the purifier
% The purifier is assumed to have a high pressure of 3 MPa. The pressure
% regulator (i.e. the valve) can open up to 90 % of its area.
R_Dry_Fix=4.457e-7*sqrt(3000000-101325)/5.4224e-4;

w_eva=100;
w_ads=1000;% Fast condesation rate of time constant of 1e-3
C_P_cell=stack_mea_rho * A_PEM*1e-4 * (L_PEM*1e4 + stack_t_gdl_A + stack_t_gdl_C)*1e-6 * stack_num_cells*stack_mea_cp;% J/K
C_p_air=1.012;%kJ/(kg.K)
C_th_Enc=0.51*0.005*5*1*1*7900;% five surfaces surround the electrolyzer.
% Considering 304L stainless steel material in the manufacturing of the
% enclosure with thickness of 5 mm and 1 x 1 x 1 dimensions.
% For heat transfer, as the upper surface is 1 x 1, then h=6.4 upward
% convective heat transfer and the vertical h =5.1 W/m*K. Then total h is
% given as follows
h_enc=6.4+5.1*4;
R_cond_enc=1/h_enc;

% R_recir_and=11; 
% We dont care about P_RECH and P_RECO as they are
% controllable by the pumps assuming perfect
% pumps, so no need to find R_recir_and and R_recir_cat
% R_recir_cat=11;



% lambda_e=11;
C_dl=A_PEM*1000*epsilon_0/L_PEM*1000;

% air density is 1.225 kg/m^3, assume 40 cm diameter of the fan and 5 m/s
% air velocity
% m_dot_fan= density*Area of fan*velocity of air
m_dot_fan=1.225*0.1256*5; 

wc_RES=1/1;
KpL=5000;
KiL=1000;
Kp_rsh=KpL;
Ki_rsh=KiL;
wc_RSH=wc_RES;

wc_RECO=1/1;
% wc_RECH=1/1;

KpRECO=0.4;
KiRECO=0.05;

% KpRECH=0.1;
% KiRECH=0.08;

wc_STRG=1/0.1;

KpSTRG=0.1e-6;
KiSTRG=0.1e-6;


%% Inverter parameters:
Tc=0.04;
wc=1/Tc;
vnom=380;
Vdc=622.26;% Vamp-ac-ph2g=311.13, modulation index 1.
Rf=0.1;
Lf=1.35e-3;
Cf=50e-6;
Kpi=10.5;
Kii=16000;
Kiv=390;
Kpv=0.05;
wnom=2*pi*50;
H=0.75;
Kdcvp=0.85*50;
Kdcvi=50;
Kdcip=0.5;
Kdcii=0.2;
Ls=0.5e-3;
fsw=8000;
Rc=0.03;
Lc=0.35e-3;
Tsim=1e-6;
m1=1.33e-4;
n1=1.33e-3;
nd1=n1*0.1;
md1=m1*0.1;
rn=10000;
C_el=1e-3;
Cs=1e-2;
R_load=10.62;
L_load=0.02032*5;



%% Initial Conditions
E_act_ic=0;
xi_ic=0;
T_stack_ic=T_ref_st;
T_SEPO_ic=273.15+env_T;
T_SEPH_ic=273.15+env_T;
T_RECO_ic=273.15+env_T;
T_COOL_ic=273.15+env_T;
T_COLD_ic=273.15+env_T;
T_RECH_ic=273.15+env_T;
T_DRY_ic=273.15+env_T;
T_Encl_ic=273.15+env_T;
m_SEPO_H2O_ic=0.070182360744771*1e6*V_SEPO*0.25/R/(env_T+273.15)*MW_H2O;
m_SEPO_O2_ic=env_p*1e6*V_SEPO*0.25/R/(env_T+273.15)*MW_O2;
m_SEPH_H2_ic=env_p*1e6*V_SEPH*0.25/R/(env_T+273.15)*MW_H2;
m_SEPH_H2O_ic=0.070182360744771*1e6*V_SEPH*0.25/R/(env_T+273.15)*MW_H2O;
m_H2_DRY_ic=env_p*1e6*V_DRY/R/(env_T+273.15)*MW_H2;
m_H2O_DRY_ic=0;
m_SEPO_ic=(L_SEPO_max+L_SEPO_min)/2*rho_H2O_avg*A_SEPO;
m_SEPH_ic=(L_SEPH_max+L_SEPH_min)/2*rho_H2O_avg*A_SEPH;
m_d_RSO_ic=0;
m_d_RSH_ic=0;
m_d_RECO_ic=0;
m_d_RECH_ic=0;
A_rest_ic=0;
gam_L_SEPO_ic=0;
gam_L_SEPH_ic=0;
gam_T_RECO_ic=0.00;
gam_T_RECH_ic=0.00;
gam_rest_ic=0.0;

% Initials=[E_act_ic;xi_ic;T_stack_ic;T_SEPO_ic;T_SEPH_ic;T_RECO_ic;T_COOL_ic;T_COLD_ic;T_RECH_ic;T_DRY_ic;T_Encl_ic;m_SEPO_H2O_ic;m_SEPO_O2_ic;m_SEPH_H2_ic;m_SEPH_H2O_ic;m_H2_DRY_ic;m_H2O_DRY_ic;m_SEPO_ic;m_SEPH_ic;m_d_RSO_ic;m_d_RSH_ic;m_d_RECO_ic;m_d_RECH_ic;A_rest_ic;gam_L_SEPO_ic;gam_L_SEPH_ic;gam_T_RECO_ic;gam_T_RECH_ic;gam_rest_ic];
%Initials=Initials(1:24);
%% Small-signal model of the electrolyzer

% ssp=struct;
ssp_I_cell=out.I_cell.Data(end);
ssp_T_stack=out.T_stack.Data(end);%T_ref_st;
J_0_ref=J_cell_0;

ssp_I_cell_0=J_0_ref*A_PEM;
% alpha=1;
ssp_R_ohm=out.R_mem.Data(end);
ssp_E_stk=out.E_stack.Data(end);
ssp_E_cell=ssp_E_stk/cells_per_stack;
ssp_E_rev=out.E_rev.Data(end);
ssp_E_act=out.E_act.Data(end);
ssp_E_ss=R/F*ssp_T_stack*log(ssp_I_cell/ssp_I_cell_0);
ssp_P_SEPH_H2=out.P_SEPH_H2.Data(end);
ssp_P_SEPH_H2O=out.P_SEPH_H2O.Data(end);
ssp_P_SEPH=ssp_P_SEPH_H2+ssp_P_SEPH_H2O;
ssp_m_H2O_SEPH=out.m_SEPH_H2O.Data(end);
ssp_m_H2_SEPH=out.m_SEPH_H2.Data(end);

ssp_sigma_PEM=out.sigma_PEM.Data(end);


ssp_P_SEPO_H2O=out.p_SEPO_H2O.Data(end);
ssp_m_H2O_SEPO=out.m_SEPO_H2O.Data(end);
ssp_m_O2_SEPO=out.m_SEPO_O2.Data(end);
ssp_P_SEPO_O2=out.p_SEPO_O2.Data(end);
ssp_P_SEPO=ssp_P_SEPO_O2+ssp_P_SEPO_H2O;
ssp_a_SEPO_H2O=out.R_HO.Data(end);

ssp_P_OXYG=env_p*1e6;

ssp_H_H2O_0=h_w_TLU(12);%kJ/kg
ssp_H_H2_0=H2_h(12);%kJ/kg
ssp_H_O2_0=O2_h(12);%kJ/kg

ssp_H_H_2_g=C_Pm_H2*(ssp_T_stack-273.15)+ssp_H_H2_0;%kJ/kg
ssp_H_O_2_g=C_Pm_O2*(ssp_T_stack-273.15)+ssp_H_O2_0;%kJ/kg
ssp_H_H_2O_g=C_Pm_H2O_gas*(ssp_T_stack-273.15)+ssp_H_H2O_0;%kJ/kg

ssp_m_d_RECH=out.m_dot_RECH.Data(end);
ssp_m_d_RECO=out.m_dot_RECO.Data(end);
ssp_m_d_OXYG=out.m_dot_OXYG.Data(end);
ssp_m_d_SPVA=out.m_dot_SPVA.Data(end);
ssp_m_d_eos=out.m_dot_H2O_eos.Data(end);
ssp_m_d_diff=out.m_dot_H2O_diff.Data(end);



ssp_T_RECH=out.T_RECH.Data(end);
ssp_T_RECO=out.T_RECO_td.Data(end);
ssp_T_ENC=out.T_Encl.Data(end);
ssp_T_RES=env_T+273.15;
ssp_T_RSH=env_T+273.15;



ssp_xi=out.xi.Data(end);




ssp_T_SEPO=out.T_SEPO.Data(end);

ssp_L_SEPO=out.L_SEPO.Data(end);


ssp_T_SEPH=out.T_SEPH.Data(end);

ssp_L_SEPH=out.L_SEPH.Data(end);


ssp_Htot_SEPO=out.Htot_SEPO_td.Data(end);
ssp_C_th_SEPO=out.C_th_SEPO_td.Data(end);

ssp_Htot_SEPH=out.Htot_SEPH_td.Data(end);
ssp_C_th_SEPH=out.C_th_SEPH_td.Data(end);

ssp_P_SPVA_H2O=out.P_SPVA_H2O.Data(end);
ssp_P_SPVA_H2=out.P_SPVA_H2.Data(end);
ssp_P_SPVA=ssp_P_SPVA_H2+ssp_P_SPVA_H2O;

P_PURI=env_p*1e6;

ssp_m_d_PURI=out.m_dot_PURI.Data(end);

ssp_m_H2_DRY=out.m_DRY_H2.Data(end);
ssp_m_H2O_DRY=out.m_DRY_H2O.Data(end);

ssp_T_DRY=out.T_DRY.Data(end);

ssp_m_d_ads=out.ssp_m_d_ads.Data(end);

ssp_A_rest=out.A_rest.Data(end);

ssp_Psat=out.P_sat_SEPO.Data(end);

ssp_m_d_evo=out.m_dot_evaO.Data(end);
ssp_m_d_evh=out.m_dot_evaH.Data(end);
ssp_Dw=out.D_H2O_m.Data(end);
ssp_CH2O_chh=out.Conc_H2O_ccl.Data(end);
ssp_CH2O_cho=out.Conc_H2O_acl.Data(end);
%C_dl=11;
lambda=out.Lambda.Data(end);

ssp_PI_L_SEPO=out.PI_L_SEPO.Data(end);
ssp_PI_L_SEPH=out.PI_L_SEPH.Data(end);
ssp_PI_T_SEPO=out.PI_T_SEPO.Data(end);
%ssp_PI_T_SEPH=out.PI_T_SEPH.Data(end);
ssp_PI_A_rest=out.PI_A_rest.Data(end);

ac_m_rso=double(ssp_PI_L_SEPO>0 && ssp_PI_L_SEPO<1);
ac_m_rsh=double(ssp_PI_L_SEPH>-1 && ssp_PI_L_SEPH<1);
ac_T_spo=double(ssp_PI_T_SEPO>0 && ssp_PI_T_SEPO<1);
%ac_T_sph=double(ssp_PI_T_SEPH>0 && ssp_PI_T_SEPH<1);
ac_A_rst=double(ssp_PI_A_rest>0 && ssp_PI_A_rest<0.9);

ssp_T_COOL=out.T_COOL_td.Data(end);
ssp_T_COLD=out.T_COLD_td.Data(end);

Dbar=ssp_E_stk/Vdc;
Iod1=out.I_od.Data(end);
Ioq1=out.I_oq.Data(end);
Vod1=out.V_od.Data(end);
Voq1=0;
Pssp=Vod1*Iod1;
Qssp=-Vod1*Ioq1;
% d01=0;
Ild1=out.I_ld.Data(end);
Ilq1=out.I_lq.Data(end);
w0=wnom-m1*Pssp;
Iloadq=out.I_loadq.Data(end);
Iloadd=out.I_loadd.Data(end);

Vid1=Vod1+Rf*Ild1-w0*Ilq1*Lf;
Viq1=Voq1+Rf*Ilq1+w0*Ild1*Lf;

Vbd1=Vod1+w0*Ioq1*Lc-Rc*Iod1;
Vbq1=Voq1-w0*Iod1*Lc-Rc*Ioq1;
Pinvbar=-ssp_E_stk*Stacks_No*ssp_I_cell;
Ib=out.i_b.Data(end);



Lin_op=[ssp_E_act;ssp_xi;ssp_T_stack;ssp_T_SEPO;ssp_T_SEPH;ssp_T_RECO;ssp_T_COOL;ssp_T_COLD;ssp_T_RECH;ssp_T_DRY;ssp_T_ENC;ssp_m_H2O_SEPO;ssp_m_O2_SEPO;ssp_m_H2_SEPH;ssp_m_H2O_SEPH;ssp_m_H2_DRY;ssp_m_H2O_DRY;ssp_L_SEPO*rho_H2O_avg*A_SEPO;ssp_L_SEPH*rho_H2O_avg*A_SEPH;ssp_PI_L_SEPO*max_flow_tank;ssp_PI_L_SEPH*max_flow_tank;ssp_m_d_RECO;ssp_A_rest;0;0;0;ssp_E_stk;Ib;Vdc;0;0;Pssp;Qssp;0;0;0;0;Ild1;Ilq1;Vod1;Voq1;Iod1;Ioq1;Iloadd;Iloadq];
%==========================================================================
% small signal subsystem of the PEM electrolysis cell

% Electrical states

z1=1/C_dl*(1-ssp_E_act/ssp_E_ss);% Icel to Eact
z2=-ssp_I_cell/C_dl/ssp_E_ss;% Eact to Eact
z3=ssp_I_cell/C_dl*ssp_E_act/ssp_E_ss^2;%E_ss to Eact

z4=R/F*log(ssp_I_cell/ssp_I_cell_0); % Tstk to E_ss
z5=R/F*ssp_T_stack/ssp_I_cell; % Icel to E_ss

z6=1/ssp_R_ohm;% Ecel to Icel
z7=-1/ssp_R_ohm;% Erev to Icel
z8=-1/ssp_R_ohm;% Eact to Icel
z9=-(ssp_E_cell-ssp_E_rev-ssp_E_act)/(ssp_R_ohm^2);% Rohm to Icel

z10=(1.5421*10^(-3)+9.523*10^(-5)*(log(ssp_T_stack)+1)+2*9.84e-8*ssp_T_stack);% Tstk to Erev (1)
z11=R/2/F*log(ssp_P_SEPH_H2*sqrt(ssp_P_SEPO_O2)/ssp_a_SEPO_H2O/(env_p*1e6)^1.5); % Tstk to Erev (2)

z12=R*ssp_T_stack/2/F/ssp_P_SEPH_H2; % P_SEPH_H2 to Erev
z13=R*ssp_T_stack/4/F/ssp_P_SEPO_O2; % P_SEPO_O2 to Erev
z14=-R*ssp_T_stack/2/F/ssp_P_SEPO_H2O; % P_SEPO_H2O to Erev

z15=0.005139*exp(1268*(1/303-1/ssp_T_stack)); %lambda to sigma
z16=((0.005139*lambda-0.00326)*exp(1268*(1/303-1/ssp_T_stack))*1268/ssp_T_stack^2); %Tstk to sigma

z17=0.5*(17.81/(ssp_Psat)-2*39.85*ssp_P_SEPO_H2O/(ssp_Psat)^2+3*36*ssp_P_SEPO_H2O^2/(ssp_Psat)^3);% P_SEPO_H2O to lambda
z18=0.5*(17.81/(ssp_Psat)-2*39.85*ssp_P_SEPH_H2O/(ssp_Psat)^2+3*36*ssp_P_SEPH_H2O^2/(ssp_Psat)^3);% P_SEPH_H2O to lambda

z19=-L_PEM/(A_PEM*ssp_sigma_PEM^2); % sigma to Rohm


a_pem_e_21=z8/(2*F*tau_e);
a_pem_e_22=-1/tau_e;%##
b_pem_e_21=1/(2*F*tau_e)*(z7*(z10+z11)+z9*(z19*z16));

b_pem_e_22=z7/(2*F*tau_e)*z12;
b_pem_e_23=z7/(2*F*tau_e)*z13;
b_pem_e_24=1/(2*F*tau_e)*(z7*z14+z9*(z19*(z15*z17)));

b_pem_e_25=1/(2*F*tau_e)*z9*(z19*(z15*z18));
b_pem_e_26=z6/(2*F*tau_e);

c_pem_e_21=z8;

d_pem_e_21=b_pem_e_21*(2*F*tau_e);%##
d_pem_e_22=b_pem_e_22*(2*F*tau_e);%##
d_pem_e_23=b_pem_e_23*(2*F*tau_e);%##
d_pem_e_24=b_pem_e_24*(2*F*tau_e);%##
d_pem_e_25=b_pem_e_25*(2*F*tau_e);%##
d_pem_e_26=b_pem_e_26*(2*F*tau_e);%##



a_pem_e_11=z2+z1*z8+z3*(z5*z8);%made_up_error /ssp_E_act is added


b_pem_e_11=z3*(z4+z5*(z7*(z10+z11)+z9*(z19*z16)))+z1*(z7*(z10+z11)+z9*(z19*z16));

b_pem_e_12=z3*z5*z7*z12+z1*z7*z12;
b_pem_e_13=z3*z5*z7*z13+z1*z7*z13;
b_pem_e_14=z3*z5*(z7*z14+z9*z19*(z15*z17))+z1*(z7*z14+z9*z19*(z15*z17));

b_pem_e_15=z3*z5*z9*(z19*(z15*z18))+z1*z9*(z19*(z15*z18));

b_pem_e_16=z3*z5*z6+z1*z6;


% thermal states:
e1=1/1000/R_th_cell; % Tstk>> >>Qstk

e2=0.0134*stack_num_cells*MW_H2O*ssp_xi;% Tstk>> >>mdeos

e3=stack_num_cells*MW_H2O*(0.0134*ssp_T_stack+0.03);% xi>> >>mdeos

e4=-(ssp_T_stack-273.15)*C_Pm_H2O_avg;% mdeos>> >>Hand
 
%-e4;% mdrco>> >>Hand

% e4 % mdH2Ocon>> >>Hand

e5= C_Pm_H2O_gas*(ssp_T_stack-273.15)+ssp_H_H2O_0;% mdH2Odif>> >>Hand

% ssp_H_O_2_g % mdO2gen>> >>Hand

e6=C_Pm_H2O_avg*(ssp_m_d_RECO-ssp_m_d_eos-ssp_xi*MW_H2O*v_H2O*stack_num_cells);% Tstk>> >>Hand (1) water

e7=C_Pm_H2O_gas*ssp_m_d_diff;% Tstk>> >>Hand (2) gas

e8=C_Pm_O2*ssp_xi*MW_O2*v_O2*stack_num_cells;% Tstk>> >>Hand (3) O2 

% +++++++++

% e5 % mdeos>> >>Hcat

% -e4 % mdrch>> >>Hcat

% e4 % mdH2Odif>> >>Hcat

% ssp_H_H_2_g % mdH2gen>> >>Hcat

e9=C_Pm_H2O_avg*(ssp_m_d_RECH-ssp_m_d_diff);% Tstk>> >>Hcat (1) water

e10=C_Pm_H2O_gas*ssp_m_d_eos;% Tstk>> >>Hcat (2) gas


e11=C_Pm_H2*ssp_xi*MW_H2*v_H2*stack_num_cells;% Tstk>> >>Hcat (3) H2


a_pem_t_11=1000/C_P_cell*(-e9-e10-e5*e2-e6-e7-e4*e2-e1-e11-e8);%

b_pem_t_18=1000/C_P_cell*C_Pm_H2O_avg*ssp_m_d_RECH;%##

b_pem_t_19=1000/C_P_cell*C_Pm_H2O_avg*ssp_m_d_RECO;%##

b_pem_t_110=1/C_P_cell/R_th_cell;%##

b_pem_t_11=1000/C_P_cell*C_Pm_H2O_avg*(ssp_T_RECH-ssp_T_stack);%##

b_pem_t_12=1000/C_P_cell*C_Pm_H2O_avg*(ssp_T_RECO-ssp_T_stack);%##

% b_pem_t_13=1000/C_P_cell*C_Pm_H2O_avg*(ssp_T_stack-273.15)-1/C_P_cell*ssp_H_H_2O_g;

b_pem_t_13=-1000/C_P_cell*e5-1000/C_P_cell*e4;

b_pem_t_15=1000*stack_num_cells/C_P_cell*(-e4*MW_H2O*v_H2O-v_H2*MW_H2*ssp_H_H_2_g-v_O2*MW_O2*ssp_H_O_2_g)-1000/C_P_cell*e3*(e4+ssp_H_H_2O_g);%made_up_error - removed

b_pem_t_16=stack_num_cells/C_P_cell*ssp_E_cell;%##

b_pem_t_17=stack_num_cells/C_P_cell*ssp_I_cell;%##


% gaseous fluidic
d_pem_g_1=-stack_num_cells*ssp_Dw*A_PEM*0.0001*MW_H2O*rho_membrane/(0.01*L_PEM)/M_membrane*z17;
d_pem_g_2=stack_num_cells*ssp_Dw*A_PEM*0.0001*MW_H2O*rho_membrane/(0.01*L_PEM)/M_membrane*z18;
d_pem_g_3=stack_num_cells*A_PEM*0.0001*MW_H2O/(0.01*L_PEM)*(ssp_CH2O_chh-ssp_CH2O_cho)*1.25e-10*exp(2416/303.15-2416/ssp_T_stack)/ssp_T_stack^2*2416;

%==========================================================================
% small signal subsystem of the Oxygen separator vessel

% chemical states:
%[m_SEPO_O2;m_SEPO_H2O]##

a_SEPO_g_11=-ssp_m_O2_SEPO*ssp_m_d_OXYG/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2;%##
a_SEPO_g_12=ssp_m_H2O_SEPO*ssp_m_d_OXYG/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2;%##


a_SEPO_g_21= ssp_m_O2_SEPO*ssp_m_d_OXYG/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2;%##
a_SEPO_g_22=-(ssp_m_H2O_SEPO)*ssp_m_d_OXYG/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2;%##


b_SEPO_g_13=-ssp_m_H2O_SEPO/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO);%##
b_SEPO_g_23=-ssp_m_O2_SEPO/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO);%##

b_SEPO_g_24=stack_num_cells*MW_O2*v_O2;%##

c_SEPO_g_11=R*ssp_T_SEPO/(MW_H2O*(V_SEPO-A_SEPO*ssp_L_SEPO));%##
c_SEPO_g_22=R*ssp_T_SEPO/(MW_O2*(V_SEPO-A_SEPO*ssp_L_SEPO));%##

d_SEPO_g_11=R*ssp_T_SEPO*ssp_m_H2O_SEPO/(MW_H2O*(V_SEPO-A_SEPO*ssp_L_SEPO)^2*rho_H2O_avg);%##
d_SEPO_g_21=R*ssp_T_SEPO*ssp_m_O2_SEPO/(MW_O2*(V_SEPO-A_SEPO*ssp_L_SEPO)^2*rho_H2O_avg);%##

d_SEPO_g_12=R*ssp_m_H2O_SEPO/(MW_H2O*(V_SEPO-A_SEPO*ssp_L_SEPO));%##
d_SEPO_g_22=R*ssp_m_O2_SEPO/(MW_O2*(V_SEPO-A_SEPO*ssp_L_SEPO));%##


c_SEPO_g_33=-w_eva;%##
d_SEPO_g_31=-w_eva*ssp_Psat*MW_H2O/R/ssp_T_SEPO/rho_H2O_avg;%##
d_SEPO_g_32=-w_eva*ssp_Psat*MW_H2O/R/(ssp_T_SEPO^2)*(V_SEPO-A_SEPO*ssp_L_SEPO);%##


% thermal state:
% T_SEPO 


b_SEPO_t_11=1/ssp_C_th_SEPO*(e6+e4*e2+e7+e8);
a_SEPO_t_11=-1/ssp_C_th_SEPO*(ssp_m_d_RECO*C_Pm_H2O_avg+0.001/R_cond_SEPO+ssp_m_d_OXYG*(ssp_m_O2_SEPO*C_Pm_O2+ssp_m_H2O_SEPO*C_Pm_H2O_gas)/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)+ssp_m_d_evo*(C_Pm_H2O_gas-C_Pm_H2O_avg));%##
b_SEPO_t_12=0.001/ssp_C_th_SEPO/R_cond_SEPO;%##

b_SEPO_t_13=1/ssp_C_th_SEPO*ssp_T_stack*C_Pm_H2O_avg-1/ssp_C_th_SEPO*ssp_T_SEPO*C_Pm_H2O_avg;%##
b_SEPO_t_14=1/ssp_C_th_SEPO*(ssp_T_RES-273.15)*C_Pm_H2O_avg;%##

b_SEPO_t_15=-1/ssp_C_th_SEPO*((ssp_m_O2_SEPO*C_Pm_O2+ssp_m_H2O_SEPO*C_Pm_H2O_gas)*(ssp_T_SEPO-273.15)+ssp_m_O2_SEPO*ssp_H_O2_0+ssp_m_H2O_SEPO*ssp_H_H2O_0)/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO);%##

b_SEPO_t_18=1/ssp_C_th_SEPO*e4*(e3+stack_num_cells*MW_H2O*v_H2O)+1/ssp_C_th_SEPO*ssp_H_O_2_g*stack_num_cells*MW_O2*v_O2;

b_SEPO_t_16=1/ssp_C_th_SEPO*e5;

b_SEPO_t_19=-ssp_Htot_SEPO/(ssp_C_th_SEPO)^2*C_Pm_H2O_gas+1/ssp_C_th_SEPO*ssp_m_d_OXYG*ssp_m_O2_SEPO/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2*(C_Pm_O2*(ssp_T_SEPO-273.15)+ssp_H_O2_0-C_Pm_H2O_gas*(ssp_T_SEPO-273.15)-ssp_H_H2O_0);%##
b_SEPO_t_110=-ssp_Htot_SEPO/(ssp_C_th_SEPO)^2*C_Pm_O2+1/ssp_C_th_SEPO*ssp_m_d_OXYG*ssp_m_H2O_SEPO/(ssp_m_O2_SEPO+ssp_m_H2O_SEPO)^2*(C_Pm_H2O_gas*(ssp_T_SEPO-273.15)+ssp_H_H2O_0-C_Pm_O2*(ssp_T_SEPO-273.15)-ssp_H_O2_0);%##

b_SEPO_t_111=-ssp_Htot_SEPO/(ssp_C_th_SEPO)^2*C_Pm_H2O_avg;%##

b_SEPO_t_17=-1/ssp_C_th_SEPO*(-C_Pm_H2O_avg*(ssp_T_SEPO-273.15)+C_Pm_H2O_gas*(ssp_T_SEPO-273.15)+ssp_H_H2O_0);%##

% Fluidic states:##


a_SEPO_f_21=-wc_RES*max_flow_tank*KpL/A_SEPO/rho_H2O_avg*ac_m_rso;  %##
a_SEPO_f_22=-wc_RES;%##
a_SEPO_f_23=wc_RES*max_flow_tank*KiL*ac_m_rso;%##

a_SEPO_f_31=-1/A_SEPO/rho_H2O_avg*ac_m_rso;  %##

b_SEPO_f_11=-MW_H2O*ssp_xi*0.0134*stack_num_cells;%##
b_SEPO_f_12=-stack_num_cells*MW_H2O*(0.03+0.0134*ssp_T_stack+v_H2O);%##

% output ##

d_SEPO_f_11=1/(2*R_SEP_AND*sqrt(ssp_P_SEPO-ssp_P_OXYG))*sign(ssp_m_d_OXYG);%##
d_SEPO_f_12=d_SEPO_f_11;%##

%==========================================================================
% small signal subsystem of the Hydrogen separator vessel

% chemical states:##

a_SEPH_g_11=-ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;%##
a_SEPH_g_12=ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;%##

a_SEPH_g_21=ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;%##
a_SEPH_g_22=-ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;%##


b_SEPH_g_11=-(ssp_m_H2_SEPH)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);%##
b_SEPH_g_21=-(ssp_m_H2O_SEPH)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);%##

b_SEPH_g_14=stack_num_cells*MW_H2*v_H2;%##
b_SEPH_g_24=stack_num_cells*MW_H2O*(0.03+0.0134*ssp_T_stack);%##

b_SEPH_g_23=stack_num_cells*MW_H2O*ssp_xi*0.0134;%##

c_SEPH_g_11=R*ssp_T_SEPH/(MW_H2*(V_SEPH-A_SEPH*ssp_L_SEPH));%##
c_SEPH_g_22=R*ssp_T_SEPH/(MW_H2O*(V_SEPH-A_SEPH*ssp_L_SEPH));%##

d_SEPH_g_11=R*ssp_T_SEPH*ssp_m_H2_SEPH/(MW_H2*(V_SEPH-A_SEPH*ssp_L_SEPH)^2*rho_H2O_avg);%##
d_SEPH_g_21=R*ssp_T_SEPH*ssp_m_H2O_SEPH/(MW_H2O*(V_SEPH-A_SEPH*ssp_L_SEPH)^2*rho_H2O_avg);%##

d_SEPH_g_12=R*ssp_m_H2_SEPH/(MW_H2*(V_SEPH-A_SEPH*ssp_L_SEPH));%##
d_SEPH_g_22=R*ssp_m_H2O_SEPH/(MW_H2O*(V_SEPH-A_SEPH*ssp_L_SEPH));%##

c_SEPH_g_33=-w_eva;%##
d_SEPH_g_31=-w_eva*ssp_Psat*MW_H2O/R/ssp_T_SEPH/rho_H2O_avg;%##
d_SEPH_g_32=-w_eva*ssp_Psat*MW_H2O/R/ssp_T_SEPH^2*(V_SEPH-A_SEPH*ssp_L_SEPH);%##

% thermal state:
% T^SEPH


a_SEPH_t_11=-1/ssp_C_th_SEPH*(ssp_m_d_RECH*C_Pm_H2O_avg+0.001/R_cond_SEPH+ssp_m_d_SPVA*(ssp_m_H2O_SEPH*C_Pm_H2O_gas+ssp_m_H2_SEPH*C_Pm_H2)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH))-1/ssp_C_th_SEPH*(C_Pm_H2O_gas-C_Pm_H2O_avg)*ssp_m_d_evh;%##

b_SEPH_t_11=1/ssp_C_th_SEPH*(e9+e10+e11+e5*e2);
b_SEPH_t_12=0.001/ssp_C_th_SEPH/R_cond_SEPH;%##

b_SEPH_t_13=1/ssp_C_th_SEPH*ssp_T_stack*C_Pm_H2O_avg-1/ssp_C_th_SEPH*ssp_T_SEPH*C_Pm_H2O_avg;%##
b_SEPH_t_14=1/ssp_C_th_SEPH*(ssp_T_RSH-273.15)*C_Pm_H2O_avg;%##

b_SEPH_t_15=-1/ssp_C_th_SEPH*((ssp_m_H2O_SEPH*C_Pm_H2O_gas+ssp_m_H2_SEPH*C_Pm_H2)*(ssp_T_SEPH-273.15)+ssp_m_H2_SEPH*ssp_H_H2_0+ssp_m_H2O_SEPH*ssp_H_H2O_0)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);%##

b_SEPH_t_16=-1/ssp_C_th_SEPH*(ssp_T_stack-273.15)*C_Pm_H2O_avg;%##

b_SEPH_t_17=-1/ssp_C_th_SEPH*(-C_Pm_H2O_avg*(ssp_T_SEPH-273.15)+C_Pm_H2O_gas*(ssp_T_SEPH-273.15)+ssp_H_H2O_0);%##

b_SEPH_t_18=1/ssp_C_th_SEPH*e5*e3+1/ssp_C_th_SEPH*(ssp_H_H_2_g)*stack_num_cells*MW_H2*v_H2;

b_SEPH_t_19=-ssp_Htot_SEPH/(ssp_C_th_SEPH)^2*C_Pm_H2+1/ssp_C_th_SEPH*ssp_m_d_SPVA*ssp_m_H2O_SEPH/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2*(C_Pm_H2O_gas*(ssp_T_SEPH-273.15)+ssp_H_H2O_0-C_Pm_H2*(ssp_T_SEPH-273.15)-ssp_H_H2_0);%##
b_SEPH_t_110=-ssp_Htot_SEPH/(ssp_C_th_SEPH)^2*C_Pm_H2O_gas+1/ssp_C_th_SEPH*ssp_m_d_SPVA*ssp_m_H2_SEPH/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2*(C_Pm_H2*(ssp_T_SEPH-273.15)+ssp_H_H2_0-C_Pm_H2O_gas*(ssp_T_SEPH-273.15)-ssp_H_H2O_0);%##

b_SEPH_t_111=-ssp_Htot_SEPH/(ssp_C_th_SEPH)^2*C_Pm_H2O_avg;%##

% Fluidic states:
% m^SEPH ##

a_SEPH_f_21=-wc_RSH*max_flow_tank*Kp_rsh/A_SEPH/rho_H2O_avg*ac_m_rsh;  %##
a_SEPH_f_22=-wc_RSH;%##
a_SEPH_f_23=wc_RSH*max_flow_tank*Ki_rsh*ac_m_rsh;%##

a_SEPH_f_31=-1/A_SEPH/rho_H2O_avg*ac_m_rsh;  %##

% output ##
d_SEPH_f_11=1/(2*sqrt(abs(ssp_P_SEPH_H2+ssp_P_SEPH_H2O-ssp_P_SPVA_H2O-ssp_P_SPVA_H2)*R_SEP_CAT))*sign(ssp_m_d_SPVA);%##
d_SEPH_f_12=d_SEPH_f_11;%##
d_SEPH_f_13=-d_SEPH_f_11;%##
d_SEPH_f_14=-d_SEPH_f_11;%##
%==========================================================================
% small-signal subsysem of the RECO

% thermal states:##

a_RECO_t_11=-1/C_th_RECO*(C_Pm_H2O_avg*ssp_m_d_RECO+0.001/R_cond_RECO+1/R_ech);%##
a_RECO_t_12=1/C_th_RECO*1/R_ech;%##

a_RECO_t_21=1/C_th_COOL/R_ech;%##
a_RECO_t_22=-1/C_th_COOL*(C_Pm_H2O_avg*cooling_flow+1/R_ech);%##
a_RECO_t_23=1/C_th_COOL*C_Pm_H2O_avg*cooling_flow;%##

a_RECO_t_32=1/C_th_COLD*C_Pm_H2O_avg*cooling_flow;%##
a_RECO_t_33=-1/C_th_COLD*C_Pm_H2O_avg*cooling_flow;%##

b_RECO_t_11=1/C_th_RECO*C_Pm_H2O_avg*ssp_m_d_RECO;%##
b_RECO_t_12=0.001/C_th_RECO/R_cond_RECO;%##

b_RECO_t_13=1/C_th_RECO*C_Pm_H2O_avg*(ssp_T_SEPO-ssp_T_RECO);%##

% fluidic states##

a_RECO_f_11=-wc_RECO;%##
a_RECO_f_12=wc_RECO*max_flow_and*KiRECO*ac_T_spo;%##
b_RECO_f_11=wc_RECO*max_flow_and*KpRECO*ac_T_spo;%##

%==========================================================================
% small-signal subsysem of the RECH
% thermal states:##

a_RECH_t_11=-1/C_th_RECH*(C_Pm_H2O_avg*ssp_m_d_RECH+0.001/R_cond_RECH);%##


b_RECH_t_11=1/C_th_RECH*C_Pm_H2O_avg*ssp_m_d_RECH;%##
b_RECH_t_12=0.001/C_th_RECH/R_cond_RECH;%##

b_RECH_t_13=1/C_th_RECH*C_Pm_H2O_avg*(ssp_T_SEPH-ssp_T_RECH);%##

% small-signal subsystem of the purification unit

% thermal states
% T^DRY ##

a_DRY_t_11=1/C_th_DRY*((C_Pm_H2O_gas-2*C_Pm_H2O_avg)*ssp_m_d_ads+...
                        -ssp_m_d_PURI*(ssp_m_H2O_DRY*C_Pm_H2O_gas+ssp_m_H2_DRY*C_Pm_H2)/(ssp_m_H2O_DRY+ssp_m_H2_DRY)+...
                       -0.001/R_cond_DRY...
                       );%##

b_DRY_t_11=1/C_th_DRY*ssp_m_d_SPVA*(ssp_m_H2O_SEPH*C_Pm_H2O_gas+ssp_m_H2_SEPH*C_Pm_H2)/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);%##
b_DRY_t_12=0.001/C_th_DRY/R_cond_DRY;%##

b_DRY_t_13=1/C_th_DRY*ssp_m_H2O_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2*(C_Pm_H2O_gas*(ssp_T_DRY-273.15)+ssp_H_H2O_0-C_Pm_H2*(ssp_T_DRY-273.15)-ssp_H_H2_0);%##
b_DRY_t_14=1/C_th_DRY*ssp_m_H2_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2*(C_Pm_H2*(ssp_T_DRY-273.15)+ssp_H_H2_0-C_Pm_H2O_gas*(ssp_T_DRY-273.15)-ssp_H_H2O_0);%##

b_DRY_t_15=-1/C_th_DRY*ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2*(C_Pm_H2O_gas*(ssp_T_SEPH-273.15)+ssp_H_H2O_0-C_Pm_H2*(ssp_T_SEPH-273.15)-ssp_H_H2_0);%##
b_DRY_t_16=-1/C_th_DRY*ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2*(C_Pm_H2*(ssp_T_SEPH-273.15)+ssp_H_H2_0-C_Pm_H2O_gas*(ssp_T_SEPH-273.15)-ssp_H_H2O_0);%##

b_DRY_t_17=-1/C_th_DRY*((ssp_m_H2O_DRY*C_Pm_H2O_gas+ssp_m_H2_DRY*C_Pm_H2)*(ssp_T_DRY-273.15)+ssp_m_H2_DRY*ssp_H_H2_0+ssp_m_H2O_DRY*ssp_H_H2O_0)/(ssp_m_H2O_DRY+ssp_m_H2_DRY);%##
b_DRY_t_18=1/C_th_DRY*b_SEPH_t_15*-ssp_C_th_SEPH;%##
b_DRY_t_19=1/C_th_DRY*((C_Pm_H2O_gas-2*C_Pm_H2O_avg)*(ssp_T_DRY-273.15)+ssp_H_H2O_0);%##


% chemical
% [m_DRY_O2;m_DRY_H2;m_DRY_H2O]##

a_DRY_g_11=-ssp_m_H2O_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2;%##
a_DRY_g_12=ssp_m_H2_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2;%##

a_DRY_g_21=ssp_m_H2O_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2;%##
a_DRY_g_22=-ssp_m_H2_DRY*ssp_m_d_PURI/(ssp_m_H2O_DRY+ssp_m_H2_DRY)^2;%##

b_DRY_g_11=ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;%##
b_DRY_g_12=-ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;%##

b_DRY_g_21=-ssp_m_H2O_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;%##
b_DRY_g_22=ssp_m_H2_SEPH*ssp_m_d_SPVA/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH)^2;%##


b_DRY_g_13=ssp_m_H2_SEPH/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);%##
b_DRY_g_23=ssp_m_H2O_SEPH/(ssp_m_H2O_SEPH+ssp_m_H2_SEPH);%##

b_DRY_g_14=-ssp_m_H2_DRY/(ssp_m_H2O_DRY+ssp_m_H2_DRY);%##
b_DRY_g_24=-ssp_m_H2O_DRY/(ssp_m_H2O_DRY+ssp_m_H2_DRY);%##

c_DRY_g_11=R*ssp_T_DRY/(MW_H2*V_DRY);%##
c_DRY_g_22=R*ssp_T_DRY/(MW_H2O*V_DRY);%##

d_DRY_g_11=R*ssp_m_H2_DRY/(MW_H2*V_DRY);%##
d_DRY_g_21=R*ssp_m_H2O_DRY/(MW_H2O*V_DRY);%##


% small-signal subsystem of the unit enclosure

% thermal states:##

a_Enc_t_11=-0.001/C_th_Enc*(1/R_th_cell+1/R_cond_RECH+1/R_cond_DRY+1/R_cond_RECO+1/R_cond_SEPH+1/R_cond_SEPO+1/R_cond_enc+1000*C_p_air*m_dot_fan);%##
b_Enc_t_11=0.001/C_th_Enc*1/R_th_cell;%##
b_Enc_t_12=0.001/C_th_Enc*1/R_cond_SEPO;%##
b_Enc_t_13=0.001/C_th_Enc*1/R_cond_SEPH;%##
b_Enc_t_14=0.001/C_th_Enc*1/R_cond_RECO;%##
b_Enc_t_15=0.001/C_th_Enc*1/R_cond_RECH;%##
b_Enc_t_16=0.001/C_th_Enc*1/R_cond_DRY;%##

b_Enc_t_17=1/C_th_Enc*(C_p_air*m_dot_fan+1/R_cond_enc);%##

% small-signal subsystem of the Hydrogen storage tank##

a_STRG_g_11=-wc_STRG;%##
a_STRG_g_12=wc_STRG*KiSTRG*ValvArea*ac_A_rst;%##


b_STRG_g_11=wc_STRG*KpSTRG*ValvArea*ac_A_rst;%##
b_STRG_g_12=wc_STRG*KpSTRG*ValvArea*ac_A_rst;%##

c_STRG_g_11=sqrt(ssp_P_SPVA-P_PURI)/R_Dry_Fix*sign(ssp_m_d_PURI);%##

d_STRG_g_11=ssp_A_rest/R_Dry_Fix/2/sqrt(ssp_P_SPVA-P_PURI)*sign(ssp_m_d_PURI);%##
d_STRG_g_12=d_STRG_g_11;%##


%% small-signal complete model


a01=[a_pem_e_11	0	        b_pem_e_11	b_pem_e_13*d_SEPO_g_22+b_pem_e_14*d_SEPO_g_12	b_pem_e_12*d_SEPH_g_12+b_pem_e_15*d_SEPH_g_22	0	0	0	0	0	0	b_pem_e_14*c_SEPO_g_11	b_pem_e_13*c_SEPO_g_22	b_pem_e_12*c_SEPH_g_11	b_pem_e_15*c_SEPH_g_22	0	0	b_pem_e_13*d_SEPO_g_21+b_pem_e_14*d_SEPO_g_11	b_pem_e_12*d_SEPH_g_11+b_pem_e_15*d_SEPH_g_21	0	0	0	0	0	0	0	0];%##
a02=[a_pem_e_21	a_pem_e_22	b_pem_e_21	b_pem_e_23*d_SEPO_g_22+b_pem_e_24*d_SEPO_g_12	b_pem_e_22*d_SEPH_g_12+b_pem_e_25*d_SEPH_g_22	0	0	0	0	0	0	b_pem_e_24*c_SEPO_g_11	b_pem_e_23*c_SEPO_g_22	b_pem_e_22*c_SEPH_g_11	b_pem_e_25*c_SEPH_g_22	0	0	b_pem_e_23*d_SEPO_g_21+b_pem_e_24*d_SEPO_g_11	b_pem_e_22*d_SEPH_g_11+b_pem_e_25*d_SEPH_g_21	0	0	0	0	0	0	0	0];%##

a03=[b_pem_t_16*c_pem_e_21	b_pem_t_15	a_pem_t_11+b_pem_t_13*d_pem_g_3+b_pem_t_16*d_pem_e_21	b_pem_t_13*d_pem_g_1*d_SEPO_g_12+b_pem_t_16*d_pem_e_23*d_SEPO_g_22+b_pem_t_16*d_pem_e_24*d_SEPO_g_12	b_pem_t_13*d_pem_g_2*d_SEPH_g_22+b_pem_t_16*d_pem_e_22*d_SEPH_g_12+b_pem_t_16*d_pem_e_25*d_SEPH_g_22	b_pem_t_19	0	0	b_pem_t_18	0	b_pem_t_110	b_pem_t_13*d_pem_g_1*c_SEPO_g_11+b_pem_t_16*d_pem_e_24*c_SEPO_g_11	b_pem_t_16*d_pem_e_23*c_SEPO_g_22	b_pem_t_16*d_pem_e_22*c_SEPH_g_11	b_pem_t_13*d_pem_g_2*c_SEPH_g_22+b_pem_t_16*d_pem_e_25*c_SEPH_g_22	0	0	b_pem_t_13*d_pem_g_1*d_SEPO_g_11+b_pem_t_16*d_pem_e_23*d_SEPO_g_21+b_pem_t_16*d_pem_e_24*d_SEPO_g_11	b_pem_t_13*d_pem_g_2*d_SEPH_g_21+b_pem_t_16*d_pem_e_22*d_SEPH_g_11+b_pem_t_16*d_pem_e_25*d_SEPH_g_21	0	0	b_pem_t_12+b_pem_t_11   0	0	0	0	0];%##

a04=[0	b_SEPO_t_18	b_SEPO_t_11+b_SEPO_t_16*d_pem_g_3	a_SEPO_t_11+b_SEPO_t_15*d_SEPO_f_11*d_SEPO_g_12+b_SEPO_t_15*d_SEPO_f_12*d_SEPO_g_22+b_SEPO_t_16*d_pem_g_1*d_SEPO_g_12+b_SEPO_t_17*d_SEPO_g_32	b_SEPO_t_16*d_pem_g_2*d_SEPH_g_22	0	0	0	0	0	b_SEPO_t_12	b_SEPO_t_19+b_SEPO_t_15*d_SEPO_f_11*c_SEPO_g_11+b_SEPO_t_16*d_pem_g_1*c_SEPO_g_11+b_SEPO_t_17*c_SEPO_g_33	b_SEPO_t_110+b_SEPO_t_15*d_SEPO_f_12*c_SEPO_g_22	0	b_SEPO_t_16*d_pem_g_2*c_SEPH_g_22	0	0	b_SEPO_t_111+b_SEPO_t_15*d_SEPO_f_11*d_SEPO_g_11+b_SEPO_t_15*d_SEPO_f_12*d_SEPO_g_21+b_SEPO_t_16*d_pem_g_1*d_SEPO_g_11+b_SEPO_t_17*d_SEPO_g_31	b_SEPO_t_16*d_pem_g_2*d_SEPH_g_21	b_SEPO_t_14	0	b_SEPO_t_13	0	0	0	0	0];%##

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

A_ele=[a01;a02;a03;a04;a05;a06;a07;a08;a09;a10;a11;a12;a13;a14;a15;a16;a17;a18;a19;a20;a21;a22;a24;a25;a26;a27;a29];


B_ele=[b_pem_e_16;b_pem_e_26;b_pem_t_17+b_pem_t_16*d_pem_e_26;zeros(size(A_ele,2)-3,1)];
%C_ele=eye(size(A_ele,2));
C_icl=[c_pem_e_21	0	d_pem_e_21	d_pem_e_23*d_SEPO_g_22+d_pem_e_24*d_SEPO_g_12	d_pem_e_22*d_SEPH_g_12+d_pem_e_25*d_SEPH_g_22	0	0	0	0	0	0	d_pem_e_24*c_SEPO_g_11	d_pem_e_23*c_SEPO_g_22	d_pem_e_22*c_SEPH_g_11	d_pem_e_25*c_SEPH_g_22	0	0	d_pem_e_23*d_SEPO_g_21+d_pem_e_24*d_SEPO_g_11	d_pem_e_22*d_SEPH_g_11+d_pem_e_25*d_SEPH_g_21	0	0	0	0	0	0	0	0];

D_icl=d_pem_e_26;

EigenEle=eig(A_ele);




%% Inverter and Buck Boost converter linearization:


A_dc_lc=[0            1/C_el                0
        -1/Ls         0                     Dbar/Ls
        0            -Dbar/Cs               Pinvbar/(Vdc^2*Cs)];

B_dc_lc_1=[-Stacks_No/C_el
           0
           0];%##


B_dc_lc_2=[0;Vdc/Ls;-Ib/Cs];%##

B_dc_lc_3=[0;0;-1/Cs/Vdc];


A_dc_ct=[0       0
         Kdcvi   0];%##

B_dc_ct_1=[0 0 1
           0 0 Kdcvp];%##

B_dc_ct_2=[0
           -Stacks_No];%##

C_dc_ct=[Kdcip*Kdcvi Kdcii];%##

D_dc_ct_1=[0 0 Kdcip*Kdcvp];%##
D_dc_ct_2=-Stacks_No*Kdcip;%##


A_DC=[A_dc_lc+B_dc_lc_2*D_dc_ct_1          B_dc_lc_2*C_dc_ct
     B_dc_ct_1                          A_dc_ct];%##

B_DC_cl=[B_dc_lc_1+B_dc_lc_2*D_dc_ct_2
      B_dc_ct_2];%##

B_DC_p=[B_dc_lc_3
        zeros(2,1)];%##

C_DC_ec=[1/cells_per_stack 0 0 0 0];

C_DC_ed=[0 0 1 0 0];%##

C_DC=[C_DC_ec
    C_DC_ed];%##
%++++++++++++++++++++++++++++++++++
Ap=[-wc 0
    0 -wc];%##
Bp_1=[0     0   wc*Iod1           wc*Ioq1         wc*Vod1           wc*Voq1
      0     0   -wc*Ioq1          wc*Iod1         wc*Voq1           -wc*Vod1];%##
Cpw_1=[-m1+md1*wc 0];%##
Cpv_1=[0 -n1+nd1*wc
       0 0];%##
Dpw_1=[0     0   -md1*wc*Iod1     -md1*wc*Ioq1     -md1*wc*Vod1  -md1*wc*Voq1];%##

Dpv_1=[ 0     0   nd1*wc*Ioq1     -nd1*wc*Iod1    -nd1*wc*Voq1    nd1*wc*Vod1
        0     0   0                 0               0                0       ];%##

Bv1=[1 0
     0 1];%##      % Voltage Controller Integral states with the references voltages 
 
Bv2=[0 0 -1 0 0 0
     0 0 0 -1 0 0];%##% Voltage Controller Integral states with the voltage and current measurements 

Cv=[Kiv 0
    0 Kiv];%##     %Inductor current references with the integral states 

Dv1=[Kpv 0
     0 Kpv];%##    %Inductor current references with the references voltages
 
Dv2=[0 0 -Kpv       -wnom*Cf    H 0
     0 0 wnom*Cf    -Kpv        0 H];%##%Inductor current references with the voltage and current measurements

Bc1=[1 0
     0 1];%##      % Current Controller Integral states with the references currents 
 
Bc2=[-1 0 0 0 0 0
     0 -1 0 0 0 0];%##% Current Controller Integral states with the voltage and current measurements 

Cc=[Kii 0
    0   Kii];%##     %Input voltage references with the integral states 

Dc1=[Kpi 0
     0   Kpi];%##    %Input voltage references with the references inductor currents 
 
Dc2=[-Kpi       -wnom*Lf 0 0 0 0
     wnom*Lf    -Kpi     0 0 0 0];%## %Input voltage references with the voltage and current measurements

ALCL=[-Rf/Lf  w0     -1/Lf  0       0               0               
      -w0     -Rf/Lf 0      -1/Lf   0               0               
      1/Cf    0      0      w0      -1/Cf           0               
      0       1/Cf   -w0    0       0               -1/Cf           
      0       0      1/Lc   0       -Rc/Lc-rn/Lc    w0              
      0       0      0      1/Lc    -w0             -Rc/Lc-rn/Lc];%##   %Voltages and Currents of LCL states

BLCL1=[1/Lf 0
        0   1/Lf
        0   0
        0   0
        0   0
        0   0]; %Voltages references of VSC %##
    
BLCL2=[0            0
       0            0
       0            0
       0            0
       -rn/Lc       0
       0            -rn/Lc];%current source of the MG%##
    
BLCL3_1=[Ilq1
        -Ild1
        Voq1
        -Vod1
        Ioq1 
        -Iod1];%## %The coupling terms  (NOT common)

BLCL4_1=[0 1/Lf*Vid1/Vdc
         0 1/Lf*Viq1/Vdc
         0 0
         0 0
         0 0
         0 0]; %##%Voltages of VSC 

BLCL5_1=[0                           0
         0                           0
         0                           0
         0                           0
         rn/Lc                       0
         0                           rn/Lc];%##

ALoad=[-rn/L_load-R_load/L_load    w0  
       -w0                         -rn/L_load-R_load/L_load];%##
Bload1=[0 0 0 0 rn/L_load    0
        0 0 0 0 0            rn/L_load];%##%Output current to load current
Bload2=[Iloadq
        -Iloadd];%##
Bload3=[rn/L_load    0
        0            rn/L_load];%##%Source current to load current
%Inverter states [P Q phid phiq gamd gamq ild ilq vod voq iod ioq igd igq]
Ainv_1=[Ap                                                                       zeros(2,2)           zeros(2,2)           Bp_1                                                     zeros(2,2)
        Bv1*Cpv_1                                                                zeros(2,2)           zeros(2,2)           Bv2+Bv1*Dpv_1                                            zeros(2,2)
        Bc1*Dv1*Cpv_1                                                            Bc1*Cv               zeros(2,2)           Bc1*Dv2+Bc2+Bc1*Dv1*Dpv_1                                zeros(2,2)
        BLCL1*Dc1*Dv1*Cpv_1+BLCL3_1*Cpw_1                                        BLCL1*Dc1*Cv         BLCL1*Cc             ALCL+BLCL1*(Dc1*Dv2+Dc2+Dc1*Dv1*Dpv_1)+BLCL3_1*Dpw_1     BLCL5_1
        Bload2*Cpw_1                                                             zeros(2,2)           zeros(2,2)           Bload1+Bload2*Dpw_1                                      ALoad];
Eig_inv=eig(Ainv_1);

Binv_is=[zeros(6,2)
        BLCL2
        Bload3];%##

Binv_ed=[zeros(2,2)
         zeros(2,2)
         zeros(2,2)
         BLCL4_1
         zeros(2,2)];%##
p_il=[Ild1 Ilq1];
p_ed=[Vid1/Vdc;Viq1/Vdc];
p_cc=[Vid1 Viq1 0 0 0 0];
C_AC_p=[p_il*Dc1*Dv1*Cpv_1         p_il*Dc1*Cv         p_il*Cc       p_il*(Dc1*Dv2+Dc2+Dc1*Dv1*Dpv_1)+p_cc     zeros(1,2)];
D_AC_p=p_il*p_ed;


A_tot=[A_ele                                    B_ele*C_DC_ec                                           zeros(size(A_ele,2),size(Ainv_1,2))
       B_DC_cl*C_icl                            A_DC+B_DC_cl*D_icl*C_DC_ec+B_DC_p*D_AC_p*C_DC_ed        B_DC_p*C_AC_p
       zeros(size(Ainv_1,2),size(A_ele,2))      Binv_ed*C_DC                                            Ainv_1];

A_tot=A_tot([1:26,28:end],[1:26,28:end]);
Eig_tot=eig(A_tot);

B_tot=[zeros(26,2)
       zeros(5,2)
       Binv_is];
C_tot=eye(45);
D_tot=zeros(size(A_tot,2),size(B_tot,2));

figure('color','white')
xlabel("Real")
ylabel("Imaginary")
grid on;
%box on
hold on
plot(Eig_tot,'k<','LineWidth',1.0)
xlim([-3000 500])
ylim([-3200 3200])
% legend('All DGs within limit','DG3 reached the limit','DG2 & DG3 reached the limit')
% axes('position', [0.63 0.65 0.25 0.25])
% box on
% hold on 
% grid minor
% plot(Eig_tot,'ko','LineWidth',1.0)
% xlim([-0.6 0.1])
% ylim([-2.2 2.2])

[V, D_eig] = eig(A_tot,"vector");  % Eigenvectors (V) and Eigenvalues (D_eig)
W = inv(V)';          % Left eigenvectors (W = inv(V)')

% Compute Participation Factors (P_ij = V_ij * W_ji)
jj=1;
% j is index on columns (modes)
% i is index on rows (states)
while jj<46
 ii=1;
 while ii<46
 P(ii,jj)=abs(real(W(ii,jj)*V(ii,jj)));
 ii=ii+1;
 end
 jj=jj+1;
end

% Rank states based on dominant participation factors
[neworder, dom_states] = sort(sum(P, 2), 'descend');

r = 38;

selected_states = dom_states(1:r);
n=size(A_tot,2);

T=eye(n);
T=T(dom_states,:);


A_new=T*A_tot*T';
B_new=T*B_tot;
C_new=C_tot*T';
D_new=D_tot;

A_11=A_new(1:r,1:r);
A_12=A_new(1:r,r+1:end);
A_21=A_new(r+1:end,1:r);
A_22=A_new(r+1:end,r+1:end);

B_1=B_new(1:r,:);
B_2=B_new(r+1:end,:);

C_11=C_new(1:r,1:r);
C_12=C_new(1:r,r+1:end);
C_21=C_new(r+1:end,1:r);
C_22=C_new(r+1:end,r+1:end);


D_1=D_new(1:r,:);
D_2=D_new(r+1:end,:);

A_22_inv=pinv(A_22);
As=A_11-A_12*A_22_inv*A_21;
eigenred=eig(As);
% eig(blkdiag(As,A_22));

Bs=B_1-A_12*A_22_inv*B_2;
Css=C_11-C_12*A_22_inv*A_21;
Cff=C_21-C_22*A_22_inv*A_21;
Ds=D_1-C_12*A_22_inv*B_2;
Df=D_2-C_22*A_22_inv*B_2;
Lin_op_red=Lin_op;
Lin_op_red=Lin_op_red(1:r);

% 



% axesHandles = findall(gcf, 'Type', 'axes');
% xlim([29 35]);
% ylim([80 83.5]);
% ylim([190 215]);
% ylim([622 623]);
% ylim([-18000 -16000]);
% ylim([5500 5650]);
% ylim([370 380]);
% ylim([-50 -40]);
% ylim([-15.1 -14.95]);
% ylim([25.8 26.0]);
% ylim([-15.2 -15.05]);
% ylim([0.224 0.23]);
% ylim([1.9e-4 2.2e-4]);
% ylim([333.15 333.25]);
% ylim([93.468 93.469]);
% ylim([-7e-3 7e-3]);
% ylim([0.015 0.025]);
