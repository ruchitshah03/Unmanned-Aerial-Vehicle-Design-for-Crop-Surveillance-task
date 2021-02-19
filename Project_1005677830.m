clear; clc;

%%
length = 500;
width = 200;
in2m = 0.0254;
rho = 1.225;            % considering at sea level
V = 15;                 
S = 0.4;                % assumption
Wbat = 0.232*9.81;      % assumption
e = 0.65;
AR = 7;
Clmax = 1.2;
n_max = 1.5;
Cd0a = 0.0072;
Cd0b = 0.0062;
Wea = 11.7672;
Web = 8.4037;
We = Wea + Web*S;
Wcam = 0.5*9.81;
Cd0 = Cd0a + Cd0b/S;
W = We + Wbat + Wcam;
q = 0.5*rho*(V^2);
Cl = W/(q*S);
if Cl > Clmax
    fprintf('Cl is greater than Clmax. Check W or S again');
end
g = 9.81;       

%% Path Comparison

TE_comp = zeros(4,2);
for i = 1:4
    path = strcat('path',num2str(i),'.m');
    run(path);
    
    % Drag
    Cd = Cd0 + ((nb^2)*(W^2))/(pi*e*AR*(q^2)*(S^2));
    
    % Thrust Required
    Tr_st = q*S*Cd0 + (W^2)/(pi*e*AR*q*S);
    Tr_tb = q*S*Cd0 + (((nb^2)*(W^2))/(pi*e*AR*q*S));
    Tr_ts = q*S*Cd0 + (((ns^2)*(W^2))/(pi*e*AR*q*S));
    
    % minimum thrust required
    Tr = max([Tr_st Tr_tb Tr_ts]);
    TE_comp(i,1) = Tr;
    
    % Energy Required
    E1 = Tr_st*ls;
    E2 = Tr_tb*arc_b;
    E3 = Tr_ts*arc_s;
    Ereq = E1+E2+E3;
    TE_comp(i,2) = Ereq;                % Energy required
    TE_comp(i,3) = (dist_tot/V)/60;     % time required
    TE_comp(i,4) = TE_comp(i,1)*V;      % Power required
    
end

%% Path Selection
% Path is selected based on the Energy Required. The path with lowest energy is selected for the design. 
[Emin, i] = min(TE_comp(:,2));
E = Emin;
Treq = TE_comp(i,1);
Preq = Treq*V;
min_dist = TE_comp(i,3)*60*V;
T_endu = TE_comp(i,3);

%% Propeller Selection
prop = ["apcsf_8x3.8"       % APC 8x3.8 SF
    "apcsf_8x6"             % APC 8x6 SF
    "apcsf_9x6"             % APC 9x6 SF
    "apcsf_9x7.5"           % APC 9x7.5 SF
    "apcsf_10x4.7"];        % APC 10x4.7 SF
propd = [8 8 9 9 10];
prop_dat = [];
for i = 1:5
    prop_nm = prop(i);
    prop_d = propd(i);
    count = 1;
    for RPM = 3000:1000:7000
        dat = read_prop_data(prop_nm,RPM,'G:\UofT\4_Fall_2020\AER1216\Project\UIUC-propDB\UIUC-propDB\volume-1\data');
        n = RPM/60;         % unit rev/sec
        D = prop_d*in2m;
        J = V/(n*D);
        const = Treq/(rho*(D^2)*(V^2));
        Ct1 = const*(J^2);
        Jmat = dat(:,1);
        Ctmat = dat(:,2);
        Cpmat = dat(:,3);
        etamat = dat(:,4);
        Ct2 = interp1(Jmat,Ctmat,J,'linear');
        Cp2 = interp1(Jmat,Cpmat,J,'linear');
        Cq2 = Cp2/(2*pi);
        eta2 = (Ct2*J)/Cp2;
        Q2 = Cq2*rho*(n^2)*(D^5);
        % plot(Jmat,Ctmat);
        prop_dat(count,1) = RPM/1000;
        prop_dat(count,2) = J;
        prop_dat(count,3) = Ct1;
        prop_dat(count,4) = Ct2;
        prop_dat(count,5) = Cp2;
        prop_dat(count,6) = eta2;
        prop_dat(count,7) = Cq2;
        prop_dat(count,8) = Q2;
        count = count + 1;
    end
    propeller_data(:,:,i) = prop_dat;
    % propeller and rpm with very close Ct1 and Ct2 are selected for required thrust (Ct2>=Ct1) 
end

%%
% Manually select the propeller from the matrix of propeller_data and then proceed further.
% Selecting Prop APC 9x7.5 SF with RPM 6000
prp_sel = "apcsf_9x7.5";
prp_dia = 9;            % inch
D = prp_dia*in2m;       % metres
RPM = 6000;
[datf, rpm] = read_prop_data(prp_sel,RPM,'G:\UofT\4_Fall_2020\AER1216\Project\UIUC-propDB\UIUC-propDB\volume-1\data');
n = RPM/60;
Jdata = datf(:,1);
Ctdata = datf(:,2);
Cpdata = datf(:,3);
etadata = datf(:,4);
J = V/(n*D);
Ct = interp1(Jdata,Ctdata,J,'linear');
Cp = interp1(Jdata,Cpdata,J,'linear');
Cq = Cp/(2*pi);
Q = Cq*rho*(n^2)*(D^5);
P = Cp*rho*(n^3)*(D^5);
T = Ct*rho*(n^2)*(D^4);
eta = (Ct*J)/Cp;

%% Motor Selection
% AXI 2217/12: Kv=1380RPM/V, rm=0.061, i0=0.7A, imax=32A
% AXI 2217/20: Kv=840RPM/V, rm=0.185, i0=0.4A, imax=20A

% ESC
re = 0.05;
Q = Q;
omega = 2*pi*n;
mKv = [1380 840];
mrm = [0.061 0.185];
mi0 = [0.7 0.4];
mimax = [32 20];
for i = 1:2
    Kv = mKv(i);              % RPM/V
    Kv = Kv*(2*pi)/60;      % rad/s/V
    Kt = 1/Kv;
    rm = mrm(i);
    i0 = mi0(i);
    imax = mimax(i);
    im(i) = (Q/Kt) + i0;
    if im(i)>imax
        fprintf('im is greater than imax. Change motor or prop');
    end
    Vmi(i) = (omega/Kv) + im(i)*rm;
    V_e0(i) = Vmi(i) + im(i)*re;
    i_e0(i) = im(i);
    P_esc(i) = V_e0(i)*i_e0(i);
    P_mi(i) = Vmi(i)*im(i);                      % Power across motor input
    eff_m_p(i) = (P/P_mi(i))*100;          
end
[eff_max, motor_number] = max(eff_m_p);

%% Battery
% Turnigy 30C Lipo Battery (2 or 3 cell)
% lets take 3 cell 1300 mAh 30C battery

E_eta = 0.72;                        % 80% draining battery and 90% efficient = 0.8*0.9 = 0.72
Eb_reqT = T_endu*E_eta*P_esc(motor_number)*60;
Eb_reqR = (min_dist*T)/E_eta;
Eb_req = max(Eb_reqT, Eb_reqR);
DCR = 30;

cell = [2 3];

for i = 1:2
    nc = cell(i);
    v_nom = 3.7;
    vb(i) = v_nom*nc;
    ke(i) = V_e0(1)/vb(i);
    if ke(i) < 1
        cell_f = cell(i);
        ib(i) = i_e0(1)*ke(i);
        Cap(i) = Eb_req*1000/(3600*vb(i));
        ib_max(i) = DCR*Cap(i)/1000;
    end
end

% So capacity of the battery required is 
% Thus, battery selected is 30C 2650 mAh 3 cell Lipo battery 
% weight is 0.232 kg. 
