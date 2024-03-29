%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part 1 Blood Flow Mathematical Model. 
%Compliance = change in volume / change in pressure
%Flow = Pressure/ Resistance
%%%% Initialization %%%%

Rate = 4 / 3;
omega = 2 * pi * Rate;

%Diaphragm
Md = 4;    %%Mass of the jug? 
g = 9.81;  %Gravitational Constant
h = 0.01;   %
Ad = 0.06; %in m2 %Area of Diaphram
Aj = 0.03; %Area of Jug?
deltaX = (Aj * h) / Ad;
Plung = (Md * g) / Aj; %Lung Pressure
Kd = (Ad * Plung) / deltaX;  %Spring Constant Equation

%Resistances:
Rvo = 2700; % Vital organs: Heart & Brain (mmHg/(L/sec))
Rp = 180; % Pulmonary arteries, capillaries, veins (mmHg/(L/sec))
Rsa = 60; % Small in line restance of aorta (mmHg/(L/sec))
Rsv = 60;  % Small in line resistance of vena cava (mmHg/(L/sec))
Rl = 5400; % Legs (mmHg/(L/sec))
Rair = 0.014;

%Compliances:
Cao = 0.00104167; % Thoracic aorta (L/mmHg)
Caa = 0.00052083; % Abdominal aorta (external pressure) (L/mmHg)
Civc =  30 * Caa; % I. vena cava (external pressure) (L/mmHg)
Crh = 30 * Cao; % Right heart: S.vena cava, RA, RV (L/mmHg)
Clung = (1 / 5) * 1.4

%Pressures:
Paa = 0; 
Pivc = 0;
Pao = 0;
Prh = 0;
deltaPmax = 100;
Cpp = Pao - Prh; %Central Perfusion Pressure or whatever it's called for abdominal CPR

%Figures 
f1 = figure('Name', 'Abdominal Aortic Pressure');
f2 = figure('Name', 'Vena Cava Pressure');
f3 = figure('Name', 'Thoracic Aorta Pressure');
f4 = figure('Name', 'Right Heart Pressure');
f5 = figure('Name', 'Central Perfusion Pressure');

%Equations:
dPaa_dt = 0;
dPivc_dt = 0;
dPao_dt = 0;
dPrh_dt = 0;
dDPext_dt = 0;
dV = 0;
dPlung = 0;
flowin = 0;
flowout = 0;
Pmax = 100;

%Time constraints:
time = 0;
deltaT = 0.000001;
endTime = 10;

time_place = 0:0.01:10;
pre_placeholder_Paa = zeros(1,1001);
pre_placeholder_Pao = zeros(1,1001);
pre_placeholder_Pivc = zeros(1,1001);
pre_placeholder_Prh = zeros(1,1001);
pre_placeholder_Cpp = zeros(1,1001);
%dV_store = zeroes(1, 1001);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = 1;
%%%% Update Pressure:
for time = 0:deltaT:endTime
    %index = index + 1;
    %plot([time, Paa], [time, Pivc], [time, Pao], [time, Prh]);
    %hold on;
    
    
	dDPext_dt = (deltaPmax / 2) * omega * sin(omega * time); %change in external pressure  %Probably want to replace 50 with deltaPmax / 2
    dV = dDPext_dt * (power(Ad, 2) / (power(Kd, 2))); 
    dPlung = (dV / Clung) * time;
    flowin = Pmax * power(Ad, 2) * ((omega * sin(omega * time)) / (power(Kd, 2)));
    flowout = dPlung / Rair;
    
    
    %%%Valve; We did not use the max function to model the action of heart valves but this code works
    Pao_seg = Prh - Pao;
    if(Pao_seg < 0)
        Pao_seg = 0;
    end
    
    Prh_seg = Prh - Pao;%checking for neg value for dPrh; 
    if(Prh_seg < 0)
        Prh_seg = 0;
    end
	  
    dPaa_dt = dDPext_dt + ( (1/Cao) * ( ((1/Rsa)*(Pao - Paa)) - ((1/Rl)*(Paa - Pivc)) ) );
    dPivc_dt = dDPext_dt + ( (1/Civc) * ( ((1/Rl)*(Paa - Pivc)) - ((1/Rsv)*(Pivc - Prh)) ) );
    dPao_dt = (1/Cao) * ( ((1/Rp) * Pao_seg) - ((1/Rsa) * (Pao - Paa)) - ((1/Rvo) * (Pao - Prh)) );
    dPrh_dt = (1/Crh) * ( ((1/Rsv) * (Pivc - Prh)) - ((1/Rp) * Prh_seg) - ((1/Rvo) * (Pao - Prh)) );
  
    Paa = Paa + (deltaT * dPaa_dt);
    Pivc = Pivc + (deltaT * dPivc_dt);
    Pao = Pao + (deltaT * dPao_dt);
    Prh = Prh + (deltaT * dPrh_dt);
    Cpp = Pao - Prh;
    
     if mod(time,.01) == 0
         pre_placeholder_Paa(index) = Paa;
         pre_placeholder_Pivc(index) = Pivc;
         pre_placeholder_Prh(index) = Prh;
         pre_placeholder_Pao(index) = Pao;
         pre_placeholder_Cpp(index) = Cpp;
         index = index + 1;
         
     end
end

CppFunct = sineFit(time_place, pre_placeholder_Cpp)
CppModel = double(CppFunct(2)) * (sin(2 * pi * double(CppFunct(3)) * x + double(CppFunct(4)))) + double(CppFunct(1));
areaCpp = double((int(CppModel, 3, 3.75)) / 0.75)


figure(f1);
%plot(time_place,pre_placeholder_Paa);
xlabel('Time');
ylabel('Pressure');
title('Abdominal Aortic Pressure vs Time');

figure(f2);
%plot(time_place, pre_placeholder_Pivc);
xlabel('Time');
ylabel('Pressure');
title('Vena Cava Pressure vs Time');

figure(f3);
%plot(time_place, pre_placeholder_Pao);
xlabel('Time');
ylabel('Pressure');
title('Thoracic Aorta Pressure vs Time');

figure(f4);
%plot(time_place, pre_placeholder_Prh);
xlabel('Time');
ylabel('Pressure');
title('Right Heart Pressure vs Time');

figure(f5);
%plot(time_place, pre_placeholder_Cpp);
xlabel('Time');
ylabel('Pressure');
title('Central Perfusion Pressure vs Time');
hold on

%plot(time_place, CppModel(time_place), 'b-');