%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Initialization %%%%

Rate = 4 / 3;
omega = 2 * pi * Rate;

%Diaphragm
Md = 4;    %%Mass of the jug? 
g = 9.81;  %Gravitational Constant
h = 0.01;   %
Ad = 0.06; %in m2
Aj = 0.03;
deltaX = (Aj * h) / Ad;
Plung = (Md * g) / Aj;

Kd = (Ad * Plung) / deltaX

%Resistances:
Rvo = 2700; % Vital organs: Heart & Brain (mmHg/(L/sec))
Rp = 180; % Pulmonary arteries, capillaries, veins (mmHg/(L/sec))
Rsa = 60; % Small in line restance of aorta (mmHg/(L/sec))
Rsv = 60;  % Small in line resistance of vena cava (mmHg/(L/sec))
Rl = 5400; % Legs (mmHg/(L/sec))
Rair = 0.014;

%Compliances:
Clung = (1 / 5) * 1.4;

%Pressures:
deltaPmax = 100;


%Figures 
f1 = figure('Name', 'Abdominal Aortic Pressure');
f2 = figure('Name', 'Vena Cava Pressure');
f3 = figure('Name', 'Thoracic Aorta Pressure');
f4 = figure('Name', 'Right Heart Pressure');
f5 = figure('Name', 'Central Perfusion Pressure');
f6 = figure('Name', 'Exhaled Volume (L) vs time during OAC-CPR');


%Equations:
dPaa_dt = 0;
dPivc_dt = 0;
dPao_dt = 0;
dPrh_dt = 0;
dDPext_dt = 0;
dVol = 0; %Volume (13)
dPlung = 0;
flowin = 0;
flowout = 0;
flowair = 0;
Pmax = 100;
Vol = 0;

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
dV_store = zeros(1, 1001);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = 1;
%%%% Update Pressure:
for time = 0:deltaT:endTime
    %index = index + 1;
    %plot([time, Paa], [time, Pivc], [time, Pao], [time, Prh]);
    %hold on;
    
    %11, 12, 13
	dDPext_dt = (Pmax / 2) * omega * sin(omega * time); %change in external pressure  %Probably want to replace 50 with deltaPmax / 2
    dVol = dDPext_dt * (power(Ad, 2) / Kd);
    flowin = (Pmax * power(Ad, 2) * omega * sin(omega * time)) / (2 * Kd);
    flowout = Plung / Rair;
    %dPlung = ((Pmax * power(Ad, 2) * omega * sin(omega * time) / (Clung * 2 * Kd)) - dPlung / Rair) * deltaT; %(12)
    %dPlung = (1 / Clung) * (((power(Ad, 2) * dDPext_dt) / Kd) - (Plung / Rair));
    dPlung = (1 / Clung) * (flowin - flowout);
    flowair = Plung / Rair;
    dVlung = dPlung * deltaT;
    %Vol = flowair * time;
    
    %part 12
    %Lung Prssure
    %dPlung / dt = (Ad ^2) * 
	Plung = Plung + (deltaT * dPlung); %%Pressure
    
    %part 13
    %Exaled Volume Volume
    %%Not sure about this part
     Vol = Vol + dVlung;
    
     if mod(time,.01) == 0
         dV_store(index) = Vol;
         index = index + 1; 
     end
    
    
   
    
     
end


%Determinations of tidal volume. 

%CppFunct = sineFit(time_place, pre_placeholder_Cpp);
%CppModel = double(CppFunct(2)) * (sin(2 * pi * double(CppFunct(3)) * x + double(CppFunct(4)))) + double(CppFunct(1));
%areaCpp = double((int(CppModel, 3, 3.75)) / 0.75);


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

figure(f6);
plot(time_place, dV_store);
xlabel('Time');
ylabel('Volume of diaphragm');
title('Diaphragm volume vs Time');
hold on
