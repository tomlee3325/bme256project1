%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialization %%%%

Rate = 4 / 3; % 80/60 but simpified the constant abdominal compression rate in Hz
omega = 2 * pi * Rate; % angular frequency
deltaPmax = 170;

%Diaphragm
Md = 4;    %%Mass of the jug? 
g = 9.81;  %Gravitational Constant
h = 0.01;  % cm
Ad = 0.06; % area of abdomen
Aj = 0.03; % area of jug
deltaX = (Aj * h) / Ad;
P_lung = (Md * g) / Aj;
Kd = (Ad * P_lung) / deltaX;

%Resistances:
Rvo = 2700; % Vital organs: Heart & Brain (mmHg/(L/sec))
Rp = 180; % Pulmonary arteries, capillaries, veins (mmHg/(L/sec))
Rsa = 60; % Small in line restance of aorta (mmHg/(L/sec))
Rsv = 60;  % Small in line resistance of vena cava (mmHg/(L/sec))
Rl = 5400; % Legs (mmHg/(L/sec))
Rair = 1.02978; %%R airway%%

%Compliances:
Cao = 0.00104167; % Thoracic aorta (L/mmHg)
Caa = 0.00052083; % Abdominal aorta (external pressure) (L/mmHg)
Civc =  30 * Caa; % I. vena cava (external pressure) (L/mmHg)
Crh = 30 * Cao; % Right heart: S.vena cava, RA, RV (L/mmHg)
Clung = (1 / 5) * 1.4;

%Pressures initialize:
Paa = 0; %abdominal aorta
Pivc = 0; %inferior vena cava
Pao = 0; %thoracic aorta
Prh = 0; % Right heart
Cpp = Pao - Prh; %Central Perfusion Pressure or whatever it's called for abdominal CPR
meanCpp = 0; % coronary perfusion pressure

%Figures 
f1 = figure('Name', 'Abdominal Aortic Pressure');
f2 = figure('Name', 'InferiorVena Cava Pressure');
f3 = figure('Name', 'Thoracic Aorta Pressure');
f4 = figure('Name', 'Right Heart Pressure');
f5 = figure('Name', 'Coronary Perfusion Pressure');
f6 = figure('Name', 'Exhaled Volume (L) vs time during OAC-CPR');


%Equations Initial:
dPaa_dt = 0; % abdominal aorta
dPivc_dt = 0; % inferior vena cava
dPao_dt = 0; % thoracic aorta
dPrh_dt = 0; % right heart
n = 0;
dVol = 0; %Volume (13)
dPlung = 0;% pressure of lung
flowin = 0; %air flow in
flowout = 0; %air flow out
flowair = 0; %rate of airflow


%Time constraints:
deltaT = 0.000001; % increments by .000001 
endTime = 10; % 10 seconds
time_place = 0:0.01:10; % array for plotting graph

%Place holder to store Pressures
pre_placeholder_Paa = zeros(1,1001); % abdominal aorta 
pre_placeholder_Pao = zeros(1,1001); % thoracic aorta
pre_placeholder_Pivc = zeros(1,1001); % inferior vena cava
pre_placeholder_Prh = zeros(1,1001); % right heart
pre_placeholder_Cpp = zeros(1,1001); % coronary perfusion pressure
dV_store = zeros(1, 1001); % change of volume
Plung = zeros(1,1001); % Pressure of lung
Vol = zeros(1, 1001); % volume in lung
timecount = 0:deltaT:endTime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = 1;
%%%% Update Pressure:
for time = 0:deltaT:endTime
    
  %change in external pressure
	dDPext_dt = (deltaPmax/2) * omega * sin(omega * time); 
	  
   %calcualting the change in pressure
    dPaa_dt = dDPext_dt + ( (1/Caa) * ( ((1/Rsa)*(Pao - Paa)) - ((1/Rl)*(Paa - Pivc)) ) );
    dPivc_dt = dDPext_dt + ( (1/Civc) * ( ((1/Rl)*(Paa - Pivc)) - ((1/Rsv)*(Pivc - Prh)) ) );
    dPao_dt = (1/Cao) * ( max(0,((Prh-Pao)/Rp)) - ((Pao-Paa)/Rsa) - ((Pao-Prh)/Rvo) );
    dPrh_dt = (1/Crh) * (((1/Rsv) * (Pivc - Prh)) - max(0,((Prh-Pao)/Rp)) + ((1/Rvo) * (Pao-Prh)));
   
  %Updates the Pressure for each iteration
    Paa = Paa + (deltaT * dPaa_dt);
    Pivc = Pivc + (deltaT * dPivc_dt);
    Pao = Pao + (deltaT * dPao_dt);
    Prh = Prh + (deltaT * dPrh_dt);
    
    Cpp = Pao - Prh; % finding the coronary perfusion pressure 
    meanCpp = meanCpp + Cpp;
    n = n+1;
    
    %updating matrix for place
     if mod(time,.01) == 0
         pre_placeholder_Paa(index) = Paa;
         pre_placeholder_Pivc(index) = Pivc;
         pre_placeholder_Prh(index ) = Prh;
         pre_placeholder_Pao(index) = Pao;
         pre_placeholder_Cpp(index) = Cpp;
         index = index + 1;
         
     end
end
meanCpp = meanCpp/n; % calcualte true mean perfusion pressure

index_2 = 1;
for time = deltaT:deltaT:endTime

    
    %11, 12, 13
	dDPext_dt = (deltaPmax / 2) * omega * sin(omega * time); %change in external pressure  %Probably want to replace 50 with deltaPmax / 2
    dPlung = ((133* 1000* dDPext_dt * power(Ad, 2))  / Kd - (Plung(index_2) / Rair)) * (1 / Clung) ; %(12)
    index_2 = index_2 + 1;
    %part 12
	Plung(index_2) = Plung(index_2 - 1) + (deltaT * dPlung); %%Pressure
    Vol(index_2) = (Plung(index_2) / Rair) * deltaT + Vol(index_2 - 1); %Volume
    %part 13
    %Exaled Volume Volume
   
end

%%%% ploting the graphs 
figure(f1);
plot(time_place,pre_placeholder_Paa);
xlabel('Time (sec)');
ylabel('Pressure (mmHg)');
title('Abdominal Aortic Pressure vs Time');

figure(f2);
plot(time_place, pre_placeholder_Pivc);
xlabel('Time (sec)');
ylabel('Pressure mmHg');
title('Inferior Vena Cava Pressure vs Time');

figure(f3);
plot(time_place, pre_placeholder_Pao);
xlabel('Time (sec)');
ylabel('Pressure (mmHg)');
title('Thoracic Aorta Pressure vs Time');

figure(f4);
plot(time_place, pre_placeholder_Prh);
xlabel('Time (sec)');
ylabel('Pressure (mmHg)');
title('Right Heart Pressure vs Time');

figure(f5);
plot(time_place, pre_placeholder_Cpp);
xlabel('Time (sec)');
ylabel('Pressure (mmHg)');
title('Coronary Perfusion Pressure vs Time');

figure(f6);
plot(timecount, Vol);
xlabel('Time(sec)');
ylabel('Volume of diaphragm (L)');
title('Diaphragm[Lung] volume(L)during OAC-CPR');
hold on

