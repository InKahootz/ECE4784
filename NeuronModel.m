%% Neuron Excitation - Hudgkin Huxley Model

% Simulation Time 
duration = 100; % millisec
stepsize = 0.01;
totalsteps = duration/stepsize; 

t = linspace(0,duration,totalsteps);%% Constants
gb_K = 36;  % mS/cm^2
gb_Na = 120; % ^^
gb_L = 0.3;  % ^^
E_K = -12;  % mV
E_Na = 115; % ^^
E_L = 10.6; % ^^
C_m = 1; %Membrane cap
I(1:totalsteps) = 10; %Membrane Current, microAmps/cm^2

% You can specify other currents here.
% A .5 ms pulse of 5 uA/cm^2 would be I(101:150) = 5; At 1ms, it would
% pulse for .5ms given the simulation is 100ms long with a step of .01 ms


%% Equations
V_m = 0;
% Gating vars
a_m = 0.1 * ((25 - V_m)/(exp((25 - V_m)/10) - 1));
b_m = 4 * exp(-V_m/18);
a_n = 0.01 * ((10 - V_m)/(exp((10 - V_m)/10) -1));
b_n = 0.125 * exp(-V_m/80);
a_h = 0.07 * exp(-V_m/20);
b_h = 1 / (exp((30 - V_m)/10) + 1);

m = a_m / (a_m + b_m);
n = a_n / (a_n + b_n);
h = a_h / (a_h + b_h);

% Currents
I_Na = m.^3 * gb_Na * h * (V_m - E_Na);
I_K = n.^ 4 * gb_K * (V_m - E_K);
I_L = gb_L * (V_m - E_L);
I_ion = I(1) - I_K - I_Na - I_L;

% Now that all variables and equations are setup we can simulate
% All of the previous var will change in each given step so we need to
% recalculate each of them using Euler's Method
% y(n+1) = y(n) + stepsize * f(t(n),y(n))
for i = 1:(totalsteps-1) %steps - 1 b/c we dont start at 0
    % Gating vars
    a_m(i) = 0.1 * ((25 - V_m(i))/(exp((25 - V_m(i))/10) - 1));
    b_m(i) = 4 * exp(-V_m(i)/18);
    a_n(i) = 0.01 * ((10 - V_m(i))/(exp((10 - V_m(i))/10) -1));
    b_n(i) = 0.125 * exp(-V_m(i)/80);
    a_h(i) = 0.07 * exp(-V_m(i)/20);
    b_h(i) = 1 / (exp((30 - V_m(i))/10) + 1);

    % Currents
    I_Na = m(i).^3 * gb_Na * h(i) * (V_m(i) - E_Na);
    I_K = n(i).^ 4 * gb_K * (V_m(i) - E_K);
    I_L = gb_L * (V_m(i) - E_L);
    I_ion = I(i) - I_K - I_Na - I_L;
    
    % Euler Method
    V_m(i+1) = V_m(i) + stepsize * I_ion / C_m;
    m(i+1) = m(i) + stepsize * (a_m(i) * (1 - m(i)) - b_m(i) * m(i));
    n(i+1) = n(i) + stepsize * (a_n(i) * (1 - n(i)) - b_n(i) * n(i));
    h(i+1) = h(i) + stepsize * (a_h(i) * (1 - h(i)) - b_h(i) * h(i));
    
end


% Plot Voltage vs time
plot(t,V_m-70)  % Plot the Y values (V_m) at times t.
axis([0 100 -100 50]) % Set the viewing min, max on x and y axis
xlabel 'Time (ms)' % Label axis
ylabel 'Voltage (mV)'
title 'Membrane Potential' % Plot title
legend 'Voltage' % Legend 

% Plotting conductances
% The equations for conductance come from the paper
g_K = gb_K .* (n .^ 4);        % Equation 6
g_Na = (m .^ 3) .* h .* gb_Na; % Equation 14
figure
hold on
plot(t,g_K, t, g_Na)
xlabel 'Time (ms)'
ylabel('Conductance (mS / cm^2)')
title 'Conductance'
legend('g_K', 'g_{Na}') %Matlab interprets TeX here.
