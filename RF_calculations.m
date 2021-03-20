%  ____  ____      _    __  __  ____ ___
% |  _ \|  _ \    / \  |  \/  |/ ___/ _ \
% | | | | |_) |  / _ \ | |\/| | |  | | | |
% | |_| |  _ <  / ___ \| |  | | |__| |_| |
% |____/|_| \_\/_/   \_\_|  |_|\____\___/
%                           research group
%                             dramco.be/
%
%  KU Leuven - Technology Campus Gent,
%  Gebroeders De Smetstraat 1,
%  B-9000 Gent, Belgium
%
%         File: RF_calculations
%      Created: 2020-03-21
%       Author: Chesney buyle
%      Version: 1.0
%
%  Description: calculation of RF WPT for Entrepreneurship
%
%  Commissiond by company C (optionally)
%
%  License L (optionally)
%

%% TRANSMITTER

% Transmit power regulations
% Link: https://afar.net/tutorials/fcc-rules/

P_transmitter = 27; % in dBm 
G_transmitter = 0; % in dBi 
G_receiver = 2; % in dBi
operating_frequency = 869.525e6; % in Hz
distance = 10; % in meters

% Converting variables to non-dB scale
P_transmitter = 10^(P_transmitter/10);
G_transmitter = 10^(G_transmitter/10);
G_receiver = 10^(G_receiver/10);

% Calculate additional variables
lambda = physconst('LightSpeed')/operating_frequency;

%% Propagation (based on Friis transmission equation)

% Ditances at which the total received power will be calculated
meas_dist = 0.1:0.1:10; % steps of 10cm

% Power received at antenna
P_delivered_to_antennas = zeros(1, numel(meas_dist));
P_delivered_to_antennas_dbm = zeros(1, numel(meas_dist));

for dist = 1:numel(meas_dist)
    
    P_received = (P_transmitter*G_transmitter*G_receiver*lambda^2)/...
    (1000*(4*pi*meas_dist(dist))^2); % in Watts
    
    % Power delivered to antenna = power received at receiver * (1 - S11
    % antenna)
    P_delivered_to_antennas(:, dist) = P_received;
    P_delivered_to_antennas_dbm(:, dist) = ...
        10.*log10(P_delivered_to_antennas(:, dist)*1000); %wrt to 1 mW - so multiply by 1000
end
close all;
figure
plot(meas_dist, P_delivered_to_antennas_dbm);
title('Received power at antenna of mobile device');
xlabel('Distance [m]');
ylabel('Received power [dBm]');

%% Energy harvester (and corresponding efficiencies)
full_eff = readmatrix('full_efficiency_868mhz.csv', 'DecimalSeparator', ',');

% Power available after the boost converter
P_out_boost = zeros(1, numel(meas_dist)); % in Watt
P_out_boost_dbm = zeros(1, numel(meas_dist)); % in dBm
check_efficiencies = zeros(1, numel(meas_dist));

for dist = 1:numel(meas_dist)
    p_in = P_delivered_to_antennas_dbm(1, dist);
    
    % Check if received power at antenna is within specifications
    if (p_in < 10 && p_in > -19.5)
        
        % Search efficiency (from datasheet) corresponding with received power
        [val, idx] = min(abs(full_eff(:, 1) - p_in));
        efficiency = full_eff(idx, 2)/100;
        
        % For debugging purposes (TODO: delete)
        check_efficiencies(1, dist) = efficiency;
           
        P_out_boost(1, dist) = P_delivered_to_antennas(1, dist) * efficiency; % in Watt
        P_out_boost_dbm(1, dist) = 10*log10(P_delivered_to_antennas(1, dist) * 1000);

    else 
        P_out_boost(1, dist) = nan;   
        P_out_boost_dbm(1, dist) = nan;
    end
end

figure
plot(meas_dist, P_out_boost*1e6); % in µW
title('Power available after boost converter')
xlabel('distance [m]');
ylabel('power [µW]')

% The boost converter delivers the power to the battery
% No information is available about the efficiency of charging 
% for example a battery
% The series resistance of a battery or for example supercap must taken
% into account to calculate this power delivery efficiency.
% For now, we suppose that this efficiency is 100%

efficiency_boost_to_energy_storage = 1;
P_in_battery = P_out_boost*efficiency_boost_to_energy_storage;
P_in_battery_dbm = 10*log10(P_out_boost*1000);

%% Consuming energy from the energy storage
% A LDO converts the battery voltage to a steady voltage (e.g. 3.3V)
% Of course, the LDO itself has a certain efficiency, which is given in the
% datasheet for both a low voltage LDO and high voltage LDO. 
% The efficiency depends on the battery voltage and the voltage it is 
% converted to (theoretically the efficiency is Vout/Vin). 
% If we assume the high voltage LDO and that Vout is 3.3V, 
% the lowest efficiency that we get is around 70% for a Li-ion battery.
efficiency_LDO = 0.7;

% Moreover, some energy is lost while discharing the battery itself but is
% highly dependent of the current consumption, battery voltage etc. So
% again for now, we suppose this efficiency is 100%
efficiency_discharge_battery = 1;

P_out_LDO = P_in_battery*efficiency_discharge_battery*efficiency_LDO;
P_out_LDO_dbm = 10*log10(P_out_LDO*1000);

figure
plot(meas_dist, P_out_LDO*1e6); % in µW
title('Power available after LDO')
xlabel('distance [m]');
ylabel('power [µW]')

%% Combine data and write to file
overall_efficiency_sim = P_out_LDO./P_delivered_to_antennas;

fid = fopen('power_calculations.txt', 'w');
fprintf(fid', 'Distance [m]; RF input power [µW]; P_out_LDO [µW]\n');

for dist = 1:numel(meas_dist)
    fprintf(fid, '%.4f; %.4f; %.4f\n', meas_dist(dist), ...
        P_delivered_to_antennas(dist)*1e6, P_out_LDO(dist)*1e6);
end

close(fid);

% Write data to file
% writematrix(data, 'power_calculations.csv');


