%% CSI simulation

% clear
% close all

%% Basic simulation for 3 antennas

% Number of Paths D=8
D = 5; 
% Set incident angle and weight for each path 

num_of_pkt = 50;
csi_simulated_pkt = zeros(num_of_pkt,3,57);
csi_trace_30 = cell(num_of_pkt,1);

true_aoa=28;
true_tof=8e-9;

for pkt=1:num_of_pkt
    
incident_angle = randi([-90 90],[1 D]);
incident_attenuation = 1*(randn([1 D])+1j*(randn([1 D])));
tof = randi([10 300],[D 1]) * 1e-9;
% incident_attenuation = [10 20 30 40 50 60 70 80];
% set 45 degree as the major direction
% set time of flight, 1ns~0.3m, 10ns~3m, 100ns~30m
incident_angle(3) = true_aoa;
incident_attenuation(3) = incident_attenuation(3)*1.5;
tof(3) = true_tof;

% % Number of Paths D=2
% D = 2; 
% M = 3;
% % Set incident angle and weight for each path 
% incident_angle = [17 20];
% incident_attenuation = (randi([2 6],[1 2])+1j*(randi([2 4],[1 2])));
% % set 45 degree as the major direction
% % incident_attenuation(2) = 30 + 30j;
% % set time of flight, 1ns~0.3m, 10ns~3m, 100ns~30m
% tof = randi([8 40],[1 2]) * 1e-9;
% % tof(3) = 20e-9;

%% Fomulate steering matrix and incident quanties, given by the number of multipath
c = 3e8; % speed of light
% f = 2.412e9; % central frequency
f = 2.437e9; % central frequency
fs = 312.5e3; % 312.5 kHz
d = 0.07; % the minimal distance in between is 0.012 
twopi = 2*pi;
deg2rad = pi/180;
phiD = exp(-1j*twopi*d*sin(incident_angle*deg2rad)*f/c);
omega = exp(-1j*twopi*fs*tof);

A = [ones(1,D)
    phiD 
    phiD.^2];

F = zeros(D,57);
FF = omega.^(0:56);

for row=1:D
   F(row,:) = incident_attenuation(row).*FF(row,:);    
end

%% Simulate CSI measurements

X1 = A*F;
snr = 40;
% X  = awgn(X1,snr,'measured');
csi_simulated = X1;

scidx20M = [-28:2:-2 -1:2:27 28]+29;
% c_scidx20M = [-27:2:-3 0 2:2:26]+29; % couterpart of csidx20M
% 
csi_simulated_30 = csi_simulated(:,scidx20M);

csi_cell.timestamp_low = 4;
csi_cell.bfee_count = 1;
csi_cell.Nrx = 3;
csi_cell.Ntx = 1;
csi_cell.rssi_a = 39;
csi_cell.rssi_b = 36;
csi_cell.rssi_c = 33;
csi_cell.noise = -78;
csi_cell.agc = 24;
csi_cell.perm = [1 2 3];
csi_cell.rate = 8454;
csi_cell.csi = zeros(1,3,30);
csi_cell.csi(1,:,:) = csi_simulated_30;

csi_trace_30(pkt) = {csi_cell};
end

save mat_csi_30_simulated.mat csi_trace_30 true_aoa true_tof
