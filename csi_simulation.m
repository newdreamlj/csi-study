%% CSI simulation

% clear
% close all

%% Basic simulation for 3 antennas

% Number of Paths D=8
D = 8; 
% Set incident angle and weight for each path 
incident_angle = [72 70 -40 1 3 5 7 9];

num_of_pkt = 1;
csi_simulated_pkt = zeros(num_of_pkt,3,57);

for pkt=1:num_of_pkt

incident_attenuation = (randn([1 D])+1j*(randn([1 D])));
% incident_attenuation = [10+0.2j 3-0.5j];
% incident_attenuation = [10 20 30 40 50 60 70 80];
% set 45 degree as the major direction
%incident_attenuation(3) = 5;
% set time of flight, 1ns~0.3m, 10ns~3m, 100ns~30m
% tof = randi([1 50],[D 1]) * 1e-9;
tof = [200 100 120 130 140 150 160 170]' .* 1e-9;
% tof(3) = 15e-9;

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
f = 2.412e9; % central frequency
fs = 312.5e3; % 312.5 kHz
d = 0.012; % the minimal distance in between is 0.012 
twopi = 2*pi;
deg2rad = pi/180;
phiD = exp(-1j*twopi*d*sin(incident_angle*deg2rad)*f/c);  % 1*D
omega = exp(-1j*twopi*fs*tof);   % D*1

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
snr = 25;
X  = awgn(X1,snr,'measured');
% X  = awgn(X1,snr);
%X = X1;
%% Traditional MUSIC algorithm
% 
% Rxx=X*X'/M;
% % InvS=inv(Rxx); %%%%
% [eVector,eValue]=eig(Rxx);%%%% 
% EVA=diag(eValue)';
% [EVA,I]=sort(EVA);
% eVector=fliplr(eVector(:,I));
% En=eVector(:,D+1:M);
% dd = (0:2)*d;
% % MUSIC
% for iang = 1:361
%         angle_x(iang)=(iang-181)/2;
%         phim=deg2rad*angle_x(iang);
%         a=exp(-1i*twopi*dd*sin(phim)*f/c).';
% %        SP(iang)=(a'*a)/(a'*En*En'*a);
%         SP(iang)=1/(a'*En*En'*a);
% end
% 
% % 
% SP=abs(SP);
% SPmax=max(SP);
% SP=10*log10(SP/SPmax);
% h=plot(angle_x,SP);
% set(h,'Linewidth',2)
% xlabel('angle (degree)')
% ylabel('magnitude (dB)')
% axis([-90 90 -60 0])
% set(gca, 'XTick',[-90:30:90])
% grid on  


%% Spoi-fi solution
csi_simulated_pkt(pkt,:,:) = X;
csi_simulated = X;
[tofs, rads, Pmu] = csi_find_aoa_spotfi(csi_cell,csi_simulated);

            figure(10);
            surf(tofs*1e9,rads*180/pi,Pmu)
            drawnow;
            pause(0.8);
            Pmu_mirror = [Pmu; flipud(Pmu)];
            maxima = find_maxima(Pmu_mirror);
            fprintf('idx %d: \n',idx);
            for k=1:size(maxima,1)
                if maxima(k,1)<=length(rads)
                    AoA = rads(maxima(k,1))*180/pi;
                    ToF = maxima(k,2);
                    fprintf('   AoA=%d  ToF=%d ns\n',ceil(AoA),ToF);
    %             else
    %                 fprintf('idx %d: miss\n',idx);
                end
            end


% scidx20M = [-28:2:-2 -1:2:27 28]+29;
% c_scidx20M = [-27:2:-3 0 2:2:26]+29; % couterpart of csidx20M
% 
% csi_simulated_30 = csi_simulated(:,scidx20M);

end


%[aoa, rads, pmu] = csi_find_aoa(csi_cell,squeeze(csi_simulated_pkt(1,:,:)));
% [tof, rads, Pmu] = csi_find_aoa_multipkt(csi_cell,csi_simulated_pkt,num_of_pkt);

%figure(10);
% figure(101),surf(tof*1e9,rads*180/pi,Pmu)
% figure(102),polar(rads,Pmu(:,15));
%figure(10);
%plot(aoas);
% idx_aoas


% M = 30;
% N = M - D;
% 
% derad = pi/180;        % deg -> rad
% radeg = 180/pi;
% twpi = 2*pi;
% kelm = 8;               % number
% dd = 0.5;               % space 
% d=0:dd:(kelm-1)*dd;     % 
% iwave = 3;              % number of DOA
% theta = [19 28 64];     % angle
% snr = 37;               % input SNR (dB)
% n = 8;                 % 
% A=exp(-1i*twpi*d'*sin(theta*derad));%%%% direction matrix
% S=randn(iwave,n);
% X=A*S;
% X1=awgn(X,snr,'measured');
% Rxx=X1*X1'/n;
% % InvS=inv(Rxx); %%%%
% [EV,D]=eig(Rxx);%%%% 
% EVA=diag(D)';
% [EVA,I]=sort(EVA);
% EVA=fliplr(EVA);
% EV=fliplr(EV(:,I));
% 
% % MUSIC
% for iang = 1:361
%         angle(iang)=(iang-181)/2;
%         phim=derad*angle(iang);
%         a=exp(-1i*twpi*d*sin(phim)).';
%         L=iwave;    
%         En=EV(:,L+1:kelm);
% %        SP(iang)=(a'*a)/(a'*En*En'*a);
%         SP(iang)=1/(a'*En*En'*a);
% end
% 
% % 
% SP=abs(SP);
% SPmax=max(SP);
% SP=10*log10(SP/SPmax);
% h=plot(angle,SP);
% set(h,'Linewidth',2)
% xlabel('angle (degree)')
% ylabel('magnitude (dB)')
% axis([-90 90 -60 0])
% set(gca, 'XTick',[-90:30:90])
% grid on  
