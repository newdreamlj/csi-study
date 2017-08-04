 %% AoA test, run both in MATLAB and Octave

% clear
%load csi_good

% csi_trace = read_bf_file('../sample_data/log.all_csi.6.7.6');
% csi_trace = read_bf_file('../../../csi-data/csi2-0609-13.dat');
% csi_trace = read_bf_file('../../../csi-data/csi-20170803-400-70--45.dat');
% csi_trace = read_bf_file('../../../csi-data/csi-20170804-320-7-60-1.dat');
%csi_trace = read_bf_file('../../../csi-data/csi-20170804-320-7-45-2.dat');
csi_trace = read_bf_file('../../../csi-data/csi-20170804-320-7-45-3.dat'); % 72000pkt in 45s
% aoas = zeros(1,length(csi_good))
% csi_trace = csi_trace_30;
% load mat_csi_30_simulated.mat
% csi_trace = csi_trace_30;

dataset = [];
countdown = 20;

for idx=1200:length(csi_trace)
    if csi_trace{idx}.Nrx == 3
        countdown = countdown - 1;
        for tx=1:1 % csi_trace{idx}.Ntx
            e_csi = csi_extend_57(csi_trace{idx}.csi(tx,:,:));
            % do estimation
            [tofs, rads, Pmu] = csi_find_aoa_spotfi(csi_trace{idx},e_csi);
%             figure(10);
%             surf(tofs*1e9,rads*180/pi,Pmu)
%     xlabel('ToF (ns)')
%     ylabel('AoA (degree)') 
%     zlabel('magnitude (dB)')
%     grid on  
%             shading interp;
%             drawnow;
%             pause(0.8);
            Pmu_mirror = [Pmu; flipud(Pmu)];
            maxima = find_maxima(Pmu_mirror);
            fprintf('idx %d: \n',idx);
            for k=1:size(maxima,1)
                if maxima(k,1)<=length(rads)
                    AoA = rads(maxima(k,1))*180/pi;
                    ToF = maxima(k,2);
                    fprintf('   AoA=%d  ToF=%d ns\n',ceil(AoA),ToF);
                    dataset(end+1,:) = [AoA,ToF,Pmu(maxima(k,1),maxima(k,2))];
    %             else
    %                 fprintf('idx %d: miss\n',idx);
                end
            end
        end
    end
    if countdown <= 0
        break;
    end
end

figure(17);
scatter(dataset(:,1),dataset(:,2));
xlabel('AoA (degree)') 
ylabel('ToF (ns)')
grid on  

figure(18);
statistic = tabulate(dataset(:,1));
bar(statistic(:,1),statistic(:,2));
xlabel('AoA (degree)') 
ylabel('')
grid on  



% csi_find_aoa(csi_good{564})
% aoas
% figure
% plot(aoas);
% idx_aoas