%% AoA test, run both in MATLAB and Octave

% clear
%load csi_good

% csi_good = read_bf_file('../sample_data/log.all_csi.6.7.6');
% csi_good = read_bf_file('../../../csi-data/csi2-0609-13.dat');
csi_trace = read_bf_file('../../../csi-data/csi-20170731-tp.dat');
% aoas = zeros(1,length(csi_good));
% csi_trace = csi_trace_30;

aoas = [];
idx_aoas = [];
countdown = 5;

num_of_pkt = 1;
cnt_pkt = 0;
csi_simulated_pkt = zeros(num_of_pkt,3,57);

for idx=1000:length(csi_trace)
    if csi_trace{idx}.Nrx == 3
        cnt_pkt = cnt_pkt + 1;
        countdown = countdown - 1;
        e_csi = csi_extend_57(csi_trace{idx});
        csi_simulated_pkt(cnt_pkt,:,:) = e_csi;
        if cnt_pkt == num_of_pkt
            cnt_pkt = 0;
            % do estimation
            [tofs, rads, Pmu] = csi_find_aoa_multipkt(csi_trace{idx},csi_simulated_pkt,num_of_pkt);
            figure(10);
            surf(tofs*1e9,rads*180/pi,Pmu)
            drawnow;
            pause(0.8);
            idx_aoas(end+1) = idx;
            last_valid_idx = idx;
            % fprintf('\nAoA=%d\n',aoa);
        end        
    end
    if countdown == 0
        break;
    end
end

% csi_find_aoa(csi_good{564})
% aoas
% figure
% plot(aoas);
% idx_aoas