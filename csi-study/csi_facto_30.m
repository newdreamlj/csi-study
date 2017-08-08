%% AoA test, run both in MATLAB and Octave

% clear
%load csi_good

% csi_good = read_bf_file('../sample_data/log.all_csi.6.7.6');
% csi_good = read_bf_file('../../../csi-data/csi2-0609-13.dat');
csi_good = read_bf_file('../../../csi-data/20170727-2.dat');
% aoas = zeros(1,length(csi_good));

last_valid_idx = 0;
aoas = [];
idx_aoas = [];
countdown = 2;

figure(10);
hold off
x_sc = 1:30;

for idx=1:length(csi_good)
    if csi_good{idx}.Nrx == 3
        countdown = countdown - 1;
        csi_good{idx}.csi = csi_good{idx}.csi(1,:,:);
        if last_valid_idx == 0
            last_valid_idx = idx;
        else
            % csi_good{idx}.csi = csi_good{idx}.csi * 0.8 + csi_good{last_valid_idx}.csi * 0.2;
            [aoa, rads, pmu] = csi_find_aoa_30(csi_good{idx});
            aoas(end+1) = aoa;
            polar(rads,pmu);
            drawnow;
            pause(0.1);
            idx_aoas(end+1) = idx;
            last_valid_idx = idx;
        end
        %fprintf('\nAoA=%d\n',aoa);
        if countdown == 0
            break;
        end
    end
end

% csi_find_aoa(csi_good{564})

% figure(11)
% plot(aoas);
% idx_aoas