%% Plot Phases for mutiple packets


% csi_trace = read_bf_file('~/csi-data/csi-20170803-400-70-45.dat');
% csi_trace = read_bf_file('../../../csi-data/csi-20170803-400-70-45.dat');
csi_trace = read_bf_file('../../../csi-data/csi-20170804-320-7-60-1.dat');
% csi_frames = csi_trace_30;

aoas = [];
idx_aoas = [];
countdown = 20;

x_sc = 1:57;

% figure(1), clf
% figure(2), clf
% figure(3), clf
figure(66),clf

for idx=3070:length(csi_trace)
    if csi_trace{idx}.Nrx == 3
        countdown = countdown - 1;
        csi_trace{idx}.csi = csi_trace{idx}.csi(1,:,:);
        e_csi = csi_extend_57(csi_trace{idx}.csi);

        csi_matrix = e_csi; % from 3 columns to 3 row
        amp = abs(csi_matrix);
        fi = angle(csi_matrix);
        fii = fi;

        % phase smooth
        for rx=1:csi_trace{idx}.Nrx
            offset = 0;
            for i=2:length(fi(rx,:))
                if abs(fi(rx,i) - fi(rx,i-1)) > pi
                    if fi(rx,i-1) > fi(rx,i)
                        offset = offset + 2 * pi;
                    else
                        offset = offset - 2 * pi;
                    end
                end
                fii(rx,i) = fi(rx,i) + offset;
            end
        end
        
%         % original phase plot
%         figure(11)
%         hold on
%         plot(x_sc,squeeze(fi(1,:)),'b-o') 
%         %plot(x_sc,squeeze(fi(2,:)),'r-o')
%         %plot(x_sc,squeeze(fi(3,:)),'y-o')
%         hold off
%         title('Original phase of CSI for each subcarrier')
%         xlabel('Subcarriers')
%         ylabel('Phase')
% 
         % linerlized phase plot
%         figure(1),
         subplot(2,2,1)
         plot(x_sc,squeeze(amp(1,:)),'b-*') 
         hold on
         plot(x_sc,squeeze(amp(2,:)),'r-*')
         plot(x_sc,squeeze(amp(3,:)),'g-*')
         hold off
         title('Smoothed amplitude of CSI for each subcarrier')
         xlabel('Subcarriers')
         ylabel('amplitude')
         
%         figure(2),clf
         subplot(2,2,2)
         plot(x_sc,squeeze(fii(1,:)),'b-o') 
         hold on
         plot(x_sc,squeeze(fii(2,:)),'r-o')
         plot(x_sc,squeeze(fii(3,:)),'g-o')
         hold off
         title('Smoothed phase of CSI for each subcarrier')
         xlabel('Subcarriers')
         ylabel('Phase')
         
         
%         figure(3),clf
        pdp = (abs(ifft(squeeze(csi_matrix(1,:))))+abs(ifft(squeeze(csi_matrix(2,:))))+abs(ifft(squeeze(csi_matrix(3,:)))))/3;
        pdp = pdp - 44 -  csi_trace{idx}.agc;
        rssi = (csi_trace{idx}.rssi_a+csi_trace{idx}.rssi_b+csi_trace{idx}.rssi_c)/3 - 44 - csi_trace{idx}.agc;
        pdp_peak = max(pdp);
        
        fprintf('rssi:%f  pdp_peak:%f\n',rssi,pdp_peak);
         subplot(2,2,3)
         plot(x_sc*0.05,pdp,'b-.')
%          plot(x_sc*0.05,abs(ifft(squeeze(csi_matrix(1,:)))),'b-.')
%          hold on
%          plot(x_sc*0.05,abs(ifft(squeeze(csi_matrix(2,:)))),'r-.')
%          plot(x_sc*0.05,abs(ifft(squeeze(csi_matrix(3,:)))),'g-.')
%          hold off
         title('PDP')
         xlabel('time (us)')
         ylabel('amplitude')
         
%         % %% Sanitizing ToF Estimates
%         a = polyfit([x_sc x_sc x_sc],[fii(1,:) fii(2,:) fii(3,:)],1);
%         figure(3)
%         hold on
%         plot(x_sc,fii(1,:) - polyval(a,x_sc,1),'b-o')
%         % plot(x_sc,fii(2,:)' - polyval(a,x_sc,1),'r-o')
%         % plot(x_sc,fii(3,:)' - polyval(a,x_sc,1),'g-o')
%         hold off
%         title('filtered CSI, 30 sub-carriers')
%                 
        drawnow;
        pause(0.5);
    end
    
    if countdown == 0
        break;
    end
end

