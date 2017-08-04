
function csi57 = csi_extend_57(csi)

    % csi_frame = csi_good{1167};

    % narginchk(1,1)
    csi30 = squeeze(csi);
    amp30 = abs(csi30);
    fi30 = angle(csi30);
    fii30 = fi30;

    %%%%
    % original plot
    % plot(1:30,squeeze(fii(:,1)),'b-o') 
    % plot(1:30,squeeze(fii(:,2)),'r-o')
    % plot(1:30,squeeze(fii(:,3)),'y-o')
    %%%%

    % phase smooth
    for rx=1:3
        offset = 0;
        for i=2:length(fi30(rx,:))
            if abs(fi30(rx,i) - fi30(rx,i-1)) > pi-1
                if fi30(rx,i-1) > fi30(rx,i)
                    offset = offset + 2 * pi;
                else
                    offset = offset - 2 * pi;
                end
            end
            fii30(rx,i) = fi30(rx,i) + offset;
        end
    end

    x_sc = 1:30;
%     figure(11), clf
%     hold on
%     plot(x_sc,squeeze(fi30(1,:)),'b-o') 
%     plot(x_sc,squeeze(fi30(2,:)),'r-o')
%     plot(x_sc,squeeze(fi30(3,:)),'y-o')
%     hold off
%     title('Phase of original CSI, 30 sub-carriers')
%     drawnow;
%     pause(1)

%     figure(2), clf
%     hold on
%     plot(x_sc,squeeze(fii30(1,:)),'b-o') 
%     plot(x_sc,squeeze(fii30(2,:)),'r-o')
%     plot(x_sc,squeeze(fii30(3,:)),'y-o')
%     hold off
%     title('Phase of smoothed CSI, 30 sub-carriers')

    % %amp30 = amp30 - 44 - csi_frame.agc; % check get_total_rss for more information
    % figure(2)
    % plot(squeeze(amp30(:,1)),'b-*')
    % hold on
    % plot(squeeze(amp30(:,2)),'r-*')
    % plot(squeeze(amp30(:,3)),'y-*')
    % hold off
    % title('Amplitude of CSI, 30 sub-carriers')

    scidx20M = [-28:2:-2 -1:2:27 28]+29;
    c_scidx20M = [-27:2:-3 0 2:2:26]+29; % couterpart of csidx20M
    amp57 = zeros(3,57);
    fii57 = zeros(3,57);

    amp57(:,scidx20M) = amp30;
    fii57(:,scidx20M) = fii30;

    fii57(:,c_scidx20M) = (fii57(:,c_scidx20M-1)+fii57(:,c_scidx20M+1))/2;
    amp57(:,c_scidx20M) = 10*log10((10.^(amp57(:,c_scidx20M-1)./10)+10.^(amp57(:,c_scidx20M+1)./10))./2);

    [csi_matrix_real,csi_matrix_imag] = pol2cart(fii57,amp57);
    csi57 = csi_matrix_real+1j*csi_matrix_imag;

    % 
    % x_sc57 = 1:57;
    % 
    % figure(3)
    % plot(x_sc57,squeeze(fii57(:,1)),'b-o') 
    % hold on
    % plot(x_sc57,squeeze(fii57(:,2)),'r-o')
    % plot(x_sc57,squeeze(fii57(:,3)),'y-o')
    % hold off
    % title('Phase of CSI, extended to 57 sub-carriers')
    % 
    % figure(4)
    % plot(squeeze(amp57(:,1)),'b-*')
    % hold on
    % plot(squeeze(amp57(:,2)),'r-*')
    % plot(squeeze(amp57(:,3)),'y-*')
    % hold off
    % title('Amplitude of CSI, extended to 57 sub-carriers')
    % 
    % a = [polyfit(x_sc,fii30(1,:)',1);polyfit(x_sc,fii30(2,:)',1);polyfit(x_sc,fii30(3,:)',1)];
    % figure(5)
    % plot(x_sc,fii30(1,:)'-polyval(a(1,:),x_sc,1))
    % hold on
    % plot(x_sc,fii30(2,:)'-polyval(a(2,:),x_sc,1))
    % plot(x_sc,fii30(3,:)'-polyval(a(3,:),x_sc,1))
    % hold off
    % title('filtered CSI, 30 sub-carriers')
    % 
    % b = [polyfit(x_sc57,fii57(:,1)',1);polyfit(x_sc57,fii57(:,2)',1);polyfit(x_sc57,fii57(:,3)',1)];
    % figure(6),
    % plot(x_sc57,fii57(:,1)'-polyval(b(1,:),x_sc57,1))
    % hold on
    % plot(x_sc57,fii57(:,2)'-polyval(b(2,:),x_sc57,1))
    % plot(x_sc57,fii57(:,3)'-polyval(b(3,:),x_sc57,1))
    % hold off
    % title('filtered CSI, extended to 57 sub-carriers')

    % title('Phase of CSI for each subcarrier')
    % xlabel('Subcarriers')
    % ylabel('Phase')

    % 

end