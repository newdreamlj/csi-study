%% Implementation of Spot-Fi
% Basic plots of csi data
%%%%%%%%%%
% csi_trace = read_bf_file('sample_data/log.all_csi.6.7.6');
% csi_trace = read_bf_file('csi-0605-2.dat');
% load csi_trace
% csi_good = csi_trace;
%%%%%%%%%%

function [tofs,rads,Pmu] = csi_find_aoa_spotfi(csi_frame, e_csi)

    narginchk(2,2)

    x_sc = 1:57;

    csi_matrix = e_csi; % 
    amp = abs(csi_matrix);
    fi = angle(csi_matrix);
    fii = fi;

    % % % for rx=1:csi_frame.Nrx
    % % %     phase_comp=0;
    % % %     if fi(2,rx)>fi(1,rx) || fi(1,rx) - fi(2,rx) > pi
    % % %         factor = 1;
    % % %     else
    % % %         factor = -1;
    % % %     end
    % % % 
    % % %     for i=2:length(fi(:,rx))
    % % %         if (fi(i,rx) - fi(i-1,rx)) * factor < 0
    % % %             phase_comp = phase_comp + 2 * factor * pi;
    % % %         end
    % % %         fii(i,rx) = fi(i,rx) + phase_comp;
    % % %     end
    % % % end

    % phase smooth
    for rx=1:csi_frame.Nrx
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


    % % original phase plot
    % figure(1)
    % plot(x_sc,squeeze(fi(1,:)),'b-o') 
    % hold on
    % plot(x_sc,squeeze(fi(2,:)),'r-o')
    % plot(x_sc,squeeze(fi(3,:)),'y-o')
    % hold off
    % title('Original phase of CSI for each subcarrier')
    % xlabel('Subcarriers')
    % ylabel('Phase')
    
    % linerlized phase plot
%     figure(2)
%     plot(x_sc,squeeze(fii(1,:)),'b-o') 
%     hold on
%     plot(x_sc,squeeze(fii(2,:)),'r-o')
%     plot(x_sc,squeeze(fii(3,:)),'y-o')
%     hold off
%     title('Smoothed phase of CSI for each subcarrier')
%     xlabel('Subcarriers')
%     ylabel('Phase')
%     drawnow;
%     
    % % % plot amplitude
    % % figure(3), bar(squeeze(amp(:,:)))
    % % axis([0 31 -inf inf])
    % % title('Amplitude of CSI for each subcarrier')
    % % xlabel('Subcarriers')
    % % ylabel('dB? Need to make sure')
    % 
    % % plot amplitude in dB scale
    % amp = amp - 44 - csi_frame.agc; % check get_total_rss for more information
    % figure(4)
    % plot(squeeze(amp(1,:)),'-*')
    % hold on
    % plot(squeeze(amp(2,:)),'-*')
    % plot(squeeze(amp(3,:)),'-*')
    % hold off
    % axis([0 31 -inf inf])
    % title('Amplitude of CSI for each subcarrier')
    % xlabel('Subcarriers')
    % ylabel('dB Need to make sure')


    %% Sanitizing ToF Estimates

    % a = [polyfit(x_sc,fii(:,1)',1);polyfit(x_sc,fii(:,2)',1);polyfit(x_sc,fii(:,3)',1)];
%     a = polyfit([x_sc x_sc x_sc],[fii(1,:) fii(2,:) fii(3,:)],1);

    % figure(5)
    % plot(x_sc,fii(:,1)' - polyval(a(1,:),x_sc,1))
    % hold on
    % plot(x_sc,fii(:,2)' - polyval(a(2,:),x_sc,1))
    % plot(x_sc,fii(:,3)' - polyval(a(3,:),x_sc,1))
    % hold off
    % title('filtered CSI, 30 sub-carriers')

%         fi = fii - [polyval(a,x_sc,1); polyval(a,x_sc,1); polyval(a,x_sc,1)];
%         [csi_matrix_real,csi_matrix_imag] = pol2cart(fi,amp);
%         csi_matrix = csi_matrix_real + 1j*csi_matrix_imag;

%     figure(9)
%     plot(x_sc,squeeze(fii(1,:)),'b-o') 
%     hold on
%     plot(x_sc,squeeze(fii(2,:)),'r-o')
%     plot(x_sc,squeeze(fii(3,:)),'y-o')
%     plot(x_sc,polyval(a,x_sc,1))
%     hold off
%     title('Phase of CSI for each subcarrier')
%     xlabel('Subcarriers')
%     ylabel('Phase')

    %% Smoothed CSI Matrix
    % csi_matrix(:,:,idx);
    smoothed_csi = zeros(30,86);

    for i=1:15
        smoothed_csi(i,:) = [csi_matrix(1,i:i+42),csi_matrix(2,i:i+42)]; 
    end

    for i=16:30
        smoothed_csi(i,:) = [csi_matrix(2,i-15:i+27),csi_matrix(3,i-15:i+27)]; 
    end

    %% Extract Steering Matrix
    % For a normal matrix
    % idx = 1;

    % X = csi_matrix(:,:)';
    % MUSIC_1 = X*X';

    MUSIC_S = smoothed_csi * smoothed_csi';

    [EigenVector1, EigenValue1] = eig(MUSIC_S);
    EigenValueList1 = diag(EigenValue1);

    [~, order] = sort(EigenValueList1);

    EigenVector = EigenVector1(:,order);
    % EigenValue = EigenValue1(:,order);

%     figure(8);
%     plot(EigenValueList1)
%     title('Eigenvalues of Matrix S')
%     ylabel('Eigenvalues')

    % EigenValueThreshold = 10;
    % NoiseVectorIdx=find(EigenValueList1<EigenValueThreshold);
    NoiseVectorIdx=1:22;

    %% M sensors, D multipath, N zero eigenvalues
    % M = N + D
    % Empirically, D=6~8 significant reflectors (SpotFi 3.1)
    % M = 30; % sensors, elements
    % D = 15;  % multipath
    % N = M - D;  % size of null space

    En = EigenVector(:,NoiseVectorIdx);
    % A = EigenVector(:,N+1:M); % A'*E=0  -> E'*A=0
    % A'*En
    %% Algorithm

    degrees = -90:1:90; % linespace
    [~,deg_tot] = size(degrees);
    tofs = 1e-9:1e-9:300e-9;
    [~,tof_tot] = size(tofs);

    c = 3e8; % speed of light
    % f = 2.412e9; % central frequency
    f = 2.437e9;
    fs = 312.5e3;
    d = 0.07;
    SP = zeros(deg_tot,tof_tot);
    % tic

    for deg_idx = 1:deg_tot
        deg = degrees(deg_idx);
        theta = deg*pi/180;
        phi = exp(-1j*2*pi*d*sin(theta)*f/c);
        for tof_idx = 1:tof_tot
            tof = tofs(tof_idx);
            omega = exp(-1j*2*pi*fs*tof);
            half = omega.^((0:14)');
            a = [half;
                half.*phi];

            % SP(deg_idx,tof_idx)=(a'*a)/(a'*En*En'*a);
            SP(deg_idx,tof_idx)=1/(a'*En*En'*a);
        end
    end
    % toc

    Pmu=abs(SP);
    Pmu_max=max(Pmu);
    % Pmu_db=10*log10(Pmu./Pmu_max);

%     figure(11);
%     surf(tofs*1e9,degrees,Pmu)
%     xlabel('ToF (ns)')
%     ylabel('AoA (degree)') 
%     zlabel('magnitude (dB)')
%     grid on  
%     drawnow;

    % figure, plot(degrees,Pmu_db(:,sig_tof));
    % xlabel('AoA')
    % ylabel('Magnitude when Tof=10 ns')



    %% One side
    % oneside = fix(deg_tot/4):fix(deg_tot*3/4);
    % rad_to_plot = 3/2*pi-degrees(oneside)'*pi/180;
    % Pmu_to_plot = Pmu(oneside,sig_tof);
    % % polar(rad_to_plot,Pmu_to_plot);
    % 
    % [~,sortk] = sort(Pmu_to_plot);
    % aoa = rad_to_plot(sortk(end))*180/pi;
    % rads = rad_to_plot;
    % pmu = Pmu_to_plot;

    rads = degrees'*pi/180;

end
