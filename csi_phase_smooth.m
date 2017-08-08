function fii = csi_phase_smooth(fi,Nrx)

fii = fi;

    for rx = 1:Nrx
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

end