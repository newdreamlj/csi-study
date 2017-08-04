
function maxima = find_maxima(Pmu)
    %% test
%     A =[1 1 1 1 1 1 1 1;
%         1 2 2 7 2 1 1 1;
%         1 2 3 9 9 3 2 1;
%         1 2 2 9 9 1 1 1;
%         1 2 3 1 1 1 1 1;
%         1 1 1 1 1 1 1 1];
% 
%     
%     [lr_c,lr_r] = find(diff(2*(sign(diff(A,1,1))==1)-1,1,1)==-2);
%     lr_c = lr_c + 1;
%     [ud_c,ud_r] = find(diff(2*(sign(diff(A,1,2))==1)-1,1,2)==-2);
%     ud_r = ud_r + 1;
%     
% %     [lr_c,lr_r] = find(diff(sign(diff(A,1,1)),1,1)==-2);
% %     lr_c = lr_c + 1;
% %     [ud_c,ud_r] = find(diff(sign(diff(A,1,2)),1,2)==-2);
% %     ud_r = ud_r + 1;
%     
%     AA = zeros(size(A));
%     for idx=1:length(lr_c)
%         AA(lr_c(idx),lr_r(idx)) = 1;    
%     end
%     for idx=1:length(ud_c)
%         AA(ud_c(idx),ud_r(idx)) = AA(ud_c(idx),ud_r(idx)) + 1;
%     end
%     
%     [maxima_c,maxima_r] = find(AA==2);
%     maxima = [maxima_c maxima_r];
%     
    %% implementation
    
    A = Pmu;
    [lr_c,lr_r] = find(diff(2*(sign(diff(A,1,1))==1)-1,1,1)==-2);
    lr_c = lr_c + 1;
    [ud_c,ud_r] = find(diff(2*(sign(diff(A,1,2))==1)-1,1,2)==-2);
    ud_r = ud_r + 1;
    
    AA = zeros(size(A));
    for idx=1:length(lr_c)
        AA(lr_c(idx),lr_r(idx)) = 1;    
    end
    for idx=1:length(ud_c)
        AA(ud_c(idx),ud_r(idx)) = AA(ud_c(idx),ud_r(idx)) + 1;
    end
    
    [maxima_c,maxima_r] = find(AA==2);
    maxima = [maxima_c maxima_r];

end