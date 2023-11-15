function [OSPA2_output] = OSPA2(TrackSet1, TrackSet2, p, c)
% The implementation of OSPA^(2)
% Reference:
% M. Beard, B. T. Vo and B. Vo, "A Solution for Large-Scale Multi-Object Tracking," in IEEE Transactions on Signal Processing, vol. 68, pp. 2754-2769, 2020, doi: 10.1109/TSP.2020.2986136.
Nt1 = length(TrackSet1);
Nt2 = length(TrackSet2);
CostMatrix = zeros(Nt1,Nt2);
for idx1 = 1:Nt1
    for idx2 = 1:Nt2
        t1 = TrackSet1{idx1};
        t2 = TrackSet2{idx2};
        time1 = t1(1,:);
        time2 = t2(1,:);
        timeUnion = union(time1, time2);
        timeUnionLen = length(timeUnion);
        d = zeros(1,timeUnionLen);
        parfor k = 1:timeUnionLen
            tempTime = timeUnion(k);
            index1 = find(time1 == tempTime);
            index2 = find(time2 == tempTime);
            if isempty(index1) && isempty(index2)
                d(k) = 0;
            elseif length(index1) ~= length(index2)
                d(k) = c;
            else
                X1 = t1(2:end, index1);
                X2 = t2(2:end, index2);
                d(k) = min(c, sum(abs(X1-X2)));
            end
        end
        CostMatrix(idx1, idx2) = sum(d) / timeUnionLen;
    end
end
CostMatrix = CostMatrix.^p;
[cost, ~] = HungarianAlgorithm(CostMatrix);
if Nt1 > Nt2
    OSPA2_output = ((cost + c^p*(Nt1-Nt2))/Nt1)^(1/p);
else
    OSPA2_output = ((cost + c^p*(Nt2-Nt1))/Nt2)^(1/p);
end
end

