function [OSPA2_output] = OSPA2Visualization(TrackSet1,TrackSet2, p, c, N, Time)
% The implementation of OSPA^(2)
% Reference:
% M. Beard, B. T. Vo and B. Vo, "A Solution for Large-Scale Multi-Object Tracking," in IEEE Transactions on Signal Processing, vol. 68, pp. 2754-2769, 2020, doi: 10.1109/TSP.2020.2986136.
TimeLength = length(Time);
OSPA2_output = zeros(1,TimeLength);
Nt1 = length(TrackSet1);
Nt2 = length(TrackSet2);
TimeInit = Time(1);
Step = Time(2) - Time(1);
parfor t = 1:TimeLength
    tempTrackSet1 = cell(0,0);
    tempTrackSet2 = cell(0,0);
    tempTime = max(Time(t)-N*Step,TimeInit):Step:Time(t);
    count = 1;
    for idx = 1:Nt1
        t1 = TrackSet1{idx};
        time1 = t1(1,:);
        tempIdx = [];
        for i = 1:length(tempTime)
            index = find(time1 == tempTime(i));
            tempIdx = [tempIdx, index];
        end
        if ~isempty(tempIdx)
            tempTrackSet1{count} = t1(:,tempIdx);
            count = count + 1;
        end
    end
    count = 1;
    for idx = 1:Nt2
        t2 = TrackSet2{idx};
        time2 = t2(1,:);
        tempIdx = [];
        for i = 1:length(tempTime)
            index = find(time2 == tempTime(i));
            tempIdx = [tempIdx, index];
        end
        if ~isempty(tempIdx)
            tempTrackSet2{count} = t2(:,tempIdx);
            count = count + 1;
        end
    end
    OSPA2_output(t) = OSPA2(tempTrackSet1, tempTrackSet2, p, c);
end
end

