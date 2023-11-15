function [distance] = GOSPAT(TrackSet1,TrackSet2, p, c, gamma, Time)
% Reference:
% A. F. Garcia-Fernandez, A. S. Rahmathullah and L. Svensson, "A Metric on
% the Space of Finite Sets of Trajectories for Evaluation of Multi-Target
% Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68,
% pp. 3917-3928, 2020, doi: 10.1109/TSP.2020.2005309.
% c: cut-off parameter
% p: power
% gamma: switching penalty
% Time: the time indices
TimeLength = length(Time);
TrackSetLen1 = length(TrackSet1);
TrackSetLen2 = length(TrackSet2);
pi = cell(1,TimeLength);
distance1 = 0;
for k = 1:TimeLength
    % obtain the assignment vector
    StateSet1 = [];
    StateIndex1 = [];
    StateSet2 = [];
    StateIndex2 = [];
    for j = 1:TrackSetLen1
        t = TrackSet1{j}(1,:);
        index = find(t == Time(k));
        if ~isempty(index)
            StateSet1 = [StateSet1, TrackSet1{j}(2:end,index)];
            StateIndex1 = [StateIndex1, j];
        end
    end
    for j = 1:TrackSetLen2
        t = TrackSet2{j}(1,:);
        index = find(t == Time(k));
        if ~isempty(index)
            StateSet2 = [StateSet2, TrackSet2{j}(2:end,index)];
            StateIndex2 = [StateIndex2, j];
        end
    end
    StateIndexLen1 = length(StateIndex1);
    StateIndexLen2 = length(StateIndex2);
    if StateIndexLen1 <= StateIndexLen2
        CostMatrix = zeros(StateIndexLen1, StateIndexLen1 + StateIndexLen2);
        for i = 1:StateIndexLen1
            for j = 1:(StateIndexLen1+StateIndexLen2)
                if j <= StateIndexLen2
                    X1 = StateSet1(:,i);
                    X2 = StateSet2(:,j);
                    CostMatrix(i,j) = sum(abs(X1-X2));
                else
                    if i == (j - StateIndexLen2)
                        CostMatrix(i,j) = c;
                    else
                        CostMatrix(i,j) = inf;
                    end
                end
            end
        end
        [~, AssignmentMatrix] = HungarianAlgorithm(CostMatrix);
        pi_new = zeros(0,2);
        theta_count = 0;
        temp = 0;
        for i = 1:StateIndexLen1
            tempIndex = find(AssignmentMatrix(i,:) == 1);
            if ~isempty(tempIndex) && tempIndex <= StateIndexLen2
                pi_new = [pi_new; StateIndex1(i), StateIndex2(tempIndex)];
                theta_count = theta_count + 1;
                temp = temp + CostMatrix(i, tempIndex)^p;
            else
                pi_new = [pi_new; StateIndex1(i), 0];
            end
        end
    else
        CostMatrix = zeros(StateIndexLen1 + StateIndexLen2, StateIndexLen2);
        for i = 1:(StateIndexLen1 + StateIndexLen2)
            for j = 1:StateIndexLen2
                if i <= StateIndexLen1
                    X1 = StateSet1(:,i);
                    X2 = StateSet2(:,j);
                    CostMatrix(i,j) = sum(abs(X1-X2));
                else
                    if j == (i - StateIndexLen1)
                        CostMatrix(i,j) = c;
                    else
                        CostMatrix(i,j) = inf;
                    end
                end
            end
        end
        [~, AssignmentMatrix] = HungarianAlgorithm(CostMatrix');
        AssignmentMatrix = AssignmentMatrix';
        pi_new = zeros(0,2);
        theta_count = 0;
        temp = 0;
        for i = 1:StateIndexLen1
            tempIndex = find(AssignmentMatrix(i,:) == 1);
            if ~isempty(tempIndex)
                pi_new = [pi_new; StateIndex1(i), StateIndex2(tempIndex)];
                theta_count = theta_count + 1;
                temp = temp + CostMatrix(i, tempIndex)^p;
            else
                pi_new = [pi_new; StateIndex1(i), 0];
            end
        end
    end
    pi_final = [];
    for i = 1:TrackSetLen1
        tempIndex = find(pi_new(:,1) == i);
        if isempty(tempIndex)
            pi_final = [pi_final; i, 0];
        else
            pi_final = [pi_final; pi_new(tempIndex,:)];
        end
    end
    distance1 = distance1 + temp + c^p/2 * (StateIndexLen1 + StateIndexLen2 - 2*theta_count);
    pi{k} = pi_final;
end
distance2 = 0;
for k = 1:(TimeLength-1)
    pi_old = pi{k};
    pi_new = pi{k+1};
    temp = 0;
    for i = 1:TrackSetLen1
        x_idx_old = find(pi_old(:,1) == i);
        x_idx_new = find(pi_new(:,1) == i);
        y_idx_old = pi_old(x_idx_old,2);
        y_idx_new = pi_new(x_idx_new,2);
        if y_idx_old == y_idx_new
            temp = temp + 0;
        elseif y_idx_old ~= y_idx_new && y_idx_old ~= 0 && y_idx_new ~= 0
            temp = temp + 1;
        else
            temp = temp + 1/2;
        end
    end
    distance2 = distance2 + temp * gamma^p;
end
distance = (distance1 + distance2)^(1/p);
end

