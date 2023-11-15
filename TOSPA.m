function [ OSPA ] = TOSPA( TrackSet1, TrackSet2, p, alpha, cc, Time )
% The implementation of algorithm on TOSPA
% REFERFENCE:
% Branko Ristic, Ba-Ngu Vo, Daniel Clark. Performance evaluation of
% multi-target tracking using the OSPA metric.

Nt1 = length(TrackSet1);
Nt2 = length(TrackSet2);
TimeLen = length(Time);
% Computing association matrix:
c = zeros(Nt1,Nt2);
for i = 1:Nt1
    for j = 1:Nt2
        tempC1 = 0;
        tempC2 = 0;
        T1 = TrackSet1{i};
        T2 = TrackSet2{j};
        for k = 1:TimeLen
            eks = find(T1(1,:) == Time(k));
            ekt = find(T2(1,:) == Time(k));
            if ~isempty(eks) && ~isempty(ekt)
                xks = T1(2:end,eks);
                xkt = T2(2:end,ekt);
                tempC1 = tempC1 + (xks - xkt)' * (xks - xkt);
            end
            eks = ~isempty(eks);
            ekt = ~isempty(ekt);
            tempC2 = tempC2 + eks * ekt;
        end
        if tempC2 > 0
            c(i,j) = tempC1 / exp(tempC2);
        else
            c(i,j) = inf;
        end
    end
end
% Obtain the best association.
[~, idx] = HungarianAlgorithm(c);
BestAssociation = zeros(1,Nt1);
parfor j = 1:Nt1
    tempIndex = find(idx(j,:)==1);
    if isempty(tempIndex)
        BestAssociation(j) = 0;
    else
        BestAssociation(j) = tempIndex(1);
    end
end
% Computing the OSPA value:
OSPA = zeros(1,length(Time));
for k = 1:length(Time)
    X = [];
    Y = [];
    for i = 1:length(TrackSet1)
        T1 = TrackSet1{i};
        tempIndex = find(T1(1,:) == Time(k));
        if ~isempty(tempIndex)
            X = [X, T1(2:end,tempIndex(1))];
        end
    end
    for i = 1:length(TrackSet2)
        T2 = TrackSet2{i};
        tempIndex = find(T2(1,:) == Time(k));
        if ~isempty(tempIndex)
            Y = [Y, T2(2:end,tempIndex(1))];
        end
    end
    m = size(X,2);
    n = size(Y,2);
    if m == 0 && n == 0
        OSPA(k) = 0;
        continue;
    end
    if isempty(X) || isempty(Y)
        if isempty(X) && isempty(Y)
            distance = 0;
            OSPA(k) = distance;
            continue;
        else
            distance = cc;
            OSPA(k) = distance;
            continue;
        end
    end
    A = zeros(m,n);
    for i = 1:m
        for j = 1:n
            try
                A(i,j) = (X(:,i) - Y(:,j))' * (X(:,i) - Y(:,j));
                A(i,j) = sqrt(A(i,j));
            catch
                pause;
            end
        end
    end
    if m > n
        [cost, idxMatrix] = HungarianAlgorithm(A);
        minIndex = zeros(1,m);
        parfor j = 1:m
            tempIndex = find(idxMatrix(j,:)==1);
            if isempty(tempIndex)
                minIndex(j) = 0;
            else
                minIndex(j) = tempIndex;
            end
        end
        minScore = 0;
        for j = 1:m
            if minIndex(j) == 0
                minScore = minScore + cc^p;
            else
                tempA = A(j,minIndex(j));
                tempA = min(tempA, cc);
                minScore = minScore + tempA^p;
            end
        end
        tempSumAlpha = 0;
        for i = 1:m
            if minIndex(i) == BestAssociation(i) && minIndex(i) ~= 0
                tempSumAlpha = tempSumAlpha + 0;
            else
                tempSumAlpha = tempSumAlpha + alpha^p;
            end
        end
        minScore = minScore + tempSumAlpha;
        distance = (minScore / m)^(1/p);
    else
        At = A';
        [cost, idxMatrix] = HungarianAlgorithm(At);
        minIndex = zeros(1,n);
        parfor j = 1:n
            tempIndex = find(idxMatrix(j,:)==1);
            if isempty(tempIndex)
                minIndex(j) = 0;
            else
                minIndex(j) = tempIndex;
            end
        end
        minScore = 0;
        for j = 1:n
            if minIndex(j) == 0
                minScore = minScore + cc^p;
            else
                tempA = At(j,minIndex(j));
                tempA = min(tempA, cc);
                minScore = minScore + tempA;
            end
        end
        tempSumAlpha = 0;
        for i = 1:n
            tempAssociationIndex = find(BestAssociation == i);
            if isempty(tempAssociationIndex)
                tempSumAlpha = tempSumAlpha + alpha^p;
            else
                tempAssociationIndex = tempAssociationIndex(1);
                if minIndex(i) ~= tempAssociationIndex
                    tempSumAlpha = tempSumAlpha + alpha^p;
                else
                    tempSumAlpha = tempSumAlpha + 0;
                end
            end
        end
        minScore = minScore + tempSumAlpha;
        distance = (minScore / n)^(1/p);
    end
    OSPA(k) = distance;
end

