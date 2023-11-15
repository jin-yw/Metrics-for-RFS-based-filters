function [ distance ] = OSPA( X, Y, p, c )
%% The implementation of algorithm on OSPA indicator 
m = size(Y, 2);         % Number of estimations
n = size(X, 2);         % Number of real states
if isempty(X) || isempty(Y)
    if isempty(X) && isempty(Y)
        distance = 0;
    else
        distance = c;
    end
else
    %% Obtaining the association matrix
    A = zeros(m,n);
    for i = 1:m
        for j = 1:n
            A(i,j) = (X(:,j) - Y(:,i))'* (X(:,j) - Y(:,i));     % Euclidean distance
            A(i,j) = sqrt(A(i,j));
        end
    end
    if m > n
        Score = [];
        Matrix = cell(0,0);
        SelectionIndex = cell(0,0);
        count = 1;
        for i = 1:m
            if i == 1
                Score(count) = c^p;
                Matrix{count} = A;
                Matrix{count}(1,:) = -1;
                SelectionIndex{count} = 0;
                count = count + 1;
                for j = 1:n
                    Score(count) = A(i,j)^p;
                    Matrix{count} = A;
                    Matrix{count}(:,j) = -1;
                    Matrix{count}(i,:) = -1;
                    SelectionIndex{count} = j;
                    count = count + 1;
                end
            else
                tempLength = length(SelectionIndex);
                for j = 1:tempLength
                    tempMatrix = Matrix{j};
                    for k = 1:n
                        if tempMatrix(i,k) == -1
                            continue;
                        else
                            Score(count) = Score(j) + tempMatrix(i,k)^p;
                            Matrix{count} = tempMatrix;
                            Matrix{count}(i,:) = -1;
                            Matrix{count}(:,k) = -1;
                            SelectionIndex{count} = [SelectionIndex{j}, k];
                            count = count + 1;
                        end
                    end
                    Score(j) = Score(j) + c^p;
                    SelectionIndex{j} = [SelectionIndex{j},0];
                    Matrix{j}(i,:) = -1;
                end
            end
        end
        minScore = min(Score);
        distance = (minScore / n)^(1/p);
    else
        Score = [];
        Matrix = cell(0,0);
        SelectionIndex = cell(0,0);
        count = 1;
        for i = 1:n
            if i == 1
                Score(count) = c^p;
                Matrix{count} = A;
                Matrix{count}(:,1) = -1;
                SelectionIndex{count} = 0;
                count = count + 1;
                for j = 1:m
                    Score(count) = A(j,i)^p;
                    Matrix{count} = A;
                    Matrix{count}(:,1) = -1;
                    Matrix{count}(j,:) = -1;
                    SelectionIndex{count} = j;
                    count = count + 1;
                end
            else
                tempLength = length(SelectionIndex);
                for j = 1:tempLength
                    tempMatrix = Matrix{j};
                    for k = 1:m
                        if tempMatrix(k,i) == -1
                            continue;
                        else
                            Score(count) = Score(j) + tempMatrix(k,i)^p;
                            Matrix{count} = tempMatrix;
                            Matrix{count}(k,:) = -1;
                            Matrix{count}(:,i) = -1;
                            SelectionIndex{count} = [SelectionIndex{j},k];
                            count = count + 1;
                        end
                    end
                    Matrix{j}(:,i) = -1;
                    Score(j) = Score(j) + c^p;
                    SelectionIndex{j} = [SelectionIndex{j},0];
                end
            end
        end
        minScore = min(Score);
        distance = (minScore / n)^(1/p);
    end
end


end

