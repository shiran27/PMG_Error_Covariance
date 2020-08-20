function cycle = milpSol(graph, numOfVisits, timeLimit, numOfApproxPoints)
        
    numOfTargets = length(graph.targets); % M
    numOfVisits; % k
    numOfApproxPoints; % P
    timeLimit; % T
    bigM = 100; % big M 
    depot = numOfTargets; % last target in the graph is to be visited last
    distMat = graph.distanceMatrix;

    % Constraint Set 1
    A1 = zeros(numOfTargets*(numOfVisits-1),numOfTargets*numOfTargets*numOfVisits);
    b1 = zeros(numOfTargets*(numOfVisits-1),1);
    count = 1;
    for i = 1:1:numOfTargets
        for v = 2:1:numOfVisits
            A = zeros(numOfTargets,numOfTargets,numOfVisits);
            A(:,i,v-1) = ones(numOfTargets,1);
            A(i,:,v) = -1*ones(1,numOfTargets);
            A1(count,:) = transpose(A(:));
            count = count + 1;
        end
    end
    A1_1 = A1;
    A1_2 = zeros(numOfTargets*(numOfVisits-1),numOfTargets*numOfVisits);
    A1_3 = zeros(numOfTargets*(numOfVisits-1),numOfTargets*numOfVisits);
    A1_4 = zeros(numOfTargets*(numOfVisits-1),1);
    A1 = [A1_1, A1_2, A1_3, A1_4]; % A1 x = b1



    % Constraint Set 2
    A2 = zeros(numOfTargets*numOfVisits,numOfTargets*numOfTargets*numOfVisits);
    b2 = zeros(numOfTargets*numOfVisits,1);
    count = 1;
    for i = 1:1:numOfTargets
        for v = 1:1:numOfVisits
            A = zeros(numOfTargets,numOfTargets,numOfVisits);
            A(i,i,v) = 1;
            A2(count,:) = transpose(A(:));
            count = count + 1;
        end
    end
    A2_1 = A2;
    A2_2 = zeros(numOfTargets*numOfVisits,numOfTargets*numOfVisits);
    A2_3 = zeros(numOfTargets*numOfVisits,numOfTargets*numOfVisits);
    A2_4 = zeros(numOfTargets*numOfVisits,1);
    A2 = [A2_1, A2_2, A2_3, A2_4]; % A2 x = b2


    % Constraint Set 3
    A3 = zeros(2,numOfTargets*numOfTargets*numOfVisits);
    b3 = ones(2,1);

    A = zeros(numOfTargets,numOfTargets,numOfVisits);
    A(depot,:,1) = ones(1,numOfTargets);
    A3(1,:) = transpose(A(:));  

    A = zeros(numOfTargets,numOfTargets,numOfVisits);
    A(:,depot,numOfVisits) = ones(numOfTargets,1);
    A3(2,:) = transpose(A(:));

    A3_1 = A3;
    A3_2 = zeros(2,numOfTargets*numOfVisits);
    A3_3 = zeros(2,numOfTargets*numOfVisits);
    A3_4 = zeros(2,1);
    A3 = [A3_1, A3_2, A3_3, A3_4]; % A3 x = b3


    % Constraint Set 4
    A4 = zeros(numOfTargets,numOfTargets*numOfTargets*numOfVisits);
    b4 = ones(numOfTargets,1);
    count = 1;
    for i = 1:1:numOfTargets
        A = zeros(numOfTargets,numOfTargets,numOfVisits);
        A(i,:,:) = ones(numOfTargets,numOfVisits);
        A4(count,:) = transpose(A(:));
        count = count + 1;
    end
    A4_1 = A4;
    A4_2 = zeros(numOfTargets,numOfTargets*numOfVisits);
    A4_3 = zeros(numOfTargets,numOfTargets*numOfVisits);
    A4_4 = zeros(numOfTargets,1);
    A4 = [A4_1, A4_2, A4_3, A4_4]; % A4 x = b4



    % Constraint Set 5
    A5 = zeros(numOfVisits,numOfTargets*numOfTargets*numOfVisits);
    b5 = zeros(numOfVisits,1);
    count = 1;
    for v = 1:1:numOfVisits
        A = zeros(numOfTargets,numOfTargets,numOfVisits);
        A(:,:,v) = 1;
        A5(count,:) = transpose(A(:));
        count = count + 1;
    end
    A5_1 = A5;
    A5_2 = zeros(numOfVisits,numOfTargets*numOfVisits);
    A5_3 = zeros(numOfVisits,numOfTargets*numOfVisits);
    A5_4 = zeros(numOfVisits,1);
    A5 = [A5_1, A5_2, A5_3, A5_4]; % A5 x = b5



    % Constraint Set 7
    A7_1 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfTargets*numOfVisits));
    A7_2 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfVisits));
    A7_3 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfVisits));
    b7 = zeros(numOfTargets*(numOfVisits-1),1);
    count = 1;
    for v = 2:1:numOfVisits
        A_1 = zeros(numOfTargets,numOfTargets,numOfVisits);
        A_1(:,:,v) = -1*ones(numOfTargets,numOfTargets).*distMat;
        for i = 1:1:numOfTargets
            A7_1(count,:) = transpose(A_1(:));

            A_2 = zeros(numOfTargets,numOfVisits);
            A_2(i,v) = 1;
            A_2(i,v-1) = -1;
            A7_2(count,:) = transpose(A_2(:));

            A_3 = zeros(numOfTargets,numOfVisits);
            A_3(i,v-1) = 1;
            A7_3(count,:) = transpose(A_3(:));

            count = count + 1;
        end
    end
    A7_4 = zeros(numOfTargets*(numOfVisits-1),1);
    A7 = [A7_1, A7_2, A7_3, A7_4]; % A7 x = b7



    % Constraint Set 8 and 9
    A8_1 = zeros(numOfTargets,(numOfTargets*numOfTargets*numOfVisits));
    A8_2 = zeros(numOfTargets,(numOfTargets*numOfVisits));
    b8 = zeros(numOfTargets,1);
    count = 1;

    A_1 = zeros(numOfTargets,numOfTargets,numOfVisits);
    A_1(depot,:,1) = -1*ones(1,numOfTargets).*distMat(depot,:);

    for i = 1:1:numOfTargets
        A8_1(count,:) = transpose(A_1(:));

        A_2 = zeros(numOfTargets,numOfVisits);
        A_2(i,1) = 1;
        if i~=depot
            A_2(i,numOfVisits) = -1;
        end
        A8_2(count,:) = transpose(A_2(:));

        count = count + 1;
    end
    A8_3 = zeros(numOfTargets,numOfTargets*numOfVisits);
    A8_4 = zeros(numOfTargets,1);
    A8 = [A8_1, A8_2, A8_3, A8_4];% A8 x = b8



    % Constraint Set 14
    A14_1 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfTargets*numOfVisits));
    A14_2 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfVisits));
    A14_3 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfVisits));
    b14 = bigM*ones(numOfTargets*(numOfVisits-1),1);
    count = 1;
    for v = 2:1:numOfVisits
        for i = 1:1:numOfTargets
            A_1 = zeros(numOfTargets,numOfTargets,numOfVisits);
            A_1(:,i,v-1) = bigM*ones(numOfTargets,1);
            A14_1(count,:) = transpose(A_1(:));

            A_2 = zeros(numOfTargets,numOfVisits);
            A_2(i,v-1) = -1;
            A14_2(count,:) = transpose(A_2(:));

            A_3 = zeros(numOfTargets,numOfVisits);
            A_3(i,v-1) = -1;
            A14_3(count,:) = transpose(A_3(:));

            count = count + 1;
        end
    end
    A14_4 = zeros(numOfTargets*(numOfVisits-1),1);
    A14 = [A14_1, A14_2, A14_3, A14_4]; % A14 x <= b14



    % Constraint Set 15
    A15_1 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfTargets*numOfVisits));
    A15_3 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfVisits));
    b15 = zeros(numOfTargets*(numOfVisits-1),1);
    count = 1;
    for v = 2:1:numOfVisits
        for i = 1:1:numOfTargets
            A_1 = zeros(numOfTargets,numOfTargets,numOfVisits);
            A_1(:,i,v-1) = -1*bigM*ones(numOfTargets,1);
            A15_1(count,:) = transpose(A_1(:));

            A_3 = zeros(numOfTargets,numOfVisits);
            A_3(i,v-1) = 1;
            A15_3(count,:) = transpose(A_3(:));

            count = count + 1;
        end
    end
    A15_2 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfVisits));
    A15_4 = zeros(numOfTargets*(numOfVisits-1),1);
    A15 = [A15_1, A15_2, A15_3, A15_4]; % A15 x <= b15



    % Constraint Set 16
    A16_2 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfVisits));
    A16_3 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfVisits));
    b16 = zeros(numOfTargets*(numOfVisits-1),1);
    count = 1;
    for v = 2:1:numOfVisits
        for i = 1:1:numOfTargets
            A_2 = zeros(numOfTargets,numOfVisits);
            A_2(i,v-1) = -1;
            A16_2(count,:) = transpose(A_2(:));

            A_3 = zeros(numOfTargets,numOfVisits);
            A_3(i,v-1) = 1;
            A16_3(count,:) = transpose(A_3(:));

            count = count + 1;
        end
    end
    A16_1 = zeros(numOfTargets*(numOfVisits-1),(numOfTargets*numOfTargets*numOfVisits));
    A16_4 = zeros(numOfTargets*(numOfVisits-1),1);
    A16 = [A16_1, A16_2, A16_3, A16_4]; % A16 x <= b16



    % Constraint Set 12
    A12_2 = zeros(numOfTargets*numOfVisits*(numOfApproxPoints+1), numOfTargets*numOfVisits);
    A12_4 = zeros(numOfTargets*numOfVisits*(numOfApproxPoints+1), 1);
    b12 = zeros(numOfTargets*numOfVisits*(numOfApproxPoints+1),1);
    timeStep = timeLimit/numOfApproxPoints;
    count = 1; 
    for i = 1:1:numOfTargets 
        for p = 0:1:numOfApproxPoints

            [m_ip, c_ip] = graph.targets(i).computeGradientAndIntersection(timeStep*p);

            for v = 1:1:numOfVisits
                A_2 = zeros(numOfTargets,numOfVisits);
                A_2(i,v) = m_ip;
                A12_2(count,:) = transpose(A_2(:));

                A_4 = -1;
                A12_4(count,:) = transpose(A_4(:));

                b12(count,1) = c_ip;

                count = count + 1;
            end
        end
    end
    A12_1 = zeros(numOfTargets*numOfVisits*(numOfApproxPoints+1),(numOfTargets*numOfTargets*numOfVisits));
    A12_3 = zeros(numOfTargets*numOfVisits*(numOfApproxPoints+1), numOfTargets*numOfVisits);
    A12 = [A12_1, A12_2, A12_3, A12_4]; % A12 x <= b12


    % objective
    numOfVars = numOfTargets*numOfTargets*numOfVisits + 2*numOfTargets*numOfVisits + 1
    f = [zeros(1,(numOfVars-1)), 1];
    intcon = 1:(numOfTargets*numOfTargets*numOfVisits);
    Aeq = [A1; A2; A2; A4; A5; A7; A8];
    size(Aeq,1)
    beq = [b1; b2; b2; b4; b5; b7; b8]; % Aeq x = beq
    A = [A14; A15; A16; A12];
    size(A,1)
    b = [b14; b15; b16; b12]; % A x <= b
    Lb = zeros(numOfVars,1);
    Ub = [ones(numOfTargets*numOfTargets*numOfVisits,1); timeLimit*ones(2*numOfTargets*numOfVisits,1); 100];

    [x1,fval1,exitflag1,output1] = intlinprog(f,intcon,A,b,Aeq,beq,Lb,Ub)
end
