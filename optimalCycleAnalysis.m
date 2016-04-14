function optimalCycleAnalysis

    clear;

    TRUE = uint8(1);
    FALSE = uint8(0);
    ROUND_ERROR = 0.0000001;

    dataFaciesCodes = zeros(1,1);
    sectionFaciesNames = cell(1,100);
    dataThick = zeros(1,1);
    sectLength = 0.0;
    markovOrderMetric = zeros(1,1);
    j = uint32(0);

    % Read in the section facies and thickness data - you might need to edit the directory names
    fileName = input('Enter succession data path and file name (without any .dat or .txt extension)\ne.g. ../outcropSections/test100 ../syntheticSections/syntheticOrdered100units10facies\nFilename is : ','s');

    % Read section data, thickness and facies
    fullFileName = strcat(fileName, '.txt');
    if exist(fullFileName, 'file')
        dataIn = load(fullFileName);
        dataThick = dataIn(:,1);
        dataFaciesCodes = dataIn(:,2);
    else
        fprintf('Succession data file %s does not exist. Please run this script again and enter a valid path and filename\n', fullFileName);
        return;
    end

    tic
    
    % Read colour map file, one colour and one text name per facies
    fullColourMapFileName = strcat(fileName, 'LithoCol.txt');
    if exist(fullColourMapFileName, 'file');
        fid = fopen(fullColourMapFileName);
        while ~feof(fid)
            file = textscan(fid, '%u8 %f %f %f %s');
        end
        fclose(fid);
    else
        fprintf('Facies colour map and name file %s does not exist. Please run this script again and enter a valid path and filename\n', fullColourMapFileName);
        return;
    end
    
    % Extract facies names, numeric codes and RGB colours from data input file as a cell array
    faciesNames = file {5}();
    faciesNumbers = file {1}; % particularly important - the original facies codes
    for k=1:4
        faciesColours(:,k)= double(file{k}());
    end

    % Calculate the basic stats on the dataFacies array
    sectionLength = max(size(dataFaciesCodes));
    maxNumbOfFacies = length(faciesColours);
    sectionMean = mean(dataThick);

    % Make a reference section that uses the facies names to ID each unit
    % rather than the numerical facies code. Need to remember that faciesNames
    % is a cell array hence the curly brackets
    % Note that as of 12.7.2015 this array is not used elsewhere in the code
    for k=1:sectionLength
        sectionFaciesNames{k} = faciesNames{dataFaciesCodes(k)};
    end

    fprintf('For %d units, total %d facies, mean unit thickness %4.3f m\n', sectionLength, maxNumbOfFacies, sectionMean);

    % Calculate and output the order metric for the entire data succession
    markovOrderMetric = calculateTPMatrixAndOrderMetric(dataFaciesCodes, 0);  % note the zero is the TP matrix plot flag, set to false to not trigger a plot
    fprintf('Markov metric for strata is %5.4f\n', markovOrderMetric);

    allFaciesCombos = perms(faciesNumbers);
    permsCount = length(allFaciesCombos);

    % allFaciesCombos now contains all the possible ways the n facies classes read in can be positioned in the TP matrix
    % In other words, a permutation 5 1 2 3 4 means that facies 1 is coded as facies 5 so that it positioned in row 5 of the TP matrix for this permutation
    % Now make a strat section in testFaciesSect using the facies coding in each row of the matrix and calculate the markov statistic for each testFaciesSect
    % Define the necessary variables to make a test section and record markov
    % results from each test section. 
    testFaciesSect = zeros(1, sectionLength); % the facies section input but will contain facies coding from current permutation
    markovOrderMetricTest = zeros(1,permsCount);
    oneOffsetDiagCheck = zeros(1,permsCount);
    markovDiagProduct = zeros(1,permsCount);
    maxMarkov = 0.0;
    minMarkov = 1000.0;
    maxDiag = 0.0;
    maxMarkovDiagProduct = 0.0;
    
    % Set variables to monitor the progress of the calculation
    progressCount = 0;
    progressIndicator = permsCount / 100.0;

    fprintf('Calculating M stats for %d permutations, * every %d iterations...\n', permsCount, progressIndicator);
    % Print a string of * symbols length permsCount as a marker for the *
    % progress indicator printed below in the main calculation loop
    for j = 1:100.0
        fprintf('*');
    end
    fprintf('Done\n');

    for j = 1:permsCount

        % Make a vertical succession with each facies occurrence coded
        % according to the facies coding in permutation j, not the original coding
        for k =1:sectionLength
            oneFacies = dataFaciesCodes(k);
            testFaciesSect(k) = allFaciesCombos(j, oneFacies); % Swap the original facies code from the data section with a substitute code from permutation j
        end

        % Calculate the markov statistic and the one offset diagonal for this permutation
        [markovOrderMetricTest(j), oneOffsetDiagCheck(j)] = calculateTPMatrixAndOrderMetric(testFaciesSect, 0);
        markovDiagProduct(j) = markovOrderMetricTest(j) * oneOffsetDiagCheck(j);
        
        % if the result is the highest markov stat yet, remember it
        if markovOrderMetricTest(j) >= maxMarkov
            maxMarkov = markovOrderMetricTest(j);
        end
        
        if oneOffsetDiagCheck(j) >= maxDiag
            maxDiag = oneOffsetDiagCheck(j);
        end
        
        if markovDiagProduct(j) >= maxMarkovDiagProduct
            maxMarkovDiagProduct = markovDiagProduct(j);
        end

        % if the result is the lowest markov stat yet, remember it 
        if markovOrderMetricTest(j) <= minMarkov
            minMarkov = markovOrderMetricTest(j);
        end

        % Output a progress indicator - useful for calculations with more permutations so that the user can see how far the calculation has got
        progressCount = progressCount + 1;
        if progressCount == progressIndicator
            fprintf('*');
            progressCount = 0;
        end
    end

    maxMarkovPerm = zeros(1,permsCount); % Record the position of each test section that gives a maximum mvalue
    minMarkovPerm = zeros(1,permsCount); % Record the position of each test section that gives a minimum mvalue
    maxDiagPerm = zeros(1,permsCount); % Record the position of each test section that gives a maximum diagonal value
    bestFaciesOrder = zeros(1, maxNumbOfFacies); % An array of the facies code permutations that give the maximum markov stat value
    worstFaciesOrder = zeros(1, maxNumbOfFacies); % An array of facies code permutations that give the minimum markov stat value
    bestDiagFaciesOrder = zeros(1, maxNumbOfFacies); % An array of facies code permutations that give the maximum one offset diagonal test value
    maxMarkovCount = 0; % The number of permutations that give the maximum markov value
    minMarkovCount = 0;  % The number of permutations that give the minimum markov value
    maxDiagCount = 0;  % The number of permutations that give the maximum one offset diagonal test value
    bestFaciesSect = zeros(1, sectionLength); % An array of the succession permutations that give the maximum markov stat value
    worstFaciesSect = zeros(1, sectionLength); % An array of the succession permutations that give the maximum markov stat value
    bestDiagFaciesSect = zeros(1, sectionLength); % An array of the succession permutations that give the maximum markov stat value

    % now loop again to find and record all of the permutations and resulting
    % strat sections that gave the maximum or minimum m and diag values found above
    % and record various information about each permutation
    for j = 1:permsCount

        % Convert to single precision for the >= test to remove possible
        % rounding error effects that can otherwise introduced small
        % differences in markovOrderMetricTest that are not "real"
        if single(markovOrderMetricTest(j)) >= single(maxMarkov)
            maxMarkovCount = maxMarkovCount + 1;
            maxMarkovPerm(maxMarkovCount) = j;
            bestFaciesOrder(maxMarkovCount, :) = allFaciesCombos(j,:); % record the facies code permutation

            % Make a vertical succession with each facies occurrence coded
            % according to the  facies coding from the jth permutation
            for k =1:sectionLength
                oneFacies = dataFaciesCodes(k);
                bestFaciesSect(maxMarkovCount, k) = allFaciesCombos(j, oneFacies); % Make a complete vertical succession with facies coding from current permutation
            end
        end

        if single(markovOrderMetricTest(j)) <= single(minMarkov)
            minMarkovCount = minMarkovCount + 1;
            minMarkovPerm(minMarkovCount) = j;
            worstFaciesOrder(minMarkovCount, :) = allFaciesCombos(j,:); % record the facies code permutation

            % Make a vertical succession with each facies occurrence coded
            % according to the  facies coding from the jth permutation
            for k =1:sectionLength
                oneFacies = dataFaciesCodes(k);
                worstFaciesSect(minMarkovCount, k) = allFaciesCombos(j, oneFacies); % Make a complete vertical succession with facies coding from current permutation
            end
        end
        
        % Select permuations with both a high m statistic and a high one-offset diagonal test value as indicated by a high markovDiagProduct
        if markovDiagProduct(j) >= maxMarkovDiagProduct - ROUND_ERROR
            maxDiagCount = maxDiagCount + 1;
            maxDiagPerm(maxDiagCount) = j;
            bestDiagFaciesOrder(maxDiagCount, :) = allFaciesCombos(j,:); % record the facies code permutation
            
            % Make a vertical succession with each facies occurrence coded
            % according to the  facies coding from the jth permutation
            for k =1:sectionLength
                oneFacies = dataFaciesCodes(k);
                bestDiagFaciesSect(maxDiagCount, k) = allFaciesCombos(j, oneFacies); % Make a complete vertical succession with facies coding from current permutation
            end
        end
    end

    % **************************************************************************************
    % Report all of the results as text output
    fullFileName = strcat(fileName, '_report.txt'); % Save the permutation results to .mat file in ascii format
    fOut = fopen(fullFileName,'w');
    fprintf('Largest markov stat value %5.4f found from %d permutations \n', maxMarkov, maxMarkovCount);
    fprintf(fOut,'Largest markov stat value %5.4f found from %d permutations \n', maxMarkov, maxMarkovCount);
    fprintf('Lowest markov stat value %5.4f found from %d permutations \n', minMarkov, minMarkovCount);
    fprintf(fOut,'Lowest markov stat value %5.4f found from %d permutations \n', minMarkov, minMarkovCount);
    fprintf('Highest m value * one-offset diagonal value %5.4f found from %d permutations \n', maxMarkovDiagProduct, maxDiagCount);
    fprintf(fOut,'Highest m value * one-offset diagonal value %5.4f found from %d permutations \n', maxMarkovDiagProduct, maxDiagCount);
    fprintf('Optimisation indicator: %5.4f From %d total permutations, %d permutations gave optimal score\n', maxDiagCount/permsCount, permsCount, maxDiagCount);
    fprintf(fOut,'Optimisation indicator: %5.4f From %d total permutations, %d permutations gave optimal score\n', maxDiagCount/permsCount, permsCount, maxDiagCount);

    k=1:maxNumbOfFacies; % NB implicit loop on k
    fprintf('Original facies: ');
    fprintf('%s ', faciesNames{k});
    fprintf('\nOptimal codings are:');
    fprintf(fOut,'Original facies: ');
    fprintf(fOut,'%s ', faciesNames{k});
    fprintf(fOut,'\nOptimal codings are:');
    for j=1:maxDiagCount
         fprintf('\n%d ',j);
         fprintf(fOut,'\n%d ',j); 
         fprintf('%s ', faciesNames{bestDiagFaciesOrder(j,k)}); % NB the implied loop on k in this - separate formatting fprints required above therefore
         fprintf(fOut,'%s ',faciesNames{bestDiagFaciesOrder(j,k)});
         fprintf('m=%5.4f diag=%5.4f m-diag product %5.4f', markovOrderMetricTest(maxDiagPerm(j)), oneOffsetDiagCheck(maxDiagPerm(j)), markovOrderMetricTest(maxDiagPerm(j)) * oneOffsetDiagCheck(maxDiagPerm(j)) );
         fprintf(fOut,'m=%5.4f diag=%5.4f m-diag product %5.4f', markovOrderMetricTest(maxDiagPerm(j)), oneOffsetDiagCheck(maxDiagPerm(j)), markovOrderMetricTest(maxDiagPerm(j)) * oneOffsetDiagCheck(maxDiagPerm(j)) );
     end

    fprintf('\n');
    fprintf(fOut,'\n');
    
    k=1:maxNumbOfFacies; % NB implicit loop on k
    fprintf('\nWorst codings are:');
    fprintf(fOut,'\nWorst codings are:');
    for j=1:minMarkovCount
        fprintf('\n%d ',j);
        fprintf(fOut,'\n%d ',j);
        fprintf('%s ', faciesNames{worstFaciesOrder(j,k)}); % NB the implied loop on k in this - separate formatting fprints required above therefore
        fprintf(fOut,'%s ', faciesNames{worstFaciesOrder(j,k)});
    end
    fprintf('\n');
    fprintf(fOut,'\n');
    fclose(fOut);

    fullFileName = strcat(fileName, '_report_perms.txt'); % Save the permutation results to .mat file in ascii format
    save(fullFileName, 'markovOrderMetricTest','-ascii');

    % **************************************************************************************
    % Plot the results graphically
    scrsz = get(0,'ScreenSize'); % screen dimensions vector
    f2 = figure('Visible','on','Position',[1 scrsz(4)/4 (scrsz(3)/3)* 2.0 (scrsz(4)/3)*2]);

    % Plot the original data vertical section
    h1 = subplot('Position',[0.04 0.1 0.075 0.85]);
    hold on
    cumThick = 0;
    for j=1:sectionLength
        fCode = dataFaciesCodes(j);
        yco = [cumThick, cumThick, cumThick + dataThick(j), cumThick + dataThick(j)];
        xco = [0, dataFaciesCodes(j), dataFaciesCodes(j), 0];
        faciesCol = [faciesColours(fCode,2) faciesColours(fCode,3) faciesColours(fCode,4)];
        patch(xco, yco, faciesCol,'EdgeColor','none');
        cumThick = cumThick + dataThick(j);
    end
    grid on;
    set(gca,'Layer','top');
    xlabel('Facies code');
    ylabel('Thickness (m)');
    set(gca,'XTick', 1:(maxNumbOfFacies),'XTickLabel', faciesNames, 'TickDir', 'out', 'XTickLabelRotation',90);

    % Plot the TP matrix for the orginal facies coding
    h2=subplot('Position',[0.2 0.55 0.25 0.4]);
    hold on
    [mStat,diag] = calculateTPMatrixAndOrderMetric(dataFaciesCodes, 1);
    axis tight;
    h2.XTickLabel = {faciesNames{1:maxNumbOfFacies}};
    h2.YTickLabel = {faciesNames{1:maxNumbOfFacies}};
    str = sprintf('TP Matrix, original facies coding, m=%5.4f diag=%4.3f', mStat, diag);
    title(str);

    % Plot the orginal facies coding order
    h3=subplot('Position',[0.47 0.55 0.05 0.4]);
    for k=1:maxNumbOfFacies
        yco = [k k k+1 k+1];
        xco = [1 1.5 1.5 1];
        faciesCol = [faciesColours(k,2) faciesColours(k,3) faciesColours(k,4)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
        labelStr = faciesNames(k);
        text(1.25, k+0.2, labelStr, 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
        axis off;
    end

    % plot the TP matrix for the first of the recorded best m value facies permutations
    if maxDiagCount > 0
        h4=subplot('Position',[0.6 0.55 0.25 0.4]);
        hold on
        [mStat,diag] = calculateTPMatrixAndOrderMetric(bestDiagFaciesSect(1,:), 1);
        axis tight;
        labelStr = cell(maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
        colourCode = zeros(maxNumbOfFacies, 3);
        for k=1:maxNumbOfFacies
            % For each permutation, one of which is recorded in each bestDiagFaciesOrder,
            % the array elements are in facies order, so element k,1 is for facies 1 
            % and the value indicates which row in the TP matrix the facies should go for this permutation
            % So facies 1 should be labelled into the label string element indicated by elmement 1 in bestDiagFaciesOrder
            labelStr{bestDiagFaciesOrder(1,k)} = faciesNames{k}; % Put the correct facies code into label array element k
            colourCode(bestDiagFaciesOrder(1,k),1:3) = faciesColours(k,2:4);
        end

        h4.XTickLabel = labelStr;
        h4.YTickLabel = labelStr;
        str = sprintf('TP Matrix, optimal facies coding, m=%4.3f diag=%4.3f m-diag product %4.3f', mStat, diag, mStat * diag);
        title(str);

        % plot the ideal section for the best m value facies combo section
        h5=subplot('Position',[0.87 0.55 0.05 0.4]);
        for k=1:maxNumbOfFacies
            yco = [k k k+1 k+1];
            xco = [8 8.5 8.5 8];
            faciesCol = [colourCode(k,1) colourCode(k,2) colourCode(k,3)];
            patch(xco, yco, faciesCol,'EdgeColor','black');
            text(8.25, k+0.2, labelStr(k), 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
            axis off;
        end
    end

     % plot the histogram of m values from all permutations
    h6=subplot('Position',[0.2 0.1 0.3 0.35]);
    histogram(markovOrderMetricTest, 'EdgeColor','none', 'FaceColor','blue');
    xlabel('m value', 'FontSize',10);
    ylabel('Frequency', 'FontSize',10);
    grid on;
    
    % plot the histogram of m scores from the n maximum score permutations
    h7=subplot('Position',[0.6 0.1 0.3 0.35]);
    histogram(markovDiagProduct, 'EdgeColor','none', 'FaceColor','blue');
    xlabel('d', 'FontSize',10);
    ylabel('Frequency', 'FontSize',10);
    grid on;
    
    

    % plot a TP matrix for all the best m value facies combo permutations
    for j=1:maxDiagCount

        h3=figure;
        h4=subplot('Position',[0.1 0.1 0.60 0.85]);
        h5=subplot('Position',[0.8 0.1 0.15 0.85]);
        hold on;

        % Calculate and draw the TP matrix
        subplot(h4);
        [mStat,diag] = calculateTPMatrixAndOrderMetric(bestDiagFaciesSect(j,:), 1);
        axis tight;
        labelStr = cell(1,maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
        for k=1:maxNumbOfFacies
                labelStr{bestDiagFaciesOrder(j,k)} = faciesNames{k}; % Because bestFaciesOrder is not in the same order as faciesNames
                colourCode(bestDiagFaciesOrder(j,k),1:3) = faciesColours(k,2:4);
        end
        h4.XTickLabel = labelStr;
        h4.YTickLabel = labelStr;
        str = sprintf('m=%5.4f diag =%3.1f m*diag %4.3f Optimal facies coding %d', mStat, diag, mStat * diag, j );
        title(str);
        
        % Draw a vertical succession for this ideal section permutation
        subplot(h5);
        for k=1:maxNumbOfFacies
            yco = [k k k+1 k+1];
            xco = [1 2 2 1];
            faciesCol = [colourCode(k,1) colourCode(k,2) colourCode(k,3)];
            patch(xco, yco, faciesCol,'EdgeColor','black');
            text(1.3, k+0.2, labelStr(k), 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',12);
            axis off;
        end

    end

    % Repeat one final time for one of the worst facies combination cases
    h3=figure;
    h4=subplot('Position',[0.1 0.1 0.60 0.85]);
    h5=subplot('Position',[0.8 0.1 0.15 0.85]);
    hold on;
    
    % Calculate and draw the TP matrix for the worst permutation
    subplot(h4);
    [mStat,diag] = calculateTPMatrixAndOrderMetric(worstFaciesSect(1,:), 1);
    axis tight;
    labelStr = cell(1,maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
    for k=1:maxNumbOfFacies   
        labelStr{worstFaciesOrder(1,k)} = faciesNames{k}; % Put the correct facies code into label array element k
        colourCode(worstFaciesOrder(1,k),1:3) = faciesColours(k,2:4);
    end
    h4.XTickLabel = labelStr;
    h4.YTickLabel = labelStr;
    str = sprintf('TP Matrix, m=%5.4f Diag=%3.1f Worst facies coding 1:', mStat, diag);
    title(str);
    
    % Draw a vertical section for this worst section permutation
    subplot(h5);
    for k=1:maxNumbOfFacies
        yco = [k k k+1 k+1];
        xco = [1 2 2 1];
        faciesCol = [colourCode(k,1) colourCode(k,2) colourCode(k,3)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
        text(1.3, k+0.2, labelStr(k), 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',12);
        axis off;
    end
    
    seconds = toc;
    fprintf('\nDone after %f2.1 seconds\n', seconds);
end

function [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(faciesSect, plotFlag)

    % find the number of elements in the succession. NB size returns both
    % matrix dimensions so max ensures nz is the biggest which should be
    % the length of the section
    nz = max(size(faciesSect));
    
    % Find the maximum facies code used in the facies succession - this is the
    % size for both dimensions of the TP matrix which can now be defined
    nFacies = max(faciesSect);
    
    if nFacies == 0
        nFacies = 0;
    end
    
    m = 0; % number of transitions, not really needed but never mind
    TFMatrix = zeros(nFacies, nFacies);
    TPMatrix = zeros(nFacies, nFacies);
    TPMultiplierMatrix = zeros(nFacies, nFacies);
    
    % populate the multiplier matrix with value 10 on the 1-offset diagonal
    % - where the max p values will be for ABCDE strata
    for i = 1:nFacies-1
        TPMultiplierMatrix(i,i+1) = 10.0;
    end
    TPMultiplierMatrix(nFacies,1) = 10.0; % Complete the diagonal with the "wraparound" position in the matrix corner
    
    % Now loop through the elements in the succession and for each different facies from-to transition,
    % increment the appropriate cell in the matrix
    for i =1 : nz-1
        fromFacies = faciesSect(i);
        toFacies = faciesSect(i+1);
        % mark transitions between different facies
        if fromFacies > 0 && toFacies > 0 && fromFacies ~= toFacies % Make sure facies codes are not zero because zero values would record an error
            TFMatrix(fromFacies, toFacies) = TFMatrix(fromFacies, toFacies) + 1; % increment the appropriate value in the tp matrix
            m = m + 1;
        end     
    end

    % Now calculate the transition probability matrix from the transition frequency matrix
    rowSums=sum(TFMatrix,2); % Calculates the sum of each row in TF matrix and stores as vector rowSums
    for k=1:nFacies
        for j=1:nFacies
            if rowSums(k) > 0 % if rowsum > 0 divide TF value by row sum to get transition probability
                TPMatrix(k,j)=TFMatrix(k,j) / rowSums(k);
            else
                TPMatrix(k,j) = 0;
            end
        end
    end
    
    % Now calculate the Markov order metrics

    % Calculate the metric for the maximum diagonal
    diagMetric = zeros(1,nFacies-1);
    diagMetricMultiplied = zeros(1,nFacies-1);
    TPMatrixMultiplied = TPMultiplierMatrix .* TPMatrix;
    
    % Now loop through each offset diagonal in the matrix and calculate the
    % mean of the cell values in the diagonal
    % Also do this for the version of the TP matrix with a *10 weighted one-offset diagonal
    for j=1:nFacies-1;
    
        diagMetric(j) = (sum(diag(TPMatrix,j)) + sum(diag(TPMatrix,-(nFacies-j))) )/ nFacies; % calculate a mean for the jth FULL diagonal in the matrix
        diagMetricMultiplied(j) = (sum(diag(TPMatrixMultiplied,j)) + sum(diag(TPMatrixMultiplied,-(nFacies-j))) )/ nFacies; % Also calculate a mean based on the weighted one-offset diagonal
    end
    
    markovOrderMetric = max(diagMetric)- min(diagMetric);
    oneOffsetMValueDiagCheck = max(diagMetricMultiplied)- min(diagMetricMultiplied);

    if plotFlag == 1

        a1 = gca;
        TPCellSize = 1.0;
        matrixTopYco = (nFacies + 1) * TPCellSize;
        matrixBottYco = 0.0;
        
        for i = 1:nFacies
            for j = 1:nFacies
                yco = [matrixBottYco+((j-0.5)*TPCellSize) matrixBottYco+((j-0.5)*TPCellSize) matrixBottYco+((j+0.5)*TPCellSize) matrixBottYco+((j+0.5)*TPCellSize)];
                xco = [(i-0.5)*TPCellSize (i+0.5)*TPCellSize (i+0.5)*TPCellSize (i-0.5)*TPCellSize];
                labelStr = sprintf('%3.2f',TPMatrix(j,i));
                
                if TPMatrix(j,i) <= 0.5
                    colVect = [1 TPMatrix(j,i)/0.5 TPMatrix(j,i)/1.1111111]; % gradation of colours from red (p=0) to orange-yellow (p=0.5)
                else
                    colVect = [1-(TPMatrix(j,i)) 1 1-(TPMatrix(j,i)/1.111111)]; % gradation of colours from orange-yellow (p=0.5) to green (p=1)
                end
               
                if i==j
                    colVect = [0.9 0.9 0.9]; % Matrix diagonal is a special case so colour light grey
                end
                
                patch(xco, yco, colVect);
                
                if i ~= j % So do not add text labels to the diagonal
                  text(double(((i-0.7)*TPCellSize))+(TPCellSize*0.33), double(matrixBottYco+(j*TPCellSize)-0.3), labelStr,'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
                end
            end
        end
        
        set(a1,'XTick',1:nFacies);
        set(a1,'YTick',1:nFacies);
       
        xlabel('Facies Code: To', 'FontSize',10);
        ylabel('Facies code: From', 'FontSize',10);
        
    end
end





