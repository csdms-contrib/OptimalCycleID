function strataWorkbench
% Code for StrataWorkbench
% Written by P Burgess, University of Liverpool
% Started 2/4/2016, incoporating various older code into a new GUI and creating basis for developing new functionality

    clear all;
    
    gui.main = 0;
    gui.f1 = 0;
    

    data.loaded = 0;    % Boolean flag to show if data loaded or not; if not, subsequent functions disabled
    data.faciesCodes = zeros(1,1);
    data.faciesThick = zeros(1,1);
    data.sectLength = 0.0;
    data.maxNumbOfFacies = 0;
    data.sectionMean = 0.0;
    
    results.markovOrderMetric = zeros(1,1);

    initializeGUI(gui, data, results);
end

function initializeGUI(gui, data, results)

    iteration = 0;
    dummyMap = zeros(1,1);
    
    %  Create and then hide the GUI window as it is being constructed.
    
    % ScreenSize is a four-element vector: [left, bottom, width, height]:
    scrsz = get(0,'ScreenSize'); % vector 
    scrWidthProportion = 0.75;
    scrHeightIncrement = scrsz(4)/20; % Use this to space controls down right side of the main window
    controlStartY = (scrsz(4) * 0.8) - (scrHeightIncrement / 2);
    controlStartX = (scrsz(3) * scrWidthProportion) - 420;
    
    % position requires left bottom width height values. screensize vector
    % is in this format 1=left 2=bottom 3=width 4=height
    gui.main = figure('Visible','off','Position',[1 scrsz(4)*scrWidthProportion scrsz(3)*scrWidthProportion scrsz(4)*0.8]);
   
    
   %  Construct the control panel components.
  
   hSectionFpathLabel = uicontrol('style','text','string','Section data file path:','Position',[controlStartX+40, controlStartY-scrHeightIncrement, 200, 15]);
   hSectionFpath = uicontrol('Style','edit','String','sectionData/','Position',[controlStartX+200, controlStartY-scrHeightIncrement, 200, 25]);
   hSectionFnameLabel = uicontrol('style','text','string','Section data filename:','Position',[controlStartX+40, (controlStartY-scrHeightIncrement*2), 200, 15]);
   hSectionFname = uicontrol('Style','edit','String','sectSwaps20recoded.txt','Position',[controlStartX+200, (controlStartY-scrHeightIncrement*2), 200, 25]);
   hSectionColourMapFnameLabel = uicontrol('style','text','string','Section colour codes filename:','Position',[controlStartX+20, controlStartY-(scrHeightIncrement*3), 200, 15]);
   hSectionColourMapFname = uicontrol('Style','edit','String','sectSwaps20recodedLithoCol.txt','Position',[controlStartX+200, controlStartY-(scrHeightIncrement*3), 200, 25]);
   
   hLoadSectionData = uicontrol('Style','pushbutton','String','Load and plot section data',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*4),200,25],...
          'BackgroundColor',[0.6 1.0 0.6],...
          'Callback',{@loadSectionButton_callback});
      
    hPlotTPMatrix = uicontrol('Style','pushbutton','String','Plot TP Matrix for observed strata',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*5),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@plotTPMatrix_callback});
      
    hPlotFaciesFrequencyHistogram = uicontrol('Style','pushbutton','String','Plot facies frequency histogram',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*6),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@plotFaciesFrequencyHistogram_callback});
      
    hCalculateAndPlotRandomModel = uicontrol('Style','pushbutton','String','Calculate and plot random models',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*7),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@CalculateAndPlotRandomModel_callback});
      
    hCalculateAndPlotIdealCycles = uicontrol('Style','pushbutton','String','Calculate and plot ideal cycles',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*8),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@CalculateAndPlotOptimisedCycles_callback});
      
    hCalculateAndPlotIdealCycles = uicontrol('Style','pushbutton','String','Close any windows and clear data',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*15),200,25],...
          'BackgroundColor',[1.0 0.4 0.4],...
          'Callback',{@resetButton_callback});
    
      
   % Assign the GUI a name to appear in the window title.
   set(gui.main,'Name','Strata Workbench')
   % Move the GUI to the center of the screen.
   movegui(gui.main,'center')
   % Make the GUI visible.
   set(gui.main,'Visible','on');

    function loadSectionButton_callback(source,eventdata)
        fprintf('\n\n===================================================== Loading Data =====================================================\n');
        
        dataFileNameAndPath = strcat(get(hSectionFpath,'String'),get(hSectionFname,'String'));
        coloursFileNameAndPath = strcat(get(hSectionFpath,'String'),get(hSectionColourMapFname,'String'));

        [data, success] = loadSectionData(data, dataFileNameAndPath, coloursFileNameAndPath);
        if success == 1
            gui = plotObservedSection(gui, data);
            data.loaded = 1; % Set flag to enable rest of the button functions
        end
    end

    function plotTPMatrix_callback(source,eventdata)
        
        if data.loaded
            fprintf('\n\n===================================================== Plot TP Matrix for Observed Strata =====================================================\n');
            gui.sp2 = subplot('Position',[0.2 0.55 0.3 0.4]);
            [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(data.faciesCodes, data.faciesNames, 1, gui.sp2); % 1 is plot flag, so plot the results

            fprintf('Section  gives m statistic %5.4f and one-offset diagonal check %5.4f\n', markovOrderMetric, oneOffsetMValueDiagCheck);
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function plotFaciesFrequencyHistogram_callback(source,eventdata)
        
        if data.loaded
            fprintf('\n\n===================================================== Plot Facies Frequency Histogram for Observed Strata =====================================================\n');
            gui = plotFaciesFrequencyHistogram(gui, data);
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function CalculateAndPlotRandomModel_callback(source, eventsdata)
        
        if data.loaded
            fprintf('\n\n===================================================== Calculate Random Models and Compare =====================================================\n');
            gui = calculateAndPlotRandomModels(gui, data); % 1 is plot flag, so plot the results
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function CalculateAndPlotOptimisedCycles_callback(source, eventsdata)
        
        
        if data.loaded
            fprintf('\n\n===================================================== Calculate Optimal Cycles =====================================================\n');
            gui = calculateAndPlotOptimisedCycles(gui, data); % 1 is plot flag, so plot the results
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function resetButton_callback(source, eventdata) 

        if isfield(gui,'f2')
           close(gui.f2);
           gui = rmfield(gui,'f2');
        end

        if isfield(gui,'f3')
           close(gui.f3);
           gui = rmfield(gui,'f3');
        end

        if isfield(gui,'sp1')
           cla(gui.sp1);
        end

        if isfield(gui,'sp2')
           cla(gui.sp2);
        end

        if isfield(gui,'sp3')
           cla(gui.sp3);
        end
        
        data.loaded = 0;    % Boolean flag to show if data loaded or not; if not, subsequent functions disabled
        data.faciesCodes = zeros(1,1);
        data.faciesThick = zeros(1,1);
        data.sectLength = 0.0;
        data.maxNumbOfFacies = 0;
        data.sectionMean = 0.0;
    end
end

function [data, success] = loadSectionData(data, dataFileNameAndPath, coloursFileNameAndPath)
    
    % if we are loading a new section, need to clear out all the old
    % data, but cant just use clear because this gets rid of all arrays, so
    % need to think of a way to do this here...

    data.faciesCodes = zeros(1,1);
    data.thick = zeros(1,1);
    data.sectLength = 0.0;
    j = uint32(0);
    success = 1; % Assume loading will go ok, and reset to false below if not...

    % Read section data, thickness and facies
    if exist(dataFileNameAndPath, 'file')
        dataIn = load(dataFileNameAndPath);
        data.faciesThick = dataIn(:,1);
        data.faciesCodes = dataIn(:,2);
    else
        message = sprintf('Succession data file %s does not exist. Please correct the filename and try again.\n', dataFileNameAndPath);
        m1 = msgbox(message,'Data not loaded');
        success = 0;
        return;
    end
    
    % Read colour map file, one colour and one text name per facies
    if exist(coloursFileNameAndPath, 'file');
        fid = fopen(coloursFileNameAndPath);
        while ~feof(fid)
            file = textscan(fid, '%u8 %f %f %f %s');
        end
        fclose(fid);
        
    else
        message = sprintf('Facies colour map and name file %s does not exist. Please correct the filename and try again.\n', coloursFileNameAndPath);
        m1 = msgbox(message,'Data not loaded');
        success = 0;
        return;
    end
    
    % Extract facies names, numeric codes and RGB colours from data input file as a cell array
    data.faciesNames = file {5}(); % Because the 5th element in each row of the litho col file should be the facies name
    data.faciesNumbers = file {1}; % And the original faceis code is in 1st elemnt in each row - particularly important for optimal cycle analysis
    data.faciesColours = zeros(length(data.faciesNumbers),4); % Dimension colours matrix based on number of facies read from file
    for k=1:4
        data.faciesColours(:,k)= double(file{k}());
    end

    % Calculate the basic stats on the dataFacies array
    data.sectLength = max(size(data.faciesCodes));
    data.maxNumbOfFacies = length(data.faciesColours);
    data.sectionMean = mean(data.faciesThick);

    % Make a reference section that uses the facies names to ID each unit
    % rather than the numerical facies code. Need to remember that faciesNames
    % is a cell array hence the curly brackets
    % Note that as of 12.7.2015 this array is not used elsewhere in the code
    for k=1:data.sectLength
        sectionFaciesNames{k} = data.faciesNames{data.faciesCodes(k)};
    end

    fprintf('Loaded vertical section data from %s\nFor %d units, total %d facies, mean unit thickness %4.3f m\n', dataFileNameAndPath, data.sectLength, data.maxNumbOfFacies, data.sectionMean);
end

function gui = plotObservedSection(gui,data)

     % Plot the original data vertical section
    gui.sp1 = subplot('Position',[0.04 0.1 0.075 0.85]);
    cla(gui.sp1); % Clears the axis of whatever has been previously plotted
    
    hold on
    cumThick = 0;
    for j=1:data.sectLength
        fCode = data.faciesCodes(j);
        yco = [cumThick, cumThick, cumThick + data.faciesThick(j), cumThick + data.faciesThick(j)];
        xco = [0, data.faciesCodes(j), data.faciesCodes(j), 0];
        faciesCol = [data.faciesColours(fCode,2) data.faciesColours(fCode,3) data.faciesColours(fCode,4)];
        patch(xco, yco, faciesCol,'EdgeColor','none');
        cumThick = cumThick + data.faciesThick(j);
    end
    grid on;
    set(gca,'Layer','top');
    xlabel('Facies code');
    ylabel('Thickness (m)');
    set(gca,'XTick', 1:(data.maxNumbOfFacies),'XTickLabel', data.faciesNames, 'TickDir', 'out', 'XTickLabelRotation',90);
    title('Observed section');
    hold off;
end

function [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(faciesCodes, faciesNames, plotFlag, g1)

    % find the number of elements in the succession. NB size returns both
    % matrix dimensions so max ensures nz is the biggest which should be
    % the length of the section
    nz = max(size(faciesCodes));
    
    % Find the maximum facies code used in the facies succession - this is the
    % size for both dimensions of the TP matrix which can now be defined
    nFacies = max(faciesCodes);
    
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
        fromFacies = faciesCodes(i);
        toFacies = faciesCodes(i+1);
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
        
        labelStr = cell(1,nFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
        for k=1:nFacies
            labelStr{k} = faciesNames{k}; % Because bestFaciesOrder is not in the same order as faciesNames
        end
        g1.XTickLabel = labelStr;
        g1.YTickLabel = labelStr;
        
        set(g1,'XTick',1:nFacies);
        set(g1,'YTick',1:nFacies);
       
        xlabel('Facies Code: To', 'FontSize',10);
        ylabel('Facies code: From', 'FontSize',10);
        axis tight;
    end
end

function gui = plotFaciesFrequencyHistogram(gui, data)
% Plot the histogram of the facies frequencies

    gui.sp3 = subplot('Position',[0.2 0.075 0.3 0.4]);
    cla(gui.sp3);
    hold on
    
    frequencyDataFacies=histc(data.faciesCodes,1:data.maxNumbOfFacies);
    for j=1:data.maxNumbOfFacies
        xco = [j-1, j-1, j, j];
        yco = [0, frequencyDataFacies(j), frequencyDataFacies(j), 0];
        faciesCol = [data.faciesColours(j,2) data.faciesColours(j,3) data.faciesColours(j,4)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
    end

    set(gca,'XTick', 0.5:(data.maxNumbOfFacies-0.5),'XTickLabel', data.faciesNames);
    xlabel('Lithofacies');
    ylabel('Frequency');
    
    hold off;
    
    % Check is all facies are present in the section, and if not, issue a warning to recode
    zeroFreqWarningFlag = 0;
    for i=1:data.maxNumbOfFacies;
        if frequencyDataFacies(i) == 0
             zeroFreqWarningFlag = 1;
        end
    end
    if zeroFreqWarningFlag
        message = sprintf('At least one facies code has no occurences in this section\nConsider re-coding facies with concurruent numbering to avoid misleading results.');
        m1 = msgbox(message,'Missing Facies');
    end 
end

function gui = calculateAndPlotRandomModels(gui, data)

    TRUE = uint8(1);
    FALSE = uint8(0);
    maxIterations = 5000;
    numberOfSwaps = data.sectLength;
    maxRun = 3;
    minRun = 0;
    runBinIncrement = 0.05;
    runRange = maxRun - minRun;

    % Calculate and output the order metric for the entire data succession
    [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(data.faciesCodes, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01]));  % note the zero is the TP matrix plot flag, set to false to not trigger a plot
    runsOrderMetric = calculateRunsOrderMetric(data.faciesThick);
    fprintf('Markov metric for strata is %4.3f\nRuns analysis metric for strata is %5.4f\n', markovOrderMetric, runsOrderMetric);

    % Now calculate the metrics for many iterations of a random model

    % Shuffle the observed section and calculate the order metric and DW stat each time
    for j = 1:maxIterations;
        [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(data.faciesCodes, data.faciesThick, numberOfSwaps);
        multiMarkovOrderMetricDataShuffled(j) = calculateTPMatrixAndOrderMetric(shuffledFacies, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01])); % note the zero is the TP matrix plot flag, set to false so no plot
        multiRunsOrderMetricDataShuffled(j) = calculateRunsOrderMetric(shuffledThick); 
    end

    % Stats on the shuffled section random model
    meanMultiMarkovDataShuffled = mean(multiMarkovOrderMetricDataShuffled);
    stdDevMultiMarkovDataShuffled = std(multiMarkovOrderMetricDataShuffled);
    bins = 0:0.02:1.00; % because 0<m<=1
    multiMarkovOrderMetricDataShuffledHistData=histc(multiMarkovOrderMetricDataShuffled, bins) / maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
    markovPValueSum = sum(multiMarkovOrderMetricDataShuffledHistData(int16(markovOrderMetric*50:length(bins)))); % area under curve from m to max m value 1

    meanMultiRunsDataShuffled = mean(multiRunsOrderMetricDataShuffled);
    stdDevMultiRunsDataShuffled = std(multiRunsOrderMetricDataShuffled);
    bins = 0:0.05:runRange; % 0 is the minimum run metric, 3 a generally maximum value (this is what runRange should be set to)
    multiRunsOrderMetricDataShuffledHistData = histc(multiRunsOrderMetricDataShuffled, bins)/maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
    runBinIndex = 1 + int16(((runsOrderMetric-minRun)/runRange)*(runRange/runBinIncrement)); % Position of runs stat in the histogram
    runsPValueSum = sum(multiRunsOrderMetricDataShuffledHistData(runBinIndex:length(multiRunsOrderMetricDataShuffledHistData))); % area under curve from r to max run value

    fprintf('For %d iterations of a SHUFFLED DATA model\n', maxIterations);
    fprintf('Markov stats mean %5.4f std dev %5.4f Markov order metric P Value %5.4f\n', meanMultiMarkovDataShuffled, stdDevMultiMarkovDataShuffled, markovPValueSum);
    fprintf('Runs stats mean %5.4f std dev %5.4f Runs analysis metric P Value %5.4f\n', meanMultiRunsDataShuffled, stdDevMultiRunsDataShuffled, runsPValueSum);

    % Plot the results
    scrsz = get(0,'ScreenSize'); % screen dimensions vector

    % Plot the original vertical section, alongside the facies and thickness
    % order plots showing the outcrop datapoint and the multiple iteration shuffled section
    % metric frequency distribution
    gui.f2 = figure('Visible','on','Position',[1 scrsz(4)/4 (scrsz(3)/2) (scrsz(4)/3)*2]);
    
    subplot('Position',[0.075 0.1 0.075 0.85]);
    plotShuffledSection(data);

    % Subplot for the Markov order analysis histogram   
    subplot('Position',[0.25 0.57 0.65 0.38]);
    hold on
    bins = 0:0.02:1.00; % Make sure bins is set correctly for Markov plots
    bar(bins, multiMarkovOrderMetricDataShuffledHistData, 'EdgeColor','none', 'BarWidth', 1, 'FaceColor',[0.2 0.4 0.7]); % Colour is dark slate blue
    maxFreq = max(multiMarkovOrderMetricDataShuffledHistData) *1.1; % This is needed to scale the plot
    x = [markovOrderMetric markovOrderMetric];
    y = [0 maxFreq];
    line(x,y, 'color', [0.80 0.00 0.00], 'linewidth', 3.0); % Colour is dark red
    grid on;
    axis([0 1 0 Inf]);
    set(gca,'Layer','top');
    xlabel('Markov Order Metric for Facies');
    ylabel('Relative Frequency');

    % Subplot for the runs analysis histogram    
    subplot('Position',[0.25 0.1 0.65 0.38]);
    hold on
    bins = minRun:runBinIncrement:maxRun; % Make sure that bins is set correctly for Runs analysis plots
    bar(bins, multiRunsOrderMetricDataShuffledHistData, 'EdgeColor','none', 'BarWidth', 1, 'FaceColor',[0.2 0.4 0.7]);
    maxFreq = max(multiRunsOrderMetricDataShuffledHistData) * 1.1; % This is needed to scale the plot

    lineColor = [0.80 0.00 0.00];
    x = [runsOrderMetric runsOrderMetric];
    y = [0 maxFreq]; % Draw the data line from y=0 to y=max frequency of the three histograms
    line(x,y, 'color', lineColor, 'linewidth', 3.0);

    grid on;
    axis([0 Inf 0 Inf]);
    set(gca,'Layer','top');
    xlabel('Runs Analysis Order Metric for Thickness ');
    ylabel('Relative Frequency');
    
    % Complete analysis with a textbox message summary
    message = sprintf('Markov metric %4.3f\nRuns analysis metric %5.4f\nMarkov order metric P Value %5.4f\nRuns analysis metric P Value %5.4f\n', markovOrderMetric, runsOrderMetric, markovPValueSum, runsPValueSum);
    m1 = msgbox(message,'Random model comparison');
end

% function plotShuffledSection(gui, data)
function plotShuffledSection(data)

    % Plot a shuffled iteration of the observed vertical section
    numberOfSwaps = data.sectLength;
    [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(data.faciesCodes, data.faciesThick, numberOfSwaps);
    cumThick = 0;
    for j=1:data.sectLength
        fCode = shuffledFacies(j);
        yco = [cumThick, cumThick, cumThick + shuffledThick(j), cumThick + shuffledThick(j)];
        xco = [0, shuffledFacies(j), shuffledFacies(j), 0];
        faciesCol = [data.faciesColours(fCode,2) data.faciesColours(fCode,3) data.faciesColours(fCode,4)];
        patch(xco, yco, faciesCol,'EdgeColor','none');
        cumThick = cumThick + shuffledThick(j);
    end
    grid on;
    set(gca,'Layer','top');
    xlabel('Facies code');
    ylabel('Thickness (m)');
    set(gca,'XTick', 1:(data.maxNumbOfFacies),'XTickLabel', data.faciesNames, 'TickDir', 'out', 'XTickLabelRotation',90);
    title('Shuffled section');
end

function [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(sectFacies, sectThick, totalSwaps)
% function to shuffle the facies succession to ensure a random configuration

    % Make copies of the original data in new arrays that will be used to
    % store the shuffled sections
    shuffledFacies = sectFacies;
    shuffledThick = sectThick;
    n = uint16(max(size(shuffledFacies)));
    j = 0;
    infiniteStopCount = 0;
    
    while j < totalSwaps && infiniteStopCount < 1000000
        
        % Select two unit numbers randomly to be swapped
        unit1 = uint16((rand * (n-1)) + 1);
        unit2 = uint16((rand * (n-1)) + 1);
        
        % Need to check above and below for both positions that swapping will not put same
        % facies adjacent to one another and cause a transition to self
        swapFacies1 = shuffledFacies(unit1);
        if unit1 > 1 swapFacies1Below = shuffledFacies(unit1-1); else swapFacies1Below = 0;end
        if unit1 < n swapFacies1Above = shuffledFacies(unit1+1); else swapFacies1Above = 0;end
        
        swapFacies2 = shuffledFacies(unit2);
        if unit2 > 1 swapFacies2Below = shuffledFacies(unit2-1); else swapFacies2Below = 0;end
        if unit2 < n swapFacies2Above = shuffledFacies(unit2+1); else swapFacies2Above = 0;end
        
        % So compare facies in their new positions with the facies above and below and
        % only swap and increment loop counter if NOT the same...
        if swapFacies1Below ~= swapFacies2 && swapFacies1Above ~= swapFacies2 && swapFacies2Below ~= swapFacies1 && swapFacies2Above ~= swapFacies1
            
            %Swap the facies
            temp = shuffledFacies(unit1);
            shuffledFacies(unit1) = shuffledFacies(unit2);
            shuffledFacies(unit2) = temp;

            %Swap the thicknesses
            temp = shuffledThick(unit1);
            shuffledThick(unit1) = shuffledThick(unit2);
            shuffledThick(unit2) = temp;

            j = j + 1;
        end
        
        infiniteStopCount = infiniteStopCount + 1;
    end
end


function runsOrderMetric = calculateRunsOrderMetric(thicknesses)

    % find the number of units in the succession and declare arrays accordingly
    nz = max(size(thicknesses));
    deltaThick = zeros(1,nz);
    runsUp = zeros(1,nz);
    runsDown = zeros(1,nz);

    % Calculate the change in thickness between successive units
    i = 1:nz-1;
    j =2:nz; % so j = i + 1 therefore thickness change is thickness(j) - thickness(i)
    deltaThick(i) = thicknesses(j) - thicknesses(i);

    if deltaThick(1) > 0 runsUp(1) = 1; end
    if deltaThick(1) < 0 runsDown(1) = 1; end

    for i=2:nz
        if deltaThick(i) > 0 runsUp(i) = runsUp(i-1)+1; end
        if deltaThick(i) < 0 runsDown(i) = runsDown(i-1)+1; end
    end

    runsUpNormSum = (sum(runsUp)/nz);
    runsDownNormSum = (sum(runsDown)/nz);
    runsOrderMetric = (runsUpNormSum + runsDownNormSum);
end

function gui = calculateAndPlotOptimisedCycles(gui, data) % 1 is plot flag, so plot the results

    TRUE = uint8(1);
    FALSE = uint8(0);
    ROUND_ERROR = 0.0000001;

    tic

    % Calculate the basic stats on the dataFacies array
    sectionLength = max(size(data.faciesCodes));
    maxNumbOfFacies = length(data.faciesColours);
    sectionMean = mean(data.faciesThick);

    % Make a reference section that uses the facies names to ID each unit
    % rather than the numerical facies code. Need to remember that faciesNames
    % is a cell array hence the curly brackets
    % Note that as of 12.7.2015 this array is not used elsewhere in the code
    for k=1:sectionLength
        sectionFaciesNames{k} = data.faciesNames{data.faciesCodes(k)};
    end

    fprintf('For %d units, total %d facies, mean unit thickness %4.3f m\n', sectionLength, maxNumbOfFacies, sectionMean);

    % Calculate and output the order metric for the entire data succession
    markovOrderMetric = calculateTPMatrixAndOrderMetric(data.faciesCodes, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01]));  % note the zero is the TP matrix plot flag, set to false to not trigger a plot
    fprintf('Markov metric for strata is %5.4f\n', markovOrderMetric);

    allFaciesCombos = perms(data.faciesNumbers);
    permsCount = length(allFaciesCombos);

    % allFaciesCombos now contains all the possible ways the n facies classes read in can be positioned in the TP matrix
    % In other words, a permutation 5 1 2 3 4 means that facies 1 is coded as facies 5 so that it positioned in row 5 of the TP matrix for this permutation
    % Now make a strat section in testFaciesSect using the facies coding in each row of the matrix and calculate the markov statistic for each testFaciesSect
    % Define the necessary variables to make a test section and record markov
    % results from each test section. 
    testFaciesSect = zeros(1, sectionLength); % the facies section input but will contain facies coding from a single permutation
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
    
    % Loop through all of the facies coding permutations and find the
    % minimum, maximum m statistics, plus the max m-diagonal product
    for j = 1:permsCount

        % Make a vertical succession with each facies occurrence coded
        % according to the facies coding in permutation j, not the original coding
        for k =1:sectionLength
            oneFacies = data.faciesCodes(k);
            testFaciesSect(k) = allFaciesCombos(j, oneFacies); % Swap the original facies code from the data section with a substitute code from permutation j
        end

        % Calculate the markov statistic and the one offset diagonal for this permutation
        [markovOrderMetricTest(j), oneOffsetDiagCheck(j)] = calculateTPMatrixAndOrderMetric(testFaciesSect, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01]));
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

    maxMarkovPerm = zeros(1,permsCount); % Records the array index position of each test section that gives a maximum mvalue
    minMarkovPerm = zeros(1,permsCount); % Records the array index position of each test section that gives a minimum mvalue
    maxDiagPerm = zeros(1,permsCount); % Records the array index position of each test section that gives a maximum diagonal value
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
                oneFacies = data.faciesCodes(k);
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
                oneFacies = data.faciesCodes(k);
                worstFaciesSect(minMarkovCount, k) = allFaciesCombos(j, oneFacies); % Make a complete vertical succession with facies coding from current permutation
            end
        end
        
        % Select permuations with both a high m statistic and a high one-offset diagonal test value as indicated by a high markovDiagProduct
        if markovDiagProduct(j) >= maxMarkovDiagProduct - ROUND_ERROR
            maxDiagCount = maxDiagCount + 1;
            maxDiagPerm(maxDiagCount) = j;
            bestDiagFaciesOrder(maxDiagCount, :) = allFaciesCombos(j,:); % record the facies code permutation, gives same results as optimalCycleAnalysis
            
            % Make a vertical succession with each facies occurrence coded
            % according to the  facies coding from the jth permutation
            for k =1:sectionLength
                oneFacies = data.faciesCodes(k);
                bestDiagFaciesSect(maxDiagCount, k) = allFaciesCombos(j, oneFacies); % Make a complete vertical succession with facies coding from current permutation
            end
        end
    end

    % **************************************************************************************
    % Report all of the results as text output
    fullFileName = strcat('optimalCyles', '_report.txt'); % Save the permutation results to .mat file in ascii format
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
    fprintf('%s ', data.faciesNames{k});
    fprintf('\nOptimal codings are:');
    fprintf(fOut,'Original facies: ');
    fprintf(fOut,'%s ', data.faciesNames{k});
    fprintf(fOut,'\nOptimal codings are:');
    for j=1:maxDiagCount
         fprintf('\n%d ',j);
         fprintf(fOut,'\n%d ',j); 
         fprintf('%s ', data.faciesNames{bestDiagFaciesOrder(j,k)}); % NB the implied loop on k in this - separate formatting fprints required above therefore
         fprintf(fOut,'%s ',data.faciesNames{bestDiagFaciesOrder(j,k)});
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
        fprintf('%s ', data.faciesNames{worstFaciesOrder(j,k)}); % NB the implied loop on k in this - separate formatting fprints required above therefore
        fprintf(fOut,'%s ', data.faciesNames{worstFaciesOrder(j,k)});
    end
    fprintf('\n');
    fprintf(fOut,'\n');
    fclose(fOut);

%    fullFileName = strcat(fileName, '_report_perms.txt'); % Save the permutation results to .mat file in ascii format
    fullFileName = strcat('optimalCycles', '_report_perms.txt'); % Save the permutation results to .mat file in ascii format
    save(fullFileName, 'markovOrderMetricTest','-ascii');
    
    
    % Now finally, compare the first optimally coded section with random
    % shuffled sections to see if the m value gives any evidence of an
    % orderd succession
    maxIterations = 5000;
    numberOfSwaps = data.sectLength;

    % Calculate and output the order metric for the first best diagonal optimally coded succession
    [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(bestDiagFaciesSect(1,:), data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01]));  % note the zero is the TP matrix plot flag, set to false to not trigger a plot
    fprintf('Markov metric for optimised facies order strata is %4.3f\n', markovOrderMetric);

    % Now calculate the metrics for many iterations of a random model
    % Shuffle the observed section and calculate the order metric and DW stat each time
    for j = 1:maxIterations;
        [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(data.faciesCodes, data.faciesThick, numberOfSwaps);
        multiMarkovOrderMetricDataShuffled(j) = calculateTPMatrixAndOrderMetric(shuffledFacies, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01])); % note the zero is the TP matrix plot flag, set to false so no plot
    end

    % Stats on the shuffled section random model
    meanMultiMarkovDataShuffled = mean(multiMarkovOrderMetricDataShuffled);
    stdDevMultiMarkovDataShuffled = std(multiMarkovOrderMetricDataShuffled);
    bins = 0:0.02:1.00; % because 0<m<=1
    multiMarkovOrderMetricDataShuffledHistData=histc(multiMarkovOrderMetricDataShuffled, bins) / maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
    markovPValueSum = sum(multiMarkovOrderMetricDataShuffledHistData(int16(markovOrderMetric*50:length(bins)))); % area under curve from m to max m value 1

    fprintf('For %d iterations of a SHUFFLED DATA model\n', maxIterations);
    fprintf('Markov stats mean %5.4f std dev %5.4f For optimised facies coding Markov order metric P Value %5.4f\n', meanMultiMarkovDataShuffled, stdDevMultiMarkovDataShuffled, markovPValueSum);


% *************************************************************************************
% Plot the results graphically
    % scrsz = left bottom width height

     scrsz = get(0,'ScreenSize'); % screen dimensions vector
     gui.f3 = figure('Visible','on','Position',[1 1 (scrsz(3)*0.75) (scrsz(4)*0.75)]);

     % TP matrix for the original facies coding
     h1 = subplot('Position', [0.05 0.55 0.20 0.3]);
    [mStat,diag] = calculateTPMatrixAndOrderMetric(data.faciesCodes, data.faciesNames, 1, h1);
    axis tight;
    labelStr = cell(maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
    h1.XTickLabel = {data.faciesNames{1:maxNumbOfFacies}}; % because this is a plot for the original order of facies row coding, just use that
    h1.YTickLabel = {data.faciesNames{1:maxNumbOfFacies}};
    str = sprintf('TP Matrix, original coding\nm=%4.3f diag=%4.3f m*diag %4.3f', mStat, diag, mStat * diag);
    title(str);

    % Plot the original facies coding order
    h2=subplot('Position',[0.27 0.55 0.03 0.3]);
    for k=1:maxNumbOfFacies
        yco = [k k k+1 k+1];
        xco = [1 2 2 1];
        faciesCol = [data.faciesColours(k,2) data.faciesColours(k,3) data.faciesColours(k,4)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
        labelStr = data.faciesNames(k);
        text(1.25, k+0.2, labelStr, 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
        axis off;
    end
    
    % TP matrix for the worst permutation
    h3=subplot('Position',[0.35 0.55 0.20 0.3]);
    [mStat,diag] = calculateTPMatrixAndOrderMetric(worstFaciesSect(1,:), data.faciesNames, 1, h3);
    axis tight;
    colourCode = zeros(maxNumbOfFacies, 3); % make sure colour code matrix is ready but empty
    labelStr = cell(1,maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
    for k=1:maxNumbOfFacies
        % For each jth permutation in each worstFaciesOrder the array elements are in original facies order, 
        % so element j,k gives the facies k row code for the jth worst permutation, in this case worst permutation 1
        labelStr{worstFaciesOrder(1,k)} = data.faciesNames{k}; 
        colourCode(worstFaciesOrder(1,k),1:3) = data.faciesColours(k,2:4);
    end
    h3.XTickLabel = labelStr;
    h3.YTickLabel = labelStr;
    str = sprintf('Worst TP Matrix, worst facies coding\nm=%4.3f m*diags %4.3f', mStat, mStat * diag);
    title(str);
    
    % Plot a worst permutation facies coding order
    h4=subplot('Position',[0.57 0.55 0.03 0.3]);
    for k=1:maxNumbOfFacies
        yco = [k k k+1 k+1];
        xco = [1 2 2 1];
        faciesCol = [colourCode(k,1) colourCode(k,2) colourCode(k,3)]; % Colour codes are set in the TP plot code above
        patch(xco, yco, faciesCol,'EdgeColor','black');
        text(1.25, k+0.2, labelStr(k), 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
        axis off;
    end
    
    % TP matrix for the first best permutation
    h5=subplot('Position',[0.65 0.55 0.20 0.3]);
    [mStat,diag] = calculateTPMatrixAndOrderMetric(bestDiagFaciesSect(1,:), data.faciesNames, 1, h5);
    axis tight;
    colourCode = zeros(maxNumbOfFacies, 3); % make sure colour code matrix is ready but empty
    labelStr = cell(1,maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
    for k=1:maxNumbOfFacies
        % For each jth permutation in each bestDiagFaciesOrder the array elements are in original facies order, 
        % so element j,k gives the facies k row code for the jth best permutation, in this case best permutation 1
        labelStr{bestDiagFaciesOrder(1,k)} = data.faciesNames{k}; 
        colourCode(bestDiagFaciesOrder(1,k),1:3) = data.faciesColours(k,2:4);
    end
    h5.XTickLabel = labelStr;
    h5.YTickLabel = labelStr;
    str = sprintf('Best TP Matrix, optimal facies coding\nm=%4.3f m*diags %4.3f', mStat, mStat * diag);
    title(str);

    % Plot a best permutation facies coding order
    h6=subplot('Position',[0.87 0.55 0.03 0.3]);
    for k=1:maxNumbOfFacies
        yco = [k k k+1 k+1];
        xco = [1 2 2 1];
        faciesCol = [colourCode(k,1) colourCode(k,2) colourCode(k,3)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
        text(1.25, k+0.2, labelStr(k), 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
        axis off;
    end

     % plot the histogram of m values from all permutations
    h4=subplot('Position',[0.05 0.1 0.18 0.35]);
    histogram(markovOrderMetricTest, 'EdgeColor','none', 'FaceColor','blue');
    xlabel('m value', 'FontSize',10);
    ylabel('Frequency', 'FontSize',10);
    grid on;
    
    % plot the histogram of diagonal scores from the n maximum score permutations
    h5=subplot('Position',[0.27 0.1 0.18 0.35]);
    histogram(markovDiagProduct, 'EdgeColor','none', 'FaceColor','blue');
    xlabel('d', 'FontSize',10);
    ylabel('Frequency', 'FontSize',10);
    grid on;

    % Subplot for the Markov order analysis histogram   
    h6 = subplot('Position',[0.50 0.1 0.40 0.35]);
    bins = 0:0.02:1.00; % Make sure bins is set correctly for Markov plots
    bar(bins, multiMarkovOrderMetricDataShuffledHistData, 'EdgeColor','none', 'BarWidth', 1, 'FaceColor',[0.2 0.4 0.7]); % Colour is dark slate blue
    maxFreq = max(multiMarkovOrderMetricDataShuffledHistData) *1.1; % This is needed to scale the plot
    x = [markovOrderMetric markovOrderMetric];
    y = [0 maxFreq];
    line(x,y, 'color', [0.80 0.00 0.00], 'linewidth', 3.0); % Colour is dark red
    grid on;
    axis([0 1 0 Inf]);
    set(gca,'Layer','top');
    xlabel('Markov Order Metric for Facies');
    ylabel('Relative Frequency');
    
     % Complete analysis with a textbox message summary
    message = sprintf('Calculated %d facies coding premutations\nLargest markov stat value %5.4f found from %d permutations \nHighest m value * one-offset diagonal value %5.4f found from %d permutations\nFor optimised facies coding Markov order metric P Value %5.4f\n',...
        permsCount, maxMarkov, maxMarkovCount, maxMarkovDiagProduct, maxDiagCount, markovPValueSum);
    m1 = msgbox(message,'Optimised cycle analysis');
end



