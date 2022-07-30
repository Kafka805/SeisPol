function  [rectilinearity, planar, azimuth, incident, segments] = SeisPol(fileEast,fileNorth,fileVert,period,window)
    %This function imports a seismic signal in the mSEED format, filters the
    %signal using Butterworth filters, then calculates the polarization of
    %the signal using methods found in Jurkevics (1988), Bostock (2012), ...
    %Huesca-Perez (2014). IT IS HEAVILY RECOMMENDED that one import and plot
    %their signal files using the included rdmseed function (or other signal import
    %tool) BEFORE running this analysis tool to determine if the signal meets
    %requirements for data collection.
    %   Inputs: 'fileEast/North/Vert' must all be singular, non-multiplexed mSEED
    %   files. 'period' must be a duration object. 'window' must be defined
    %   as a timerange object with the format "YYYY-MM-DD HH:MM:SS", existing within
    %   the period of time contained by the mSEED signal. It is recommended to
    %   define these inputs prior to calling the function, so that they may be
    %   changed as needed in the workspace.
%% NECESSARY FUNCTIONS: 
%             - rdmseed.m (reads the miniSEED file and imports to MATLAB)
%             - filtermerge.m (automatically pads the signal data and applies
%                    the Butterworth filter. Edit this function to change the 
%                    nature of the filter.)
%             - buttern_filter.m (the butterworth filter function)
%             - polarityAngles.m (performs the polarization analysis)
%             - order.m (nested within polarityAngles. Organizes eigenvalues
%                     and vectors)
%% Polarization Analysis of Tremor Signals
% Austin Abreu, UCSC 2021, Under the direction of Susan Schwartz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load datasets from files
isMat = strfind(fileEast,'.mat');

if isempty(isMat) == 1
    resp = input('Files found in miniSeed. Would you like to convert to mat for faster processing? Y/N [Y]: ','s');
    disp('Importing miniSeed files...' + newline)
    
    East_raw = rdmseed(fileEast);
        BHE = [cat(1,East_raw.d) cat(1,East_raw.t)];
        nameEast = string(append(East_raw(1).NetworkCode, '.',...
            strtrim(East_raw(1).StationIdentifierCode),'.',...
            East_raw(1).ChannelIdentifier,'.',num2str(East_raw(1).RecordStartTime(2))));
        
    North_raw = rdmseed(fileNorth);
        BHN = [cat(1,North_raw.d) cat(1,North_raw.t)];
        nameNorth = string(append(North_raw(1).NetworkCode, '.',...
            strtrim(North_raw(1).StationIdentifierCode),'.',...
            North_raw(1).ChannelIdentifier,'.',num2str(North_raw(1).RecordStartTime(2))));
        
    Vert_raw = rdmseed(fileVert);
        BHZ = [cat(1,Vert_raw.d) cat(1,Vert_raw.t)];
        nameVert = string(append(Vert_raw(1).NetworkCode, '.',...
            strtrim(Vert_raw(1).StationIdentifierCode),'.',...
            Vert_raw(1).ChannelIdentifier,'.',num2str(Vert_raw(1).RecordStartTime(2))));
        
        
    if matches(resp,'Y') == 1 || matches(resp,'y') == 1
        disp('Converting Files...' + newline)
        
        fE = append(nameEast,'.mat');
        save(fE,'East_raw')
        disp('Saved '+fE+' in current working directory.')
            
        fN = append(nameNorth,'.mat');
        save(fN, 'North_raw')
        disp('Saved '+fN+' in current working directory.')
        
        fV = append(nameVert,'.mat');
        save(fV, 'Vert_raw')
        disp('Saved '+fV+' in current working directory.')
        disp('Conversion done, continuing analysis' + newline)
        
    elseif matches(resp,'abort') == 1
        return
        
    end
            
    elseif isempty(isMat) == 0
        disp('Importing .mat files...' + newline)

        load(fileEast);
        load(fileNorth);
        load(fileVert);

        BHE = [cat(1,East_raw.d) cat(1,East_raw.t)];
        nameEast = string(append(East_raw(1).NetworkCode, '.',...
            strtrim(East_raw(1).StationIdentifierCode),'.',...
            East_raw(1).ChannelIdentifier,'.',num2str(East_raw(1).RecordStartTime(2))));

        BHN = [cat(1,North_raw.d) cat(1,North_raw.t)];
        nameNorth = string(append(North_raw(1).NetworkCode, '.',...
            strtrim(North_raw(1).StationIdentifierCode),'.',...
            North_raw(1).ChannelIdentifier,'.',num2str(North_raw(1).RecordStartTime(2))));

        BHZ = [cat(1,Vert_raw.d) cat(1,Vert_raw.t)];
        nameVert = string(append(Vert_raw(1).NetworkCode, '.',...
            strtrim(Vert_raw(1).StationIdentifierCode),'.',...
            Vert_raw(1).ChannelIdentifier,'.',num2str(Vert_raw(1).RecordStartTime(2))));

    end
    
stationName = string(strtrim(East_raw(1).StationIdentifierCode));
sr = East_raw(1).SampleRate
%% Apply the filter and construction script to consolidate the data into arrays of equal size   
[signal] = filtermerge(BHE, BHN, BHZ, sr);

%% Applies the investigation window to the filtered signal.
    inspect = signal(window,:);
%% Create a moving average filter to smooth out the signal (we will use this for viewing the signal)
     %filteredBHE = movmean(abs(inspect.E),240);

%% Plot the inspection area
    SignalPlot = figure(1);
    % East Component Signal Plot
    subplot(3,1,1)
    plot(inspect.time,inspect.E,'Green')
    titleE = append(nameEast, ' with 2nd Order Butterworth Bandpass (2-8Hz)');
    title(titleE)
    xlabel('Time')
    ylabel('Counts')
    datetick('x','HH:MM:SS','keeplimits')
    axis tight

    subplot(3,1,2)
    % North Component Signal Plot
    plot(inspect.time,inspect.N,'Blue')
    titleN = append(nameNorth, ' with 2nd Order Butterworth Bandpass (2-8Hz)');
    title(titleN)
    xlabel('Time')
    ylabel('Counts')
    datetick('x','HH:MM:SS','keeplimits')
    axis tight

    subplot(3,1,3)
    % Vertical Component Signal Plot
    plot(inspect.time,inspect.Z,'Red')
    titleV = append(nameVert, ' with 2nd Order Butterworth Bandpass (2-8Hz)');
    title(titleV)
    xlabel('Time')
    ylabel('Counts')
    datetick('x','HH:MM:SS','keeplimits')
    axis tight

%% Applies covariance and SVD methods to resolve eigenvalues
%This method parses a window of time (Outer selection) and shifts the frame
%(Inner selection) within that window of time.

    %The frame length. This frame shifts half of a window length
    frameShift = period/2;
    
    %Create the window using timetable methods
    analysisFrame = [inspect.time(1) inspect.time(1)+period];
    
    %calculate the number of outputs for variable initialization:
    numOuts = ceil((inspect.time(end) - inspect.time(1)) / frameShift);
    
    %initialize variables for the loop
    rectilinearity = zeros(numOuts,1);
    planar = zeros(numOuts,1);
    azimuth = zeros(numOuts,1);
    incident = zeros(numOuts,1);

    %The tapered cosine window method. We'll blend the edges with 0.5
    %   amplitude scaling.
    pointWindow = (seconds(period)*sr);
    cTW = tukeywin(pointWindow, 0.5);
    
    div = pointWindow/2;
   
% The loop performs all the relevant analytic operations and compiles the results.
for idx = 1:numOuts
    %wndw = isbetween(inspect.time,analysisFrame(1),analysisFrame(2));
    if length(cTW) ~= length(inspect.E(div))
        cTW = tukeywin(length(inspect.E(div)), 0.5);
    end
        [rectilinearity(idx),planar(idx),azimuth(idx),incident(idx)] = ...
            eigenAngles(inspect.E(((idx-1)*div+1):idx*div).*cTW,inspect.N(((idx-1)*div+1):idx*div).*cTW,inspect.Z(((idx-1)*div+1):idx*div).*cTW);
end



%create an array of timing alignments for the X-axis of figures
segments = inspect.time(1)+frameShift:frameShift:inspect.time(end)+frameShift;

%calculates difference in eigenvalues
%eigenDiff = abs(rectilinearity-planar);

%logical array defining values that fall w/in criteria
criterion = [rectilinearity >= 0.4]; 
%% Correlation Plot
correlationPlot = figure(8);

a1 = subplot(4,1,1)
    plot(inspect.time,inspect.N,'Green')
        titleE = append(nameNorth);
        title(titleE), xlabel('Time'), ylabel('Counts'), datetick('x','HH:MM:SS','keeplimits')
        ylim([-2e-8 2e-8])
        axis tight
        
% a2 = subplot(5,1,2)
%     plot(inspect.time,inspect.N,'Blue')
%         titleE = append(nameNorth);
%         title(titleE), xlabel('Time'), ylabel('Counts'), datetick('x','HH:MM:SS','keeplimits')
%         ylim([-2e-8 2e-8])
%         axis tight

a2 = subplot(4,1,2)
    plot(segments,rectilinearity, 'b+')
        hold on
    plot(segments(criterion),rectilinearity(criterion), 'r+')
    plot(segments,movmean(rectilinearity,(seconds(360)/period) + 1))
        hold off
         title('Rectilinearity of ' + stationName + ' Signal, ' + string(period) + ' Analysis Window')
    title('Rectilinear Polarization Ratio, DLP')
    xlabel('Time'), ylabel('DLP'), axis tight, ylim([0 1])
        
a3 = subplot(4,1,3)
    %scatter(segments, incident, 'bo')
        %hold on
    scatter(segments(criterion), incident(criterion),'bx')
        %hold off
    title('Ray Angle of ' + stationName + ' Signal, ' + string(period)+ ' Analysis Windows')
    xlabel('Time'), ylabel('Ray Angle (degrees)'), axis tight

a4 = subplot(4,1,4)
    %scatter(segments, azimuth, 'mo')
        %hold on
    scatter(segments, azimuth,'bx')
       % hold off
    title('Back-Azimuth of ' + stationName + ' Signal with ' + string(period)+ ' Analysis Windows')
    xlabel('Time'), ylabel('Back-Azimuth (degrees)'), axis tight
    
    %linkaxes([a1,a2,a3,a4],'x')
%% File writiing syntax
% If you would like to save more variables, follow the format under the save command.
% timeBounds = [datest(inspect.time(1)) datest(datetime(inspect.time(end)+seconds(1)))];
% currentPath = pwd;
% fileName = append(extract(fileEast,'Analysis_10_18_25\'),replace(datestr(timeOrigin),":","-"));
% oldPath = cd(fileName);
% data = struct('Linearity',rectilinearity,'Polarization',planar,'Azimuth',...
%     azimuth','Incident',incident,'Time_Segments',segments);
% save data.mat
% save SignalPlot.fig
% save(fileName, 'data'); %declare additional after fileName

%% Polarization ratio calculator
function [dlp, dpp, azimuthal, incident] = eigenAngles(east,north,vert)
    matrix = [east north vert];
    CM = (1/(max(size(matrix))-1))* (matrix'*matrix); %Covariance Matrix
    [eVecs,eVals] = eig(CM);
    
    [v,d] = order(eVecs,eVals); %Calls the order function to organize the eigs
    dlp = [1 - ((d(2) + d(3)) / (d(1)))]; %Rectilinear polarization factor
    dpp = [1 - ((d(3))/(d(1) + d(2)))]; %Planar polarization facator
    alpha = v(1,1)/norm(v(:,1)); %direction cosine of the largest eigenvalue
    beta = v(1,2)/norm(v(:,1)); %direction cosine of the 2nd eigenvalue
    gamma = v(1,3)/norm(v(:,1)); %direction cosine of the 3rd eigenvalue
    incident = atan2d(sqrt(alpha^2 + beta^2),gamma); %incident angle
    azimuthal = atan2d((alpha*sign(gamma)),(beta*sign(gamma))); %azimuthal angle
    
    %constrain the azimuthal angle
    if azimuthal < 0
        azimuthal = azimuthal+360;
    elseif  azimuthal > 360
        azimuthal = azimuthal-360;
    end
    %constrain the incident angle
    if incident < 0
        incident = incident+180;
    elseif  incident > 180
        incident = incident-180;
    end
end

function [filteredsignal] = filtermerge(east, north, vert, dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check if each input is the same size before they are merged
% into one array. 
east = [detrend(east(:,1))-mean(east(:,1)) east(:,2)];
north = [detrend(north(:,1))-mean(north(:,1)) north(:,2)];
vert = [detrend(vert(:,1))-mean(vert(:,1)) vert(:,2)];

if ~isequal(length(east), length(north), length(vert))
    %determines whether the signals are equal in length. Note this is determining a negative statement.
    u = length(east);
    o = length(north);
    t = length(vert);
    arraySize = max([u o t]);
    if arraySize == u  %Now we need to equalize the timeseries to match the signal length.
        timeHandOff = east(:,2);
    elseif arraySize == o
        timeHandOff = north(:,2);
    else
        timeHandOff = vert(:,2);
    end
    %and we pad the array to match the size of the smallest array
    n = padarray(east, arraySize - u,0,'post');
    k = padarray(north, arraySize - o,0,'post');
    l = padarray(vert, arraySize - t,0,'post');
else
    %if the primary statement is false, we don't run the interior of the loop and simply hand off the arrays.
    n = east;
    k = north;
    l = vert;
    timeHandOff = east(:,2);
end


time = datetime(timeHandOff,'ConvertFrom','datenum'); %converts date numbers to datetime.
%Calls the butterworth filter function to filter the signal. See the
%function for parameter explanations.
% if filterOn == True || 1

cornerL = 2;
cornerH = 8;
order = 4;
    E = filbutt(n(:,1),dt,cornerL,cornerH);
    N = filbutt(k(:,1),dt,cornerL,cornerH);
    Z = filbutt(l(:,1),dt,cornerL,cornerH);
%     E = buttern_filter(n(:,1),order,cornerL,cornerH,dt);
%     N = buttern_filter(k(:,1),order,cornerL,cornerH,dt);
%     Z = buttern_filter(l(:,1),order,cornerL,cornerH,dt);
% else
%     BHE = n(:,1);
%     BHN = k(:,1);
%     BHZ = l(:,1);
%Collect all arrays into a timetable and export. Timetables allow for easy
%analysis using time & duration categorized programming methods.
filteredsignal = timetable(time,E,N,Z);
end
end