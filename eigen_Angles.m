%% Austin Abreu | UCSC 2021 | Under the direction of Susan Schwartz, UCSC Seismology
% Calculates the covariance matrix, polarization numbers, and angles of
% signal approach for the input signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Changelog:
%        - 1.0: Added and confirmed angle calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dlp, dpp, azimuthal, incident] = eigen_Angles(east,north,vert)
    matrix = [vert north east];
    CM = (1/(max(size(matrix))-1))* (matrix'*matrix); %Covariance Matrix
    [eVecs,eVals] = eig(CM);
    
    [v,d] = order(eVecs,eVals); %Calls the order function to organize the eigs
    dlp = [1 - ((d(2) + d(3)) / (2*d(1)))]; %Rectilinear polarization factor
    dpp = [1 - ((2*d(3))/(d(1) + d(2)))]; %Planar polarization facator
    alpha = v(1,1)/norm(v(:,1)); %direction cosine of the largest eigenvalue
    beta = v(2,1)/norm(v(:,1)); %direction cosine of the 2nd eigenvalue
    gamma = v(3,1)/norm(v(:,1)); %direction cosine of the 3rd eigenvalue
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