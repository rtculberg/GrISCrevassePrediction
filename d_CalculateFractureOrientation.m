%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 2024-10-11
%
% Applies a Hough transform to binary crevasse maps from Culberg et al.
% (2022) "Shallow fracture buffers high elevation runoff in Northwest
% Greenland" to estimate the dominant fracture orientation on some grid. 
%
% Inputs:
% Binary crevasse images from Northwest Greenland
% Strain rate maps over Northwest Greenland
% Annual surface temperature map over Northwest Greenland
%
% Outputs:
% csv file with two columns; column 1 is the orientation of observed
% fractures in degrees as measured counterclockwise from the positive x
% axis, column 2 is the orientation of the first principle stress as
% measured counterclockwise from the positive x axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load Data

clear;

% Directory where binary crevasse maps are saved
in_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\NWIceSlabs\no_border\";

% List of image files to analyze
images = ["WV01_20120803164856" "QB02_20120729152314" "QB02_20120731154958" ...
    "QB02_20120731155001" "QB02_20120731155004" "WV01_20120713164417" ...
    "WV01_20120713164418" "WV01_20120713164419" "WV01_20120803164853" ...
    "WV01_20120803164854" "WV01_20120803164855" ...
    "WV01_20120802153817" "WV01_20120802153816" "WV01_20120802153815" ...     
    "WV01_20120802153814" "WV01_20120802153813" "WV01_20120713005153" ...
    "WV01_20120713005152" "WV01_20120713005151"];

% Directory where Greenland strain rate maps are save
file_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\IceSlabCrevasseData\StrainRate\";
% Load strain rate maps for Northwest Greenland from Experiment 18
% (MEASURES multi-decadal velocities, Gaussian filter with a 1000 m
% bandwidth)
[e_xy, e_R] = readgeoraster(strcat(file_dir, "Exp18_exy_MAvg.tif"));
[e_xx, ~] = readgeoraster(strcat(file_dir, "Exp18_exx_MAvg.tif"));
[e_yy, ~] = readgeoraster(strcat(file_dir, "Exp18_eyy_MAvg.tif"));

% Directory where the RACMO Surface temperatures are saved
temp_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\IceSlabCrevasseData\RACMOTraining\";
% Load the RACMO annual average surface temperature interpolated onto the
% same grid at the velocity data
[temp, ~] = readgeoraster(strcat(temp_dir,"Clip_RACMOAnnual_250m.tif"));

%% Set up Pixel Coordinate Reference for the Strain Maps

e_X = e_R.XWorldLimits(1):e_R.CellExtentInWorldX:e_R.XWorldLimits(2);
e_Y = fliplr(e_R.YWorldLimits(1):e_R.CellExtentInWorldY:e_R.YWorldLimits(2));

%% Calculate Stresses

% Creep rate factor constants from Cuffey and Patterson (2010)
A_star = 3.5e-25;           % Pa^-3 s^-1
T = temp + 273.15;   % K
T_start = 263;       % K
R = 8.314;                  % J mol^-1 K^-1
Qc = 6e4;                   % J mol^-1
Qc2 = 115e3;                % J mol^-1

% Calculate creep rate factor for every grid cell from surface temperature
A = zeros(size(T));
for p = 1:size(T,1)
    for m = 1:size(T,2)
        if T(p,m) > 273.15
            T(p,m) = 273.15;
            A(p,m) = A_star.*exp((-Qc2/R).*((1./T(p,m)) - (1./T_start)));
        elseif T(p,m) > 263.15
            A(p,m) = A_star.*exp((-Qc2/R).*((1./T(p,m)) - (1./T_start)));
        else
            A(p,m) = A_star.*exp((-Qc/R).*((1./T(p,m)) - (1./T_start)));
        end
    end
end

% Convert all strain rates from 1/yr to 1/s 
e_xy_s = e_xy*(1/3.154e7);
e_xx_s = e_xx*(1/3.154e7);
e_yy_s = e_yy*(1/3.154e7);

% Calculate effective strain rate
e_E = sqrt(0.5.*(e_xx_s.^2 + e_yy_s.^2) + e_xy_s.^2);

% Calculate deviatoric stresses from Glen's flow law
n = 3;
tau_xx = (A.^(-1/n)).*(e_E.^((1-n)/n)).*e_xx_s;
tau_yy = (A.^(-1/n)).*(e_E.^((1-n)/n)).*e_yy_s;
tau_xy = (A.^(-1/n)).*(e_E.^((1-n)/n)).*e_xy_s;

% Calculate Cauchy stresses
sigma_xx = 2*tau_xx + tau_yy;
sigma_yy = 2*tau_yy + tau_xx;

%%

f = [];    % Data array to save fracture orientation
s = [];    % Data array to save principal stress orientation 
for k =1   %:length(images)
    fprintf("Processing %s\n", images(k));

    % Load binary fracture maps
    file = strcat(in_dir, "FracMap_", images(k), ".tif");
    [fractures, frac_R] = readgeoraster(file);

    % Load mask files that mask edge effect errors in the binary fracture
    % maps
    mask_file = strcat(in_dir, "DataMask_", images(k), ".tif");
    [mask, ~] = readgeoraster(mask_file);

    % Set up pixel coordinate reference for the fracture maps
    frac_X = frac_R.XWorldLimits(1):frac_R.CellExtentInWorldX:frac_R.XWorldLimits(2);
    frac_Y = fliplr(frac_R.YWorldLimits(1):frac_R.CellExtentInWorldY:frac_R.YWorldLimits(2));

    % Segment strain maps to match fracture image extent
    [~, x1] = min(abs(e_X - frac_X(1)));
    [~, x2] = min(abs(e_X - frac_X(end)));
    [~, y1] = min(abs(e_Y - frac_Y(1)));
    [~, y2] = min(abs(e_Y - frac_Y(end)));
    
    seg_X = e_X(x1:x2);
    seg_Y = e_Y(y1:y2);

    x2 = x2 - 1;
    y2 = y2 - 1;
    e_xy_seg = tau_xy(y1:y2, x1:x2);
    e_xx_seg = sigma_xx(y1:y2, x1:x2);
    e_yy_seg = sigma_yy(y1:y2, x1:x2);
    
    % Calculate fracture orientation and principal stress orientation for
    % each 250 m x 250 m strain pixel
    frac_angle = zeros(size(e_xy_seg));
    stress_angle = zeros(size(e_xy_seg));
    for m = 1:size(e_xy_seg,1)
        for p = 1:size(e_xy_seg,2)
            % Subset the fracture maps
            [~, x3] = min(abs(frac_X - seg_X(p)));
            [~, x4] = min(abs(frac_X - seg_X(p+1)));
            [~, y3] = min(abs(frac_Y - seg_Y(m)));
            [~, y4] = min(abs(frac_Y - seg_Y(m+1)));
            x4 = x4 - 1;
            y4 = y4 - 1;
            frac_img = fractures(y3:y4, x3:x4);
            mask_image = mask(y3:y4, x3:x4);

            % Calculate principal stress orientation in the grid cell
            stress_angle(m,p) = rad2deg(0.5*atan2(2*e_xy_seg(m,p),(e_xx_seg(m,p) - e_yy_seg(m,p))));

            % If we are looking at a valid grid cell and the fracture
            % density in that grid cell exceeds 0.001, then calculate the
            % observed fracture orientation
            if mean(mask_image(:)) >= 0.75 && mean(frac_img(:)) > 0.001

                % Run the Hough transform and extract the top 50 peaks
                [ht,theta,rho] = hough(frac_img);
                peaks = houghpeaks(ht,50,"Threshold", 0.5*max(ht(:)));

                % Remove noise peaks at around 45 degrees
                for w = 1:size(peaks,1)
                    if abs(theta(peaks(w,2))) == 45
                        peaks(w,:) = 0;
                    end
                end

                % Delete any peaks located in row 0 (noise)
                ind = find(peaks(:,1) == 0);
                peaks(ind,:) = [];

                % Take the most dominant (first)peak as the fracture orientation
                % and rotate so it is measured relative to the same axis as
                % the stress orientation
                if size(peaks,1) >= 1
                    frac_angle(m,p) = -theta(peaks(1,2)) + 90;
                else
                    frac_angle(m,p) = NaN;
                end
            else
                frac_angle(m,p) = NaN;
            end
        end
    end

    % Remove any NaN values from the arrays
    a = frac_angle(:);

    b = stress_angle(:) + 90;
    ind = find(isnan(a));

    a(ind) = [];
    b(ind) = [];

    ind = find(isnan(b));

    a(ind) = [];
    b(ind) = [];

    % Store results
    f = [f; a];
    s = [s; b];

end

%% Save Data

writematix([f s], "FractureOrientation.csv");


