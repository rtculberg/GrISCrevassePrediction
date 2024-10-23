%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 2024-10-11
%
% Script to generate strain rate maps from remotely sensed velocities over
% Greenland training regions. This code relies on and repackages the
% logarithmic strain code first published in Alley et al (2018) 
% "Continent-wide estimates of Antarctic strain rates from Landsat 8
% -derived velocity grids (doi: 10.1017/jog.2018.23). This code was used to
% generate strain rate maps for the suite of optimization experiments
% described in Section 2.3 of Culberg et al (202X) "von Mises stress a
% robust predictor of ice slab fracture in Greenalnd". 
%
% Inputs:
% v_x and v_y velocity maps for the area of interest
% Bedmachine Greenland ice thickness interpolated onto the same grid as the
% velocity maps
%
% Outputs:
% Nx3 GeoTIFF file for each grid-aligned strain rates over the region of
% interest (e_xx, e_yy, and e_xy). 
% N = number of optimization experiments
% Output files are saved in the `out_dir` location specified at the top of
% this script and following the naming convention: 
% Exp[experiment number]_eXX_[velocity file code].tiff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Add utilities directory to path to access strain rate codes
addpath("./utilities");

% Directory in which to save the output strain rate maps
out_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\IceSlabCrevasseData\StrainRateNE\";

% Directory where the velocity maps are located
vel_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\IceSlabCrevasseData\Velocity\";
% List of absolute paths to the v_x velocity maps for each data source
% under consideration (MEASURES multi-decadal, MEASURES 2012-2013, and
% ITS_LIVE multi-decadal)
vx_files = [strcat(vel_dir, "greenland_vel_mosaic250_vx_NE.tif") ...
    strcat(vel_dir, "greenland_vel_mosaic500_2012_2013_vx_NE.tif") ...
    strcat(vel_dir, "ITSLIVE_average_240m_vx_NE.tif")];
% List of absolute paths to the v_y velocity maps for each data source
vy_files = [strcat(vel_dir, "greenland_vel_mosaic250_vy_NE.tif") ...
    strcat(vel_dir, "greenland_vel_mosaic500_2012_2013_vy_NE.tif") ...
    strcat(vel_dir, "ITSLIVE_average_240m_vy_NE.tif")];
% No data value for each of the three velocity grids
vxNoData = [-2e9 -2e9 -32767];
vyNoData = [-2e9 -2e9 -32767];
% Shorthand code for each of the velocity data sources:
% MAvg: MEASURES multi-decadal average velocities
% M12: MEASURES 2012-2013 winter velocities
% ILAvg: ITS_LIVE multi-decadal average velocities
file_code = ["MAvg" "M12" "ILAvg"];

% Directory where interpolated BedMachine ice thickness maps are stored
thick_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\IceSlabCrevasseData\IceThickness\";
% Absolute paths to the interpolated ice thickness map matching each of the
% three velocity grids
thick_files = [strcat(thick_dir, "BM_Thickness_250m_NE.tif") ...
    strcat(thick_dir, "BM_Thickness_500m_NE.tif") ...
    strcat(thick_dir, "BM_Thickness_240m_NE.tif")];
% No data value for the ice thickness grid
thickNoData = 0;

% Velocity grid resolution in meters for each velocity data source
pixel_size = [250 500 240];
% Tolerance for the adaptive time stepping scheme - should not need to
% change this
tol = 10^-4;
% Set to 1 if the positive y-direction is in the upwards direction on the 
% screen, and -1 if the positive y-direction is downwards
ydir = 1;

% Velocity smoothing filter to use for each experiment
% Box: box filter
% Tri: triangular filter
% Gauss: gaussian filter
% None: use ice-thickness dependent smoothing
filt_name = ["Box" "Box" "Box" "Box" "Box" "Box" "Box" ...
            "Tri" "Tri" "Tri" "Tri" "Tri" "Tri" "Tri" ...
            "Gauss" "Gauss" "Gauss" "Gauss" "Gauss" "Gauss" "Gauss" ...
            "None" "None" "None" "None" "None" "None" "None"];
% Smoothing half length scale when using ice-thickness dependent smoothing,
% set to 0 for all experiments that use kernel-based filters
thick_multiplier = [0 0 0 0 0 0 0 ...
                    0 0 0 0 0 0 0 ...
                    0 0 0 0 0 0 0 ...
                    0.5 1 1.5 2 2.5 3 3.5];
% Set to 1 if using an ice thickness grid for the smoothing, otherwise set
% to 0
thick_grid = [0 0 0 0 0 0 0 ...
              0 0 0 0 0 0 0 ...
              0 0 0 0 0 0 0 ...
              1 1 1 1 1 1 1];
% Kernel width in meters for box and triangular filters, kernel bandwidth
% in meters for the gaussian filter, set to 0 for ice thickness dependent
% smoothing
smooth_length = [1500 2500 3500 4500 5500 6500 7500 ...
                1500 2500 3500 4500 5500 6500 7500 ...
                1500 2500 3500 4500 5500 6500 7500 ...
                0 0 0 0 0 0 0];
% Generate ID numbers for each experiment
experiment = 1:1:length(thick_grid);

%% Generate Strain Maps for Each Experiment

for i = 1:length(vx_files)
    fprintf("Processing %s.\n", vx_files(i));
    for j = 1:2

        fprintf("Filter Type: %s, Filter Length: %d, Ice Thickness Multiplier: %f\n", filt_name(j), smooth_length(j), thick_multiplier(j));
        
        % Generate smoothing kernels for each scenario
        if filt_name(j) == "Box"
            % Find kernel size in number of pixels that best approximates
            % the desired kernel size in meters
            kernel_size = round(smooth_length(j)/pixel_size(i));
            % Ensure that the filter size is odd to get centered estimates
            if mod(kernel_size,2) ~= 1
                kernel_size = kernel_size - 1;
            end
            % Generate the kernel
            kernel = (1/(kernel_size*kernel_size))*ones(kernel_size,kernel_size);
            % Set the virtual strain stake length scale
            length_scale = pixel_size(i);
            % Set the max strain diamond size (should always be length
            % scale divided by pixel size)
            maxR = ceil(length_scale/pixel_size(i));
        elseif filt_name(j) == "Tri"
            % Find kernel size in number of pixels that best approximates
            % the desired kernel size in meters
            kernel_size = round(smooth_length(j)/pixel_size(i));
            % Ensure that the filter size is odd to get centered estimates
            if mod(kernel_size,2) ~= 1
                kernel_size = kernel_size - 1;
            end
            % Generate the kernel by convolving two box kernels and
            % normalizing the output
            starter = (1/(kernel_size*kernel_size))*ones(kernel_size,kernel_size);
            kernel = conv2(starter,starter,"same");
            kernel = kernel/sum(kernel(:));
            % Set the virtual strain stake length scale
            length_scale = pixel_size(i);
            % Set the max strain diamond size (should always be length
            % scale divided by pixel size)
            maxR = ceil(length_scale/pixel_size(i));
        elseif filt_name(j) == "Gauss"
            % Generate a Gaussian kernel with a bandwidth equal to
            % 0.25*(L-500) where L is the kernel width in meters
            kernel = round(0.25*(smooth_length(j)-500)/pixel_size(i));
            % Set the virtual strain stake length scale
            length_scale = pixel_size(i);
            % Set the max strain diamond size (should always be length
            % scale divided by pixel size)
            maxR = ceil(length_scale/pixel_size(i));
        else
            % No kernel needed for ice thickness dependent smoothing
            length_scale = 0;
            kernel = 0;
            maxR = ceil((thick_multiplier(j)*2500)/pixel_size(i));
        end

        % Pass input values to logarithmic strain code from Alley et al.
        % (2018)
        calculateStrain(out_dir, vx_files(i), vy_files(i), filt_name(j), kernel, thick_files(i),...
        vxNoData(i), vyNoData(i), thickNoData, pixel_size(i), tol, ydir, thick_grid(j),...
        thick_multiplier(j), length_scale, maxR, file_code(i), experiment(j));

    end
end

