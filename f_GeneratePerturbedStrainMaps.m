%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 2024-10-11
%
% Script to generate strain rate maps from preturbed velocities over
% the training regions. This code relies on and repackages the
% logarithmic strain code first published in Alley et al (2018) 
% "Continent-wide estimates of Antarctic strain rates from Landsat 8
% -derived velocity grids (doi: 10.1017/jog.2018.23). This code was used to
% generate strain rate maps for 240 model ensemble used to quantify
% the tensile strength of ice slabs
%
% Inputs:
% v_x and v_y velocity maps for the GrIS
% Gaussian random fields for velocity perturbation (10 x Gaussian
% spatial covariance function, 10 x Exponential spatial covariance
% function)
%
% Outputs:
% Nx3 GeoTIFF file for each grid-aligned strain rates over the training 
% regions (e_xx, e_yy, and e_xy). 
% N = number of velocity perturbation experiments
% Output files are saved in the `out_dir` location specified at the top of
% this script and following the naming convention: 
% Exp[experiment number]_eXX_[covariance file code].tiff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Add utilities directory to path to access strain rate codes
addpath("./utilities");

% Directory in which to save the output strain rate maps
out_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\IceSlabCrevasseData\ErrorMaps";

% Directory where the velocity maps are located
vel_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\IceSlabCrevasseData\Velocity\";
% Absolute paths to v_x and v_y velocity fields
vx_file = strcat(vel_dir, "greenland_vel_mosaic250_vx_training.tif");
vy_file = strcat(vel_dir, "greenland_vel_mosaic250_vy_training.tif");
% No data value for the velocity grids
vxNoData = -2e9;
vyNoData = -2e9;
% Shorthand code for each of the velocity data sources:
% MAvg: MEASURES multi-decadal average velocities
% M12: MEASURES 2012-2013 winter velocities
% ILAvg: ITS_LIVE multi-decadal average velocities
file_code = "MAvg";

% Velocity grid resolution in meters for each velocity data source
pixel_size = 250;
% Tolerance for the adaptive time stepping scheme - should not need to
% change this
tol = 10^-4;
% Set to 1 if the positive y-direction is in the upwards direction on the 
% screen, and -1 if the positive y-direction is downwards
ydir = 1;

filt_name = "Gauss";                   % Filter type (gaussian)
thick_multiplier = 0;                  % 0 since using a kernel filter
thick_grid = 0;                        % 0 since using a kernel filter
smooth_length = 4500;                  % Kernel width in meters

% Directory where the error maps for the velocity grids are saved
error_dir = "C:\Users\rtc84\Documents\Data\Greenland\Crevasses\IceSlabCrevasseData\ErrorMaps\";
% Absolute paths to v_x error and v_y error fields
ex_file = strcat(error_dir, "greenland_vel_mosaic250_ex_training.tif");
ey_file = strcat(error_dir, "greenland_vel_mosaic250_ey_training.tif");

% Array of ID numbers for each velocity perturbation iteration
err_files = 0:1:9;
% Type of correlation function for GRFs
err_type = ["Gaussian" "Exponential"];

%% Load in the velocity files

% Read in velocity geotiffs
[vx_in, vx_info] = readgeoraster(vx_file); % x-component of velocity
[vy_in, ~] = readgeoraster(vy_file); % y-component of velocity
tiffinfo=geotiffinfo(vx_file); % info from either geotiff

% Read in error fields for the velocity maps
[ex_in, ~] = readgeoraster(ex_file); % x-component of velocity
[ey_in, ~] = readgeoraster(ey_file); % y-component of velocity

%% Process the strain maps

for i = 1:length(err_files)
    fprintf("Processing %d.\n", err_files(i));
    for j = 1:length(err_type)
        fprintf("Processing %s.\n", err_type(j));

        % Read in the Gaussian Random Field (GRF) with the appropriate
        % correlation structure
        if err_type(j) == "Gaussian"
            rgf = readmatrix(strcat(error_dir, "Gaussian\field_", num2str(err_files(i)), ".txt"));
        else
            rgf = readmatrix(strcat(error_dir, "Exponential\field_", num2str(err_files(i)), ".txt"));
        end
        % Create perturbed velocity maps by adding the error-scaled GRF to the
        % original velocity measurements
        vx_new = vx_in + ex_in.*rgf(1:size(ex_in,1),1:size(ex_in,2));
        vy_new = vy_in + ey_in.*rgf(1:size(ex_in,1),1:size(ex_in,2));

        % Create the appropriate smoothing kernel (see GenerateStrainMaps.m
        % for detailed comments on this processing)
        if filt_name  == "Box"
            kernel_size = round(smooth_length/pixel_size);
            if mod(kernel_size,2) ~= 1
                kernel_size = kernel_size - 1;
            end
            kernel = (1/(kernel_size*kernel_size))*ones(kernel_size,kernel_size);
            length_scale = pixel_size;
            maxR = ceil(length_scale/pixel_size);
        elseif filt_name == "Tri"
            kernel_size = round(smooth_length/pixel_size);
            if mod(kernel_size,2) ~= 1
                kernel_size = kernel_size - 1;
            end
            starter = (1/(kernel_size*kernel_size))*ones(kernel_size,kernel_size);
            kernel = conv2(starter,starter,"same");
            kernel = kernel/sum(kernel(:));
            length_scale = pixel_size(i);
            maxR = ceil(length_scale/pixel_size);
        elseif filt_name == "Gauss"
            kernel = round(0.25*(smooth_length-500)/pixel_size);
            length_scale = pixel_size;
            maxR = ceil(length_scale/pixel_size);
        else
            length_scale = 0;
            kernel = 0;
            maxR = ceil((thick_multiplier*2500)/pixel_size);
        end

        % Pass input values to logarithmic strain code from Alley et al.
        % (2018). The function has been modified slightly to save outputs
        % under the appropriate file names.
        calculateErrorStrain(out_dir, vx_new, vy_new, tiffinfo, filt_name, kernel, thick_files,...
        vxNoData, vyNoData, thickNoData, pixel_size, tol, ydir, thick_grid,...
        thick_multiplier, length_scale, maxR, err_type(j), err_files(i));

    end
end

