function piv_3d_10(piv_parameters);
% This function is designed to process two or three dimensional PIV data using
% the parameters specified by the 'piv_parameters' data structure.
%
% This code is based upon the code 'basic_3d_rpc_processing_05.m'.
%
% Updates on previous versions:
%
%   Version 08
%
%       Saves the PIV processing parameters data structure to the PIV
%       vector field write directory for later reference.
%
%   Version 09
%
%       Added the ability to use several different outlier vector
%       replacement methods including a local mean interpolation (default
%       in previous versions), a Laplacian interpolation, and a Delaunay
%       interpolation.
%
%   Version 10
%
%       This updated the order of the spatial and temporal dimensions to be
%       consistent with other PIV software versions and to minimize the
%       number of overhead operations (since MATLAB prefers some orders
%       over others).
%
%       The ability to import other data formats was also included.  The
%       code can now import 2D images in one of the following formats: 
%       'bmp', 'jpg', 'jp2', 'png', 'ppm', 'tif', 'mat', 'dcm', 'cdf', or 
%       'h5' while 3D images may be one of the following formats: 'mat', 
%       'dcm', 'cdf', or 'h5'.
%
%       Any NaN velocity vector values that are produced during the 
%       correlations (likely due to simulated data not being seeded within
%       a window) are considered outliers and replaced accordingly.
%
% Things to work on, modify, fix, et cetera:
%
%   Change the code that is optimized for C to code that is optimized for
%   matlab.
%
%   Check on the speed of the if-then statements added to make the code run
%   for 2D images as well.
%
%   Add peak fitting algorithms besides the 3-point Gaussian.
%
%   Add the ability to record multiple correlation peaks, specifically to
%   be used during validation.  Possibly also record the peak ratio.
%
%   Add the ability to do multiple iterations of validation and replacement
%   for each pass.
%
%   Add in the ability to perform deform correlations.
%
% Authors:		Rod La Foy
% First Written On:	31 October 2014
% Last Modified On:	21 November 2014

% This extracts the vector field saving directory from the data processing structure
vector_field_write_directory=piv_parameters.general_parameters.vector_field_write_directory;
% This saves a copy of the processing parameters to the vector field write
% directory for future reference
save([vector_field_write_directory,'piv_processing_parameters.mat'],'piv_parameters');

% This extracts the number of passes to perform
pass_number=piv_parameters.general_parameters.pass_number;

% This is the initial image frame to load
initial_image_frame=piv_parameters.general_parameters.initial_image_frame;
% This is the final image frame to load
final_image_frame=piv_parameters.general_parameters.final_image_frame;
% This is the number of frames to step each iteration
image_frame_step=piv_parameters.general_parameters.image_frame_step;
% This is the number of frames to step between correlation frame pairs
image_correlation_step=piv_parameters.general_parameters.image_correlation_step;

% This is the Boolean value stating whether to perform window deformation
perform_window_deformation=piv_parameters.general_parameters.perform_window_deformation;
% If window deformation is to be performed, this loads the window
% deformation parameters
if perform_window_deformation;
    % This is the minimum number of window deformation operations to
    % perform
    window_deformation_iteration_min=piv_parameters.general_parameters.window_deformation_iteration_min;
    % This is the maximum number of window deformation operations to
    % perform
    window_deformation_iteration_max=piv_parameters.general_parameters.window_deformation_iteration_max;
    % This is the convergence threshhold of the window deformation process
    window_deformation_threshhold=piv_parameters.general_parameters.window_deformation_threshhold;
end;

% This is a Boolean value stating whether to perform pyramid correlations
perform_pyramid_correlations=piv_parameters.general_parameters.perform_pyramid_correlations;

% This is the image filename extension
image_filename_extension=piv_parameters.general_parameters.image_filename_extension;
% This is the image variable name (if the variable is saved within a data
% object that contain multiple variables, ie a 'mat' file)
image_variable_name=piv_parameters.general_parameters.image_variable_name;

% This is a list of the images within the image read directory
image_filename_list=dir([piv_parameters.general_parameters.image_read_directory,piv_parameters.general_parameters.image_filename_prefix,'*',image_filename_extension]);

% This is the filename of the first image file (for the purpose of determing the image
% file size)
first_image_filename=[piv_parameters.general_parameters.image_read_directory,image_filename_list(1).name];

% This finds the size of the images being loaded which can be 'bmp', 'jpg',
% 'jp2', 'png', 'ppm', or 'tif' images in 2D or 'mat' or 'dcm' images in 3D
if any(strcmpi(image_filename_extension,{'bmp','jpg','jp2','png','ppm','tif'}));
    
    % This loads in the first image
    I_First_Image=imread(first_image_filename);
    % This measures the image size
    image_size=size(I_First_Image);
    % This clears the image from memory
    clear('I_First_Image');
    % This adds the third dimension as 1 if it does not exist
    if length(image_size)==2;
        % This adds the third dimension of the image size as 1
        image_size(3)=1;
    end;
    % This clears the first image from memory
    clear('I_First_Image');
    
elseif strcmpi(image_filename_extension,'mat');
    
    % This creates a matlab mat object containing information about the image
    first_image_object=matfile(first_image_filename);
    % This extracts the image size information
    image_size=size(first_image_object,image_variable_name);
    % This adds the third dimension as 1 if it does not exist
    if length(image_size)==2;
        % This adds the third dimension of the image size as 1
        image_size(3)=1;
    end;
    
elseif strcmpi(image_filename_extension,'dcm');
    
    % This reads int the DICOM image
    I_First_Image=dicomread(first_image_filename);
    % This extracts the image size information
    image_size=size(I_First_Image);
    % This adds the third dimension as 1 if it does not exist and removes
    % the third dimension if the image is 3D since the 3rd dimension will
    % correspond to the colormap
    if length(image_size)==2;
        % This adds the third dimension of the image size as 1
        image_size(3)=1;
    elseif length(image_size)==4;
        % This removes the 3rd dimension since this corresponds to the
        % colormap fod DICOM files which is not used for PIV data
        image_size(3)=[];
    end;
    % This clears the first image from memory
    clear('I_First_Image');
    
elseif strcmpi(image_filename_extension,'cdf');
    
    % This reads in the Common Data Format image into a cell array
    I_First_Image=cdfread(first_image_filename,'Variables',image_variable_name);
    % This extracts the image size information from the cell array
    image_size=size(I_First_Image{1});
    % This adds the third dimension as 1 if it does not exist and removes
    % the third dimension if the image is 3D since the 3rd dimension will
    % correspond to the colormap
    if length(image_size)==2;
        % This adds the third dimension of the image size as 1
        image_size(3)=1;
    end;
    % This clears the first image from memory
    clear('I_First_Image');
    
elseif strcmpi(image_filename_extension,'h5');
    
    % This reads in the HDF5 image
    I_First_Image=hdf5read(first_image_filename,image_variable_name);
    % This extracts the image size information
    image_size=size(I_First_Image);
    % This adds the third dimension as 1 if it does not exist and removes
    % the third dimension if the image is 3D since the 3rd dimension will
    % correspond to the colormap
    if length(image_size)==2;
        % This adds the third dimension of the image size as 1
        image_size(3)=1;
    end;
    % This clears the first image from memory
    clear('I_First_Image');
    
end;

% These are the coordinates of the estimated velocities (initially nothing is assummed
% about the velocity field so a zero displacement is used and since the extrapolation value
% is set to 0, only one arbitrary coordinate must be specified)
if image_size(3)==1;
    % This calculates the initial interpolation coordinates with only one
    % kk coordinate value
    [jj_position,ii_position,kk_position]=meshgrid([1,image_size(1)],[1,image_size(2)],[1]);
    % This is the estimated (initially zero) velocity for the 2D array
    ii_velocity=zeros(2,2,1);
    jj_velocity=zeros(2,2,1);
    kk_velocity=zeros(2,2,1);
elseif image_size(3)>1;
    % This calculates the initial interpolation coordinates with the full
    % 3D array of coordinates
    [jj_position,ii_position,kk_position]=meshgrid([1,image_size(1)],[1,image_size(2)],[1,image_size(3)]);
    % This is the estimated (initially zero) velocity for the 3D array
    ii_velocity=zeros(2,2,2);
    jj_velocity=zeros(2,2,2);
    kk_velocity=zeros(2,2,2);
end;

% This iterates through the passes
for pass_index=1:pass_number;
    
    % This displays that the current pass is being processed
    fprintf('\n\nCompleting pass %d of %d passes.\n',pass_index,pass_number);
    
    % This extracts the effective current window resolution
    window_resolution=piv_parameters.pass_parameters(pass_index).window_resolution;
    % This extracts the full window size
    window_size=piv_parameters.pass_parameters(pass_index).window_size;
    % This extracts the window gridding method for the current pass
    window_gridding_method=piv_parameters.pass_parameters(pass_index).window_gridding_method;
    
    % This is the Boolean value stating whether to perform validation on the vector field
    validate_vector_field=piv_parameters.pass_parameters(pass_index).validate_vector_field;
    
    % This is the Boolean value stating whether to perform smoothing of the velocity vector field
    smooth_vector_field=piv_parameters.pass_parameters(pass_index).smooth_vector_field;
    
    % If specified by the user to perform smoothing of the velocity field, this loads the
    % smoothing specific parameters
    if smooth_vector_field;
        
        % This loads the Gaussian smoothing standard deviation
        gaussian_smoothing_kernal_std=piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_std;
        % This loads the Gaussian smoothing kernal size
        gaussian_smoothing_kernal_size=piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_size;
        
    end;
    
    % If specified by the user to perform validation, this loads the loads the validation
    % specific parameters
    if validate_vector_field;
        
        % This loads the vector outlier threshhold limits
        minimum_outlier_threshhold=piv_parameters.pass_parameters(pass_index).minimum_outlier_threshhold;
        maximum_outlier_threshhold=piv_parameters.pass_parameters(pass_index).maximum_outlier_threshhold;
        
        % This loads the UOD kernal size
        uod_kernal_size=piv_parameters.pass_parameters(pass_index).uod_kernal_size;
        % This loads the UOD expected error
        uod_epsilon=piv_parameters.pass_parameters(pass_index).uod_epsilon;
        % This loads the UOD residual threshhold value
        uod_residual_threshhold=piv_parameters.pass_parameters(pass_index).uod_residual_threshhold;
        
        % This loads the outlier vector replacement method
        vector_replacement_method=piv_parameters.pass_parameters(pass_index).vector_replacement_method;
        % This loads the relevant vector replacement parameters based upon
        % the interpolation method
        if strcmp(vector_replacement_method,'local_mean');
            % This loads the minimum number of local vectors used to
            % calculate the local mean
            local_mean_minimum_valid_vector_number=piv_parameters.pass_parameters(pass_index).local_mean_minimum_valid_vector_number;
        elseif strcmp(vector_replacement_method,'laplacian_interpolation');
            % This loads the number of adjacent points to use in 
            % calculating the interpolated values (this can be either 4 or 
            % 8 in 2D and 6, 18, or 26 in 3D)
            laplacian_interpolation_adjacent_connectivity=piv_parameters.pass_parameters(pass_index).laplacian_interpolation_adjacent_connectivity;
        elseif strcmp(vector_replacement_method,'delaunay_interpolation');
            % This loads the Delaunay interpolation weighting method to be
            % used in interpolating values
            delaunay_interpolation_weighting_method=piv_parameters.pass_parameters(pass_index).delaunay_interpolation_weighting_method;
        end;
        
    end;
    
    % This calculates the grid spacing for the current pass using either
    % the amount of window overlap or the window spacing
    if strcmp(window_gridding_method,'window_overlap');
        % This is the ratio of the window overlap distance
        window_overlap=piv_parameters.pass_parameters(pass_index).window_overlap;
        % This is the spacing between the centers of the windows
        grid_spacing=round(window_resolution.*(1-window_overlap));
        % This checks if the overlap percentage is so high that grid spacing is zero (although
        % this should actually never be done)
        if grid_spacing==0;
            % This sets the grid spacing to 1 voxel
            grid_spacing=ones(size(grid_spacing));
        end;
    elseif strcmp(window_gridding_method,'window_spacing');
        % This extracts the window grid spacing from the parameters structure
        grid_spacing=piv_parameters.pass_parameters(pass_index).window_spacing;
    end;
    
    % This calculates the window domains
    [window_min,window_max,window_center,window_number]=calculate_window_domains(image_size,window_resolution,window_size,grid_spacing);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This iterates through the frames calculating the velocity fields.                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This iterates through the image frames calculating the vector fields
    for frame_index=initial_image_frame:image_frame_step:final_image_frame;
        
        % This displays that the current frame is being processed
        fprintf('\nProcessing frame %d of %d frames.\n',frame_index,final_image_frame);
        
        % If this is the second or higher pass, this loads the
        % previously calculated vector field which is used to deform
        % the images in the pyramid correlations
        if pass_index>1;
            
            % This is the filename to load the data from
            data_filename_read=[vector_field_write_directory,'pass_',sprintf('%02.0f',pass_index-1),'_frame_',sprintf('%04.0f',frame_index),'_x_',sprintf('%04.0f',frame_index+image_correlation_step),'.mat'];
            % This loads in the data file
            load(data_filename_read);
            
            % This renames the data loaded from memory to be conistent 
            % with the variables names within this code (ie Prana uses X,
            % Y, Z, U, V, W for the position and velocity values and this
            % code uses the index dimensions ii, jj, kk)
            ii_position=Y;
            jj_position=X;
            kk_position=Z;
            ii_velocity=V;
            jj_velocity=U;
            kk_velocity=W;
            
        end;
        
        % This checks whether standard correlations are supposed to be
        % performed
        if not(perform_pyramid_correlations);
            
            % This performs the standard correlations of the current frame
            [ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity]=standard_correlations(piv_parameters,pass_index,frame_index,window_number,window_min,window_max,window_center,ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity);
            
        elseif perform_pyramid_correlations;
            
            % This performs the pyramid correlations of the current frame
            [ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity]=pyramid_correlations(piv_parameters,pass_index,frame_index,window_number,window_min,window_max,window_center,ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity);
            
        end;
        
        % This performs validation of the vector field if specified by the user
        if validate_vector_field;
            
            % This displays that the vector field is being validated
            disp('Validating the vector field . . . ');
            
            % This determines the outliers of the vector field based upon several metrics
            outlier_vector_array=locate_vector_outliers(ii_velocity,jj_velocity,kk_velocity,minimum_outlier_threshhold,maximum_outlier_threshhold,uod_kernal_size,uod_epsilon,uod_residual_threshhold);
            
            % This replaces the outlier vectors using the specified method
            if strcmp(vector_replacement_method,'local_mean');
                % This replaces the outlier vectors using a local mean of valid vectors
                [ii_velocity,jj_velocity,kk_velocity]=local_mean_replacement(ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity,outlier_vector_array,local_mean_minimum_valid_vector_number);
            elseif strcmp(vector_replacement_method,'laplacian_interpolation');
                % This replaces the outlier vectors using Laplacian interpolation
                [ii_velocity,jj_velocity,kk_velocity]=laplacian_replacement(ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity,outlier_vector_array,laplacian_interpolation_adjacent_connectivity);
            elseif strcmp(vector_replacement_method,'delaunay_interpolation');
                % This replaces the outlier vectors using a weighted sum based upon the
                % Delaunay triangulation of the valid vectors
                [ii_velocity,jj_velocity,kk_velocity]=delaunay_replacement(ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity,outlier_vector_array,delaunay_interpolation_weighting_method);
            end;
            
        end;
        
        % This performs smoothing of the vector field if specified by the users
        if smooth_vector_field;
            
            % This displays that the vector field is being smoothed
            disp('Smoothing the vector field . . . ');
            
            % This smooths the vector fields
            [ii_velocity,jj_velocity,kk_velocity]=smooth_velocity_field(ii_velocity,jj_velocity,kk_velocity,gaussian_smoothing_kernal_size,gaussian_smoothing_kernal_std*ones(1,3));
            
        end;
        
        % This is the filename to save the data as
        data_filename_write=[vector_field_write_directory,'pass_',sprintf('%02.0f',pass_index),'_frame_',sprintf('%04.0f',frame_index),'_x_',sprintf('%04.0f',frame_index+image_correlation_step),'.mat'];

        % This renames the data from the current vector field to be 
        % conistent with the variables names used by other vector field 
        % code (ie Prana uses X, Y, Z, U, V, W for the position and 
        % velocity values and this code uses the index dimensions ii, jj, 
        % kk)
        Y=ii_position;
        X=jj_position;
        Z=kk_position;
        V=ii_velocity;
        U=jj_velocity;
        W=kk_velocity;
        
        % This renames the validation array to be consistent with the
        % variables used by other vector field code (ie Prana uses Valid as
        % the name for the valid vectors and this code uses
        % outlier_vector_array)
        if validate_vector_field;
            % This creates the binary array variable 'Valid' which gives
            % the locations of the vectors that were considered valid
            Valid=not(outlier_vector_array);
        else;
            % This creates an array of NaNs for the variable 'Valid' since
            % validation was not performed on the vector field
            Valid=NaN(size(X));
        end;
        
        % This saves the vector field data to memory
        save(data_filename_write,'X','Y','Z','U','V','W','Valid','-v7.3');
        
    end;
    
end;



function [ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity]=standard_correlations(piv_parameters,pass_index,frame_index,window_number,window_min,window_max,window_center,ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity);
% This function calculates standard pair correlations of the data set listed
% within the parameters variable 'piv_parameters'.  The standard correlation
% consists of taking pair-wise correlations from a time series of images.
% The input parameters of this function are given by
%
%   piv_parameters      A data structure containing all the necessary
%                       parameters required to perform the PIV processing
%
%   pass_index          A scalar integer giving the current pass of the PIV
%                       processing
%
%   frame_index         A scalar integer giving the index of the current
%                       frame being processed from the directory of images.
%                       For pyramid correlations, the velocity field will
%                       always be calculated for the frame_index and
%                       frame_index + 1 froms, ie the changing the
%                       correlation step does nothing.
%
%   window_number       A 1 x 3 vector giving the number of PIV windows to
%                       process in the first, second, and third dimensions
%                       respectively.  If the images are 2D then
%                       window_number(3) should always equal 1.
%
%   window_min          A L x 3 vector giving the minimum indices of the
%                       window domains in each respective dimension where L
%                       is the total number of windows to process.
%
%   window_max          A L x 3 vector giving the maximum indices of the
%                       window domains in each respective dimension where L
%                       is the total number of windows to process.
%
%   window_center       A L x 3 vector giving the coordinates of the center
%                       of the window domains in each respective dimension
%                       where L is the total number of windows to process.
%
%   ii_position         A M x N x P array giving the 1st dimension
%                       coordinates of the previously measured vector field
%                       if this is at least the second pass.  If this is
%                       the first pass then the coordinates correspond to
%                       an all-zero vector field.
%
%   jj_position         A M x N x P array giving the 2nd dimension
%                       coordinates of the previously measured vector field
%                       if this is at least the second pass.  If this is
%                       the first pass then the coordinates correspond to
%                       an all-zero vector field.
%
%   kk_position         A M x N x P array giving the 3rd dimension
%                       coordinates of the previously measured vector field
%                       if this is at least the second pass.  If this is
%                       the first pass then the coordinates correspond to
%                       an all-zero vector field.
%
%   ii_velocity         A M x N x P array giving the 1st dimension velocity
%                       of the previously measured vector field if this is
%                       at least the second pass.  If this is the first
%                       pass then the velocity will be set to all zeros.
%
%   jj_velocity         A M x N x P array giving the 2nd dimension velocity
%                       of the previously measured vector field if this is
%                       at least the second pass.  If this is the first
%                       pass then the velocity will be set to all zeros.
%
%   kk_velocity         A M x N x P array giving the 3rd dimension velocity
%                       of the previously measured vector field if this is
%                       at least the second pass.  If this is the first
%                       pass then the velocity will be set to all zeros.
%
% The output parameters give the measured (un-validated, un-smoothed)
% velocity field of the current frame pair and are described by
%
%   ii_position         A Q x R x S array giving the 1st dimension
%                       coordinates of the measured vector field, ie the
%                       center of the correlation window in the 1st
%                       dimension.
%
%   jj_position         A Q x R x S array giving the 2nd dimension
%                       coordinates of the measured vector field, ie the
%                       center of the correlation window in the 1st
%                       dimension.
%
%   kk_position         A Q x R x S array giving the 3rd dimension
%                       coordinates of the measured vector field, ie the
%                       center of the correlation window in the 1st
%                       dimension.
%
%   ii_velocity         A Q x R x S array giving the 1st dimension velocity
%                       of the current vector field measured at each
%                       correlation window.
%
%   jj_velocity         A Q x R x S array giving the 1st dimension velocity
%                       of the current vector field measured at each
%                       correlation window.
%
%   kk_velocity         A Q x R x S array giving the 1st dimension velocity
%                       of the current vector field measured at each
%                       correlation window.
%
% Authors: Rod La Foy
% First Created On: 14 March 2014
% Last Modified On: 14 March 2014

% This is the image filename extension
image_filename_extension=piv_parameters.general_parameters.image_filename_extension;
% This is the image variable name (if the variable is saved within a data
% object that contain multiple variables, ie a 'mat' file)
image_variable_name=piv_parameters.general_parameters.image_variable_name;

% This is a list of the images within the image read directory
image_filename_list=dir([piv_parameters.general_parameters.image_read_directory,piv_parameters.general_parameters.image_filename_prefix,'*',image_filename_extension]);

% This extracts the effective current window resolution
window_resolution=piv_parameters.pass_parameters(pass_index).window_resolution;
% This extracts the full window size
window_size=piv_parameters.pass_parameters(pass_index).window_size;

% This is the number of frames to step between correlation frame pairs
image_correlation_step=piv_parameters.general_parameters.image_correlation_step;

% This is the bulk window offset distance
bulk_window_offset=piv_parameters.pass_parameters(pass_index).bulk_window_offset;

% This is the correlation method to use for the current pass
correlation_method=piv_parameters.pass_parameters(pass_index).correlation_method;
% If the correlation method is the RPC method, this creates the RPC filter
if strcmp(correlation_method,'RPC');
    % This extracts the size of the particles in the images
    particle_diameter=piv_parameters.general_parameters.particle_diameter;
    % This creates the spectral energy filter
    spectral_filter=single(fftshift(create_spectral_filter(window_size,particle_diameter)));
elseif strcmp(correlation_method,'SCC');;
    % This sets the spectral filter to a null value
    spectral_filter=[];
end;

% This is a Boolean value stating whether to zero-mean the correlation
% windows
zero_mean_windows=piv_parameters.pass_parameters(pass_index).zero_mean_windows;

% This sets the fft library for optimal speed calculation
fftw('planner','measure');

% This determines the size of the Gaussian function to use for the window masking
gaussian_width=determine_gaussian_size_2(window_size,window_resolution);

% This creates a Gaussian windowing function
gaussian_filter=single(create_gaussian_filter(window_size,gaussian_width));

% These are the coordinate vectors
ii_position_current=reshape(window_center(:,1),window_number);
jj_position_current=reshape(window_center(:,2),window_number);
kk_position_current=reshape(window_center(:,3),window_number);

% This initializes the velocity vectors
ii_velocity_current=zeros(size(ii_position_current));
jj_velocity_current=zeros(size(ii_position_current));
kk_velocity_current=zeros(size(ii_position_current));

% This determines the domain over which the known velocity field should be
% extrapolated
ii_extrap_min=min(ii_position_current(:))-ceil(window_size(1)/2);
ii_extrap_max=max(ii_position_current(:))+ceil(window_size(1)/2);
jj_extrap_min=min(jj_position_current(:))-ceil(window_size(2)/2);
jj_extrap_max=max(jj_position_current(:))+ceil(window_size(2)/2);
kk_extrap_min=min(kk_position_current(:))-ceil(window_size(3)/2);
kk_extrap_max=max(kk_position_current(:))+ceil(window_size(3)/2);
% This adds extrapolated points to the previously measured vector field so
% that if the vector field is interpolated at points near the edge of the
% image, the interpolated points will have approximately correct values (as
% opposed to zero or NaNs)
[~,~,~,ii_velocity]=extrapolate_velocity_field(ii_position,jj_position,kk_position,ii_velocity,ii_extrap_min,ii_extrap_max,jj_extrap_min,jj_extrap_max,kk_extrap_min,kk_extrap_max);
[~,~,~,jj_velocity]=extrapolate_velocity_field(ii_position,jj_position,kk_position,jj_velocity,ii_extrap_min,ii_extrap_max,jj_extrap_min,jj_extrap_max,kk_extrap_min,kk_extrap_max);
[ii_position,jj_position,kk_position,kk_velocity]=extrapolate_velocity_field(ii_position,jj_position,kk_position,kk_velocity,ii_extrap_min,ii_extrap_max,jj_extrap_min,jj_extrap_max,kk_extrap_min,kk_extrap_max);

% If the bulk window offset is non-zero, this adds it to the previous pass
% velocity field estimate
if any(bulk_window_offset);
    % This adds the bulk window offset to the previously measured velocity
    % field
    ii_velocity=ii_velocity+bulk_window_offset(1);
    jj_velocity=jj_velocity+bulk_window_offset(2);
    kk_velocity=kk_velocity+bulk_window_offset(3);
    % This creates a Boolean value stating to apply the bulk window offset
    perform_bulk_window_offset=true;
else;
    % This creates a Boolean value stating to apply the bulk window offset
    perform_bulk_window_offset=false;
end;

% If this is the second pass or higher, this calculates the discrete window
% offset, otherwise the offset is just set to zero
if (pass_index==1)&&(not(perform_bulk_window_offset));
    % This sets the first frame minimum and maximum values to the
    % non-offset values
    frame_1_window_min=window_min;
    frame_1_window_max=window_max;
    % This sets the second frame minimum and maximum values to the
    % non-offset values
    frame_2_window_min=window_min;
    frame_2_window_max=window_max;
    % This sets the displacement offset array to all zeros
    displacement_offset=zeros(size(window_min));
else;
    % This function adds the known velocity field to the window domains, ie
    % accounts for the discrete window offset
    [frame_1_window_min,frame_1_window_max,frame_2_window_min,frame_2_window_max,displacement_offset]=dwo_window_domains(window_min,window_max,window_center,image_correlation_step,ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity);
end;

% This is the filename of the first image to load
image_1_filename=[piv_parameters.general_parameters.image_read_directory,image_filename_list(frame_index).name];
% This is the filename of the first image to load
image_2_filename=[piv_parameters.general_parameters.image_read_directory,image_filename_list(frame_index+image_correlation_step).name];

% This finds the size of the images being loaded which can be 'bmp', 'jpg',
% 'jp2', 'png', 'ppm', or 'tif' images in 2D or 'mat images in 3D
if any(strcmpi(image_filename_extension,{'bmp','jpg','jp2','png','ppm','tif'}));
    
    % This loads the first image
    I1=double(imread(image_1_filename));
    % This loads the second image
    I2=double(imread(image_2_filename));
    
elseif strcmpi(image_filename_extension,'mat');
    
    % This creates a matlab file object of the first image of the correlation pairs
    image_1_object=matfile(image_1_filename);
    % This loads the first image
    eval(['I1=double(image_1_object.',image_variable_name,');']);
    % This creates a matlab file object of the second image of the correlation pairs
    image_2_object=matfile(image_2_filename);
    % This loads the second image
    eval(['I2=double(image_2_object.',image_variable_name,');']);
    
elseif strcmpi(image_filename_extension,'dcm');
    
    % This loads the first image
    I1=double(squeeze(dicomread(image_1_filename)));
    % This loads the second image
    I2=double(squeeze(dicomread(image_2_filename)));
    
elseif strcmpi(image_filename_extension,'cdf');
    
    % This loads the first image cell array
    I1_Temp=double(cdfread(image_1_filename,'Variables',image_variable_name));
    % This loads the image from the cell array
    I1=I1_Temp{1};
    % This clears the cell array from memory
    clear('I1_Temp');
    % This loads the second image
    I2_Temp=double(cdfread(image_2_filename,'Variables',image_variable_name));
    % This loads the image from the cell array
    I2=I2_Temp{1};
    % This clears the cell array from memory
    clear('I2_Temp');
    
elseif strcmpi(image_filename_extension,'h5');
    
    % This loads the first image
    I1=double(hdf5read(image_1_filename,image_variable_name));
    % This loads the second image
    I2=double(hdf5read(image_2_filename,image_variable_name));
    
end;

% This determines the image size
image_size=zeros(1,3);
image_size(1)=size(I1,1);
image_size(2)=size(I1,2);
image_size(3)=size(I1,3);

% This iterates through the windows performing the cross-correlations
for window_index=1:prod(window_number);
    
    % This displays the calculation progress
    display_calculation_progress(window_index,1:size(window_min,1));
    
    % This extracts the correlation window from the first frame
    I1_Window=extract_correlation_window(I1,window_size,frame_1_window_min(window_index,:),frame_1_window_max(window_index,:),image_size);
    
    % This extracts the correlation window from the second frame
    I2_Window=extract_correlation_window(I2,window_size,frame_2_window_min(window_index,:),frame_2_window_max(window_index,:),image_size);
    
    % This subtracts the mean value from the windows if specified by the
    % user
    if zero_mean_windows;
        % This subtracts the mean value from the first image window
        I1_Window=I1_Window-mean(I1_Window(:));
        % This subtracts the mean value from the second image window
        I2_Window=I2_Window-mean(I2_Window(:));
    end;
    
    % This performs the cross-correlation measurement
    displacement_vector=calculate_window_correlation(I1_Window,I2_Window,gaussian_filter,spectral_filter,window_size,correlation_method);
    
    % These are the velocities of the current correlation volume
    ii_velocity_current(window_index)=displacement_vector(1)+displacement_offset(window_index,1);
    jj_velocity_current(window_index)=displacement_vector(2)+displacement_offset(window_index,2);
    kk_velocity_current(window_index)=displacement_vector(3)+displacement_offset(window_index,3);

end;

% This overwrites the input velocity field coordinates with the current
% pass's coordinates
ii_position=ii_position_current;
jj_position=jj_position_current;
kk_position=kk_position_current;
% This overwrites the input velocity field vectors with the current
% pass's vectors
ii_velocity=ii_velocity_current/image_correlation_step;
jj_velocity=jj_velocity_current/image_correlation_step;
kk_velocity=kk_velocity_current/image_correlation_step;



function [ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity]=pyramid_correlations(piv_parameters,pass_index,frame_index,window_number,window_min,window_max,window_center,ii_position,jj_position,kk_position,ii_velocity,jj_velocity,kk_velocity);
% This function calculates the pyramid correlations of the data set listed
% within the parameters variable 'piv_parameters'.  The pyramid correlation
% consist of taking a series of correlations of adjacent frames within a
% time resolved image sequence and then averaging the correlations from the
% same time seperation.  This should only be done with over-resolved data
% (in time).  The input parameters of this function are given by
%
%   piv_parameters      A data structure containing all the necessary
%                       parameters required to perform the PIV processing
%
%   pass_index          A scalar integer giving the current pass of the PIV
%                       processing
%
%   frame_index         A scalar integer giving the index of the current
%                       frame being processed from the directory of images.
%                       For pyramid correlations, the velocity field will
%                       always be calculated for the frame_index and
%                       frame_index + 1 froms, ie the changing the
%                       correlation step does nothing.
%
%   window_number       A 1 x 3 vector giving the number of PIV windows to
%                       process in the first, second, and third dimensions
%                       respectively.  If the images are 2D then
%                       window_number(3) should always equal 1.
%
%   window_min          A L x 3 vector giving the minimum indices of the
%                       window domains in each respective dimension where L
%                       is the total number of windows to process.
%
%   window_max          A L x 3 vector giving the maximum indices of the
%                       window domains in each respective dimension where L
%                       is the total number of windows to process.
%
%   window_center       A L x 3 vector giving the coordinates of the center
%                       of the window domains in each respective dimension 
%                       where L is the total number of windows to process.
%
%   ii_position         A M x N x P array giving the 1st dimension
%                       coordinates of the previously measured vector field
%                       if this is at least the second pass.  If this is
%                       the first pass then the coordinates correspond to
%                       an all-zero vector field.
%
%   jj_position         A M x N x P array giving the 2nd dimension
%                       coordinates of the previously measured vector field
%                       if this is at least the second pass.  If this is
%                       the first pass then the coordinates correspond to
%                       an all-zero vector field.
%
%   kk_position         A M x N x P array giving the 3rd dimension
%                       coordinates of the previously measured vector field
%                       if this is at least the second pass.  If this is
%                       the first pass then the coordinates correspond to
%                       an all-zero vector field.
%
%   ii_velocity         A M x N x P array giving the 1st dimension velocity
%                       of the previously measured vector field if this is
%                       at least the second pass.  If this is the first
%                       pass then the velocity will be set to all zeros.
%
%   jj_velocity         A M x N x P array giving the 2nd dimension velocity
%                       of the previously measured vector field if this is
%                       at least the second pass.  If this is the first
%                       pass then the velocity will be set to all zeros.
%
%   kk_velocity         A M x N x P array giving the 3rd dimension velocity
%                       of the previously measured vector field if this is
%                       at least the second pass.  If this is the first
%                       pass then the velocity will be set to all zeros.
%
% The output parameters give the measured (un-validated, un-smoothed)
% velocity field of the current frame pair and are described by
%
%   ii_position         A Q x R x S array giving the 1st dimension
%                       coordinates of the measured vector field, ie the
%                       center of the correlation window in the 1st
%                       dimension.
%
%   jj_position         A Q x R x S array giving the 2nd dimension
%                       coordinates of the measured vector field, ie the
%                       center of the correlation window in the 1st
%                       dimension.
%
%   kk_position         A Q x R x S array giving the 3rd dimension
%                       coordinates of the measured vector field, ie the
%                       center of the correlation window in the 1st
%                       dimension.
%
%   ii_velocity         A Q x R x S array giving the 1st dimension velocity
%                       of the current vector field measured at each 
%                       correlation window.
%
%   jj_velocity         A Q x R x S array giving the 1st dimension velocity
%                       of the current vector field measured at each 
%                       correlation window.
%
%   kk_velocity         A Q x R x S array giving the 1st dimension velocity
%                       of the current vector field measured at each 
%                       correlation window.
%
% This function currently converts the images to single precision to
% increase the computational speed.  Additionally the image deformation
% interpolation is performed using a cubic interpolation and the matlab
% function 'griddedInterpolant'.  This should ideally be changed to a sinc
% interpolation using the Blackman kernal, however the sinc interpolation
% is currently very slow.
%
% Currently ~50% of the calculation time is taken by calculating the
% cross-correlation, which cannot be sped up significantly, but other 
% issues that should be addressed to speed the code up, make it more 
% accurate, or clearer to read include:
%
%   Forming an interpolant function for the previously measured velocity
%   field to calculate 'dii_shift' et cetera instead of using the
%   functions 'interp2' and 'interp3'.
%
%   Currently the 'griddedInterpolant' function returns NaN values for
%   out-of-range inputs - this is supposedly fixed (or will be fixed) in
%   newer versions of MATLAB.  This would get around calculating the
%   variables 'out_of_range_indices_1' and 'out_of_range_indices_2'.
%
%   The interpolation operation required to calculate 'CI' can probably be
%   formulated so that the 'fftshift' function doesn't need to be performed
%   first by suitably transforming the coordinates.
%
%   The interpolation operation required to calculate 'CI' can probably be
%   rewritten as a simple matrix multiplication since the coordiantes being
%   interpolated to will always be the same.
%
%   The cases where the code is testing whether the data is 2D or 3D can
%   probably be combined/streamlined in several places - this issue is
%   primarilly due to 'interp3' not working on 2D data.
%
% Authors: Rod La Foy
% First Created On: 10 March 2014
% Last Modified On: 18 March 2014

% This is the image filename extension
image_filename_extension=piv_parameters.general_parameters.image_filename_extension;
% This is the image variable name (if the variable is saved within a data
% object that contain multiple variables, ie a 'mat' file)
image_variable_name=piv_parameters.general_parameters.image_variable_name;

% This is a list of the images within the image read directory
image_filename_list=dir([piv_parameters.general_parameters.image_read_directory,piv_parameters.general_parameters.image_filename_prefix,'*',image_filename_extension]);

% This is the optimal number of correlation frames to use in the
% pyramid correlation
pyramid_optimal_frame_number=piv_parameters.general_parameters.pyramid_optimal_frame_number;
% This is the number of levels to use in the pyramid correlation
pyramid_level_number=piv_parameters.general_parameters.pyramid_level_number;
% This is the temporary directory to store the correlation
% planes/volumes used in calculating the pyramid correlation
pyramid_correlation_write_directory=piv_parameters.general_parameters.pyramid_correlation_write_directory;

% This extracts the effective current window resolution
window_resolution=piv_parameters.pass_parameters(pass_index).window_resolution;
% This extracts the full window size
window_size=piv_parameters.pass_parameters(pass_index).window_size;

% This is the bulk window offset distance
bulk_window_offset=piv_parameters.pass_parameters(pass_index).bulk_window_offset;

% This is the correlation method to use for the current pass
correlation_method=piv_parameters.pass_parameters(pass_index).correlation_method;
% If the correlation method is the RPC method, this creates the RPC filter
if strcmp(correlation_method,'RPC');
    % This extracts the size of the particles in the images
    particle_diameter=piv_parameters.general_parameters.particle_diameter;
    % This creates the spectral energy filter
    spectral_filter=single(fftshift(create_spectral_filter(window_size,particle_diameter)));
end;

% This is a Boolean value stating whether to zero-mean the correlation
% windows
zero_mean_windows=piv_parameters.pass_parameters(pass_index).zero_mean_windows;

% This sets the fft library for optimal speed calculation
fftw('planner','measure');

% This determines the size of the Gaussian function to use for the window masking
gaussian_width=determine_gaussian_size_2(window_size,window_resolution);

% This creates a Gaussian windowing function
gaussian_filter=single(create_gaussian_filter(window_size,gaussian_width));

% These are the coordinate vectors
ii_position_current=reshape(window_center(:,1),window_number);
jj_position_current=reshape(window_center(:,2),window_number);
kk_position_current=reshape(window_center(:,3),window_number);

% This initializes the velocity vectors
ii_velocity_current=zeros(size(ii_position_current));
jj_velocity_current=zeros(size(ii_position_current));
kk_velocity_current=zeros(size(ii_position_current));

% This determines the domain over which the known velocity field should be
% extrapolated
ii_extrap_min=min(ii_position_current(:))-ceil(window_size(1)/2);
ii_extrap_max=max(ii_position_current(:))+ceil(window_size(1)/2);
jj_extrap_min=min(jj_position_current(:))-ceil(window_size(2)/2);
jj_extrap_max=max(jj_position_current(:))+ceil(window_size(2)/2);
kk_extrap_min=min(kk_position_current(:))-ceil(window_size(3)/2);
kk_extrap_max=max(kk_position_current(:))+ceil(window_size(3)/2);
% This adds extrapolated points to the previously measured vector field so
% that if the vector field is interpolated at points near the edge of the
% image, the interpolated points will have approximately correct values (as
% opposed to zero or NaNs)
[~,~,~,ii_velocity]=extrapolate_velocity_field(ii_position,jj_position,kk_position,ii_velocity,ii_extrap_min,ii_extrap_max,jj_extrap_min,jj_extrap_max,kk_extrap_min,kk_extrap_max);
[~,~,~,jj_velocity]=extrapolate_velocity_field(ii_position,jj_position,kk_position,jj_velocity,ii_extrap_min,ii_extrap_max,jj_extrap_min,jj_extrap_max,kk_extrap_min,kk_extrap_max);
[ii_position,jj_position,kk_position,kk_velocity]=extrapolate_velocity_field(ii_position,jj_position,kk_position,kk_velocity,ii_extrap_min,ii_extrap_max,jj_extrap_min,jj_extrap_max,kk_extrap_min,kk_extrap_max);

% If the bulk window offset is non-zero, this adds it to the previous pass
% velocity field estimate
if any(bulk_window_offset);
    % This adds the bulk window offset to the previously measured velocity
    % field
    ii_velocity=ii_velocity+bulk_window_offset(1);
    jj_velocity=jj_velocity+bulk_window_offset(2);
    kk_velocity=kk_velocity+bulk_window_offset(3);
    % This creates a Boolean value stating to apply the bulk window offset
    perform_bulk_window_offset=true;
else;
    % This creates a Boolean value stating to apply the bulk window offset
    perform_bulk_window_offset=false;
end;

% This initializes a cell array to store the image interpolants within
image_interpolant_array=cell(pyramid_optimal_frame_number+1,1);

% This iterates through the images that will be needed for the current
% correlation and loads them into RAM (since having to load the images
% multiple times is very slow)
for pyramid_frame_index=(frame_index-(pyramid_optimal_frame_number-1)/2):(frame_index+(pyramid_optimal_frame_number-1)/2+1);
    
    % This is the filename of the first image to load
    image_filename=[piv_parameters.general_parameters.image_read_directory,image_filename_list(pyramid_frame_index).name];
    
    % This finds the size of the images being loaded which can be 'bmp', 'jpg',
    % 'jp2', 'png', 'ppm', or 'tif' images in 2D or 'mat images in 3D
    if any(strcmpi(image_filename_extension,{'bmp','jpg','jp2','png','ppm','tif'}));
        % This loads the current image
        I=single(imread(image_filename));
    elseif strcmpi(image_filename_extension,'mat');
        % This creates a matlab file object of the current image of the correlation pairs
        image_object=matfile(image_filename);
        % This loads the current image
        eval(['I=single(image_object.',image_variable_name,');']);
    elseif strcmpi(image_filename_extension,'dcm');
        % This loads the current image
        I=single(squeeze(dicomread(image_filename)));
    elseif strcmpi(image_filename_extension,'cdf');
        % This loads the current image
        I_Temp=double(cdfread(image_filename,'Variables',image_variable_name));
        % This loads the image from the cell array
        I=I_Temp{1};
        % This clears the cell array from memory
        clear('I_Temp');
    elseif strcmpi(image_filename_extension,'h5');
        % This loads the current image
        I=double(hdf5read(image_filename,image_variable_name));
    end;
    
    % This is the index into the 4th dimension in which to save the array
    image_index=pyramid_frame_index-(frame_index-(pyramid_optimal_frame_number-1)/2)+1;
    
    % This creates an interpolant for the current image and stores it
    % within the image interpolant array
    image_interpolant_array{image_index}=griddedInterpolant(I,'cubic');
    
    % If this is the first iteration, this measures the image size
    if image_index==1;
        % This is the image size
        image_size=zeros(1,3);
        image_size(1)=size(I,1);
        image_size(2)=size(I,2);
        image_size(3)=size(I,3);
    end;
    
    % This clears the current image from memory
    clear('I');
    
end;

% This is the vector across the window domain
ii_window_coordinates_vector=(-window_size(1)/2):(window_size(1)/2-1);
jj_window_coordinates_vector=(-window_size(2)/2):(window_size(2)/2-1);
kk_window_coordinates_vector=(-window_size(3)/2):(window_size(3)/2-1);

% These are the arrays of the original window coordiantes
[ii_window_coordinates,jj_window_coordinates,kk_window_coordinates]=meshgrid(ii_window_coordinates_vector,jj_window_coordinates_vector,kk_window_coordinates_vector);

% This iterates through the windows performing the cross-correlations
for window_index=1:prod(window_number);
    
    % This displays the calculation progress
    display_calculation_progress(window_index,1:size(window_min,1));
    
    % This creates the coordinate arrays of the original undeformed image
    % windows in either 2D or 3D
    if window_size(3)==1;
        
        % These are the coordinate vectors of the window
        ii_window_vector=window_min(window_index,1):window_max(window_index,1);
        jj_window_vector=window_min(window_index,2):window_max(window_index,2);
        % These are the arrays of coordinate vectors
        [ii_window_array,jj_window_array]=ndgrid(ii_window_vector,jj_window_vector);
        
        % If this is the first pass and bulk window offsets are not being 
        % applied, this sets the pixel shift values of the current window 
        % to all zeros, otherwise the pixel shift values are interpolated
        % from the previous velocity field
        if (pass_index==1)&&(not(perform_bulk_window_offset));
            
            % These are the interpolation coordinates of the first image
            ii_window_array_interp_1=ii_window_array;
            jj_window_array_interp_1=jj_window_array;
            % These are the interpolation coordinates of the second image
            ii_window_array_interp_2=ii_window_array;
            jj_window_array_interp_2=jj_window_array;
            
        elseif pass_index>1;
        
            % This calculates the shifted pixel positions in
            % the ii and jj dimensions using the previously
            % calculated velocity field
            dii_shift=0.5*interp2(ii_position,jj_position,ii_velocity,ii_window_array,jj_window_array,'cubic',0);
            djj_shift=0.5*interp2(ii_position,jj_position,jj_velocity,ii_window_array,jj_window_array,'cubic',0);
            
        end;
  
    else;
        
        % These are the coordinate vectors of the window
        ii_window_vector=window_min(window_index,1):window_max(window_index,1);
        jj_window_vector=window_min(window_index,2):window_max(window_index,2);
        kk_window_vector=window_min(window_index,3):window_max(window_index,3);
        % These are the arrays of coordinate vectors
        [ii_window_array,jj_window_array,kk_window_array]=ndgrid(ii_window_vector,jj_window_vector,kk_window_vector);
        
        % If this is the first pass and bulk window offsets are not being 
        % applied, this sets the voxel shift values of the current window 
        % to all zeros, otherwise the voxel shift values are interpolated
        % from the previous velocity field
        if (pass_index==1)&&(not(perform_bulk_window_offset));
            
            % These are the interpolation coordinates of the first image
            ii_window_array_interp_1=ii_window_array;
            jj_window_array_interp_1=jj_window_array;
            kk_window_array_interp_1=kk_window_array;
            % These are the interpolation coordinates of the second image
            ii_window_array_interp_2=ii_window_array;
            jj_window_array_interp_2=jj_window_array;
            kk_window_array_interp_2=kk_window_array;
            
        elseif pass_index>1;
            
            % This calculates the shifted voxel positions in
            % the ii, jj, and kk dimensions using the previously
            % calculated velocity field
            dii_shift=0.5*interp3(ii_position,jj_position,kk_position,ii_velocity,ii_window_array,jj_window_array,kk_window_array,'cubic',0);
            djj_shift=0.5*interp3(ii_position,jj_position,kk_position,jj_velocity,ii_window_array,jj_window_array,kk_window_array,'cubic',0);
            dkk_shift=0.5*interp3(ii_position,jj_position,kk_position,kk_velocity,ii_window_array,jj_window_array,kk_window_array,'cubic',0);
            
        end;
        
    end;
  
    % This iterates through the different levels of the pyramid calculating the
    % correlation planes
    for pyramid_level_index=1:pyramid_level_number;
        
        % This iterates through the frames to process
        for pyramid_frame_index=(frame_index-(pyramid_optimal_frame_number-1)/2):(frame_index+(pyramid_optimal_frame_number-1)/2);
            
            % This is the index of the first and second frames
            frame_1_index=pyramid_frame_index;
            frame_2_index=pyramid_frame_index+pyramid_level_index;
            
            % This checks whether the second frame beyond the edge of the
            % pyramid correlation (which can result since the last frame +
            % pyramid_level_number may be beyond the range of pyramid)
            if (frame_index+(pyramid_optimal_frame_number-1)/2+1)<frame_2_index;
                % This skips to the next iteration
                continue;
            end;
            
            % These are the indices of the images into the full image
            % memory array to load
            image_index_1=frame_1_index-(frame_index-(pyramid_optimal_frame_number-1)/2)+1;
            image_index_2=frame_2_index-(frame_index-(pyramid_optimal_frame_number-1)/2)+1;

            % This calculates the displacement of the pixels/voxels within
            % the current window based upon the previously estimated
            % displacement
            if window_size(3)==1;
                
                % If this is the first pass and bulk window offsets are
                % being applied, then the pixel shift values are linearly
                % scaled from the bulk window offset
                if (pass_index==1)&&(perform_bulk_window_offset);
                    
                    % These are the interpolation coordinates of the first image
                    ii_window_array_interp_1=ii_window_array+0.5*bulk_window_offset(1)*pyramid_level_index;
                    jj_window_array_interp_1=jj_window_array+0.5*bulk_window_offset(2)*pyramid_level_index;
                    % These are the interpolation coordinates of the second image
                    ii_window_array_interp_2=ii_window_array-0.5*bulk_window_offset(1)*pyramid_level_index;
                    jj_window_array_interp_2=jj_window_array-0.5*bulk_window_offset(2)*pyramid_level_index;
                    
                end;  
                
                % If this is the second pass or higher, then the pixel 
                % shift values are linearly scaled from the interpolated 
                % positions from the previous velocity field
               if pass_index>1;

                    % These are the interpolation coordinates of the first image
                    ii_window_array_interp_1=ii_window_array+dii_shift*pyramid_level_index;
                    jj_window_array_interp_1=jj_window_array+djj_shift*pyramid_level_index;
                    % These are the interpolation coordinates of the second image
                    ii_window_array_interp_2=ii_window_array-dii_shift*pyramid_level_index;
                    jj_window_array_interp_2=jj_window_array-djj_shift*pyramid_level_index;

                end;
                
                % These are the indices of the first image interpolation points
                % that are out-of-range of the image
                out_of_range_indices_1=(ii_window_array_interp_1(:)<1)|(image_size(1)<ii_window_array_interp_1(:));
                out_of_range_indices_1=(out_of_range_indices_1)|(jj_window_array_interp_1(:)<1)|(image_size(2)<jj_window_array_interp_1(:));
                % These are the indices of the second image interpolation points
                % that are out-of-range of the image
                out_of_range_indices_2=(ii_window_array_interp_2(:)<1)|(image_size(1)<ii_window_array_interp_2(:));
                out_of_range_indices_2=(out_of_range_indices_2)|(jj_window_array_interp_2(:)<1)|(image_size(2)<jj_window_array_interp_2(:));
                
                % This extracts the window from the first 2D image
                I1_Window=image_interpolant_array{image_index_1}(ii_window_array_interp_1,jj_window_array_interp_1);
                % This extracts the window from the second 2D image
                I2_Window=image_interpolant_array{image_index_2}(ii_window_array_interp_2,jj_window_array_interp_2);
                
                % This sets the portion of the interpolated windows that is
                % out-of-range to zeros (in some versions of matlab the
                % interpolation will yield NaN values)
                I1_Window(out_of_range_indices_1)=0;
                I2_Window(out_of_range_indices_2)=0;
                
                % This clears the out of range indices from memory
                clear('out_of_range_indices_1','out_of_range_indices_2');
                
            else;
                
                % If this is the first pass and bulk window offsets are
                % being applied, then the voxel shift values are linearly
                % scaled from the bulk window offset
                if (pass_index==1)&&(perform_bulk_window_offset);
                    
                    % These are the interpolation coordinates of the first image
                    ii_window_array_interp_1=ii_window_array+0.5*bulk_window_offset(1)*pyramid_level_index;
                    jj_window_array_interp_1=jj_window_array+0.5*bulk_window_offset(2)*pyramid_level_index;
                    kk_window_array_interp_1=kk_window_array+0.5*bulk_window_offset(3)*pyramid_level_index;
                    % These are the interpolation coordinates of the second image
                    ii_window_array_interp_2=ii_window_array-0.5*bulk_window_offset(1)*pyramid_level_index;
                    jj_window_array_interp_2=jj_window_array-0.5*bulk_window_offset(2)*pyramid_level_index;
                    kk_window_array_interp_2=kk_window_array-0.5*bulk_window_offset(3)*pyramid_level_index;
                    
                end; 
                
                % If this is the second pass or higher, then the voxel 
                % shift values are linearly scaled from the interpolated 
                % positions from the previous velocity field
                if pass_index>1;    
                    
                    % These are the interpolation coordinates of the first image
                    ii_window_array_interp_1=ii_window_array+dii_shift*pyramid_level_index;
                    jj_window_array_interp_1=jj_window_array+djj_shift*pyramid_level_index;
                    kk_window_array_interp_1=kk_window_array+dkk_shift*pyramid_level_index;
                    % These are the interpolation coordinates of the second image
                    ii_window_array_interp_2=ii_window_array-dii_shift*pyramid_level_index;
                    jj_window_array_interp_2=jj_window_array-djj_shift*pyramid_level_index;
                    kk_window_array_interp_2=kk_window_array-dkk_shift*pyramid_level_index;
                    
                end;
                
                % These are the indices of the first image interpolation points
                % that are out-of-range of the image
                out_of_range_indices_1=(ii_window_array_interp_1(:)<1)|(image_size(1)<ii_window_array_interp_1(:));
                out_of_range_indices_1=(out_of_range_indices_1)|(jj_window_array_interp_1(:)<1)|(image_size(2)<jj_window_array_interp_1(:));
                out_of_range_indices_1=(out_of_range_indices_1)|(kk_window_array_interp_1(:)<1)|(image_size(3)<kk_window_array_interp_1(:));
                % These are the indices of the second image interpolation points
                % that are out-of-range of the image
                out_of_range_indices_2=(ii_window_array_interp_2(:)<1)|(image_size(1)<ii_window_array_interp_2(:));
                out_of_range_indices_2=(out_of_range_indices_2)|(jj_window_array_interp_2(:)<1)|(image_size(2)<jj_window_array_interp_2(:));
                out_of_range_indices_2=(out_of_range_indices_2)|(kk_window_array_interp_2(:)<1)|(image_size(3)<kk_window_array_interp_2(:));
                
                % This extracts the window from the first 3D image
                I1_Window=image_interpolant_array{image_index_1}(ii_window_array_interp_1,jj_window_array_interp_1,kk_window_array_interp_1);
                % This extracts the window from the second 3D image
                I2_Window=image_interpolant_array{image_index_2}(ii_window_array_interp_2,jj_window_array_interp_2,kk_window_array_interp_2);
                
                % This sets the portion of the interpolated windows that is
                % out-of-range to zeros (in some versions of matlab the
                % interpolation will yield NaN values)
                I1_Window(out_of_range_indices_1)=0;
                I2_Window(out_of_range_indices_2)=0;
                
                % This clears the out of range indices from memory
                clear('out_of_range_indices_1','out_of_range_indices_2');
                
            end;
            
            % This subtracts the mean value from the windows if specified by the
            % user
            if zero_mean_windows;
                % This subtracts the mean value from the first image window
                I1_Window=I1_Window-mean(I1_Window(:));
                % This subtracts the mean value from the second image window
                I2_Window=I2_Window-mean(I2_Window(:));
            end;
            
            % This calculates the correlation volumes using either the SCC or RPC
            % correlations
            if strcmp(correlation_method,'SCC');
                
                % This performs the SCC correlation on the two windows
                C=scc_correlation(I1_Window,I2_Window,gaussian_filter);
                
            elseif strcmp(correlation_method,'RPC');
                
                % This performs the RPC correlation on the two windows
                C=rpc_correlation(I1_Window,I2_Window,gaussian_filter,spectral_filter);
                
            end;
    
            % This clears the image windows from memory
            clear('I1_Window','I2_Window');
            
            % This shifts the correlation plane/volume for the
            % interpolation step (this can probably be removed if the
            % interpolation coordinates are appropriately transformed)
            C=fftshift(C);
            
            % This interpolates the correlation plane onto the new coordiantes
            if window_size(3)==1;
                
                % This is the coordinate vector to interpolate over
                ii_correlation_interp_vector=(pyramid_level_index/pyramid_optimal_frame_number)*ii_window_coordinates_vector;
                jj_correlation_interp_vector=(pyramid_level_index/pyramid_optimal_frame_number)*jj_window_coordinates_vector;
                
                % This is the array of positions to interpolate over
                [ii_correlation_interp,jj_correlation_interp]=meshgrid(ii_correlation_interp_vector,jj_correlation_interp_vector);

                % This interpolates the 2D correlation plane onto the new
                % coordinates
                CI=interp2(ii_window_coordinates,jj_window_coordinates,C,ii_correlation_interp,jj_correlation_interp,'cubic');
                
            else;
                
                % This is the coordinate vector to interpolate over
                ii_correlation_interp_vector=(pyramid_level_index/pyramid_optimal_frame_number)*ii_window_coordinates_vector;
                jj_correlation_interp_vector=(pyramid_level_index/pyramid_optimal_frame_number)*jj_window_coordinates_vector;
                kk_correlation_interp_vector=(pyramid_level_index/pyramid_optimal_frame_number)*kk_window_coordinates_vector;
                
                % This is the array of positions to interpolate over
                [ii_correlation_interp,jj_correlation_interp,kk_correlation_interp]=meshgrid(ii_correlation_interp_vector,jj_correlation_interp_vector,kk_correlation_interp_vector);

                % This interpolates the 2D correlation plane onto the new
                % coordinates
                CI=interp3(ii_window_coordinates,jj_window_coordinates,kk_window_coordinates,C,ii_correlation_interp,jj_correlation_interp,kk_correlation_interp,'cubic');
                
            end;
            
            % This clears the original correlation plane/volume from memory
            clear('C');
            
            % This shifts the interpolated correlation plane/volume back to
            % the original coordinates
            CI=fftshift(CI);
            
            % This is the filename to write the correlation plane to
            correlation_filename_write=['correlation_frame_',sprintf('%06.0f',frame_1_index),'_x_',sprintf('%06.0f',frame_2_index),'.dat'];
            
            % This saves the correlation plane/volume
            fid=fopen([pyramid_correlation_write_directory,correlation_filename_write],'w');
            fwrite(fid,CI(:),'single');
            fclose(fid);
            
        end;
        
    end;
    
    % This initializes an array to store the average correlation plane of all
    % the levels of the pyramid
    C_Pyramid_Mean=zeros(window_size,'single');
    
    % This iterates through the different levels of the pyramid averaging the
    % different levels of the pyramid
    for pyramid_level_index=1:pyramid_level_number;
        
        % This initializes an array to store the average correlation plane in each
        % level of the pyramid
        C_Level_Mean=zeros(window_size);
        
        % This iteratively loads the correlation planes from the the
        % current level of the pyramid to calculate the average correlation
        % plane
        for correlation_index=(frame_index-(pyramid_optimal_frame_number-1)/2):(frame_index+(pyramid_optimal_frame_number-1)/2);
            
            % This is the index of the first and second frames
            frame_1_index=correlation_index;
            frame_2_index=correlation_index+pyramid_level_index;
            
            % This checks whether the second frame beyond the edge of the
            % pyramid correlation (which can result since the last frame +
            % pyramid_level_number may be beyond the range of pyramid)
            if (frame_index+(pyramid_optimal_frame_number-1)/2+1)<frame_2_index;
                % This skips to the next iteration
                continue;
            end;
            
            % This is the filename to write the correlation plane to
            correlation_filename_read=['correlation_frame_',sprintf('%06.0f',frame_1_index),'_x_',sprintf('%06.0f',frame_2_index),'.dat'];
            
            % This loads the correlation plane/volume from
            % memory
            fid=fopen([pyramid_correlation_write_directory,correlation_filename_read],'r');
            CI=fread(fid,prod(window_size),'single');
            fclose(fid);
            CI=reshape(CI,window_size);
            
            % This adds the current correlation plane to the full
            % correlation mean of the current level
            C_Level_Mean=C_Level_Mean+CI;
            
        end;
        
        % This normalizes the correlation plane by the number of
        % constituent planes
        C_Level_Mean=C_Level_Mean/(pyramid_optimal_frame_number-pyramid_level_index+1);
        
        % This adds the current pyramid level mean to the total mean
        C_Pyramid_Mean=C_Pyramid_Mean+C_Level_Mean;
        
    end;
    
    % This normalizes the correlation plane by the number of computed
    % levels in the pyramid (by pyramid_level_number)
    C_Pyramid_Mean=C_Pyramid_Mean/pyramid_level_number;
 
    % This calculates the displacement vector from the correlation volumes
    displacement_vector=subpixel_displacement_vector(C_Pyramid_Mean,window_size)/pyramid_optimal_frame_number;
    
    % This calculates the interpolated velocity from the previously
    % measured velocity field
    if window_size(3)==1;
        
        % If this is the first pass and no bulk window offset is applied, 
        % this sets the interpolated velocity field to all zeros, if the
        % bulk window offset is applied, this sets the velocity field to
        % the bulk window offset, otherwise the velocity is interpolated 
        % from the previously measured velocity field
        if (pass_index==1)&&(not(perform_bulk_window_offset));
            
            % This sets the interpolated velocity vector at the
            % current window to zero
            ii_velocity_interp=0;
            jj_velocity_interp=0;
            kk_velocity_interp=0;
            
        elseif (pass_index==1)&&(perform_bulk_window_offset);
            
            % This sets the interpolated velocity vector at the
            % current window to equal the bulk window offset
            ii_velocity_interp=bulk_window_offset(1);
            jj_velocity_interp=bulk_window_offset(2);
            kk_velocity_interp=bulk_window_offset(3);
            
        elseif pass_index>1;
            
            % This interpolates the velocity vector in the ii and
            % jj dimensions at the current window from the
            % previously measured velocity field
            ii_velocity_interp=interp2(ii_position,jj_position,ii_velocity,ii_position_current(window_index),jj_position_current(window_index),'cubic',0);
            jj_velocity_interp=interp2(ii_position,jj_position,jj_velocity,ii_position_current(window_index),jj_position_current(window_index),'cubic',0);
            % This sets the velocity vector in the kk dimension to
            % a value of zero
            kk_velocity_interp=0;
            
        end;
        
    else;
        
        % If this is the first pass and no bulk window offset is applied, 
        % this sets the interpolated velocity field to all zeros, if the
        % bulk window offset is applied, this sets the velocity field to
        % the bulk window offset, otherwise the velocity is interpolated 
        % from the previously measured velocity field
        if (pass_index==1)&&(not(perform_bulk_window_offset));
            
            % This sets the interpolated velocity vector at the
            % current window to zero
            ii_velocity_interp=0;
            jj_velocity_interp=0;
            kk_velocity_interp=0;
            
        elseif (pass_index==1)&&(perform_bulk_window_offset);
            
            % This sets the interpolated velocity vector at the
            % current window to equal the bulk window offset
            ii_velocity_interp=bulk_window_offset(1);
            jj_velocity_interp=bulk_window_offset(2);
            kk_velocity_interp=bulk_window_offset(3);
            
        elseif pass_index>1;
            
            % This is the interpolated velocity measurement from the
            % previous pass at the current window position
            ii_velocity_interp=interp3(ii_position,jj_position,kk_position,ii_velocity,ii_position_current(window_index),jj_position_current(window_index),kk_position_current(window_index),'cubic',0);
            jj_velocity_interp=interp3(ii_position,jj_position,kk_position,jj_velocity,ii_position_current(window_index),jj_position_current(window_index),kk_position_current(window_index),'cubic',0);
            kk_velocity_interp=interp3(ii_position,jj_position,kk_position,kk_velocity,ii_position_current(window_index),jj_position_current(window_index),kk_position_current(window_index),'cubic',0);
            
        end;
        
    end;
    
    % This is the measured velocity at the current window
    ii_velocity_current(window_index)=displacement_vector(1)+ii_velocity_interp;
    jj_velocity_current(window_index)=displacement_vector(2)+jj_velocity_interp;
    kk_velocity_current(window_index)=displacement_vector(3)+kk_velocity_interp;

end;

% This overwrites the input velocity field coordinates with the current
% pass's coordinates
ii_position=ii_position_current;
jj_position=jj_position_current;
kk_position=kk_position_current;
% This overwrites the input velocity field vectors with the current
% pass's vectors
ii_velocity=ii_velocity_current;
jj_velocity=jj_velocity_current;
kk_velocity=kk_velocity_current;



function displacement_vector=calculate_window_correlation(I1,I2,gaussian_window_filter,spectral_filter,window_size,method);
% This performs the volumetric correlation between windows I1 and I2.  The
% gaussian filter is applied to the windows, the correlation is performed,
% and a sub-pixel fit is made to the correlation volume.

% This calculates the correlation volumes using either the SCC or RPC
% correlations
if strcmp(method,'SCC');

	% This performs the SCC correlation on the two windows
	C=scc_correlation(I1,I2,gaussian_window_filter);
    
elseif strcmp(method,'RPC');

	% This performs the RPC correlation on the two windows
	C=rpc_correlation(I1,I2,gaussian_window_filter,spectral_filter);
    
end;

% This calculates the displacement vector from the correlation volumes
displacement_vector=subpixel_displacement_vector(C,window_size);



function C=scc_correlation(I1,I2,gaussian_window_filter);
% This function performs the standard cross-correlation on the two windows 'I1' and 'I2'
% after masking the windows with a Gaussian function given by 'gaussian_window_filter'.

% This applies the Gaussian filter to the first window
I1=I1.*gaussian_window_filter;

% This applies the Gaussian filter to the second window
I2=I2.*gaussian_window_filter;

% This calculates the FFT of the first window
Y1=fftn(I1);

% This calculates the FFT of the second window
Y2=fftn(I2);

% This performs the cross-correlation in the spectral domain
Z=Y2.*conj(Y1);

% This performs the inverse FFT to return the correlation volume
C=ifftn(Z,'symmetric');



function C=rpc_correlation(I1,I2,gaussian_window_filter,spectral_filter);
% This function performs the robust phase correlation on the two windows 'I1' and 'I2'
% after masking the windows with a Gaussian function given by 'gaussian_window_filter' and
% using the RPC spectral filter defined by 'spectral_filter'.

% This applies the Gaussian filter to the first window
I1=I1.*gaussian_window_filter;

% This applies the Gaussian filter to the second window
I2=I2.*gaussian_window_filter;

% This calculates the FFT of the first window
Y1=fftn(I1);

% This calculates the FFT of the second window
Y2=fftn(I2);

% This performs the cross-correlation in the spectral domain
YC=Y2.*conj(Y1);

% This is the magnitude of the cross-correlation in the spectral domain
YC_Magnitude=sqrt(YC.*conj(YC));

% These are the indices of the non-zero magnitude components of the
% cross-correlation
non_zero_indices=(YC_Magnitude(:)~=0);

% This initializes the phase correlation array
R=YC;

% THis is the phase correlation in the spectral domain
R(non_zero_indices)=YC(non_zero_indices)./YC_Magnitude(non_zero_indices);

% This applies the spectral filter
R=R.*spectral_filter;

% This performs the inverse FFT to return the correlation volume
C=ifftn(R,'symmetric');

% This takes the absolute value of the correlation volume to ensure positive values
C=abs(C);



function displacement_vector=subpixel_displacement_vector(C,window_size);
% This function finds the subpixel fit to the correlation volume given by C.

% This is the linear index of the correlation maximum (this only takes the first maximum
% value if multiple are found - this is unlikely, but should be addressed)
max_linear_index=find(C==max(C(:)),1);

% This converts the linear index to a set of subscripted indices
[ii_max,jj_max,kk_max]=ind2sub(window_size,max_linear_index);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This performs the three-point Gaussian fitting in the ii-direction.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are the three points to fit in the ii-direction
if (ii_max>1)&&(ii_max<window_size(1));
	% These are the three points to fit if the maximum is not on an edge
	C_ii_n1=C(ii_max-1,jj_max,kk_max);
	C_ii_00=C(ii_max,jj_max,kk_max);
	C_ii_p1=C(ii_max+1,jj_max,kk_max);
elseif (ii_max==1);
	% These are the three points to fit if the maximum is on the lower edge
	C_ii_n1=C(window_size(1),jj_max,kk_max);
	C_ii_00=C(ii_max,jj_max,kk_max);
	C_ii_p1=C(ii_max+1,jj_max,kk_max);
elseif (ii_max==window_size(1));
	% These are the three points to fit if the maximum is on the upper edge
	C_ii_n1=C(ii_max-1,jj_max,kk_max);
	C_ii_00=C(ii_max,jj_max,kk_max);
	C_ii_p1=C(1,jj_max,kk_max);
end;

% This is the numerator of the subpixel three-point Gaussian fit in the ii direction
ii_num_fit=log(C_ii_n1)-log(C_ii_p1);
% This is the denominator of the subpixel three-point Gaussian fit in the ii direction
ii_den_fit=2*log(C_ii_n1)-4*log(C_ii_00)+2*log(C_ii_p1);
% This is the subpixel three-point Gaussian fit in the ii direction
ii_fit=ii_max+ii_num_fit/ii_den_fit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This performs the three-point Gaussian fitting in the jj-direction.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are the three points to fit in the jj-direction
if (jj_max>1)&&(jj_max<window_size(2));
	% These are the three points to fit if the maximum is not on an edge
	C_jj_n1=C(ii_max,jj_max-1,kk_max);
	C_jj_00=C(ii_max,jj_max,kk_max);
	C_jj_p1=C(ii_max,jj_max+1,kk_max);
elseif (jj_max==1);
	% These are the three points to fit if the maximum is on the lower edge
	C_jj_n1=C(ii_max,window_size(2),kk_max);
	C_jj_00=C(ii_max,jj_max,kk_max);
	C_jj_p1=C(ii_max,jj_max+1,kk_max);
elseif (jj_max==window_size(2));
	% These are the three points to fit if the maximum is on the upper edge
	C_jj_n1=C(ii_max,jj_max-1,kk_max);
	C_jj_00=C(ii_max,jj_max,kk_max);
	C_jj_p1=C(ii_max,1,kk_max);
end;

% This is the numerator of the subpixel three-point Gaussian fit in the jj direction
jj_num_fit=log(C_jj_n1)-log(C_jj_p1);
% This is the denominator of the subpixel three-point Gaussian fit in the jj direction
jj_den_fit=2*log(C_jj_n1)-4*log(C_jj_00)+2*log(C_jj_p1);
% This is the subpixel three-point Gaussian fit in the jj direction
jj_fit=jj_max+jj_num_fit/jj_den_fit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This performs the three-point Gaussian fitting in the kk-direction.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This checks whether the array is 3D and if not, sets the sub-pixel fit
% value to 1 (which transforms to zero velocity)
if size(C,3)==1;
    
    % This sets the subpixel fit value in the 3rd dimension to one
    kk_fit=1;
    
else;

    % These are the three points to fit in the kk-direction
    if (kk_max>1)&&(kk_max<window_size(3));
        % These are the three points to fit if the maximum is not on an edge
        C_kk_n1=C(ii_max,jj_max,kk_max-1);
        C_kk_00=C(ii_max,jj_max,kk_max);
        C_kk_p1=C(ii_max,jj_max,kk_max+1);
    elseif (kk_max==1);
        % These are the three points to fit if the maximum is on the lower edge
        C_kk_n1=C(ii_max,jj_max,window_size(3));
        C_kk_00=C(ii_max,jj_max,kk_max);
        C_kk_p1=C(ii_max,jj_max,kk_max+1);
    elseif (kk_max==window_size(3));
        % These are the three points to fit if the maximum is on the upper edge
        C_kk_n1=C(ii_max,jj_max,kk_max-1);
        C_kk_00=C(ii_max,jj_max,kk_max);
        C_kk_p1=C(ii_max,jj_max,1);
    end;
    
    % This is the numerator of the subpixel three-point Gaussian fit in the kk direction
    kk_num_fit=log(C_kk_n1)-log(C_kk_p1);
    % This is the denominator of the subpixel three-point Gaussian fit in the kk direction
    kk_den_fit=2*log(C_kk_n1)-4*log(C_kk_00)+2*log(C_kk_p1);
    % This is the subpixel three-point Gaussian fit in the kk direction
    kk_fit=kk_max+kk_num_fit/kk_den_fit;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This calculates the displacement vector from the sub-pixel fits.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This accounts for the origin being located at the (1,1,1) index
ii_fit=ii_fit-1;
jj_fit=jj_fit-1;
kk_fit=kk_fit-1;

% This accounts for the periodicity of the correlation volume for the ii displacement
if ii_fit>(window_size(1)/2);
	ii_fit=ii_fit-window_size(1);
end;

% This accounts for the periodicity of the correlation volume for the jj displacement
if jj_fit>(window_size(2)/2);
	jj_fit=jj_fit-window_size(2);
end;

% This accounts for the periodicity of the correlation volume for the kk displacement
if kk_fit>(window_size(3)/2);
	kk_fit=kk_fit-window_size(3);
end;

% This is the displacement vector of the current window
displacement_vector=[ii_fit,jj_fit,kk_fit];



function I_Window=extract_correlation_window(I,window_size,window_min,window_max,image_size);
% This function extracts the current windows from the frame 'I' where 
% 'window_size' is the size of the correlation window, 'window_min' and 
% 'window_max' are the 1 x 3 vectors of the window minimum and maximum 
% indices, and 'image_size' is the the size of the frame 'I'.

% These are the coordinate domains to extract
ii_min=window_min(1);
ii_max=window_max(1);
jj_min=window_min(2);
jj_max=window_max(2);
kk_min=window_min(3);
kk_max=window_max(3);

% This checks if the ii indices are outside the domain of the image and sets them to the
% domain of the image if true
if ii_min<1;
	ii_min=1;
end;
if ii_max>image_size(1);
	ii_max=image_size(1);
end;
% This checks if the jj indices are outside the domain of the image and sets them to the
% domain of the image if true
if jj_min<1;
	jj_min=1;
end;
if jj_max>image_size(2);
	jj_max=image_size(2);
end;
% This checks if the kk indices are outside the domain of the image and sets them to the
% domain of the image if true
if kk_min<1;
	kk_min=1;
end;
if kk_max>image_size(3);
	kk_max=image_size(3);
end;

% This extracts the current window from the frame
I_Window=I(ii_min:ii_max,jj_min:jj_max,kk_min:kk_max);

% This initilizes the image size vector
I_Window_Size=zeros(1,3);
% This is the size of the extracted image
I_Window_Size(1)=size(I_Window,1);
I_Window_Size(2)=size(I_Window,2);
I_Window_Size(3)=size(I_Window,3);

% If the extracted image is smaller than the full window size (due to the window being
% on the edge of the full image) this pads the image with zeros
if any(I_Window_Size~=window_size);
	% This initializes an array of zeros the size of the full window
	I_Window_Temp=zeros(window_size);
	% These are the indices to insert the image into the zero padded image (the
	% indices are centered in the window so that information is not lost during the
	% masking process)
	window_insert_min=floor((window_size-I_Window_Size)/2)+1;   	
	window_insert_max=window_insert_min+I_Window_Size-1;
	% This copies the current image from the frame into the window
	I_Window_Temp(window_insert_min(1):window_insert_max(1),window_insert_min(2):window_insert_max(2),window_insert_min(3):window_insert_max(3))=I_Window;
	% This renames the variable so that the zero-padded image is used
	I_Window=I_Window_Temp;
end;



function outlier_vector_array=locate_vector_outliers(U,V,W,minimum_outlier_threshhold,maximum_outlier_threshhold,kernal_size,epsilon,residual_threshhold);
% This function uses several metrics to determine the location of outlier vectors from the
% vector field defined by 'U', 'V', and 'W'.

% This determines the outliers of the vector field based upon the velocity threshholds
threshhold_outlier_vector_array=locate_threshhold_outliers(U,V,W,minimum_outlier_threshhold,maximum_outlier_threshhold);

% This determines the outliers of the vector field based upon the UOD
uod_outlier_vector_array=locate_uod_outliers(U,V,W,kernal_size,epsilon,residual_threshhold);

% This determines the NaN values of the vector field (which may be
% considered to be equivalent to outliers)
nan_outlier_vector_array=locate_nan_outliers(U,V,W);

% This is the array of outliers that failed any of the outlier tests
outlier_vector_array=(threshhold_outlier_vector_array)|(uod_outlier_vector_array)|(nan_outlier_vector_array);



function outlier_vector_array=locate_threshhold_outliers(U,V,W,minimum_outlier_threshhold,maximum_outlier_threshhold);
% This function identifies outlier vectors from the vector fields given by 'U', 'V', and
% 'W' using the threshhold values defined by 'minimum_outlier_threshhold' and 
% 'maximum_outlier_threshhold' where theses are 1 x 3 vectors.

% This initializes the array of outlier vectors
outlier_vector_array=false(size(U));

% These are the vectors that are less then the minimum outlier threshhold
outlier_vector_array=(outlier_vector_array)|(U<minimum_outlier_threshhold(1));
outlier_vector_array=(outlier_vector_array)|(V<minimum_outlier_threshhold(2));
outlier_vector_array=(outlier_vector_array)|(W<minimum_outlier_threshhold(3));

% These are the vectors that are greater then the maximum outlier threshhold
outlier_vector_array=(outlier_vector_array)|(U>maximum_outlier_threshhold(1));
outlier_vector_array=(outlier_vector_array)|(V>maximum_outlier_threshhold(2));
outlier_vector_array=(outlier_vector_array)|(W>maximum_outlier_threshhold(3));



function outlier_vector_array=locate_uod_outliers(U,V,W,kernal_size,epsilon,residual_threshhold);
% This function identifies outlier vectors of the vector fields given by 'U', 'V', and
% 'W' using a kernal size (ie the size of the block that the median is calculated over)
% given by 'kernal_size' and using a normalization level (ie the expected uncertainty in
% the PIV measurement - typically about 0.1 pixel for 2D PIV) of 'epsilon'.  The
% 'kernal_size' is a 1 x 3 vector which must have all odd integer values.  The argument
% 'residual_threshhold' is the cutoff on the UOD threshhold above which a vector is
% determined to be an outlier.

% This is the size of the vector field
[ii_vector_num,jj_vector_num,kk_vector_num]=size(U);

% This initializes the array of outlier vectors
outlier_vector_array=false(ii_vector_num,jj_vector_num,kk_vector_num);

% This iterates through the U vector extracting the kernals over which to calculate the
% median test
for ii_vector=1:ii_vector_num;
	for jj_vector=1:jj_vector_num;
		for kk_vector=1:kk_vector_num;
		
			% This is the range of indices to extract
			ii_min=ii_vector-(kernal_size(1)-1)/2;
			ii_max=ii_vector+(kernal_size(1)-1)/2;
			jj_min=jj_vector-(kernal_size(2)-1)/2;
			jj_max=jj_vector+(kernal_size(2)-1)/2;
			kk_min=kk_vector-(kernal_size(3)-1)/2;
			kk_max=kk_vector+(kernal_size(3)-1)/2;
			
			% This checks if the ii indices are outside the domain of the vector field 
			% and sets them to the domain of the vector field if true
			if ii_min<1;
				ii_min=1;
			end;
			if ii_max>ii_vector_num;
				ii_max=ii_vector_num;
			end;
			% This checks if the jj indices are outside the domain of the vector field 
			% and sets them to the domain of the vector field if true
			if jj_min<1;
				jj_min=1;
			end;
			if jj_max>jj_vector_num;
				jj_max=jj_vector_num;
			end;
			% This checks if the kk indices are outside the domain of the vector field 
			% and sets them to the domain of the vector field if true
			if kk_min<1;
				kk_min=1;
			end;
			if kk_max>kk_vector_num;
				kk_max=kk_vector_num;
			end;
			
			% This extracts the U vectors of the current block
			U_Block=U(ii_min:ii_max,jj_min:jj_max,kk_min:kk_max);
			% This extracts the V vectors of the current block
			V_Block=V(ii_min:ii_max,jj_min:jj_max,kk_min:kk_max);
			% This extracts the W vectors of the current block
			W_Block=W(ii_min:ii_max,jj_min:jj_max,kk_min:kk_max);
			
			% This is the linear index of the center vector
			center_index=sub2ind(size(U_Block),ii_vector-ii_min+1,jj_vector-jj_min+1,kk_vector-kk_min+1);
			
			% This linearizes the U vectors of the current block
			U_Block=U_Block(:);
			% This linearizes the V vectors of the current block
			V_Block=V_Block(:);
			% This linearizes the W vectors of the current block
			W_Block=W_Block(:);
			
			% This is the current U vector
			U_Vector=U_Block(center_index);
			% This is the current V vector
			V_Vector=V_Block(center_index);
			% This is the current W vector
			W_Vector=W_Block(center_index);
			
			% This removes the center vector from the U block
			U_Block(center_index)=[];
			% This removes the center vector from the V block
			V_Block(center_index)=[];
			% This removes the center vector from the W block
			W_Block(center_index)=[];
			
			% This is the median of the U block
			U_Median=median(U_Block);
			% This is the median of the V block
			V_Median=median(V_Block);
			% This is the median of the W block
			W_Median=median(W_Block);
			
			% This is the residual of the U block
			U_Residual=abs(U_Block-U_Median);
			% This is the residual of the V block
			V_Residual=abs(V_Block-V_Median);
			% This is the residual of the W block
			W_Residual=abs(W_Block-W_Median);
			
			% This is the median of the residual of the U block
			U_Residual_Median=median(U_Residual);
			% This is the median of the residual of the V block
			V_Residual_Median=median(V_Residual);
			% This is the median of the residual of the W block
			W_Residual_Median=median(W_Residual);
			
			% This is the UOD residual of the U vector
			U_UOD=abs(U_Vector-U_Median)/(U_Residual_Median+epsilon);
			% This is the UOD residual of the V vector
			V_UOD=abs(V_Vector-V_Median)/(V_Residual_Median+epsilon);
			% This is the UOD residual of the W vector
			W_UOD=abs(W_Vector-W_Median)/(W_Residual_Median+epsilon);
			
			% This adds the current vector to the total outlier vector array
			outlier_vector_array(ii_vector,jj_vector,kk_vector)=(U_UOD>residual_threshhold)|(V_UOD>residual_threshhold)|(W_UOD>residual_threshhold);
			
		end;
	end;
end;



function outlier_vector_array=locate_nan_outliers(U,V,W);
% This function identifies outlier vectors from the vector fields given by
% 'U', 'V', and 'W' using the criterion that any NaN vectors are 
% essentially equivalent to outliers.  This criterion should rarely occur, 
% but may happen in sparsely seeded synthetic images in regions of all 
% zeros.

% This initializes the array of outlier vectors
outlier_vector_array=false(size(U));

% These are the vectors that have NaN values
outlier_vector_array=(outlier_vector_array)|(isnan(U));
outlier_vector_array=(outlier_vector_array)|(isnan(V));
outlier_vector_array=(outlier_vector_array)|(isnan(W));
			
			

function [U,V,W]=local_mean_replacement(X,Y,Z,U,V,W,outlier_vector_array,minimum_valid_vector_number);
% This function replaces the bad vectors flagged by the Boolean array
% 'outlier_vector_array' using the good vectors described by the
% coordinate arrays 'X', 'Y', and 'Z' and the displacement arrays 'U',
% 'V', and 'W'.  The outlier vectors are replaced using weighted average
% of the surrounding data.  The scalar integer
% 'minimum_valid_vector_number' is the minimum number of adjacent valid
% vectors required to perform the local mean calculation.

% This is the number of windows (or vectors) in each dimension
window_number=zeros(1,3);
window_number(1)=size(X,1);
window_number(2)=size(X,2);
window_number(3)=size(X,3);

% This is a vector of linear indices corresponding to the locations of the outlier vectors
outlier_indices=find(outlier_vector_array);

% This initializes the vector fields that will contain the replaced outlier vectors
U_Interp=U;
V_Interp=V;
W_Interp=W;

% These are the the median values of the vector fields
U_Median=median(U(:));
V_Median=median(V(:));
W_Median=median(W(:));

% This iterates through the outliers identifying those that are surrounded by good vectors
% and replacing them using trilinear interpolation
for n=1:length(outlier_indices);
    
    % These are the subscripted indices of the current outlier vector
    [ii_outlier,jj_outlier,kk_outlier]=ind2sub(window_number,outlier_indices(n));
    
    % This is the radius of the block (ignoring the central point, ie a block radius of 1
    % will return a 3 x 3 x 3 block)
    block_radius=1;
    
    % This incrementally increases the block size until a sufficient number of valid
    % vectors are located
    while true;
        
        % These are the indices of the block of vectors immediately surrounding the outlier
        ii_block_min=max([1,ii_outlier-block_radius]);
        ii_block_max=min([ii_outlier+block_radius,window_number(1)]);
        jj_block_min=max([1,jj_outlier-block_radius]);
        jj_block_max=min([jj_outlier+block_radius,window_number(2)]);
        kk_block_min=max([1,kk_outlier-block_radius]);
        kk_block_max=min([kk_outlier+block_radius,window_number(3)]);
        
        % This is the block of the outlier vector identification array
        outlier_vector_block=outlier_vector_array(ii_block_min:ii_block_max,jj_block_min:jj_block_max,kk_block_min:kk_block_max);
        % This is the number of valid vectors in the current block
        valid_vector_number=sum(not(outlier_vector_block(:)));
        
        % If at least the minimum number of valid vectors is identified, then the outlier is
        % replaced using a weighted average
        if valid_vector_number>=minimum_valid_vector_number;
            
            % This extracts the current block of vector values
            U_Block=U(ii_block_min:ii_block_max,jj_block_min:jj_block_max,kk_block_min:kk_block_max);
            V_Block=V(ii_block_min:ii_block_max,jj_block_min:jj_block_max,kk_block_min:kk_block_max);
            W_Block=W(ii_block_min:ii_block_max,jj_block_min:jj_block_max,kk_block_min:kk_block_max);
            
            % This sets the outliers to NaN values in the vector field blocks
            U_Block(outlier_vector_block(:))=NaN;
            V_Block(outlier_vector_block(:))=NaN;
            W_Block(outlier_vector_block(:))=NaN;
            
            % This extracts the current block of coordinates
            X_Block=X(ii_block_min:ii_block_max,jj_block_min:jj_block_max,kk_block_min:kk_block_max);
            Y_Block=Y(ii_block_min:ii_block_max,jj_block_min:jj_block_max,kk_block_min:kk_block_max);
            Z_Block=Z(ii_block_min:ii_block_max,jj_block_min:jj_block_max,kk_block_min:kk_block_max);
            
            % This is the center coordinate
            X0=X(ii_outlier);
            Y0=Y(jj_outlier);
            Z0=Z(kk_outlier);
            
            % This is the distance vector to the outlier point
            rho=sqrt((X_Block-X0).^2+(Y_Block-Y0).^2+(Z_Block-Z0).^2);
            
            % This calculates the interpolated vector values
            U_Interp_Temp=nansum(U_Block(:).*rho(:))/sum(rho(:));
            V_Interp_Temp=nansum(V_Block(:).*rho(:))/sum(rho(:));
            W_Interp_Temp=nansum(W_Block(:).*rho(:))/sum(rho(:));
            
            % This breaks the loop
            break;
            
        else;
            
            % This increments the block radius
            block_radius=block_radius+1;
            
            % This checks if the block radius is greater than half the size of the vector field
            % and if so just assigns a value equal to the mean of the vector field
            if block_radius>=max(window_number/2);
                
                % This sets the interpolated value to the median value of the vector field
                U_Interp_Temp=U_Median;
                V_Interp_Temp=V_Median;
                W_Interp_Temp=W_Median;
                
                % This breaks the loop
                break;
                
            end;
            
        end;
        
    end;
    
    % This sets the current outlier to the interpolated value
    U_Interp(ii_outlier,jj_outlier,kk_outlier)=U_Interp_Temp;
    V_Interp(ii_outlier,jj_outlier,kk_outlier)=V_Interp_Temp;
    W_Interp(ii_outlier,jj_outlier,kk_outlier)=W_Interp_Temp;
    
end;

% This copies the interpolated arrays back into the original variables
% (primarily for consistancy with the other vector replacement functions)
U=U_Interp;
V=V_Interp;
W=W_Interp;



function [U,V,W]=laplacian_replacement(X,Y,Z,U,V,W,outlier_vector_array,adjacent_connectivity);
% This function uses Laplacian interploation to replace outlier vectors
% within a measured velocity field.  The input floating point velocity
% vector arrays 'U', 'V', and 'W' may be either 2D or 3D in size, but must
% all have the same dimensions.  The logical array 'outlier_vector_array'
% is of the same size as 'U', 'V', and 'W' and contains true values (id est
% 1s) at the location of the outlier vectors.  The scalar argument
% 'adjacent_connectivity' refers to the number of adjacent points that are
% used to calculate the interpolated velocity values.  If the input arrays
% are 2D, then 'adjacent_connectivity' may be equal to either 4 or 8
% referring to the nearest adjacent points in a 2D rectilinear array.  If
% the input arrays are 3D, then 'adjacent_connectivity' may be equal to 6,
% 18, or 26 referring to the nearest adjacent points in a 3D rectilinear
% array.

% This calculates the spacing of the grid points in the X dimension
dX_Spacing=abs(X(2,1,1)-X(1,1,1));
% This calculates the spacing of the grid points in the Y dimension
dY_Spacing=abs(Y(1,2,1)-Y(1,1,1));
% This calculates the spacing of the grid points in the Z dimension if the
% array is 3D
if size(X,3)>1;
    % This calculates the spacing of the grid points in the Z dimension
    dZ_Spacing=abs(Z(1,1,2)-Z(1,1,1));
else;
    % If the array is 2D, this sets the spacing to 1 (although the value
    % effectively is never used)
    dZ_Spacing=1;
end;

% This calculates the weighting matrices for the unknown (outlier) and
% known points to be later used in performing the Laplacian interpolation
[A,B]=calculate_laplacian_matrices(outlier_vector_array,dY_Spacing,dX_Spacing,dZ_Spacing,adjacent_connectivity);

% This extracts a vector of the valid U vector components
U_Valid=U(not(outlier_vector_array(:)));
% This calculates the interpolated values of the outliers of U
U_Interp=A\(B*U_Valid);
% This replaces the outlier U components
U(outlier_vector_array(:))=U_Interp;

% This extracts a vector of the valid V vector components
V_Valid=V(not(outlier_vector_array(:)));
% This calculates the interpolated values of the outliers of V
V_Interp=A\(B*V_Valid);
% This replaces the outlier V components
V(outlier_vector_array(:))=V_Interp;

% This extracts a vector of the valid W vector components
W_Valid=W(not(outlier_vector_array(:)));
% This calculates the interpolated values of the outliers of W
W_Interp=A\(B*W_Valid);
% This replaces the outlier W components
W(outlier_vector_array(:))=W_Interp;



function [A,B]=calculate_laplacian_matrices(replacement_array,dii_spacing,djj_spacing,dkk_spacing,adjacent_connectivity);
% This function uses Laplacian interpolation to replace the components in
% the floating point precision array 'function_array' which can have at most three
% dimensions where the component values to be replaced are specified by
% the true values of the logical binary array 'replacement_array' which
% must have the same size as the array 'function_array'.  The argument
% 'adjacent_connectivity' is a scalar giving the connected region
% surrounding each outlier point over which the outlier value is
% calculated.  For 2D arrays, 'adjacent_connectivity' can be equal to
% either 4 or 8.  For 3D arrays, 'adjacent_connectivity' can be equal to 6,
% 18, or 26.

% This is the size of the scalar function 'function_array'
[ii_array_size,jj_array_size,kk_array_size]=size(replacement_array);

% This tests whether a valid input for 'adjacent_connectivity' was input
if kk_array_size==1;
    % This tests whether 'adjacent_connectivity' equals either 4 or 8 and
    % if not displays an error
    if (adjacent_connectivity~=4)&&(adjacent_connectivity~=8);
        % This displays an error stating that the value of
        % 'adjacent_connectivity' is incorrect for the current array size
        error('The ''function_array'' is 2D, however the ''adjacent_connectivity'' variable is not set to either a 4-connected or an 8-connected region.');
    end;
elseif kk_array_size>1;
    % This tests whether 'adjacent_connectivity' equals either 6, 18, or
    % 26 and if not displays an error
    if (adjacent_connectivity~=6)&&(adjacent_connectivity~=18)&&(adjacent_connectivity~=26);
        % This displays an error stating that the value of
        % 'adjacent_connectivity' is incorrect for the current array size
        error('The ''function_array'' is 3D, however the ''adjacent_connectivity'' variable is not set to a 6-connected, an 18-connected, or a 26-connected region.');
    end;
end;

% This is the number of the components of the function 'f' that need to be
% replaced
replacement_number=sum(replacement_array(:));
% This is the number of the components of the function 'f' that do not need
% to be replaced
retainment_number=numel(replacement_array)-replacement_number;

% This creates an array that gives the number of components to be replaced
% at or below the current linear index
linear_replacement_cumulative_sum=reshape(cumsum(replacement_array(:)),[ii_array_size,jj_array_size,kk_array_size]);
% This creates an array that gives the number of components not to be
% replaced at or below the current linear index
linear_retainment_cumulative_sum=reshape(cumsum(not(replacement_array(:))),[ii_array_size,jj_array_size,kk_array_size]);

% These are the indices of the outlier vectors
[ii_replace_indices,jj_replace_indices,kk_replace_indices]=ind2sub([ii_array_size,jj_array_size,kk_array_size],find(replacement_array));

% This creates the list of adjacent coordinate indices based upon the
% defined point connectivity
if adjacent_connectivity==4;
    % These are the differences in the distance to the 4 adjacent points in
    % two dimensions
    dii_adjacent=[-1,+1,00,00];
    djj_adjacent=[00,00,-1,+1];
    dkk_adjacent=[00,00,00,00];
elseif adjacent_connectivity==8;
    % These are the differences in the distance to the 8 adjacent points in
    % two dimensions
    dii_adjacent=[-1,-1,-1,00,00,+1,+1,+1];
    djj_adjacent=[-1,00,+1,-1,+1,-1,00,+1];
    dkk_adjacent=[00,00,00,00,00,00,00,00];
elseif adjacent_connectivity==6;
    % These are the differences in the distance to the 6 adjacent points in
    % three dimensions
    dii_adjacent=[-1,+1,00,00,00,00];
    djj_adjacent=[00,00,-1,+1,00,00];
    dkk_adjacent=[00,00,00,00,-1,+1];
elseif adjacent_connectivity==18;
    % These are the differences in the distance to the 18 adjacent points in
    % three dimensions
    dii_adjacent=[-1,00,00,00,+1,-1,-1,-1,00,00,+1,+1,+1,-1,00,00,00,+1];
    djj_adjacent=[00,-1,00,+1,00,-1,00,+1,-1,+1,-1,00,+1,00,-1,00,+1,00];
    dkk_adjacent=[-1,-1,-1,-1,-1,00,00,00,00,00,00,00,00,+1,+1,+1,+1,+1];
elseif adjacent_connectivity==26;
    % These are the differences in the distance to the 26 adjacent points in
    % three dimensions
    dii_adjacent=[-1,-1,-1,00,00,00,+1,+1,+1,-1,-1,-1,00,00,+1,+1,+1,-1,-1,-1,00,00,00,+1,+1,+1];
    djj_adjacent=[-1,00,+1,-1,00,+1,-1,00,+1,-1,00,+1,-1,+1,-1,00,+1,-1,00,+1,-1,00,+1,-1,00,+1];
    dkk_adjacent=[-1,-1,-1,-1,-1,-1,-1,-1,-1,00,00,00,00,00,00,00,00,+1,+1,+1,+1,+1,+1,+1,+1,+1];
end;

% This calculates the weighting vector for the adjacent points
weighting_vector=(1./sqrt((dii_spacing*dii_adjacent).^2+(djj_spacing*djj_adjacent).^2+(dkk_spacing*dkk_adjacent).^2))';
% This normalizes the weighting vector so that it sums to 1
weighting_vector=weighting_vector/sum(weighting_vector);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates the Weighting Matrices Relating the Known and Unknown Values    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the ratio of value to be replaced to the total number of function
% values
replacement_ratio=replacement_number/numel(replacement_array);
% This is the approximate number of non-zero values in the replacement
% weighting matrix
approximate_replacement_nonzero_number=round(adjacent_connectivity*replacement_ratio*replacement_number);
% This ensures that the replacement matrix has at least enough entries for
% each functional value to be replaced
if approximate_replacement_nonzero_number<replacement_number;
    % This sets the approximate nonzero replacement number equal to the
    % total number of functional values to be replaced
    approximate_replacement_nonzero_number=replacement_number;
end;

% This is the size of the replacement weighting matrix
A_ii_size=replacement_number;
A_jj_size=replacement_number;
% This initializes the replacement weighting matrix index vectors
A_ii=(1:replacement_number)';
A_jj=(1:replacement_number)';
% This initializes the replacement weighting matrix values
A_value=zeros(approximate_replacement_nonzero_number,1);
% This sets the values corresponding with the functional values to be
% replaced equal to one
A_value(1:replacement_number)=1;

% This initializes an index into the replacement vectors corresponding with
% the last component calculated (in this case none have been calculated so
% it is zero)
A_Last_Index=replacement_number;

% This is the approximate number of non-zero values in the retainment
% weighting matrix
approximate_retainment_nonzero_number=round(adjacent_connectivity*(1-replacement_ratio)*replacement_number);

% This is the size of the retainment weighting matrix
B_ii_size=replacement_number;
B_jj_size=retainment_number;
% This initializes the retainment weighting matrix index vectors
B_ii=zeros(approximate_retainment_nonzero_number,1);
B_jj=zeros(approximate_retainment_nonzero_number,1);
% This initializes the replacement weighting matrix values
B_value=zeros(approximate_retainment_nonzero_number,1);

% This initializes an index into the retainment vectors corresponding with
% the last component calculated (in this case none have been calculated so
% it is zero)
B_Last_Index=0;

% This iterates through the list of components of the function 'f' to be
% replaced creating sparse weighting arrays to solve for the replaced
% components
for replace_component_index=1:replacement_number;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculates the Adjacent Components in the Array                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % These are the subscripted indices of the current outlier
    ii_replace_index_current=ii_replace_indices(replace_component_index);
    jj_replace_index_current=jj_replace_indices(replace_component_index);
    kk_replace_index_current=kk_replace_indices(replace_component_index);
    
    % These are the subscipted indices of the nearest points to the
    % current component to be replaced
    ii_adjacent_indices=ii_replace_index_current+dii_adjacent;
    jj_adjacent_indices=jj_replace_index_current+djj_adjacent;
    kk_adjacent_indices=kk_replace_index_current+dkk_adjacent;
    
    % These are the indices that are in range of the current array
    in_range_indices=(1<=ii_adjacent_indices)&(ii_adjacent_indices<=ii_array_size)&(1<=jj_adjacent_indices)&(jj_adjacent_indices<=jj_array_size)&(1<=kk_adjacent_indices)&(kk_adjacent_indices<=kk_array_size);
    % This converts the subscripted indices that are in range into linear
    % indices (the 'sub2ind' function is not used since it is slower
    % and the previous command checks for out of range indices)
    adjacent_linear_indices=ii_array_size*jj_array_size*(kk_adjacent_indices(in_range_indices)-1)+ii_array_size*(jj_adjacent_indices(in_range_indices)-1)+ii_adjacent_indices(in_range_indices);
    
    % Since some of the adjacent ponts may be out of range, this is a
    % factor to rescale the weighting vector by to ensure that it still
    % sums to 1
    weight_scaling_factor=length(weighting_vector)/sum(in_range_indices);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculates Replacement Matrix Components                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This is a logical vector giving the members of the possible adjacent
    % points that need to be replaced
    current_replacement_points=replacement_array(adjacent_linear_indices);
    
    % These are the linear indices of the adjacent points that need to be
    % replaced (some of which may be zero - these are adjacent points that
    % are to be retained)
    %adjacent_replacement_linear_indices=replacement_array(adjacent_linear_indices).*adjacent_linear_indices;
    adjacent_replacement_linear_indices=current_replacement_points.*adjacent_linear_indices;
    % This removes the zero-valued components of the adjacent replacement
    % linear index vector
    adjacent_replacement_linear_indices(adjacent_replacement_linear_indices==0)=[];
    % This is the number of adjacent componets that need to be replaced
    adjacent_replacement_number=length(adjacent_replacement_linear_indices);
    
    % These are the weights of the adjacent replacement values
    replacement_weights=-weight_scaling_factor*weighting_vector(current_replacement_points);
    
    % This is the vector of indices to modify in the replacement vectors
    A_Index_Vector=(A_Last_Index+1):(A_Last_Index+adjacent_replacement_number);
    % This updates the last index into A that was modified
    A_Last_Index=A_Last_Index+adjacent_replacement_number;
    
    % This adds current adjacent component indices to the full replacement
    % index vectors
    A_ii(A_Index_Vector)=replace_component_index*ones(adjacent_replacement_number,1);
    A_jj(A_Index_Vector)=linear_replacement_cumulative_sum(adjacent_replacement_linear_indices)';
    % This adds the current adjacent component values to the full
    % replacement value vector
    A_value(A_Index_Vector)=replacement_weights;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculates Retainment Matrix Components                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % These are the linear indices of the adjacent points that do not need
    % to be replaced (some of which may be zero - these are adjacent points that
    % are to be retained)
    adjacent_retainment_linear_indices=not(current_replacement_points).*adjacent_linear_indices;
    % This removes the zero-valued components of the adjacent retainment
    % linear index vector
    adjacent_retainment_linear_indices(adjacent_retainment_linear_indices==0)=[];
    % This is the number of adjacent componets that do not need to be
    % replaced
    adjacent_retainment_number=length(adjacent_retainment_linear_indices);
    
    % These are the weights of the adjacent retainment values
    retainment_weights=weight_scaling_factor*weighting_vector(not(current_replacement_points));
    
    % This is the vector of indices to modify in the retainment vectors
    B_Index_Vector=(B_Last_Index+1):(B_Last_Index+adjacent_retainment_number);
    % This updates the last index into B that was modified
    B_Last_Index=B_Last_Index+adjacent_retainment_number;
    % This adds current adjacent component indices to the full retainment
    % index vectors
    B_ii(B_Index_Vector)=replace_component_index*ones(adjacent_retainment_number,1);
    B_jj(B_Index_Vector)=linear_retainment_cumulative_sum(adjacent_retainment_linear_indices)';
    % This adds the current adjacent component values to the full
    % retainment value vector
    B_value(B_Index_Vector)=retainment_weights;
    
end;

% If the last index of A is less than the length of the vectors in A, this
% cuts off the end of the vectors
if (A_Last_Index<length(A_ii))||(A_Last_Index<length(A_value));
    % This removes the remaining portions of the replacement array vectors
    A_ii((A_Last_Index+1):end)=[];
    A_jj((A_Last_Index+1):end)=[];
    A_value((A_Last_Index+1):end)=[];
end;

% If the last index of B is less than the length of the vectors in B, this
% cuts off the end of the vectors
if (B_Last_Index<length(B_ii))||(B_Last_Index<length(B_value));
    % This removes the remaining portions of the retainment array vectors
    B_ii((B_Last_Index+1):end)=[];
    B_jj((B_Last_Index+1):end)=[];
    B_value((B_Last_Index+1):end)=[];
end;

% This creates the sparse replacement weighting matrix
A=sparse(A_ii,A_jj,A_value,A_ii_size,A_jj_size);
% This creates the sparse retainment weighting matrix
B=sparse(B_ii,B_jj,B_value,B_ii_size,B_jj_size);



function [U,V,W]=delaunay_replacement(X,Y,Z,U,V,W,outlier_vector_array,weighting_method);
% This function uses a Delaunay triangulation of the valid vector data to
% interpolate the vector data onto the outlier data points.  The input
% floating point velocity  vector arrays 'U', 'V', and 'W' may be either
% 2D or 3D in size, but must all have the same dimensions.  The logical
% array 'outlier_vector_array' is of the same size as 'U', 'V', and 'W'
% and contains true values (id est 1s) at the location of the outlier
% vectors.

% This is the number of windows (or vectors) in each dimension
window_number=zeros(1,3);
window_number(1)=size(X,1);
window_number(2)=size(X,2);
window_number(3)=size(X,3);

% This extracts the valid coordinate points
X_Valid=X(not(outlier_vector_array(:)));
Y_Valid=Y(not(outlier_vector_array(:)));
Z_Valid=Z(not(outlier_vector_array(:)));

% This extracts the outlier coordinate poitns
X_Interp=X(outlier_vector_array(:));
Y_Interp=Y(outlier_vector_array(:));
Z_Interp=Z(outlier_vector_array(:));

% This checks whether the 'scatteredInterpolant' function exists in the
% current version of matlab and if not, calls the 'TriScatteredInterp'
% function which was used in older versions of matlab
if exist('scatteredInterpolant','file')==2;
    
    % This switches between performing 2D and 3D interpolation depending
    % upon the size of the vector field
    if window_number(3)>1;
        
        % This extracts a vector of the valid U vectors
        U_Valid=U(not(outlier_vector_array(:)));
        % This creates an interpolant function based upon the valid U 
        % vector component values
        Interp_Function=scatteredInterpolant(X_Valid,Y_Valid,Z_Valid,U_Valid,weighting_method,'linear');
        % This calculates the interpolated values of the outliers of U
        U_Interp=Interp_Function(X_Interp,Y_Interp,Z_Interp);
        % This replaces the outlier U components
        U(outlier_vector_array(:))=U_Interp;
        
        % This extracts a vector of the valid V vectors
        V_Valid=V(not(outlier_vector_array(:)));
        % This modifies the interpolant function to use the valid V vector
        % components
        Interp_Function.Values=V_Valid;
        % This calculates the interpolated values of the outliers of V
        V_Interp=Interp_Function(X_Interp,Y_Interp,Z_Interp);
        % This replaces the outlier V components
        V(outlier_vector_array(:))=V_Interp;
        
        % This extracts a vector of the valid W vectors
        W_Valid=W(not(outlier_vector_array(:)));
        % This modifies the interpolant function to use the valid W vector
        % components
        Interp_Function.Values=W_Valid;
        % This calculates the interpolated values of the outliers of W
        W_Interp=Interp_Function(X_Interp,Y_Interp,Z_Interp);
        % This replaces the outlier W components
        W(outlier_vector_array(:))=W_Interp;
        
    elseif window_number(3)==1;
        
        % This extracts a vector of the valid U vectors
        U_Valid=U(not(outlier_vector_array(:)));
        % This creates an interpolant function based upon the valid U 
        % vector component values
        Interp_Function=scatteredInterpolant(X_Valid,Y_Valid,U_Valid,weighting_method,'linear');
        % This calculates the interpolated values of the outliers of U
        U_Interp=Interp_Function(X_Interp,Y_Interp);
        % This replaces the outlier U components
        U(outlier_vector_array(:))=U_Interp;
        
        % This extracts a vector of the valid V vectors
        V_Valid=V(not(outlier_vector_array(:)));
        % This modifies the interpolant function to use the valid V vector
        % components
        Interp_Function.Values=V_Valid;
        % This calculates the interpolated values of the outliers of V
        V_Interp=Interp_Function(X_Interp,Y_Interp);
        % This replaces the outlier V components
        V(outlier_vector_array(:))=V_Interp;
        
        % This extracts a vector of the valid W vectors
        W_Valid=W(not(outlier_vector_array(:)));
        % This modifies the interpolant function to use the valid W vector
        % components
        Interp_Function.Values=W_Valid;
        % This calculates the interpolated values of the outliers of W
        W_Interp=Interp_Function(X_Interp,Y_Interp);
        % This replaces the outlier W components
        W(outlier_vector_array(:))=W_Interp;
        
    end;
    
else;
    
    % This switches between performing 2D and 3D interpolation depending
    % upon the size of the vector field
    if window_number(3)>1;
        
        % This uses the most recent version of the Delaunay triangulation
        % algorithm
        if exist('delaunayTriangulation','file')==2;
            % This computes the Delaunay triangulation of the valid points
            Delaunay_Tri=delaunayTriangulation(X_Valid,Y_Valid,Z_Valid);
        else;
            % This computes the Delaunay triangulation of the valid points
            Delaunay_Tri=DelaunayTri(X_Valid,Y_Valid,Z_Valid);
        end;
        
        % This extracts a vector of the valid U vectors
        U_Valid=U(not(outlier_vector_array(:)));
        % This creates an interpolant function based upon the valid vector
        % component values
        U_Interp_Function=TriScatteredInterp(Delaunay_Tri,U_Valid,weighting_method);
        % This calculates the interpolated values of the outliers of U
        U_Interp=U_Interp_Function(X_Interp,Y_Interp,Z_Interp);
        % This replaces the outlier U components
        U(outlier_vector_array(:))=U_Interp;
        % This identifies any points in the interpolated function that have
        % NaN values due to extrapolation
        nan_indices=isnan(U(:));
        % If there are any NaN values, the interpolation is repeated for
        % these values, using a nearest neighbor estimate (which will not
        % produce NaN values for extrapolated points)
        if any(nan_indices);
            % This creates an interpolant function based upon the valid vector
            % component values
            U_Interp_Function=TriScatteredInterp(Delaunay_Tri,U_Valid,'nearest');
            % This calculates the interpolated values of the outliers of U
            U_Interp=U_Interp_Function(X(nan_indices),Y(nan_indices),Z(nan_indices));
            % This replaces the outlier U components
            U(nan_indices)=U_Interp;
        end;
        
        % This extracts a vector of the valid V vectors
        V_Valid=V(not(outlier_vector_array(:)));
        % This creates an interpolant function based upon the valid vector
        % component values
        V_Interp_Function=TriScatteredInterp(Delaunay_Tri,V_Valid,weighting_method);
        % This calculates the interpolated values of the outliers of V
        V_Interp=V_Interp_Function(X_Interp,Y_Interp,Z_Interp);
        % This replaces the outlier V components
        V(outlier_vector_array(:))=V_Interp;
        % This identifies any points in the interpolated function that have
        % NaN values due to extrapolation
        nan_indices=isnan(V(:));
        % If there are any NaN values, the interpolation is repeated for
        % these values, using a nearest neighbor estimate (which will not
        % produce NaN values for extrapolated points)
        if any(nan_indices);
            % This creates an interpolant function based upon the valid vector
            % component values
            V_Interp_Function=TriScatteredInterp(Delaunay_Tri,V_Valid,'nearest');
            % This calculates the interpolated values of the outliers of V
            V_Interp=V_Interp_Function(X(nan_indices),Y(nan_indices),Z(nan_indices));
            % This replaces the outlier U components
            V(nan_indices)=V_Interp;
        end;
        
        % This extracts a vector of the valid W vectors
        W_Valid=W(not(outlier_vector_array(:)));
        % This creates an interpolant function based upon the valid vector
        % component values
        W_Interp_Function=TriScatteredInterp(Delaunay_Tri,W_Valid,weighting_method);
        % This calculates the interpolated values of the outliers of W
        W_Interp=W_Interp_Function(X_Interp,Y_Interp,Z_Interp);
        % This replaces the outlier W components
        W(outlier_vector_array(:))=W_Interp;
        % This identifies any points in the interpolated function that have
        % NaN values due to extrapolation
        nan_indices=isnan(W(:));
        % If there are any NaN values, the interpolation is repeated for
        % these values, using a nearest neighbor estimate (which will not
        % produce NaN values for extrapolated points)
        if any(nan_indices);
            % This creates an interpolant function based upon the valid vector
            % component values
            W_Interp_Function=TriScatteredInterp(Delaunay_Tri,W_Valid,'nearest');
            % This calculates the interpolated values of the outliers of W
            W_Interp=W_Interp_Function(X(nan_indices),Y(nan_indices),Z(nan_indices));
            % This replaces the outlier W components
            W(nan_indices)=W_Interp;
        end;
        
    elseif window_number(3)==1;
        
        % This uses the most recent version of the Delaunay triangulation
        % algorithm
        if exist('delaunayTriangulation','file')==2;
            % This computes the Delaunay triangulation of the valid points
            Delaunay_Tri=delaunayTriangulation(X_Valid,Y_Valid);
        else;
            % This computes the Delaunay triangulation of the valid points
            Delaunay_Tri=DelaunayTri(X_Valid,Y_Valid);
        end;
        
        % This extracts a vector of the valid U vectors
        U_Valid=U(not(outlier_vector_array(:)));
        % This creates an interpolant function based upon the valid vector
        % component values
        U_Interp_Function=TriScatteredInterp(Delaunay_Tri,U_Valid,weighting_method);
        % This calculates the interpolated values of the outliers of U
        U_Interp=U_Interp_Function(X_Interp,Y_Interp);
        % This replaces the outlier U components
        U(outlier_vector_array(:))=U_Interp;
        % This identifies any points in the interpolated function that have
        % NaN values due to extrapolation
        nan_indices=isnan(U(:));
        % If there are any NaN values, the interpolation is repeated for
        % these values, using a nearest neighbor estimate (which will not
        % produce NaN values for extrapolated points)
        if any(nan_indices);
            % This creates an interpolant function based upon the valid vector
            % component values
            U_Interp_Function=TriScatteredInterp(Delaunay_Tri,U_Valid,'nearest');
            % This calculates the interpolated values of the outliers of U
            U_Interp=U_Interp_Function(X(nan_indices),Y(nan_indices));
            % This replaces the outlier U components
            U(nan_indices)=U_Interp;
        end;
        
        % This extracts a vector of the valid V vectors
        V_Valid=V(not(outlier_vector_array(:)));
        % This creates an interpolant function based upon the valid vector
        % component values
        V_Interp_Function=TriScatteredInterp(Delaunay_Tri,V_Valid,weighting_method);
        % This calculates the interpolated values of the outliers of V
        V_Interp=V_Interp_Function(X_Interp,Y_Interp);
        % This replaces the outlier V components
        V(outlier_vector_array(:))=V_Interp;
        % This identifies any points in the interpolated function that have
        % NaN values due to extrapolation
        nan_indices=isnan(V(:));
        % If there are any NaN values, the interpolation is repeated for
        % these values, using a nearest neighbor estimate (which will not
        % produce NaN values for extrapolated points)
        if any(nan_indices);
            % This creates an interpolant function based upon the valid vector
            % component values
            V_Interp_Function=TriScatteredInterp(Delaunay_Tri,V_Valid,'nearest');
            % This calculates the interpolated values of the outliers of V
            V_Interp=V_Interp_Function(X(nan_indices),Y(nan_indices));
            % This replaces the outlier U components
            V(nan_indices)=V_Interp;
        end;
        
        % This extracts a vector of the valid W vectors
        W_Valid=W(not(outlier_vector_array(:)));
        % This creates an interpolant function based upon the valid vector
        % component values
        W_Interp_Function=TriScatteredInterp(Delaunay_Tri,W_Valid,weighting_method);
        % This calculates the interpolated values of the outliers of W
        W_Interp=W_Interp_Function(X_Interp,Y_Interp);
        % This replaces the outlier W components
        W(outlier_vector_array(:))=W_Interp;
        % This identifies any points in the interpolated function that have
        % NaN values due to extrapolation
        nan_indices=isnan(W(:));
        % If there are any NaN values, the interpolation is repeated for
        % these values, using a nearest neighbor estimate (which will not
        % produce NaN values for extrapolated points)
        if any(nan_indices);
            % This creates an interpolant function based upon the valid vector
            % component values
            W_Interp_Function=TriScatteredInterp(Delaunay_Tri,W_Valid,'nearest');
            % This calculates the interpolated values of the outliers of W
            W_Interp=W_Interp_Function(X(nan_indices),Y(nan_indices));
            % This replaces the outlier W components
            W(nan_indices)=W_Interp;
        end;
        
    end;
    
end;



function [U_Smooth,V_Smooth,W_Smooth]=smooth_velocity_field(U,V,W,Kernal_Size,sigma);
% This function smooths the velocity fields given by 'U', 'V', and 'W' using a Gaussian
% convolution kernal of dimensions given by 'Kernal_Size' with a diameter specified by
% the Gaussian standard deviation 'sigma'.

% These are the vectors to create the Gaussian kernal over
ii_vector=-(Kernal_Size(1)-1)/2:1:(Kernal_Size(1)-1)/2;
jj_vector=-(Kernal_Size(2)-1)/2:1:(Kernal_Size(2)-1)/2;
kk_vector=-(Kernal_Size(3)-1)/2:1:(Kernal_Size(3)-1)/2;

% This creates a coordinate array for the Gaussian kernal
[ii_array,jj_array,kk_array]=ndgrid(ii_vector,jj_vector,kk_vector);

% This creates the Gaussian kernal
h=exp(-(ii_array.^2)/(2*sigma(1)^2)-(jj_array.^2)/(2*sigma(2)^2)-(kk_array.^2)/(2*sigma(3)^2));

% This rescales the kernal so that the vectors magnitude isn't changed
h=h/sum(h(:));

% This smooths the vector fields
U_Smooth=imfilter(U,h,'circular','same');
V_Smooth=imfilter(V,h,'circular','same');
W_Smooth=imfilter(W,h,'circular','same');



function [ii_extrap_position,jj_extrap_position,kk_extrap_position,extrap_scalar_field]=extrapolate_velocity_field(ii_position,jj_position,kk_position,scalar_field,ii_extrap_min,ii_extrap_max,jj_extrap_min,jj_extrap_max,kk_extrap_min,kk_extrap_max);
% This function extrapolates the velocity field given by
% 'measured_velocity' which is known at the coordinates given by
% 'ii_position', 'jj_position', and 'kk_position' to coordinates outside of
% the known domain such that the velocity field may later be interpolated
% correctly.  The scalars 'ii_min', 'ii_max', 'jj_min', 'jj_max', 'kk_min',
% and 'kk_max' give the domain over which the output velocity field should
% be calculated.

% This calculates the current domain of the scalar field
ii_scalar_min=min(ii_position(:));
ii_scalar_max=max(ii_position(:));
jj_scalar_min=min(jj_position(:));
jj_scalar_max=max(jj_position(:));
kk_scalar_min=min(kk_position(:));
kk_scalar_max=max(kk_position(:));

% This calculates the spacing between the ii values
h_ii=unique(diff(ii_position,1,1));
% This calculates the spacing between the jj values
h_jj=unique(diff(jj_position,1,2));
% If the array is 2D in size, the spacing is set to one, otherwise if the
% array is 3D, this calculates the kk spacing
if (kk_scalar_min==kk_scalar_max);
    % This sets the spacing in the kk dimension to one
    h_kk=1;
else;
    % This calculates the spacing between the kk values
    h_kk=unique(diff(kk_position,1,3));
end;
    
% This calculates the minimum coordinate value in the ii dimension to
% extrapolate the scalar field to
if ii_extrap_min<ii_scalar_min;
    % This calculates the minimum extrapolation coordinate in the ii
    % dimension
    ii_extrap_min=ii_scalar_min-h_ii*ceil((ii_scalar_min-ii_extrap_min)/h_ii);
else;
    % This sets the minimum extrapolation coordinate in the ii
    % dimension to the minimum current ii coordinate
    ii_extrap_min=ii_scalar_min;
end;
% This calculates the maximum coordinate value in the ii dimension to
% extrapolate the scalar field to
if ii_extrap_max>ii_scalar_max;
    % This calculates the maximum extrapolation coordinate in the ii
    % dimension
    ii_extrap_max=ii_scalar_max+h_ii*ceil((ii_extrap_max-ii_scalar_max)/h_ii);
else;
    % This sets the maximum extrapolation coordinate in the ii
    % dimension to the maximum current ii coordinate
    ii_extrap_max=ii_scalar_max;
end;

% This calculates the minimum coordinate value in the jj dimension to
% extrapolate the scalar field to
if jj_extrap_min<jj_scalar_min;
    % This calculates the minimum extrapolation coordinate in the jj
    % dimension
    jj_extrap_min=jj_scalar_min-h_jj*ceil((jj_scalar_min-jj_extrap_min)/h_jj);
else;
    % This sets the minimum extrapolation coordinate in the jj
    % dimension to the minimum current jj coordinate
    jj_extrap_min=jj_scalar_min;
end;
% This calculates the maximum coordinate value in the jj dimension to
% extrapolate the scalar field to
if jj_extrap_max>jj_scalar_max;
    % This calculates the maximum extrapolation coordinate in the jj
    % dimension
    jj_extrap_max=jj_scalar_max+h_jj*ceil((jj_extrap_max-jj_scalar_max)/h_jj);
else;
    % This sets the maximum extrapolation coordinate in the jj
    % dimension to the maximum current jj coordinate
    jj_extrap_max=jj_scalar_max;
end;

% If ths array is 2D, then this sets the extrapolation values equal to the
% current minimum and maximum, otherwise if the array is 3D, this
% calculates the new extrapolation minimum and maximum values
if (kk_scalar_min==kk_scalar_max);
    % This sets the minimum extrapolation coordinate in the kk
    % dimension to the minimum current kk coordinate
    kk_extrap_min=kk_scalar_min;
    % This sets the maximum extrapolation coordinate in the kk
    % dimension to the maximum current kk coordinate
    kk_extrap_max=kk_scalar_max;
else;
    % This calculates the minimum coordinate value in the kk dimension to
    % extrapolate the scalar field to
    if kk_extrap_min<kk_scalar_min;
        % This calculates the minimum extrapolation coordinate in the kk
        % dimension
        kk_extrap_min=kk_scalar_min-h_kk*ceil((kk_scalar_min-kk_extrap_min)/h_kk);
    else;
        % This sets the minimum extrapolation coordinate in the kk
        % dimension to the minimum current kk coordinate
        kk_extrap_min=kk_scalar_min;
    end;
    % This calculates the maximum coordinate value in the kk dimension to
    % extrapolate the scalar field to
    if kk_extrap_max>kk_scalar_max;
        % This calculates the maximum extrapolation coordinate in the kk
        % dimension
        kk_extrap_max=kk_scalar_max+h_kk*ceil((kk_extrap_max-kk_scalar_max)/h_kk);
    else;
        % This sets the maximum extrapolation coordinate in the kk
        % dimension to the maximum current kk coordinate
        kk_extrap_max=kk_scalar_max;
    end;
end;

% This creates the coordinate vectors to extrapolate the scalar field to
ii_extrap_vector=ii_extrap_min:h_ii:ii_extrap_max;
jj_extrap_vector=jj_extrap_min:h_jj:jj_extrap_max;
kk_extrap_vector=kk_extrap_min:h_kk:kk_extrap_max;

% These are the coordinate arrays to extrapolate the vector field to
[ii_extrap_position,jj_extrap_position,kk_extrap_position]=meshgrid(ii_extrap_vector,jj_extrap_vector,kk_extrap_vector);

% This determines the indices of the points that need to be extrapolated in
% the ii dimension
extrap_indices=((ii_extrap_position(:)<ii_scalar_min)|(ii_scalar_max<ii_extrap_position(:)));
% This determines the indices of the points that need to be extrapolated in
% the ii and jj dimensions
extrap_indices=((extrap_indices)|(jj_extrap_position(:)<jj_scalar_min)|(jj_scalar_max<jj_extrap_position(:)));
% If the kk dimension scalar minimum and maximum are different, then this 
% determines the indices of the points that need to be extrapolated in the 
% ii, jj, or kk dimensions
if kk_scalar_min~=kk_scalar_max;
    % This determines the indices of the points that need to be 
    % extrapolated in the ii, jj, or kk dimensions
    extrap_indices=((extrap_indices)|(kk_extrap_position(:)<kk_scalar_min)|(kk_scalar_max<kk_extrap_position(:)));
end;

% This initializes an array to store the known scalar field and the
% extrapolated field
extrap_scalar_field=zeros(size(ii_extrap_position));

% This sets the values of this field within the known field equal to the
% known field
extrap_scalar_field(not(extrap_indices))=scalar_field(:);

% This sets the values of the field outside of the known field to the
% extrapolated values
extrap_scalar_field(extrap_indices)=taylor_series_extrapolation(ii_position,jj_position,kk_position,scalar_field,ii_extrap_position(extrap_indices),jj_extrap_position(extrap_indices),kk_extrap_position(extrap_indices));



function WI=taylor_series_extrapolation(X,Y,Z,W,XI,YI,ZI);
% This function is used to perform interpolation and extrapolation of a
% scalar field.  The scalar field is defined by the function
%
%   W(X,Y,Z)
%
% where 'X', 'Y', 'Z', and 'W' are N x M x P arrays (and P may be equal to
% one for 2D arrays).  The function is interpolated at the points 'XI',
% 'YI', and 'ZI' which may have any dimensions.  The interpolated or
% extrapolated output array 'WI' has the same dimensions as the input
% interpolation coordinates 'XI', 'YI', and 'ZI'.
%
% The interpolation/extrapolation operation is performed by finding the
% nearest point in the coordinates 'X', 'Y', and 'Z' and performing a first
% order Taylor series expansion about this point.  This produces a linear
% interpolation or extrapolation from the known data.
%
% To decrease sensitivity to noise in the input data array 'W', the
% interpolated/extrapolated value returned in 'WI' is taken as the Gaussian
% weighted sum of the Taylor series of the surrounding points.
%
% This function is named as an extrapolation function since it will have
% relatively low performance as an interpolation function (due to being
% linear and only using a limited number of adjacent points), but will 
% produce a smooth extrapolation beyond the edge of the known scalar field.
%
% The diameter of the Gaussian weighting function is designed such that the
% nearest point contributes approximately half of the signal to the output
% 'WI' while the adjacent points contribute the other half of the signal.
% This weighting was chosen with the assumption that the input data will 
% be a scalar field where the adjacent components of 'W' will have some
% correlation (but the correlation falls off beyond a distance of 
% approximately 1 unit).
%
% The size of the Gaussian kernal was taken such that along it's edges the
% weighting is approximately 1% or less of the maximum weighting.
%
% Authors: Rod La Foy
% First Created On: 26 March 2014
% Last Modified On: 26 March 2014

% This is the number of points to perform the interpolation at
interpolation_point_number=numel(XI);

% This is the standard deviation of the Gaussian function (which is chosen
% so that the nearest vector contributes approximately 50% of the
% information to the interpolated value while nearby vectors contribute the
% remaining 50% of the information)
gaussian_sigma=0.8;

% This is the kernal size for smoothing the scalar field about the nearest
% point (the size of this is chosen such that vectors on the edge of the
% kernal contribute approximately 1% or less of the total signal)
kernal_size=[5,5,5];

% This is the "radius" of the smoothing kernal
kernal_radius=(kernal_size-1)/2;

% These are the coordinate vectors for calculating the Gaussian function
ii_vector=-kernal_radius(1):kernal_radius(1);
jj_vector=-kernal_radius(2):kernal_radius(2);
kk_vector=-kernal_radius(3):kernal_radius(3);
% These are the coordinate arrays for calculating the Gaussian function
[ii_array,jj_array,kk_array]=ndgrid(ii_vector,jj_vector,kk_vector);
% This is the radius to the origin of the coordinates
gaussian_radius_sqrd=ii_array.^2+jj_array.^2+kk_array.^2;
% This calculates the Gaussian function
gaussian_kernal=exp(-gaussian_radius_sqrd/(2*gaussian_sigma^2));

% This is the size of the scalar field
scalar_field_size(1)=size(X,1);
scalar_field_size(2)=size(X,2);
scalar_field_size(3)=size(X,3);

% This calculates the spacing between the X values
h_x=unique(diff(X,1,1));
% This calculates the spacing between the Y values
h_y=unique(diff(Y,1,2));
% This calculates the spacing between the Z values
h_z=unique(diff(Z,1,3));

% If the array is 2D, this sets the spacing in the third dimension 1 so
% that the function is compatible with 2D and 3D data
if scalar_field_size(3)==1;
    % This sets the spacing in the Z direction to 1
    h_z=1;
end;

% If there is more than one spacing an error is returned
if (numel(h_x)>1)||(numel(h_y)>1)||(numel(h_z)>1);
     % This displays an error
     error(' The gradient cannot be calculated with this coordinate grid: the spacing between the coordinates is non-uniform.');
end;
% If the spacing is zero in any dimension an error is returned
if (abs(h_x)<eps(1))||(abs(h_y)<eps(1))||(abs(h_z)<eps(1));
    % This displays an error
    error(' The gradient cannot be calculated with this coordinate grid: the spacing between the coordinates is non-uniform.');
end;

% If the data is 3D, this calculates the gradient in all three
% dimensions, otherwise the gradient is only calculated in the first
% two dimensions and set equal to zero in the third
if scalar_field_size(3)>1;
    % This calculates the gradient with respect to W using the matlab code
    [dWdX,dWdY,dWdZ]=gradient(W,h_y,h_x,h_z);
else;
    % This calculates the gradient with respect to W using the matlab code
    [dWdX,dWdY]=gradient(W,h_y,h_x);
    % This sets the gradient in the third dimension to zero
    dWdZ=zeros(scalar_field_size);
end;

% This initializes the output array
WI=zeros(size(XI));

% This iterates through the array of interpolation points calculated the
% interpolated values
for interpolation_point_index=1:interpolation_point_number;
    
    % This is the distance to the known points within the array
    rho_sqrd=(X(:)-XI(interpolation_point_index)).^2+(Y(:)-YI(interpolation_point_index)).^2+(Z(:)-ZI(interpolation_point_index)).^2;
    
    % This is the index of the nearest point to the current interpolation
    % point
    [~,nearest_point_index]=min(rho_sqrd);
    
    % This converts the index of the nearest point to subscripted indices
    [ii_nearest,jj_nearest,kk_nearest]=ind2sub(scalar_field_size,nearest_point_index);
    
    % This is the domain of the scalar field to extract
    ii_window_min=ii_nearest-kernal_radius(1);
    ii_window_max=ii_nearest+kernal_radius(1);
    jj_window_min=jj_nearest-kernal_radius(2);
    jj_window_max=jj_nearest+kernal_radius(2);
    kk_window_min=kk_nearest-kernal_radius(3);
    kk_window_max=kk_nearest+kernal_radius(3);
    
    % This checks if the ii indices are outside the domain of the scalar field and sets them to the
    % domain of the scalar field if true
    if ii_window_min<1;
        ii_window_min=1;
    end;
    if ii_window_max>scalar_field_size(1);
        ii_window_max=scalar_field_size(1);
    end;
    % This checks if the jj indices are outside the domain of the scalar field and sets them to the
    % domain of the scalar field if true
    if jj_window_min<1;
        jj_window_min=1;
    end;
    if jj_window_max>scalar_field_size(2);
        jj_window_max=scalar_field_size(2);
    end;
    % This checks if the jj indices are outside the domain of the scalar field and sets them to the
    % domain of the scalar field if true
    if kk_window_min<1;
        kk_window_min=1;
    end;
    if kk_window_max>scalar_field_size(3);
        kk_window_max=scalar_field_size(3);
    end;
    
    % This is the domain of the smoothing kernal to extract
    ii_kernal_min=ii_window_min-ii_nearest+kernal_radius(1)+1;
    ii_kernal_max=ii_window_max-ii_nearest+kernal_radius(1)+1;
    jj_kernal_min=jj_window_min-jj_nearest+kernal_radius(2)+1;
    jj_kernal_max=jj_window_max-jj_nearest+kernal_radius(2)+1;
    kk_kernal_min=kk_window_min-kk_nearest+kernal_radius(3)+1;
    kk_kernal_max=kk_window_max-kk_nearest+kernal_radius(3)+1;
    
    % This extracts the current X coordinates
    X_Window=X(ii_window_min:ii_window_max,jj_window_min:jj_window_max,kk_window_min:kk_window_max);
    % This extracts the current Y coordinates
    Y_Window=Y(ii_window_min:ii_window_max,jj_window_min:jj_window_max,kk_window_min:kk_window_max);
    % This extracts the current Z coordinates
    Z_Window=Z(ii_window_min:ii_window_max,jj_window_min:jj_window_max,kk_window_min:kk_window_max);
    
    % This extracts the current window from the scalar field
    W_Window=W(ii_window_min:ii_window_max,jj_window_min:jj_window_max,kk_window_min:kk_window_max);
    
    % This extracts the current window from the X derivative
    dWdX_Window=dWdX(ii_window_min:ii_window_max,jj_window_min:jj_window_max,kk_window_min:kk_window_max);
    % This extracts the current window from the Y derivative
    dWdY_Window=dWdY(ii_window_min:ii_window_max,jj_window_min:jj_window_max,kk_window_min:kk_window_max);
    % This extracts the current window from the Z derivative
    dWdZ_Window=dWdZ(ii_window_min:ii_window_max,jj_window_min:jj_window_max,kk_window_min:kk_window_max);
    
    % This extracts the same portion of the Gaussian smoothin kernal
    gaussian_kernal_temp=gaussian_kernal(ii_kernal_min:ii_kernal_max,jj_kernal_min:jj_kernal_max,kk_kernal_min:kk_kernal_max);
    
    % This normalizes the Gaussian kernal so that the sum is equal to one
    gaussian_kernal_temp=gaussian_kernal_temp/sum(gaussian_kernal_temp(:));
    
    % This calculates the interpolated point as a function of each of the
    % elements within the current window
    WI_Window=W_Window+0*dWdX_Window.*(XI(interpolation_point_index)-X_Window)+...
        0*dWdY_Window.*(YI(interpolation_point_index)-Y_Window)+...
        0*dWdZ_Window.*(ZI(interpolation_point_index)-Z_Window);
    
    % This calculates the interpolated point as a weighted sum of the
    % adjacent points to the nearest point (with the weighting scaled by
    % the Gaussian function)
    WI(interpolation_point_index)=sum(WI_Window(:).*gaussian_kernal_temp(:));
    
end;



function p=determine_gaussian_size_2(window_size,resolution_size);
% This function is used to determine the diameter of a Gaussian function
% used to mask PIV windows to prevent aliasing during the Fourier
% transform.
%
% The input argument 'window_size' refers to the full size of
% the PIV window prior to performing the masking operation and may be a
% scaler or a vector.  The input argument 'resolution_size' refers to the
% effective resolution that the user wants to attain from the masked PIV
% window, ie the window may be 64 x 64 in size (window_size=[64,64]) but
% the user wants the effective resolution after masking to be equal to 32 x
% 32 (resolution_size=[32,32]).
%
% The output argument 'p' is used to contruct the Gaussian masking function
% of the following form
%
%  G(x) = exp(-((p*x)^2)/2)                                             (1)
%
% where the independent variable lies in the domain
%
%  -(window_size-1)/2 <= x <= (window_size-1)/2                         (2)
%
% The variable 'p' is calculated by solving the integral equation
%
%  resolution_size = int(G(x),-(window_size-1)/2,(window_size-1)/2)     (3)
%
% or the equivalent expression
%
%  resolution_size = sqrt(2*pi)*erf(window_size*p/(2*sqrt(2))/p         (4)
%
% given by simplifying the integral.  These equations are taken from
% "Assessment of advanced windowing techniques for digital particle image
% velocimetry (DPIV)" Eckstein 2009.
%
% This function is written to replace the 'findwidth' function in prana.
% The equation is solved by re-writing equation (4) in the form
% 
%  0 = f(x) = erf(b*x)/x-a                                              (5)
% 
% and solving for the independent variable using Halley's method.
%
% Written By: Rod La Foy
% Written On: 6 June 2013

% These are the ratios of the resolution size to window size
size_ratio=resolution_size./window_size;

% This initializes the output vector p
p=zeros(size(window_size));

% This pre-computes the square root of pi for later use
pi_sqrt=sqrt(pi);

% This is the tolerance with which to find the root
tolerance=1e-4;

% This iterates through the number of dimensions to solve for (probably
% either 2 or 3)
for dimension_number=1:length(window_size);

	% This checks if the size ratio is equal to one or greater
	if size_ratio(dimension_number)>=1;
		
		% This sets the Gaussian width variable to 0
		p(dimension_number)=0;
		
	else;
    
		% This is the first constant term in the equation
		a=resolution_size(dimension_number)/sqrt(2*pi);
		% This is the second constant term in the equation
		b=window_size(dimension_number)/(2*sqrt(2));

		% This is the initial guess of value of the Gaussian width term
		p_temp=1;
	
		% This iterates to the root
		while true;
		
			% This is the error function of the b*p term
			bp_erf=erf(b*p_temp);
			% This is the exponential term of the (b*p)^2 term
			bp_exp=exp((b*p_temp)^2);
		
			% This is the numerator of the difference term
			dp_numerator=a*p_temp-bp_erf;
		
			% This initializes the denominator term
			dp_denominator=0;
			% This adds the first term to the denominator
			dp_denominator=dp_denominator+a;
			% This adds the second term to the denominator
			dp_denominator=dp_denominator+2*b/(pi_sqrt*bp_exp);
			% This adds the third term to the denominator
			dp_denominator=dp_denominator+(2*b^3*p_temp^2*dp_numerator)/(2*b*p_temp-bp_exp*pi_sqrt*bp_erf);
		
			% This is the expected difference to the root
			dp=dp_numerator/dp_denominator;
		
			% This is the new root approximation
			p_temp=p_temp-dp;
		
			% This breaks if the error is less then the tolerance
			if abs(erf(b*p_temp)/p_temp-a)<tolerance;
				% This breaks the loop
				break;
			end;
		
		end;
	
		% This is the Gaussian function diameter term
		p(dimension_number)=p_temp;
		
	end;
    
end;



function gaussian_window_filter=create_gaussian_filter(window_size,gaussian_width);
% This function creates a Gaussian filter to apply in the spatial domain to
% the correlation windows to remove the effects of aliasing

% This is the standard deviation of the filter
sigma=1./gaussian_width;

% These are the vectors of coordinates used to create the Gaussian function
ii_index_vector=linspace(-(window_size(1)-1)/2,(window_size(1)-1)/2,window_size(1));
jj_index_vector=linspace(-(window_size(2)-1)/2,(window_size(2)-1)/2,window_size(2));
kk_index_vector=linspace(-(window_size(3)-1)/2,(window_size(3)-1)/2,window_size(3));

% These are the full array of coordinates
[ii_index_array,jj_index_array,kk_index_array]=ndgrid(ii_index_vector,jj_index_vector,kk_index_vector);

% This creates the Gaussian function
gaussian_window_filter=exp(-(ii_index_array.^2)/(2*sigma(1)^2)-(jj_index_array.^2)/(2*sigma(2)^2)-(kk_index_array.^2)/(2*sigma(3)^2));



function [frame_1_window_min,frame_1_window_max,frame_2_window_min,frame_2_window_max,displacement_offset]=dwo_window_domains(window_min,window_max,window_center,image_correlation_step,ii_coordinates,jj_coordinates,kk_coordinates,ii_displacement,jj_displacement,kk_displacement);
% This function adds the discrete window offset to the window positions
% based upon the known velocity field.  The non-offset window positions are
% given by the variables
%
%  'window_min'
%  'window_max'
%  'window_center'
%
% where these variables are N x 3 arrays where N is the total number of
% windows and column vectors correspond to the 1st, 2nd, and 3rd dimension
% coordinates.  The known velocity field is given at the coordinates
%
%  'ii_coordinates'
%  'jj_coordinates'
%  'kk_coordinates'
%
% which are I x J x K arrays where I is the number of known vector
% coordinates in the 1st dimension, J is the number of known vector
% coordinates in the 2nd dimension, and K is the number of known vector
% coordinates in teh 3rd dimension.  The known velocity field values are
% given by the arrays
%
%  'ii_displacement'
%  'jj_displacement'
%  'kk_displacement'
%
% which are also I x J x K in size.  The output arrays
%
%  'frame_1_window_min'
%  'frame_1_window_max'
%  'frame_2_window_min'
%  'frame_2_window_max'
%
% are N x 3 vectors in size and given the minimum and maximum range of the
% first and second frames respectively.  The output array
%
%  'displacement_offset'
%
% is also N x 3 in size and corresponds to the offset that must be added
% back to the measured velocity field due to the discrete window offset.

% This scales the velocity estimates by the number of number of frames
% between the correlations
ii_displacement=image_correlation_step*ii_displacement;
jj_displacement=image_correlation_step*jj_displacement;
kk_displacement=image_correlation_step*kk_displacement;

% This checks whether the array of velocity values is 2D or 3D and chooses
% the appropriate interpolotion function
if size(ii_coordinates,3)==1;
    
    % If the vector field is greater then 3 x 3, then a cubic interpolation
    % is performed, otherwise a linear interpolation is performed
    if (size(ii_coordinates,1)>2)&&(size(ii_coordinates,2)>2);
        
        % This is the estimated ii-displacement at the center of the windows
        estimated_window_dii=interp2(ii_coordinates,jj_coordinates,ii_displacement,window_center(:,1),window_center(:,2),'cubic',0);
        % This is the estimated jj-displacement at the center of the windows
        estimated_window_djj=interp2(ii_coordinates,jj_coordinates,jj_displacement,window_center(:,1),window_center(:,2),'cubic',0);
        % This is the estimated kk-displacement at the center of the windows
        estimated_window_dkk=interp2(ii_coordinates,jj_coordinates,kk_displacement,window_center(:,1),window_center(:,2),'cubic',0);
        
    else;
        
        % This is the estimated ii-displacement at the center of the windows
        estimated_window_dii=interp2(ii_coordinates,jj_coordinates,ii_displacement,window_center(:,1),window_center(:,2),'linear',0);
        % This is the estimated jj-displacement at the center of the windows
        estimated_window_djj=interp2(ii_coordinates,jj_coordinates,jj_displacement,window_center(:,1),window_center(:,2),'linear',0);
        % This is the estimated kk-displacement at the center of the windows
        estimated_window_dkk=interp2(ii_coordinates,jj_coordinates,kk_displacement,window_center(:,1),window_center(:,2),'linear',0);
        
    end;
    
else;
    
    % If the vector field is greater then 3 x 3, then a cubic interpolation
    % is performed, otherwise a linear interpolation is performed
    if (size(ii_coordinates,1)>2)&&(size(ii_coordinates,2)>2)&&(size(ii_coordinates,3)>2);
        
        % This is the estimated ii-displacement at the center of the windows
        estimated_window_dii=interp3(ii_coordinates,jj_coordinates,kk_coordinates,ii_displacement,window_center(:,1),window_center(:,2),window_center(:,3),'cubic',0);
        % This is the estimated jj-displacement at the center of the windows
        estimated_window_djj=interp3(ii_coordinates,jj_coordinates,kk_coordinates,jj_displacement,window_center(:,1),window_center(:,2),window_center(:,3),'cubic',0);
        % This is the estimated kk-displacement at the center of the windows
        estimated_window_dkk=interp3(ii_coordinates,jj_coordinates,kk_coordinates,kk_displacement,window_center(:,1),window_center(:,2),window_center(:,3),'cubic',0);
        
    else;
        
        % This is the estimated ii-displacement at the center of the windows
        estimated_window_dii=interp3(ii_coordinates,jj_coordinates,kk_coordinates,ii_displacement,window_center(:,1),window_center(:,2),window_center(:,3),'linear',0);
        % This is the estimated jj-displacement at the center of the windows
        estimated_window_djj=interp3(ii_coordinates,jj_coordinates,kk_coordinates,jj_displacement,window_center(:,1),window_center(:,2),window_center(:,3),'linear',0);
        % This is the estimated kk-displacement at the center of the windows
        estimated_window_dkk=interp3(ii_coordinates,jj_coordinates,kk_coordinates,kk_displacement,window_center(:,1),window_center(:,2),window_center(:,3),'linear',0);
        
    end;

end;
    
% This is half of the measured velocity at the window coordinates
rounded_half_estimated_displacement=round([estimated_window_dii,estimated_window_djj,estimated_window_dkk]/2);

% This calculates the first frame minimum domain values
frame_1_window_min=window_min-floor(rounded_half_estimated_displacement);
% This calculates the first frame maximum domain values
frame_1_window_max=window_max-floor(rounded_half_estimated_displacement);

% This calculates the second frame minimum domain values
frame_2_window_min=window_min+ceil(rounded_half_estimated_displacement);
% This calculates the second frame maximum domain values
frame_2_window_max=window_max+ceil(rounded_half_estimated_displacement);

% This is the offset that must be added back onto the measured velocity
% field of the current pass to account for the dicrete window offset
displacement_offset=ceil(rounded_half_estimated_displacement)+floor(rounded_half_estimated_displacement);



function [window_min,window_max,window_center,window_number]=calculate_window_domains(image_size,window_resolution,window_size,grid_spacing);
% This function calculates the minimum and maximum values of the window
% indices to extract from the image.

% These are the number of windows in each dimension
window_number=ceil((image_size-window_resolution)./grid_spacing)+1;

% This initializes the window indices
window_min=zeros(prod(window_number),3);
window_max=zeros(prod(window_number),3);

% This initializes the window centers
window_center=zeros(prod(window_number),3);

% This is the full resolution covered by the PIV windows (which will be equal to or larger
% than the resolution of the full image)
full_resolution=grid_spacing.*(window_number-1)+window_resolution;

% This is the distance in indices below 1 to start the windows
negative_index_distance=floor((full_resolution-image_size)/2);

% This is a counting variable for indexing the window domains
count=0;
	
% This calculates the window domains
for kk=1:window_number(3);
	for jj=1:window_number(2);
		for ii=1:window_number(1);
		
			% This increments the counting variable
			count=count+1;
			
			% These are the minimum window indices
			window_min(count,:)=([ii,jj,kk]-1).*grid_spacing-ceil((window_size-window_resolution)/2)-negative_index_distance+1;
			% These are the maximum window indices
			window_max(count,:)=window_min(count,:)+window_size-1;

			% These are the window center indices
            window_center(count,:)=(window_min(count,:)+window_max(count,:))/2;

		end;
	end;
end;



function spectral_filter=create_spectral_filter(window_size,particle_diameter);
% This function creates the RPC spectral energy filter with a size of 'window_size' using
% a particle diameter given by 'particle_diameter' where 'bit_depth' is the bit depth of
% the three-dimensional image to be processed.

% These are the wavenumber vectors
k_ii_vector=-pi:2*pi/window_size(1):pi-2*pi/window_size(1);
k_jj_vector=-pi:2*pi/window_size(2):pi-2*pi/window_size(2);
k_kk_vector=-pi:2*pi/window_size(3):pi-2*pi/window_size(3);

% These are the wavenumber coordinates
[k_ii,k_jj,k_kk]=ndgrid(k_ii_vector,k_jj_vector,k_kk_vector);

% This is the squared radius of the wavenumber
rho=k_ii.^2+k_jj.^2+k_kk.^2;

% None of the coefficients actually matter . . . so RPC filter is just a low pass filter . . .
spectral_filter=exp(-(rho*particle_diameter^2)/16);



function display_calculation_progress(current_value,value_vector);
% This function displays a progress bar showing the percent complete of the
% currently running calculation where the calculation is being iterated
% over the vector 'value_vector' and the current iteration's value is equal
% to 'current_value'.  For example in the for loop
%
%  for ii=1:10;
%       commands . . .
%  end;
%
% the arguments would be equal to
%
%  current_value=ii;
%  value_vector=1:10;
%
% This convention holds for non-integer or non-monotonic vectors, although
% to work correctly all values of value_vector must be unique.

% This is the number of characters to display in the pogress bar (denoted
% by either '#' or '-' characters)
progress_bar_character_number=50;

% This is the index of the value_vector that corresponds to the
% current_value
[null_variable,value_index]=min(abs(current_value-value_vector));

% This is the percentage of the calculation that is completed
current_progress_decimal=value_index/length(value_vector);
% This creates the character string showing the numerical and graphical
% progress for the current iteration
current_text_string=generate_progress_string(current_progress_decimal,progress_bar_character_number);

% If this is the first iteration, then a new line is added, the progress
% bar is displayed, and the function is exited
if value_index==1;
    % This displays the portion of the string showing the percentage of 
    % the calculation that is complete that is new to this iteration
    fprintf(current_text_string);
    % This ends the function
    return;
end;

% This is the percentage of the calculation that was completed during the
% last iteration
previous_progress_decimal=(value_index-1)/length(value_vector);
% This creates the character string showing the numerical and graphical
% progress for the previous iteration
previous_text_string=generate_progress_string(previous_progress_decimal,progress_bar_character_number);

% This compares the current progress string with the previous string and if
% they are the same, then the function exits.  If they are different, then
% only the text after the difference is displayed.
if strcmp(current_text_string,previous_text_string);
    
    % If this is the last time that the progress bar will be displayed, this
    % prints a new line
    if value_index==length(value_vector);
        % This prints a new line after the progress bar
        fprintf('\n');
    end;
    
    % This exits the function without changing the progress bar
    return;
    
else;
    % This is the total number of charcters to be displayed
    string_character_number=length(current_text_string);
    
    % This is the index into the string where the strings first differ
    first_difference_index=find(not(current_text_string==previous_text_string),1,'first');
    
    % These are the locations of the double percent signs '%%'
    double_percent_indices=strfind(current_text_string,'%%');
    % This removes the double percent indices that are before the first
    % difference index
    double_percent_indices(double_percent_indices<first_difference_index)=[];

    % This is the number of characters of the previous line to delete
    delete_character_number=string_character_number-first_difference_index+1-length(double_percent_indices);
    % If this is the first iteration, then no characters are deleted
    if value_index==1;
        % This sets the number of characters to be deleted to zero
        delete_character_number=0;
        % This sets the first difference character to one
        first_difference_index=1;
    end;
    
    % This deletes the previously displayed characters back to the first
    % differing character (by displaying the 'backspace' character)
    fprintf(1,repmat('\b',1,delete_character_number));
    
    % This displays the portion of the string showing the percentage of 
    % the calculation that is complete that is new to this iteration
    fprintf(current_text_string(first_difference_index:end));
    
    % If this is the last time that the progress bar will be displayed, this
    % prints a new line
    if value_index==length(value_vector);
        % This prints a new line after the progress bar
        fprintf('\n');
    end;
    
end;



function text_string=generate_progress_string(progress_decimal,progress_bar_character_number);
% This function generates the progress bar text that will be displayed for
% the current percentage of the calculation that is completed.

% This is a string giving a numerical value of the percentage of the
% calculation completed
numerical_progress_string=sprintf('% 4.0f',round(100*progress_decimal));

% This is the prefix to the progress bar showing the numerical value of the
% percent complete
string_prefix=['Calculation',numerical_progress_string,'%% Complete:     0%% ['];

% This is the suffix to the progress bar
string_suffix='] 100%%';

% This is the number of '#' signs to display corresponding to the graphical
% representation of the percent of the calculation that is complete
completed_character_number=round(progress_decimal*progress_bar_character_number);
% This is the number of '-' signs to display corresponding to the graphical
% representation of the percent of the calculation that remains
remaining_character_number=progress_bar_character_number-completed_character_number;

% This creates the string of characters representing the graphical
% percentage of the calculation that is complete
progress_bar_string=[repmat('#',1,completed_character_number),repmat('-',1,remaining_character_number)];

% This is the total text string to be displayed
text_string=[string_prefix,progress_bar_string,string_suffix];



