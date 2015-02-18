function piv_3d_parallel_runner_02;
% This function is designed to create a parameters structure to pass into
% the piv_3d code and to then call the piv_3d code in parallel across
% multiple computers.
%
% This code is based upon the code 'basic_3d_rpc_processing_05.m'.

% This is the temporary directory to save the list of frames being
% processed and the frames completed
processing_directory='/home/rlafoy/Documents/AEThER/Random/Vortex_Shedding_03_Test_Processing/';

% This creates the PIV processing parameters
piv_parameters=create_piv_processing_parameters;

% This tests whether the piv processing parameters were correctly
% created (so that an error is at least less likely to be thrown during
% processing)
validate_processing_parameters(piv_parameters);

% This runs the parallel processing of the data set specified in the
% 'piv_parameters' data structure
run_piv_3d_parallel_processing(piv_parameters,processing_directory);



function run_piv_3d_parallel_processing(piv_parameters,processing_directory);
% This function calls 'piv_3d_07' and runs the processing one frame at a
% time using the parameters specified by 'piv_parameters'.  A list of
% frames that are processing and have been completed are saved into the
% 'processing_directory' so that multiple computers can process the data
% set simulataneously.

% This is the filename of the data file for communicating between multiple
% computers that contains a list of frames that have begun processing.
% The frames listed in this file are not necessarily completed.  All
% parallel computers must be able to read/write this file.
processing_com_filename=[processing_directory,'processing_com_file.dat'];

% This is the filename of the data file for communicating between multiple
% computers that contains a list of frames that have completed processing.
% The frames listed in this file have been fully processed.  All parallel 
% computers must be able to read/write this file.
completed_com_filename=[processing_directory,'completed_com_file.dat'];

% This checks whether the communication files exist and if not creates them
if not(exist(processing_com_filename,'file'));
    % This opens the file for writing
    fid=fopen(processing_com_filename,'w');
    % This closes the file
    fclose(fid);
end;
if not(exist(completed_com_filename,'file'));
    % This opens the file for writing
    fid=fopen(completed_com_filename,'w');
    % This closes the file
    fclose(fid);
end;

% This loads the initial image frame from the data processing structure
initial_image_frame=piv_parameters.general_parameters.initial_image_frame;
% This loads the final image frame from the data processing structure
final_image_frame=piv_parameters.general_parameters.final_image_frame;
% This loads the image frame step from the data processing structure
image_frame_step=piv_parameters.general_parameters.image_frame_step;

% This is a vector of frames to process
frame_pros_vect=(initial_image_frame:image_frame_step:final_image_frame)';

% This iterates through the frame processing vector processing the frames
% in parallel on the several computers.
while true;

    % This opens the completed file for writing
    fid=fopen(completed_com_filename,'a+');
    % This reads the file data
    completed_com_data=fread(fid,inf,'uint16');
    % This checks whether all frames have been processed and if so exits
    % the loop, but first checks whether any frames have been written.
    if not(isempty(completed_com_data));
        % This sorts the completed data vector for comparison to the vector of
        % frames that should be processed
        completed_com_data=sort(completed_com_data,1,'ascend');
        % This removes any non-unique values in the vector. (These shouldn't
        % really exist, but there is a very small probability - less than 1
        % in 10 million - that two computers will simultaneously load the
        % vector and then write that they are processing the same vector back
        % to the file; in this case the second computer will likely overwrite
        % the data file of the first so the output data shouldn't be effected.)
        completed_com_data=unique(completed_com_data);
        % This checks if the length of this vector equals the length of the
        % frame processing vector
        if length(completed_com_data)==length(frame_pros_vect);
            % This is the sum of the differences between the vectors
            frame_vect_residual=sum((double(completed_com_data)-frame_pros_vect).^2);
            % This checks whether the vectors are equal (they should be if
            % they have the same length - but weird stuff happens)
            if frame_vect_residual<1e-14;
                % This closes the communication file
                fclose(fid);
                % This ends the processing loop
                break;
            else;
                % This closes the communication file
                fclose(fid);
                % This displays an error saying something weird is going
                % on in the processing
                error(['The processing file frame list has the same ',...
                    'length as the processing frame vector, but they are ',...
                    'not equal. This likely means that a bug exists ',...
                    'somewhere in the processing code or that there was a ',...
                    'write error.']);
            end;
        else;
            % This closes the communications file
            fclose(fid);
        end;
    else;
        % This closes the communications file
        fclose(fid);
    end;

    % Since the loop has gotten this far, it means that not all frames have
    % been processed.  In this case, this loads the list of frames being
    % processed and based upon this list calculates the next frame number
    % needing to be processed.
    
    % This opens the processing file for writing
    fid=fopen(processing_com_filename,'a+');
    % This reads the file data
    processing_com_data=fread(fid,inf,'uint16');
    % This checks which frames have begun processing, but first checks
    % whether any have begun processing
    if not(isempty(processing_com_data));
        % This is the maximum frame that has begun processing
        max_processing_frame=max(processing_com_data);
        % If the maximum frame does not equal the last frame then the
        % current frame is incremented, otherwise it is set as the last
        % frame
        if max_processing_frame<final_image_frame;
            % This increments the current frame to process
            current_frame=max_processing_frame+image_frame_step;
        else;
            % This sets the current processing frame to the last frame.
            % (This will result in at least one frame being redundantly
            % processed, which should be fixed, but there isn't an easy,
            % obvious way around this problem while still ensuring all
            % frames are processed.)
            current_frame=final_image_frame;
        end;
    else;
        % Otherwise this sets the current frame to process to the first
        % frame to be processed
        current_frame=initial_image_frame;
    end;
    % This writes the current processing frame to the file
    fwrite(fid,current_frame,'uint16');
    % This closes the processing file
    fclose(fid);
    
    % This rewrites the values of 'piv_parameters' for calculating the 
    % current frame
    piv_parameters.general_parameters.initial_image_frame=current_frame;
    piv_parameters.general_parameters.final_image_frame=current_frame;

    % This processes the current frame
    piv_3d_10(piv_parameters);
    
    % This opens the completed com file for writing that the current frame
    % has finished processing
    fid=fopen(completed_com_filename,'a+');
    % This writes the currently completed frame to the file
    fwrite(fid,current_frame,'uint16');
    % This closes the completed com file
    fclose(fid);
    
end;



function piv_parameters=create_piv_processing_parameters;
% This function creates the data structure 'piv_parameters' which contains
% all the data required to perform PIV on the specified data set

% This initializes the parameters data structure
piv_parameters=struct;
% This creates the general processing (non pass specific parameters)
piv_parameters.general_parameters=struct;
% This creates the first pass processing parameters
piv_parameters.pass_parameters=struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the general processing parameters                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a string specifying the directory from which to read the images
image_read_directory='/home/rlafoy/Documents/AEThER/Random/Vortex_Shedding_03/';
% This is a string specifying the filename prefix of the images to be loaded
image_filename_prefix='frame';
% This is a string specifying the filename extension of the images to
% be loaded (2D images may be one of the following formats: 'bmp', 'jpg', 
% 'jp2', 'png', 'ppm', 'tif', 'mat', 'dcm', 'cdf', or 'h5' while 3D images 
% may be one of the following formats: 'mat', 'dcm', 'cdf', or 'h5')
image_filename_extension='tif';
% This is a string specifying the name of the variable to load from the
% data file if the file format can contain multiple variables, ie a 'mat',
% 'cdf', or a 'h5' file - for other formats the variable may be left null
image_variable_name='I';

% This is a string specifying the directory to save the vector field data
% within
vector_field_write_directory='/home/rlafoy/Documents/AEThER/Random/Vortex_Shedding_03_Test_Processing/';

% This is a scalar double specifying the index (into the list of images in
% the image directory) of the first frame to be loaded (ie this image will
% be the first image in the correlation pair)
initial_image_frame=5;
% This is a scalar double specifying the index (into the list of images in
% the image directory) of the last frame to be loaded (ie this image will
% be the first image in the correlation pair)
final_image_frame=5;
% This is a scalar double specifying the number of frames to step between
% each measured vector field
image_frame_step=1;
% This is a scalar double specifying the number of frames to step between
% the first and second images in the correlation pair
image_correlation_step=1;

% This is a scalar double specifying the diameter of the particles in
% voxels (used for creating the RPC filter)
particle_diameter=2.83;

% This is a scalar double specifying the number of PIV passes to perform
% on the images
pass_number=3;

% This is a Boolean value stating whether to perform the window
% correlations in parallel
perform_parallel_processing=false;
% This is the number of processors to use in the parallel processing
parallel_processing_worker_number=1;

% This is a Boolean value stating whether to perform deformation PIV
perform_window_deformation=false;
% This is a scalar double specifying the minimum number of window
% deformation iterations to perform for each PIV pass
window_deformation_iteration_min=1;
% This is a scalar double specifying the maximum number of window
% deformation iterations to perform for each PIV pass
window_deformation_iteration_max=5;
% This is a scalar double specifying the threshhold magnitude of the
% difference of velocities in voxels between each iteration of window
% deformation below which the iteration is considered to be converged
% (ie once all the vectors change by less then this value between each
% iteration, the iterations will stop)
window_deformation_threshhold=0.5;

% This is a Boolean value stating whether to perform pyramid correlations
% on the images
perform_pyramid_correlations=false;
% This is a scalar double specifying the optimal number of frames used to
% calculate the velocity field using the pyramid correlation
pyramid_optimal_frame_number=3;
% This is a scalar double specifying the number of steps within the pyramid
% used to calculate the velocity field
pyramid_level_number=3;
% This is a string specifying the directory to save the temporary
% correlation volumes within
pyramid_correlation_write_directory='/home/rlafoy/Documents/Temp/';

% This saves the image read directory to the data processing structure
piv_parameters.general_parameters.image_read_directory=image_read_directory;
% This saves the image filename prefix to the data processing structure
piv_parameters.general_parameters.image_filename_prefix=image_filename_prefix;
% This saves the image filename extension to the data processing structure
piv_parameters.general_parameters.image_filename_extension=image_filename_extension;
% This saves the image variable name to the data processing structure
piv_parameters.general_parameters.image_variable_name=image_variable_name;
% This saves the vector field saving directory to the data processing structure
piv_parameters.general_parameters.vector_field_write_directory=vector_field_write_directory;
% This saves the initial image frame to the data processing structure
piv_parameters.general_parameters.initial_image_frame=initial_image_frame;
% This saves the final image frame to the data processing structure
piv_parameters.general_parameters.final_image_frame=final_image_frame;
% This saves the image frame step to the data processing structure
piv_parameters.general_parameters.image_frame_step=image_frame_step;
% This saves the image correlation step to the data processing structure
piv_parameters.general_parameters.image_correlation_step=image_correlation_step;
% This saves the particle image diameter to the data processing structure
piv_parameters.general_parameters.particle_diameter=particle_diameter;
% This saves the number of PIV passes to the data processing structure
piv_parameters.general_parameters.pass_number=pass_number;
% This saves the value stating whether to perform the processing in
% parallel to the data processing structure
piv_parameters.general_parameters.perform_parallel_processing=perform_parallel_processing;
% This saves the number of parallel workers used in the parallel processing
% to the data processing structure
piv_parameters.general_parameters.parallel_processing_worker_number=parallel_processing_worker_number;
% This saves the value stating whether to perform window deformation to the
% data processing structure
piv_parameters.general_parameters.perform_window_deformation=perform_window_deformation;
% This saves the minimum number of window deformations to the data processing
% structure
piv_parameters.general_parameters.window_deformation_iteration_min=window_deformation_iteration_min;
% This saves the maximum number of window deformation to the data processing
% structure
piv_parameters.general_parameters.window_deformation_iteration_max=window_deformation_iteration_max;
% This saves the image deformation threshhold criterion to the data
% processing structure
piv_parameters.general_parameters.window_deformation_threshhold=window_deformation_threshhold;
% This saves the value stating whether to perform pyramid correlations to the
% data processing structure
piv_parameters.general_parameters.perform_pyramid_correlations=perform_pyramid_correlations;
% This saves the optimal number or pyramid correlation frames to the data
% processing structure
piv_parameters.general_parameters.pyramid_optimal_frame_number=pyramid_optimal_frame_number;
% This saves the number of pyramid correlation levels to the data processing
% structure
piv_parameters.general_parameters.pyramid_level_number=pyramid_level_number;
% This saves the pyramid correlation write directory to the data processing
% structure
piv_parameters.general_parameters.pyramid_correlation_write_directory=pyramid_correlation_write_directory;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the 1st pass PIV processing parameters.                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a scalar double specifying the index of the current pass (this is here
% so that additional passes can be easily added without having to change the
% index of the 'pass_parameters' variable)
pass_index=1;

% This is a 1 x 3 vector of doubles specifying the effective size of the
% window to process after applying the Gaussian mask
window_resolution=[48,48,1];
% This is a 1 x 3 vector of doubles specifying the full size of the window
% including the mask
window_size=[96,96,1];
% This is a string that specifies how to define the location of the
% correlation windows and is equal to either 'window_overlap' or
% 'window_spacing'.
window_gridding_method='window_overlap';
% This is a 1 x 3 vector of doubles specifying the percent overlap of the
% correlation windows (which is redundant to specify along with the
% 'window_spacing' - the variable 'window_gridding_method' specifies which
% method to use)
window_overlap=[0.50,0.50,0.00];
% This is a 1 x 3 vector of doubles specifying the spacing between the
% correlation windows (which is redundant to specify along with the
% 'percent_overlap' - the variable 'window_gridding_method' specifies which
% method to use)
window_spacing=[64,64,64];

% This is a 1 x 3 vector of doubles that specify the initial bulk window
% offset (in voxels) to apply to the correlation windows (this would be used
% for cases where the user already new that that flow had a known zero
% mean component, ie in a pipe flow or a jet flow)
bulk_window_offset=[0,0,0];

% This is a string that specifies the method used to process the correlation
% volume data and can be a string equal to either 'SCC' or 'RPC' to specify
% using Standard Cross Correlation or Robust Phase Correlation
% respectively
correlation_method='RPC';

% This is a Boolean value stating whether to zero-mean the correlation windows
% prior to performing the correlation
zero_mean_windows=true;

% This is a Boolean value stating whether to validate the the vectors within
% the vector field of the current pass
validate_vector_field=true;

% These are 1 x 3 vectors of doubles specifying the minimum and maximum
% velocity values in voxels per frame below or above which the vectors are
% considered to be outliers
minimum_outlier_threshhold=[-2,-2,-2];
maximum_outlier_threshhold=[+2,+2,+2];

% This is a 1 x 3 vector specifying the size of the kernal to use for the UOD
% outlier detection (which currently must be equal to odd integers in each
% dimension))
uod_kernal_size=[3,3,1];
% This is a scalar double that specifies the minimum normalization level for
% the UOD (which corresponds to the expected uncertainty in the PIV
% measurement - typically about 0.1 for 2D PIV)
uod_epsilon=0.1;
% This is a scalar double that specifies the residual threshhold for the UOD
% measurement above which vectors are considered to be outliers
uod_residual_threshhold=1.5;

% This is a string specifying the method used to replace the identified
% outlier vectors.  The string may be equal to one of the following method
% strings:
%
%  Method                       2D Speed        3D Speed        Error
%   'local_mean'                 Slow            Moderate        High
%   'laplacian_interpolation'    Moderate        Fast            Low
%   'delaunay_interpolation'     Fast            Very Slow       Low
%
% where the 'Speed' refers to computation speed and the 'Error' refers to
% the interpolation error.
vector_replacement_method='laplacian_interpolation';
% This is the minimum number of valid vectors used to calculate the local 
% mean value used for vector replacement by 'local_mean'.  Depending upon
% replacement method chosen, this variable may not be used.
local_mean_minimum_valid_vector_number=8;
% This is the number of adjacent points used to calculate the vector
% replacement used by 'laplacian_interpolation'.  This number is based upon
% the connectivity with the outlier vector and can be either 4 or 8 in 2D 
% and 6, 18, or 26 in 3D.  Depending upon replacement method chosen, 
% this variable may not be used.
laplacian_interpolation_adjacent_connectivity=4;
% This is the interpolation method used during the vector replacement used
% by 'delaunay_interpolation'.  The method can be equal to 'natural',
% 'linear', or 'nearest' (for a description, open the 'TriScatteredInterp'
% or 'scatteredInterpolant' function references).  Depending upon 
% replacement method chosen, this variable may not be used.
delaunay_interpolation_weighting_method='natural';

% This is a Boolean value stating whether to smooth the vector field
smooth_vector_field=true;
% This is scalar double value that specifies the standard deviation of the
% Gaussian function (in units of vectors or windows equivalently) used to
% smooth the velocity field
gaussian_smoothing_kernal_std=1;
% This is a 1 x 3 vector of doubles that specify the size of the kernal, which
% should be odd in each dimension, used to smooth the velocity field (in units
% of vectors or windows equivalently)
gaussian_smoothing_kernal_size=[7,7,1];

% This adds the effective window resolution to the data processing structure
piv_parameters.pass_parameters(pass_index).window_resolution=window_resolution;
% This adds the full window size to the data processing structure
piv_parameters.pass_parameters(pass_index).window_size=window_size;
% This adds the method used to specify the window grid locations to the
% data processing structure
piv_parameters.pass_parameters(pass_index).window_gridding_method=window_gridding_method;
% This adds the ratio of overlap of the correlation windows to the data
% processing structure
piv_parameters.pass_parameters(pass_index).window_overlap=window_overlap;
% This adds the window spacing distance of the correlation windows to
% the data processing structure
piv_parameters.pass_parameters(pass_index).window_spacing=window_spacing;
% This adds the bulk window offset distance to the data processing structure
piv_parameters.pass_parameters(pass_index).bulk_window_offset=bulk_window_offset;
% This adds the correlation method to the data processing structure
piv_parameters.pass_parameters(pass_index).correlation_method=correlation_method;
% This adds the value stating whether to zero-mean the correlation windows to
% the data processing structure
piv_parameters.pass_parameters(pass_index).zero_mean_windows=zero_mean_windows;
% This adds the value stating whether to validate the vector field to the data
% processing structure
piv_parameters.pass_parameters(pass_index).validate_vector_field=validate_vector_field;
% This adds the minimum vector value outlier threshhold to the data processing
% structure
piv_parameters.pass_parameters(pass_index).minimum_outlier_threshhold=minimum_outlier_threshhold;
% This adds the maximum vector value outlier threshhold to the data processing
% structure
piv_parameters.pass_parameters(pass_index).maximum_outlier_threshhold=maximum_outlier_threshhold;
% This adds the UOD kernal size to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_kernal_size=uod_kernal_size;
% This adds the UOD expected error value to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_epsilon=uod_epsilon;
% This adds the UOD residual threshhold value to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_residual_threshhold=uod_residual_threshhold;
% This adds the outlier vector replacement method to the data processing
% structure
piv_parameters.pass_parameters(pass_index).vector_replacement_method=vector_replacement_method;
% This adds the local mean minimum valid vector number to the data
% processing structure
piv_parameters.pass_parameters(pass_index).local_mean_minimum_valid_vector_number=local_mean_minimum_valid_vector_number;
% This adds the Laplacian interpolation method adjacent point number to the
% data processing structure
piv_parameters.pass_parameters(pass_index).laplacian_interpolation_adjacent_connectivity=laplacian_interpolation_adjacent_connectivity;
% This adds the Delaunay interpolation method to the data processing
% structure
piv_parameters.pass_parameters(pass_index).delaunay_interpolation_weighting_method=delaunay_interpolation_weighting_method;
% This adds the value stating whether to smooth the vector field to the data
% processing structure
piv_parameters.pass_parameters(pass_index).smooth_vector_field=smooth_vector_field;
% This adds the Gaussian smoothing kernal size to the data processing structure
piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_std=gaussian_smoothing_kernal_std;
% This adds the Gaussian smoothing standard deviation value to the data
% processing structure
piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_size=gaussian_smoothing_kernal_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the 2nd pass PIV processing parameters.                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a scalar double specifying the index of the current pass (this is here
% so that additional passes can be easily added without having to change the
% index of the 'pass_parameters' variable)
pass_index=2;

% This is a 1 x 3 vector of doubles specifying the effective size of the
% window to process after applying the Gaussian mask
window_resolution=[48,48,1];
% This is a 1 x 3 vector of doubles specifying the full size of the window
% including the mask
window_size=[96,96,1];
% This is a string that specifies how to define the location of the
% correlation windows and is equal to either 'window_overlap' or
% 'window_spacing'.
window_gridding_method='window_overlap';
% This is a 1 x 3 vector of doubles specifying the percent overlap of the
% correlation windows (which is redundant to specify along with the
% 'window_spacing' - the variable 'window_gridding_method' specifies which
% method to use)
window_overlap=[0.50,0.50,0.00];
% This is a 1 x 3 vector of doubles specifying the spacing between the
% correlation windows (which is redundant to specify along with the
% 'percent_overlap' - the variable 'window_gridding_method' specifies which
% method to use)
window_spacing=[32,32,32];

% This is a 1 x 3 vector of doubles that specify the initial bulk window
% offset (in voxels) to apply to the correlation windows (this would be used
% for cases where the user already new that that flow had a known zero
% mean component, ie in a pipe flow or a jet flow)
bulk_window_offset=[0,0,0];

% This is a string that specifies the method used to process the correlation
% volume data and can be a string equal to either 'SCC' or 'RPC' to specify
% using Standard Cross Correlation or Robust Phase Correlation
% respectively
correlation_method='RPC';

% This is a Boolean value stating whether to zero-mean the correlation windows
% prior to performing the correlation
zero_mean_windows=true;

% This is a Boolean value stating whether to validate the the vectors within
% the vector field of the current pass
validate_vector_field=true;

% These are 1 x 3 vectors of doubles specifying the minimum and maximum
% velocity values in voxels per frame below or above which the vectors are
% considered to be outliers
minimum_outlier_threshhold=[-2,-2,-2];
maximum_outlier_threshhold=[+2,+2,+2];

% This is a 1 x 3 vector specifying the size of the kernal to use for the UOD
% outlier detection (which currently must be equal to odd integers in each
% dimension))
uod_kernal_size=[3,3,1];
% This is a scalar double that specifies the minimum normalization level for
% the UOD (which corresponds to the expected uncertainty in the PIV
% measurement - typically about 0.1 for 2D PIV)
uod_epsilon=0.1;
% This is a scalar double that specifies the residual threshhold for the UOD
% measurement above which vectors are considered to be outliers
uod_residual_threshhold=1.0;

% This is a string specifying the method used to replace the identified
% outlier vectors.  The string may be equal to one of the following method
% strings:
%
%  Method                       2D Speed        3D Speed        Error
%   'local_mean'                 Slow            Moderate        High
%   'laplacian_interpolation'    Moderate        Fast            Low
%   'delaunay_interpolation'     Fast            Very Slow       Low
%
% where the 'Speed' refers to computation speed and the 'Error' refers to
% the interpolation error.
vector_replacement_method='laplacian_interpolation';
% This is the minimum number of valid vectors used to calculate the local 
% mean value used for vector replacement by 'local_mean'.  Depending upon
% replacement method chosen, this variable may not be used.
local_mean_minimum_valid_vector_number=8;
% This is the number of adjacent points used to calculate the vector
% replacement used by 'laplacian_interpolation'.  This number is based upon
% the connectivity with the outlier vector and can be either 4 or 8 in 2D 
% and 6, 18, or 26 in 3D.  Depending upon replacement method chosen, 
% this variable may not be used.
laplacian_interpolation_adjacent_connectivity=4;
% This is the interpolation method used during the vector replacement used
% by 'delaunay_interpolation'.  The method can be equal to 'natural',
% 'linear', or 'nearest' (for a description, open the 'TriScatteredInterp'
% or 'scatteredInterpolant' function references).  Depending upon 
% replacement method chosen, this variable may not be used.
delaunay_interpolation_weighting_method='natural';

% This is a Boolean value stating whether to smooth the vector field
smooth_vector_field=true;
% This is scalar double value that specifies the standard deviation of the
% Gaussian function (in units of vectors or windows equivalently) used to
% smooth the velocity field
gaussian_smoothing_kernal_std=0.75;
% This is a 1 x 3 vector of doubles that specify the size of the kernal, which
% should be odd in each dimension, used to smooth the velocity field (in units
% of vectors or windows equivalently)
gaussian_smoothing_kernal_size=[7,7,7];

% This adds the effective window resolution to the data processing structure
piv_parameters.pass_parameters(pass_index).window_resolution=window_resolution;
% This adds the full window size to the data processing structure
piv_parameters.pass_parameters(pass_index).window_size=window_size;
% This adds the method used to specify the window grid locations to the
% data processing structure
piv_parameters.pass_parameters(pass_index).window_gridding_method=window_gridding_method;
% This adds the ratio of overlap of the correlation windows to the data
% processing structure
piv_parameters.pass_parameters(pass_index).window_overlap=window_overlap;
% This adds the window spacing distance of the correlation windows to
% the data processing structure
piv_parameters.pass_parameters(pass_index).window_spacing=window_spacing;
% This adds the bulk window offset distance to the data processing structure
piv_parameters.pass_parameters(pass_index).bulk_window_offset=bulk_window_offset;
% This adds the correlation method to the data processing structure
piv_parameters.pass_parameters(pass_index).correlation_method=correlation_method;
% This adds the value stating whether to zero-mean the correlation windows to
% the data processing structure
piv_parameters.pass_parameters(pass_index).zero_mean_windows=zero_mean_windows;
% This adds the value stating whether to validate the vector field to the data
% processing structure
piv_parameters.pass_parameters(pass_index).validate_vector_field=validate_vector_field;
% This adds the minimum vector value outlier threshhold to the data processing
% structure
piv_parameters.pass_parameters(pass_index).minimum_outlier_threshhold=minimum_outlier_threshhold;
% This adds the maximum vector value outlier threshhold to the data processing
% structure
piv_parameters.pass_parameters(pass_index).maximum_outlier_threshhold=maximum_outlier_threshhold;
% This adds the UOD kernal size to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_kernal_size=uod_kernal_size;
% This adds the UOD expected error value to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_epsilon=uod_epsilon;
% This adds the UOD residual threshhold value to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_residual_threshhold=uod_residual_threshhold;
% This adds the outlier vector replacement method to the data processing
% structure
piv_parameters.pass_parameters(pass_index).vector_replacement_method=vector_replacement_method;
% This adds the local mean minimum valid vector number to the data
% processing structure
piv_parameters.pass_parameters(pass_index).local_mean_minimum_valid_vector_number=local_mean_minimum_valid_vector_number;
% This adds the Laplacian interpolation method adjacent point number to the
% data processing structure
piv_parameters.pass_parameters(pass_index).laplacian_interpolation_adjacent_connectivity=laplacian_interpolation_adjacent_connectivity;
% This adds the Delaunay interpolation method to the data processing
% structure
piv_parameters.pass_parameters(pass_index).delaunay_interpolation_weighting_method=delaunay_interpolation_weighting_method;
% This adds the value stating whether to smooth the vector field to the data
% processing structure
piv_parameters.pass_parameters(pass_index).smooth_vector_field=smooth_vector_field;
% This adds the Gaussian smoothing kernal size to the data processing structure
piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_std=gaussian_smoothing_kernal_std;
% This adds the Gaussian smoothing standard deviation value to the data
% processing structure
piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_size=gaussian_smoothing_kernal_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the 3rd pass PIV processing parameters.                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a scalar double specifying the index of the current pass (this is here
% so that additional passes can be easily added without having to change the
% index of the 'pass_parameters' variable)
pass_index=3;

% This is a 1 x 3 vector of doubles specifying the effective size of the
% window to process after applying the Gaussian mask
window_resolution=[48,48,1];
% This is a 1 x 3 vector of doubles specifying the full size of the window
% including the mask
window_size=[96,96,1];
% This is a string that specifies how to define the location of the
% correlation windows and is equal to either 'window_overlap' or
% 'window_spacing'.
window_gridding_method='window_overlap';
% This is a 1 x 3 vector of doubles specifying the percent overlap of the
% correlation windows (which is redundant to specify along with the
% 'window_spacing' - the variable 'window_gridding_method' specifies which
% method to use)
window_overlap=[0.75,0.75,0.00];
% This is a 1 x 3 vector of doubles specifying the spacing between the
% correlation windows (which is redundant to specify along with the
% 'percent_overlap' - the variable 'window_gridding_method' specifies which
% method to use)
window_spacing=[16,16,16];

% This is a 1 x 3 vector of doubles that specify the initial bulk window
% offset (in voxels) to apply to the correlation windows (this would be used
% for cases where the user already new that that flow had a known zero
% mean component, ie in a pipe flow or a jet flow)
bulk_window_offset=[0,0,0];

% This is a string that specifies the method used to process the correlation
% volume data and can be a string equal to either 'SCC' or 'RPC' to specify
% using Standard Cross Correlation or Robust Phase Correlation
% respectively
correlation_method='RPC';

% This is a Boolean value stating whether to zero-mean the correlation windows
% prior to performing the correlation
zero_mean_windows=true;

% This is a Boolean value stating whether to validate the the vectors within
% the vector field of the current pass
validate_vector_field=true;

% These are 1 x 3 vectors of doubles specifying the minimum and maximum
% velocity values in voxels per frame below or above which the vectors are
% considered to be outliers
minimum_outlier_threshhold=[-2,-2,-2];
maximum_outlier_threshhold=[+2,+2,+2];

% This is a 1 x 3 vector specifying the size of the kernal to use for the UOD
% outlier detection (which currently must be equal to odd integers in each
% dimension))
uod_kernal_size=[3,3,1];
% This is a scalar double that specifies the minimum normalization level for
% the UOD (which corresponds to the expected uncertainty in the PIV
% measurement - typically about 0.1 for 2D PIV)
uod_epsilon=0.1;
% This is a scalar double that specifies the residual threshhold for the UOD
% measurement above which vectors are considered to be outliers
uod_residual_threshhold=1.0;

% This is a string specifying the method used to replace the identified
% outlier vectors.  The string may be equal to one of the following method
% strings:
%
%  Method                       2D Speed        3D Speed        Error
%   'local_mean'                 Slow            Moderate        High
%   'laplacian_interpolation'    Moderate        Fast            Low
%   'delaunay_interpolation'     Fast            Very Slow       Low
%
% where the 'Speed' refers to computation speed and the 'Error' refers to
% the interpolation error.
vector_replacement_method='laplacian_interpolation';
% This is the minimum number of valid vectors used to calculate the local 
% mean value used for vector replacement by 'local_mean'.  Depending upon
% replacement method chosen, this variable may not be used.
local_mean_minimum_valid_vector_number=8;
% This is the number of adjacent points used to calculate the vector
% replacement used by 'laplacian_interpolation'.  This number is based upon
% the connectivity with the outlier vector and can be either 4 or 8 in 2D 
% and 6, 18, or 26 in 3D.  Depending upon replacement method chosen, 
% this variable may not be used.
laplacian_interpolation_adjacent_connectivity=4;
% This is the interpolation method used during the vector replacement used
% by 'delaunay_interpolation'.  The method can be equal to 'natural',
% 'linear', or 'nearest' (for a description, open the 'TriScatteredInterp'
% or 'scatteredInterpolant' function references).  Depending upon 
% replacement method chosen, this variable may not be used.
delaunay_interpolation_weighting_method='natural';

% This is a Boolean value stating whether to smooth the vector field
smooth_vector_field=false;
% This is scalar double value that specifies the standard deviation of the
% Gaussian function (in units of vectors or windows equivalently) used to
% smooth the velocity field
gaussian_smoothing_kernal_std=1;
% This is a 1 x 3 vector of doubles that specify the size of the kernal, which
% should be odd in each dimension, used to smooth the velocity field (in units
% of vectors or windows equivalently)
gaussian_smoothing_kernal_size=[7,7,1];

% This adds the effective window resolution to the data processing structure
piv_parameters.pass_parameters(pass_index).window_resolution=window_resolution;
% This adds the full window size to the data processing structure
piv_parameters.pass_parameters(pass_index).window_size=window_size;
% This adds the method used to specify the window grid locations to the
% data processing structure
piv_parameters.pass_parameters(pass_index).window_gridding_method=window_gridding_method;
% This adds the ratio of overlap of the correlation windows to the data
% processing structure
piv_parameters.pass_parameters(pass_index).window_overlap=window_overlap;
% This adds the window spacing distance of the correlation windows to
% the data processing structure
piv_parameters.pass_parameters(pass_index).window_spacing=window_spacing;
% This adds the bulk window offset distance to the data processing structure
piv_parameters.pass_parameters(pass_index).bulk_window_offset=bulk_window_offset;
% This adds the correlation method to the data processing structure
piv_parameters.pass_parameters(pass_index).correlation_method=correlation_method;
% This adds the value stating whether to zero-mean the correlation windows to
% the data processing structure
piv_parameters.pass_parameters(pass_index).zero_mean_windows=zero_mean_windows;
% This adds the value stating whether to validate the vector field to the data
% processing structure
piv_parameters.pass_parameters(pass_index).validate_vector_field=validate_vector_field;
% This adds the minimum vector value outlier threshhold to the data processing
% structure
piv_parameters.pass_parameters(pass_index).minimum_outlier_threshhold=minimum_outlier_threshhold;
% This adds the maximum vector value outlier threshhold to the data processing
% structure
piv_parameters.pass_parameters(pass_index).maximum_outlier_threshhold=maximum_outlier_threshhold;
% This adds the UOD kernal size to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_kernal_size=uod_kernal_size;
% This adds the UOD expected error value to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_epsilon=uod_epsilon;
% This adds the UOD residual threshhold value to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_residual_threshhold=uod_residual_threshhold;
% This adds the outlier vector replacement method to the data processing
% structure
piv_parameters.pass_parameters(pass_index).vector_replacement_method=vector_replacement_method;
% This adds the local mean minimum valid vector number to the data
% processing structure
piv_parameters.pass_parameters(pass_index).local_mean_minimum_valid_vector_number=local_mean_minimum_valid_vector_number;
% This adds the Laplacian interpolation method adjacent point number to the
% data processing structure
piv_parameters.pass_parameters(pass_index).laplacian_interpolation_adjacent_connectivity=laplacian_interpolation_adjacent_connectivity;
% This adds the Delaunay interpolation method to the data processing
% structure
piv_parameters.pass_parameters(pass_index).delaunay_interpolation_weighting_method=delaunay_interpolation_weighting_method;
% This adds the value stating whether to smooth the vector field to the data
% processing structure
piv_parameters.pass_parameters(pass_index).smooth_vector_field=smooth_vector_field;
% This adds the Gaussian smoothing kernal size to the data processing structure
piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_std=gaussian_smoothing_kernal_std;
% This adds the Gaussian smoothing standard deviation value to the data
% processing structure
piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_size=gaussian_smoothing_kernal_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the 4th pass PIV processing parameters.                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a scalar double specifying the index of the current pass (this is here
% so that additional passes can be easily added without having to change the
% index of the 'pass_parameters' variable)
pass_index=4;

% This is a 1 x 3 vector of doubles specifying the effective size of the
% window to process after applying the Gaussian mask
window_resolution=[16,16,1];
% This is a 1 x 3 vector of doubles specifying the full size of the window
% including the mask
window_size=[32,32,1];
% This is a string that specifies how to define the location of the
% correlation windows and is equal to either 'window_overlap' or
% 'window_spacing'.
window_gridding_method='window_overlap';
% This is a 1 x 3 vector of doubles specifying the percent overlap of the
% correlation windows (which is redundant to specify along with the
% 'window_spacing' - the variable 'window_gridding_method' specifies which
% method to use)
window_overlap=[0.75,0.75,0.00];
% This is a 1 x 3 vector of doubles specifying the spacing between the
% correlation windows (which is redundant to specify along with the
% 'percent_overlap' - the variable 'window_gridding_method' specifies which
% method to use)
window_spacing=[16,16,16];

% This is a 1 x 3 vector of doubles that specify the initial bulk window
% offset (in voxels) to apply to the correlation windows (this would be used
% for cases where the user already new that that flow had a known zero
% mean component, ie in a pipe flow or a jet flow)
bulk_window_offset=[0,0,0];

% This is a string that specifies the method used to process the correlation
% volume data and can be a string equal to either 'SCC' or 'RPC' to specify
% using Standard Cross Correlation or Robust Phase Correlation
% respectively
correlation_method='RPC';

% This is a Boolean value stating whether to zero-mean the correlation windows
% prior to performing the correlation
zero_mean_windows=true;

% This is a Boolean value stating whether to validate the the vectors within
% the vector field of the current pass
validate_vector_field=true;

% These are 1 x 3 vectors of doubles specifying the minimum and maximum
% velocity values in voxels per frame below or above which the vectors are
% considered to be outliers
minimum_outlier_threshhold=[-8,-8,-8];
maximum_outlier_threshhold=[+8,+8,+8];

% This is a 1 x 3 vector specifying the size of the kernal to use for the UOD
% outlier detection (which currently must be equal to odd integers in each
% dimension))
uod_kernal_size=[3,3,3];
% This is a scalar double that specifies the minimum normalization level for
% the UOD (which corresponds to the expected uncertainty in the PIV
% measurement - typically about 0.1 for 2D PIV)
uod_epsilon=0.1;
% This is a scalar double that specifies the residual threshhold for the UOD
% measurement above which vectors are considered to be outliers
uod_residual_threshhold=2;

% This is a string specifying the method used to replace the identified
% outlier vectors.  The string may be equal to one of the following method
% strings:
%
%  Method                       2D Speed        3D Speed        Error
%   'local_mean'                 Slow            Moderate        High
%   'laplacian_interpolation'    Moderate        Fast            Low
%   'delaunay_interpolation'     Fast            Very Slow       Low
%
% where the 'Speed' refers to computation speed and the 'Error' refers to
% the interpolation error.
vector_replacement_method='laplacian_interpolation';
% This is the minimum number of valid vectors used to calculate the local 
% mean value used for vector replacement by 'local_mean'.  Depending upon
% replacement method chosen, this variable may not be used.
local_mean_minimum_valid_vector_number=8;
% This is the number of adjacent points used to calculate the vector
% replacement used by 'laplacian_interpolation'.  This number is based upon
% the connectivity with the outlier vector and can be either 4 or 8 in 2D 
% and 6, 18, or 26 in 3D.  Depending upon replacement method chosen, 
% this variable may not be used.
laplacian_interpolation_adjacent_connectivity=6;
% This is the interpolation method used during the vector replacement used
% by 'delaunay_interpolation'.  The method can be equal to 'natural',
% 'linear', or 'nearest' (for a description, open the 'TriScatteredInterp'
% or 'scatteredInterpolant' function references).  Depending upon 
% replacement method chosen, this variable may not be used.
delaunay_interpolation_weighting_method='natural';

% This is a Boolean value stating whether to smooth the vector field
smooth_vector_field=false;
% This is scalar double value that specifies the standard deviation of the
% Gaussian function (in units of vectors or windows equivalently) used to
% smooth the velocity field
gaussian_smoothing_kernal_std=0.5;
% This is a 1 x 3 vector of doubles that specify the size of the kernal, which
% should be odd in each dimension, used to smooth the velocity field (in units
% of vectors or windows equivalently)
gaussian_smoothing_kernal_size=[7,7,7];

% This adds the effective window resolution to the data processing structure
piv_parameters.pass_parameters(pass_index).window_resolution=window_resolution;
% This adds the full window size to the data processing structure
piv_parameters.pass_parameters(pass_index).window_size=window_size;
% This adds the method used to specify the window grid locations to the
% data processing structure
piv_parameters.pass_parameters(pass_index).window_gridding_method=window_gridding_method;
% This adds the ratio of overlap of the correlation windows to the data
% processing structure
piv_parameters.pass_parameters(pass_index).window_overlap=window_overlap;
% This adds the window spacing distance of the correlation windows to
% the data processing structure
piv_parameters.pass_parameters(pass_index).window_spacing=window_spacing;
% This adds the bulk window offset distance to the data processing structure
piv_parameters.pass_parameters(pass_index).bulk_window_offset=bulk_window_offset;
% This adds the correlation method to the data processing structure
piv_parameters.pass_parameters(pass_index).correlation_method=correlation_method;
% This adds the value stating whether to zero-mean the correlation windows to
% the data processing structure
piv_parameters.pass_parameters(pass_index).zero_mean_windows=zero_mean_windows;
% This adds the value stating whether to validate the vector field to the data
% processing structure
piv_parameters.pass_parameters(pass_index).validate_vector_field=validate_vector_field;
% This adds the minimum vector value outlier threshhold to the data processing
% structure
piv_parameters.pass_parameters(pass_index).minimum_outlier_threshhold=minimum_outlier_threshhold;
% This adds the maximum vector value outlier threshhold to the data processing
% structure
piv_parameters.pass_parameters(pass_index).maximum_outlier_threshhold=maximum_outlier_threshhold;
% This adds the UOD kernal size to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_kernal_size=uod_kernal_size;
% This adds the UOD expected error value to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_epsilon=uod_epsilon;
% This adds the UOD residual threshhold value to the data processing structure
piv_parameters.pass_parameters(pass_index).uod_residual_threshhold=uod_residual_threshhold;
% This adds the outlier vector replacement method to the data processing
% structure
piv_parameters.pass_parameters(pass_index).vector_replacement_method=vector_replacement_method;
% This adds the local mean minimum valid vector number to the data
% processing structure
piv_parameters.pass_parameters(pass_index).local_mean_minimum_valid_vector_number=local_mean_minimum_valid_vector_number;
% This adds the Laplacian interpolation method adjacent point number to the
% data processing structure
piv_parameters.pass_parameters(pass_index).laplacian_interpolation_adjacent_connectivity=laplacian_interpolation_adjacent_connectivity;
% This adds the Delaunay interpolation method to the data processing
% structure
piv_parameters.pass_parameters(pass_index).delaunay_interpolation_weighting_method=delaunay_interpolation_weighting_method;
% This adds the value stating whether to smooth the vector field to the data
% processing structure
piv_parameters.pass_parameters(pass_index).smooth_vector_field=smooth_vector_field;
% This adds the Gaussian smoothing kernal size to the data processing structure
piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_std=gaussian_smoothing_kernal_std;
% This adds the Gaussian smoothing standard deviation value to the data
% processing structure
piv_parameters.pass_parameters(pass_index).gaussian_smoothing_kernal_size=gaussian_smoothing_kernal_size;



function validate_processing_parameters(piv_parameters);
% This function is used to test whether the data structure 'piv_parameters'
% was set up correctly and to initialize any variables that will be needed
% during the processing.

% This tests whether the image directory exists and displays an error if not
if exist(piv_parameters.general_parameters.image_read_directory,'dir');
    
    % This is a list of the images within the directory that have the user
    % specified prefix and extension
    image_filename_list=dir([piv_parameters.general_parameters.image_read_directory,piv_parameters.general_parameters.image_filename_prefix,'*',piv_parameters.general_parameters.image_filename_extension]);
    
    % This is the number of images within the directory
    image_filename_number=length(image_filename_list);
    
    % This checks if the range of images contains the images specified by
    % the user to be processed and if not displays an error
    if (image_filename_number<(piv_parameters.general_parameters.final_image_frame+piv_parameters.general_parameters.image_correlation_step));
        % This displays an error stating that all the images cannot be found
        error('The full range of images to be processed were not found within the image read directory.');
    end;
    
else;
    
    % This displays a warning stating that the image directory does not exist
    warning('The image read directory does not exist.  This may be due to the ''piv_parameters'' data structure being created on a different computer than the one on which the data will be processed.');
    
end;

% This checks whether the vector field write directory exists and if not displays an error
if not(exist(piv_parameters.general_parameters.vector_field_write_directory,'dir'));
    
    % This displays a warning stating that the vector field write directory does
    % not exist
    warning('The vector field write directory does not exist.  This may be due to the ''piv_parameters'' data structure being created on a different computer than the one on which the data will be processed.');
    
end;

% If the user specified that window deformation should be performed, this checks whether the
% iteration parameters are logical and whether the threshhold is well defined
if piv_parameters.general_parameters.perform_window_deformation;
    
    % This displays an error stating that window deformation has yet to be
    % properly programmed into the PIV code
    error('Window deformation has not yet been implemented.  This feature will be added soon . . . ');
    
    % This checks at least one window deformation iteration will be applied
    if (piv_parameters.general_parameters.window_deformation_iteration_min<1);
        % This displays a warning stating that at least one window deformation must
        % be performed
        warning('Setting the window deformation iteration minimum to less then one will result in standard correlation operations being applied.');
    end;
    
    % This checks that the maximum number of iterations is equal to or greater than the
    % minimum number of iterations
    if (piv_parameters.general_parameters.window_deformation_iteration_min>piv_parameters.general_parameters.window_deformation_iteration_max);
        % This displays an error stating that the minimum number of window deformations
        % must be less than or equal to the maximum number of deformations
        error('The minimum number of window deformations must be less then or equal to the maximum number of window deformations');
    end;
    
end;

% If the user specifies that pyramid correlations should be performed, this checks whether the
% pyramid level is less then or equal to pyramid correlation frame number and whether the
% temporary write directory exists
if piv_parameters.general_parameters.perform_pyramid_correlations;
    
    % This checks that the pyramid optimal frame number is at least equal to one
    if (piv_parameters.general_parameters.pyramid_optimal_frame_number<1);
        % This displays an error stating that the optimal frame number must be equal
        % to at least one
        error('The optimal frame number use in the pyramid correlations must be equal to or greater than one.');
    elseif (piv_parameters.general_parameters.pyramid_optimal_frame_number==1);
        % This displays a warning stating that the pyramid correlation will be equivalent
        % to a standard correlation
        warning('Setting the optimal frame number to one in the pyramid correlations will result in standard pair correlations.');
    end;
    
    % This checks whether the pyramid level is less then or equal to the pyramid correlation
    % frame number and if not displays an error
    if piv_parameters.general_parameters.pyramid_level_number>piv_parameters.general_parameters.pyramid_optimal_frame_number;
        % This displays an error stating that the level number must be less then or
        % equal to the pyramid corraltion frame number
        error('The pyramid level number must be less then or equal to the optimal frame number.');
    end;
    
    % This checks whether the temporary pyramid write directory exists and if not, displays
    % that an error occured
    if not(exist(piv_parameters.general_parameters.pyramid_correlation_write_directory,'dir'));
        % This displays a warning that the temporary pyramid correlation write directory
        % does not exist
        warning('The pyramid correlation temporary write directory does not exist.  This may be due to the ''piv_parameters'' data structure being created on a different computer than the one on which the data will be processed.');
    end;
    
end;

% This is the number of passes to perform during the PIV processing
pass_number=piv_parameters.general_parameters.pass_number;

% This iterates through the PIV passes, checking that proper processing
% parameters were entered
for pass_index=1:pass_number;
    
    % This is a Boolean value stating whether to perform vector validation
    % during the current pass
    validate_vector_field=piv_parameters.pass_parameters(pass_index).validate_vector_field;

    % If the user specifies that vector validation should be performed, this
    % checks to make sure that valid options are specified
    if validate_vector_field;
        
        % This extracts the outlier vector replacement method
        vector_replacement_method=piv_parameters.pass_parameters(pass_index).vector_replacement_method;
        
        % If the vector replacement method is the Laplacian interpolation,
        % this checks that a valid number of adjacent vectors is entered
        if strcmp(vector_replacement_method,'laplacian_interpolation');
            % This extracts the number of adjacent vectors to use in
            % calculating the interpolated vectors
            laplacian_interpolation_adjacent_connectivity=piv_parameters.pass_parameters(pass_index).laplacian_interpolation_adjacent_connectivity;
            % This checks whether the adjacent connectivity number is one
            % of the valid options
            if not(any(laplacian_interpolation_adjacent_connectivity==[4,6,8,18,26]));
                % This displays an error stating that an incorrect number
                % of adjacent points was specified
                error('The specified number of adjacent points used in the Laplacian interpolation for vector replacement is not a valid value.  The number of adjacent points must be 4, 6, 8, 18, or 26.');
            end;
        end;
        
        % If the vector replacement method is the Delaunay interpolation,
        % this checks that a valid weighting method is entered
        if strcmp(vector_replacement_method,'delaunay_interpolation');
            % This extracts the Delaunay interpolation weighting method
            delaunay_interpolation_weighting_method=piv_parameters.pass_parameters(pass_index).delaunay_interpolation_weighting_method;
            % This checks whether the weighting method is one of the valid
            % methods
            if not(strcmp(delaunay_interpolation_weighting_method,'natural'))&&not(strcmp(delaunay_interpolation_weighting_method,'linear'))&&not(strcmp(delaunay_interpolation_weighting_method,'nearest'));
                % This displays an error stating that an incorrect Delaunay
                % interpolation weighting method was specified
                error('The specified weighting method used in the Delaunay interpolation for vector replacement is not a valid value.  The method must be one of the options ''natural'', ''linear'', or ''nearest''.');
            end;
        end;
        
    end;
    
end;
        
        


