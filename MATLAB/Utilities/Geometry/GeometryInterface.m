classdef (Abstract) GeometryInterface
    
    % Interface for geometry class to work with all TIGRE functionality
    
   properties (Abstract) %These must be defined before use
        % VARIABLE                                   DESCRIPTION                    UNITS
        %-------------------------------------------------------------------------------------
        % Distances
        DSD                                     % Distance Source Detector      (mm)
        DSO                                     % Distance Source Origin        (mm)
        % Detector parameters
        nDetector            					% number of pixels              (px)
        dDetector           					% size of each pixel            (mm)
        sDetector                               % total size of the detector    (mm)
        % Image parameters
        nVoxel                                  % number of voxels              (vx)
        sVoxel                                  % total size of the image       (mm)
        dVoxel                                  % size of each voxel            (mm)
        % Offsets
        offOrigin                               % Offset of image from origin   (mm)              
        offDetector                             % Offset of Detector            (mm)
                                                % These two can be also defined
                                                % per angle
   end
    
    methods
        function geo = checkGeo(obj, angles)
            % must convert properties into a struct before passing
            % to a mex function
            geo = struct(obj);
            geo = checkGeo(geo, angles);
        end
    end
    
    methods (Access = protected)
       function geo = struct(obj)
            fields = fieldnames(obj);
            for f = 1:length(fields)
                geo.(fields{f}) = obj.(fields{f});
            end
        end 
    end
end

