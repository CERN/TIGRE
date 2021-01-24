classdef GpuIds
    %GpuIds
    %Usage:
    %  To select GPUIDs of 'GeForce RTX 2080 Ti',
    %    gpuids = GpuIds('GeForce RTX 2080 Ti');
    %  To select GPUIDs of all installed GPUs,
    %    gpuids = GpuIds();
    
    properties
        name;
        devices;
    end
    methods
        function obj = GpuIds(varargin)
            obj.name = '';
            if nargin > 0
                obj.name = varargin{1};
            end
            if isempty(obj.name)
                %disp('GpuIds(none)');
                deviceCount = getGpuCount_mex();
                %fprintf('devCnt = %d\n', deviceCount);
                obj.devices = int32(0:(deviceCount-1));
            else
                obj.devices = getGpuIds(obj.name);
            end
        end
        
        function len=length(obj)
            %Returns the number of selected GPUs.
            len = length(obj.devices); 
        end
    end
end
