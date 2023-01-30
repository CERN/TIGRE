function bool = apply_wang_weights(geo)
    if (size(geo.offDetector,2) > 1) && length(unique(geo.offDetector(1,:)))>1
        warning('Wang weights: varying offDetector detected, Wang weights not being applied');
        bool = false;
        return
    end
    
    if geo.offDetector(1) == 0
        bool = false;
        return
    end
    
    if (numel(geo.DSO) > 1) && (length(unique(geo.DSO))>1)
        warning('Wang weights: varying DSO detected, Wang weights not being applied');
        bool = false;
        return
    end

    percent_offset = abs(geo.offDetector(1)/geo.sDetector(1)) * 100;    
    if percent_offset > 30
        warning('Wang weights: Detector offset percent: %0.2f is greater than 30 which may result in image artifacts, consider rebinning 360 degree projections to 180 degrees', percent_offset)
    end
    
    bool = true;
end