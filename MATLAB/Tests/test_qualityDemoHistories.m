function test_qualityDemoHistories()
% Execute only demo plotting with histories of unequal lengths in any order.
root=fileparts(fileparts(mfilename('fullpath')));
source=fileread(fullfile(root,'Demos','d10_QualityMeasurements.m'));
first=strfind(source,'%%RMSE plot');
last=strfind(source,'%% Compute the parameters for the target image');
plotting=source(first(1):last(1)-1);
previousVisibility=get(groot,'DefaultFigureVisible');
cleanup=onCleanup(@() restoreFigures(previousVisibility));
set(groot,'DefaultFigureVisible','off');
names={'SIRT','SART','OS-SART'};
metrics={'RMSE','CC','MSSIM','UQI'};
for lengths={[8 3 11],[11 8 3],[3 11 8],[0 3 1]}
    close('all');
    n=lengths{1};
    qualitySIRT=reshape(1:4*n(1),4,n(1));
    qualitySART=reshape(1:4*n(2),4,n(2))+100;
    qualityOSSART=reshape(1:4*n(3),4,n(3))+200;
    histories={qualitySIRT,qualitySART,qualityOSSART};
    eval(plotting);
    figures=findall(groot,'Type','figure');
    assert(numel(figures)==4,'Expected one plot per metric');
    for figureHandle=figures'
        ax=findobj(figureHandle,'Type','axes');
        row=find(cellfun(@(name) ~isempty(strfind(ax.Title.String,name)),metrics));
        assert(isscalar(row),'Expected a unique metric title');
        lines=findobj(ax,'Type','line');
            assert(numel(lines)==nnz(n),'Expected one line per nonempty history');
        for line=lines'
            index=find(strcmp(line.DisplayName,names));
            assert(isscalar(index),'Expected the correct algorithm label');
            assert(isequal(line.XData,1:n(index)));
            assert(isequal(line.YData,histories{index}(row,:)));
        end
    end
end
fprintf('PASS quality demo: unequal and empty histories retain metric rows and iteration counts\n');
end

function restoreFigures(visibility)
close('all');
set(groot,'DefaultFigureVisible',visibility);
end
