function [ordered_alpha,index_alpha]=order_subsets(angles,blocksize, mode)
index_alpha=1:length(angles);

block_alpha=mat2cell(angles      ,1,[repmat(blocksize,1,floor(length(angles)/blocksize)) mod(length(angles),blocksize)]);
index_alpha=mat2cell(index_alpha,1,[repmat(blocksize,1,floor(length(angles)/blocksize)) mod(length(angles),blocksize)]);

block_alpha=block_alpha(~cellfun('isempty',block_alpha));
index_alpha=index_alpha(~cellfun('isempty',index_alpha));


if strcmp(mode,'ordered')
    ordered_alpha=block_alpha;
    return;
end
if strcmp(mode,'random')
    neworder=randperm(length(block_alpha));
    ordered_alpha=block_alpha(neworder);
    index_alpha=index_alpha(neworder);
    return;
end
%% not finished


if strcmp(mode,'angularDistance') 
    % we need them sorted, so we need to recompute the blocks, but sorted
    [angles,sortindex]=sort(angles);
    index_alpha=1:length(angles);
    index_alpha=index_alpha(sortindex);
    
    block_alpha=mat2cell(angles      ,1,[repmat(blocksize,1,floor(length(angles)/blocksize)) mod(length(angles),blocksize)]);
    index_alpha=mat2cell(index_alpha,1,[repmat(blocksize,1,floor(length(angles)/blocksize)) mod(length(angles),blocksize)]);
    
    block_alpha=block_alpha(~cellfun('isempty',block_alpha));
    index_alpha=index_alpha(~cellfun('isempty',index_alpha));
    
    
    
    avrg=cellfun(@mean,block_alpha);
    used_avrg=[];
    % start from the beggining
    ordered_alpha{1}=block_alpha{1};
    auxindex_alpha=index_alpha;
    index_alpha{1}=auxindex_alpha{1};
    used_avrg(end+1)=avrg(1);
    for ii=2:length(block_alpha)
        dist=[];
        for jj=1:length(used_avrg)
            dist(jj,:)=abs(mod((avrg- used_avrg(jj))+pi,2*pi)-pi);
        end
        dist=bsxfun(@times,dist,all(dist,1));
        [~,midx]=max(dist(:));
        [~,avrgindx]=ind2sub(size(dist),midx);
        index_alpha{ii}=auxindex_alpha{avrgindx};
        ordered_alpha{ii}=block_alpha{avrgindx};
        used_avrg(end+1)=avrg(avrgindx);
    end
end

end