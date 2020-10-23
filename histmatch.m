% HISTMATCH Histogram matching for images
%
%   tgt_adj = HISTMATCH(tgt,ref,...)
%   Return a histogram-matched version tgt_adj of target image tgt based on
%   a provided reference image ref.
%
%   Which method to use is provided by name/value pair argument:
%
%   HISTMATCH(...,'method',method)
%    'full' - Full histogram match based on the cumulative distribution
%    function of the images
%    'analytic' - Full histogram match to a an analytic reference
%    cumulative distribution function
%    'adaptive' - Full histogram match using multiple image regions, using
%    interpolation between the regions
%    'adaptive_all' - Full histogram matching using a neighborhood around
%    every pixel
%    'point' - Point-wise rank-order histogram match
%    'partial' - Partial histogram match by affine transform
%   (default: 'partial')
%
%   Additional name/value pair arguments depend on the selected method
%   'full':
%       HISTMATCH(...,'roi', roi) or
%       HISTMATCH(...,'roi_ref',roi_ref,'roi_tgt',roi_tgt)
%       Determine the matching function using only a region of interest within
%       each image. If only roi is provided, the same is used for both images.
%       The matching function is still applied to the entire image and values
%       outside the reference range are extrapolated by default.
%       
%       HISTMATCH(...,'clip',clip)
%       'none' - Do not clip the output (default)
%       'roi' - Clip output to range in the reference ROI
%       'all' - Clip output to range in the entire reference image
%
%       HISTMATCH(...,'bins',bins)
%       Set the number of bins used in calculation of the cumulative
%       distribution function. (default: 256 bins)
%
%       HISTMATCH(...,'limit',limit)
%       Limit the maximum contribution in a given histogram bin, normalized
%       (0,1). (default: 1)
%
%   'analytic'
%       Note: ref should be provided as a structure with fields 'cdf'
%       and 'x' to specify the reference distribution function and
%       its bin centers
%
%       HISTMATCH(...,'roi', roi) or
%       HISTMATCH(...,'roi_tgt',roi_tgt)
%       (see above)
%
%       HISTMATCH(...,'clip',clip)
%       'none' - Do not clip the output (default)
%       'all' - Clip output to range in the entire reference image
%
%       HISTMATCH(...,'limit',limit)
%       (see above)
%
%   'adaptive'
%       HISTMATCH(...,'regions',regions)
%       2-element vector to specify the number of axial and lateral regions
%       in the matching grid (default: [5,5])
%
%       HISTMATCH(...,'clip',clip)
%       'none' - Do not clip the output
%       'roi' - Clip output to range in each reference ROI (default)
%       'all' - Clip output to range in the entire reference image 
%
%       HISTMATCH(...,'bins',bins)
%       Set the number of bins used in calculation of the cumulative
%       distribution function. (default: 256 bins)
%
%       HISTMATCH(...,'limit',limit)
%       (see above)
%
%   'adaptive_all'
%       HISTMATCH(...,'neighborhood',neighborhood)
%       2-element vector to specify the number of axial and lateral pixels
%       in each neighborhood around a pixel (default: [11,11])
%
%       HISTMATCH(...,'bins',bins)
%       Set the number of bins used in calculation of the cumulative
%       distribution function. (default: 10 bins)
%
%       HISTMATCH(...,'limit',limit)
%       (see above)
%
%   'point'
%       (no parameters)
%
%   'partial'
%       HISTMATCH(...,'version', version)
%       Apply an affine transformation to optimize one of the following:
%       'mu_sigma' - Match the distribution mean and variance (default)
%       'L2' - Minimize the L2 norm of the images
%
%       HISTMATCH(...,'roi', roi) or
%       HISTMATCH(...,'roi_ref',roi_ref,'roi_tgt',roi_tgt)
%       (see above)
%       
%       HISTMATCH(...,'clip',clip)
%       'none' - Do not clip the output (default)
%       'roi' - Clip output to range in the reference ROI
%       'all' - Clip output to range in the entire reference image
%   
function [tgt_adj] = histmatch(tgt,ref,varargin)

args = parse_inputs(tgt, ref, varargin); 
switch args.method
    
    % METHOD 1: Full histogram match
    case 'full' 
        [cdf_ref, centers_ref] = calculate_cdf(ref(args.roi_ref), args.bins, args.limit);
        [cdf_tgt, centers_tgt] = calculate_cdf(tgt(args.roi_tgt), args.bins, args.limit);
        M = get_transform(cdf_tgt, cdf_ref, centers_ref);
        tgt_adj = interp1(centers_tgt, M, tgt, 'linear', 'extrap');
        if(strcmp(args.clip, 'all')), tgt_adj = clip_data(tgt_adj, ref);
        elseif(strcmp(args.clip, 'roi')), tgt_adj = clip_data(tgt_adj, ref(args.roi_ref));
        end
    
    % METHOD 2: Analytic histogram match
    case 'analytic'
        centers_ref=linspace(ref.x(1),ref.x(end),length(ref.x));
        cdf_ref=interp1(ref.x,ref.cdf,centers_ref)/max(ref.cdf);
        [cdf_tgt, centers_tgt] = calculate_cdf(tgt(args.roi_tgt), args.bins, args.limit);
        M = get_transform(cdf_tgt, cdf_ref, centers_ref);
        tgt_adj = interp1(centers_tgt, M, tgt, 'linear', 'extrap');
        if(strcmp(args.clip, 'all')), tgt_adj = clip_data(tgt_adj, ref.x);
        end
        
    % METHOD 3: Adaptive histogram match using interpolation
    case 'adaptive'
        [rois, s] = grid_rois(size(tgt), args.regions);
        tgt_all = zeros([size(tgt), args.regions], 'like', tgt);
        for i=1:prod(args.regions)
            [cdf_ref, centers_ref] = calculate_cdf(ref(rois(:,:,i)), args.bins, args.limit);
            [cdf_tgt, centers_tgt] = calculate_cdf(tgt(rois(:,:,i)), args.bins, args.limit);
            M = get_transform(cdf_tgt, cdf_ref, centers_ref);
            tgt_all(:,:,i) = interp1(centers_tgt, M, tgt, 'linear', 'extrap');
            if(strcmp(args.clip, 'all')), tgt_all(:,:,i) = clip_data(tgt_all(:,:,i), ref);
            elseif(strcmp(args.clip, 'roi')), tgt_all(:,:,i) = clip_data(tgt_all(:,:,i), ref(rois(:,:,i)));
            end
        end
        tgt_adj = interp_regions(tgt_all, s);
        
    % METHOD 4: Adaptive histogram match for every point in the image
    case 'adaptive_all'
        tgt_adj=zeros(size(tgt),'like',tgt);
        for i=1:size(tgt,1)
            for j=1:size(tgt,2)
                roi = get_neighborhood(i,j,args.neighborhood,size(tgt));
                [cdf_ref, centers_ref] = calculate_cdf(ref(roi), args.bins, args.limit);
                [cdf_tgt, centers_tgt] = calculate_cdf(tgt(roi), args.bins, args.limit);
                M = get_transform(cdf_tgt, cdf_ref, centers_ref);
                tgt_adj(i,j) = interp1(centers_tgt, M, tgt(i,j),'linear','extrap');
            end
        end
        
    % METHOD 5: Pointwise (rank-order) histogram match    
    case 'point'
        [~,Iref]=sort(ref(:));
        [~,Itgt]=sort(tgt(:));
        tgt_adj=zeros(size(tgt),'like',tgt);
        tgt_adj(Itgt)=ref(Iref);
        
    % METHOD 6: Partial (mean and variance) match
    case 'partial'
        x=tgt(args.roi_tgt);
        y=ref(args.roi_ref);
        switch args.version
            case 'mu_sigma'
                a=std(y)/std(x);
                b=mean(y)-a*mean(x);
            case 'L2'
                a=(mean(x.*y)-mean(x)*mean(y))/var(x);
                b=mean(y)-a*mean(x);
        end
        tgt_adj=a*tgt+b;
        if(strcmp(args.clip, 'all')), tgt_adj = clip_data(tgt_adj, ref);
        elseif(strcmp(args.clip, 'roi')), tgt_adj = clip_data(tgt_adj, ref(args.roi_ref));
        end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the CDF and bin centers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cdf, centers] = calculate_cdf(V, bins, limit)
centers = linspace(min(V), max(V), bins);
dx = diff(centers(1:2));
edges = linspace(centers(1)-dx/2, centers(end)+dx/2, bins+1);
pdf = histcounts(V, edges, 'Normalization', 'probability');
overflow = sum(pdf(pdf>limit)-limit);
pdf(pdf>limit) = limit;
pdf = pdf+overflow/bins;
cdf = cumsum(pdf);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find the transformation from each target value to the reference value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = get_transform(cdf_tgt, cdf_ref, centers_ref)
M=zeros(1,length(cdf_tgt));
for i=1:length(M)
    ind=find(cdf_ref>=cdf_tgt(i),1,'first');
    if(isempty(ind))  % Use last point
        M(i)=centers_ref(end);
    elseif(ind==1)  % Use first point
        M(i)=centers_ref(1);
    else  % Linear interpolation
        dx_ref = diff(centers_ref(1:2));
        M(i)=(cdf_tgt(i)-cdf_ref(ind-1))/...
             (cdf_ref(ind)-cdf_ref(ind-1))*dx_ref+centers_ref(ind-1);
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find the transformation from each target value to the reference value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tgt_clip = clip_data(tgt, ref)
    range = [min(ref(:)), max(ref(:))];
    tgt_clip = max(min(tgt, range(2)), range(1));
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Produce a grid of ROI masks and the desired region size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rois, s] = grid_rois(sz, regions)
rois=false(sz(1),sz(2),regions(1),regions(2));
s = floor(sz./regions);
for i=1:regions(1)
    for j=1:regions(2)
        rois(s(1)*(i-1)+1:min(s(1)*i,sz(1)),s(2)*(j-1)+1:min(s(2)*j,sz(2)),i,j)=1;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate gridded results onto a single image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tgt_adj = interp_regions(tgt_all, s)
sz=size(tgt_all);
roi_x=floor((s(2)+1)/2:s(2):sz(2));
roi_z=floor((s(1)+1)/2:s(1):sz(1));
tgt_adj=zeros(sz(1:2),'like',tgt_all);

% Map corners
tgt_adj(1:roi_z(1),1:roi_x(1))=tgt_all(1:roi_z(1),1:roi_x(1),1,1);
tgt_adj(1:roi_z(1),roi_x(end):end)=tgt_all(1:roi_z(1),roi_x(end):end,1,end);
tgt_adj(roi_z(end):end,1:roi_x(1))=tgt_all(roi_z(end):end,1:roi_x(1),end,1);
tgt_adj(roi_z(end):end,roi_x(end):end)=tgt_all(roi_z(end):end,roi_x(end):end,end,end);

% Map edges (linear interpolation)
for c=1:length(roi_x)-1
    w=repmat(linspace(1,0,roi_x(c+1)-roi_x(c)+1),roi_z(1),1);
    tgt_adj(1:roi_z(1),roi_x(c):roi_x(c+1))=w.*tgt_all(1:roi_z(1),roi_x(c):roi_x(c+1),1,c)+...
                                        (1-w).*tgt_all(1:roi_z(1),roi_x(c):roi_x(c+1),1,c+1);   
    w=repmat(linspace(1,0,roi_x(c+1)-roi_x(c)+1),sz(1)-roi_z(end)+1,1);                          
    tgt_adj(roi_z(end):end,roi_x(c):roi_x(c+1))=w.*tgt_all(roi_z(end):end,roi_x(c):roi_x(c+1),end,c)+...
                                        (1-w).*tgt_all(roi_z(end):end,roi_x(c):roi_x(c+1),end,c+1);
end
for r=1:length(roi_z)-1
    w=repmat(linspace(1,0,roi_z(r+1)-roi_z(r)+1)',1,roi_x(1));
    tgt_adj(roi_z(r):roi_z(r+1),1:roi_x(1))=w.*tgt_all(roi_z(r):roi_z(r+1),1:roi_x(1),r,1)+...
                                        (1-w).*tgt_all(roi_z(r):roi_z(r+1),1:roi_x(1),r+1,1);   
    w=repmat(linspace(1,0,roi_z(r+1)-roi_z(r)+1)',1,sz(2)-roi_x(end)+1);
    tgt_adj(roi_z(r):roi_z(r+1),roi_x(end):end)=w.*tgt_all(roi_z(r):roi_z(r+1),roi_x(end):end,r,end)+...
                                        (1-w).*tgt_all(roi_z(r):roi_z(r+1),roi_x(end):end,r+1,end);  
end

% Map everything else (bilinear interpolation)
if(length(sz)==4 && all(sz(3:4)>1))
    wz=repmat(linspace(1,0,roi_z(2)-roi_z(1)+1)',1,roi_x(2)-roi_x(1)+1);
    wx=repmat(linspace(1,0,roi_x(2)-roi_x(1)+1),roi_z(2)-roi_z(1)+1,1);
    for c=1:length(roi_x)-1
        for r=1:length(roi_z)-1
            tgt_adj(roi_z(r):roi_z(r+1),roi_x(c):roi_x(c+1))=...
                (wx).*(wz).*tgt_all(roi_z(r):roi_z(r+1),roi_x(c):roi_x(c+1),r,c)+...
                (1-wx).*(wz).*tgt_all(roi_z(r):roi_z(r+1),roi_x(c):roi_x(c+1),r,c+1)+...
                (wx).*(1-wz).*tgt_all(roi_z(r):roi_z(r+1),roi_x(c):roi_x(c+1),r+1,c)+...
                (1-wx).*(1-wz).*tgt_all(roi_z(r):roi_z(r+1),roi_x(c):roi_x(c+1),r+1,c+1);
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get a neighborhood surrounding a point, dealing with edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi = get_neighborhood(r,c,neighborhood,sz)
nhalf=floor(neighborhood/2);
bz=[r-nhalf(1),r+nhalf(1)]; 
bz=bz-min(1,bz(1))+min(1,sz(1)-bz(2)+1); % Deal with edges
bx=[c-nhalf(2),c+nhalf(2)];
bx=bx-min(1,bx(1))+min(1,sz(2)-bx(2)+1); % Deal with edges
roi=false(sz);
roi(bz(1):bz(2),bx(1):bx(2))=1;
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function args = parse_inputs(tgt, ref, args_in)
p=inputParser;
p.KeepUnmatched=1;
addOptional(p,'method','partial')
p.parse(args_in{:});

switch p.Results.method
    case 'full'
        addOptional(p,'bins',256)
        addOptional(p,'roi',[])
        addOptional(p,'roi_ref',[])
        addOptional(p,'roi_tgt',[])
        addOptional(p,'clip','none')
        addOptional(p,'limit',1)
    case 'analytic'
        addOptional(p,'roi',[])
        addOptional(p,'roi_tgt',[])
        addOptional(p,'clip','none')
        addOptional(p,'limit',1)
    case 'adaptive'
        addOptional(p,'bins',256)
        addOptional(p,'regions',[5,5])
        addOptional(p,'clip','roi')
        addOptional(p,'limit',1)
    case 'adaptive_all'
        addOptional(p,'bins',10)
        addOptional(p,'neighborhood',[11,11])
        addOptional(p,'limit',1)
    case 'point'
    case 'partial'
        addOptional(p,'version','mu_sigma')
        addOptional(p,'roi',[])
        addOptional(p,'roi_ref',[])
        addOptional(p,'roi_tgt',[])
        addOptional(p,'clip','none')
    otherwise
        error('Unknown method %s',p.Results.method)
end
p.parse(args_in{:});
args=p.Results;

switch args.method
    case 'full'
        if(~isempty(p.Results.roi)), args.roi_ref = p.Results.roi;
        elseif(~isempty(p.Results.roi_ref)), args.roi_ref = p.Results.roi_ref;
        else, args.roi_ref = isfinite(ref);
        end
        assert(all(size(args.roi_ref)==size(ref)),'Reference ROI size must match image')
        
        if(~isempty(p.Results.roi)), args.roi_tgt = p.Results.roi;
        elseif(~isempty(p.Results.roi_tgt)), args.roi_tgt=p.Results.roi_tgt;
        else, args.roi_tgt=isfinite(tgt);
        end
        assert(all(size(args.roi_tgt)==size(tgt)),'Target ROI size must match image')
        
        assert(any(strcmp(p.Results.clip,{'none','all','roi'})),'Invalid clipping option')
    case 'analytic'
        if(~isempty(p.Results.roi)), args.roi_tgt = p.Results.roi;
        elseif(~isempty(p.Results.roi_tgt)), args.roi_tgt=p.Results.roi_tgt;
        else, args.roi_tgt=isfinite(tgt);
        end
        assert(all(size(args.roi_tgt)==size(tgt)),'Target ROI size must match image')
        
        assert(isfield(ref,'x'),'Must provide reference bin centers: ref.x')
        assert(isfield(ref,'cdf'),'Must provide reference cdf: ref.cdf')
        args.bins = length(ref.x);
        assert(any(strcmp(p.Results.clip,{'none','all'})),'Invalid clipping option')
    case 'adaptive'
        assert(any(strcmp(p.Results.clip,{'none','all','roi'})),'Invalid clipping option')
        assert(all(size(tgt)==size(ref)),'Image sizes must match')
    case 'adaptive_all'
        assert(all(size(tgt)==size(ref)),'Image sizes must match')
    case 'point'
        assert(all(size(tgt)==size(ref)),'Image sizes must match')
    case 'partial'
        if(~isempty(p.Results.roi)), args.roi_ref = p.Results.roi;
        elseif(~isempty(p.Results.roi_ref)), args.roi_ref = p.Results.roi_ref;
        else, args.roi_ref = isfinite(ref);
        end
        assert(all(size(args.roi_ref)==size(ref)),'Reference ROI size must match image')
        
        if(~isempty(p.Results.roi)), args.roi_tgt = p.Results.roi;
        elseif(~isempty(p.Results.roi_tgt)), args.roi_tgt = p.Results.roi_tgt;
        else, args.roi_tgt = isfinite(tgt);
        end
        assert(all(size(args.roi_tgt)==size(tgt)),'Target ROI size must match image')
        
        assert(any(strcmp(p.Results.clip,{'none','all','roi'})),'Invalid clipping option')
        assert(any(strcmp(p.Results.version,{'mu_sigma','L2'})),'Invalid version option')
end

end