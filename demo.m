%% Load data
load demo/sample_data
labels={'Original','Reference','Partial','Full','Point','Adaptive'};
colors=lines(4);
colors=[0 1 1;1 0 1;colors];

%% Demo 1 - Arbitrary transform
figure(1);clf
set(gcf,'Name','Power law detector')

clim=[-60 0];

subplot(2,4,1)
orig=db(env.^2);
imagesc(x*1e3,z*1e3,orig,clim);axis image;colormap gray
title('Original (squared)')
xlabel('Lateral (mm)')
ylabel('Axial (mm)')

subplot(2,4,2)
ref=db(env);
imagesc(x*1e3,z*1e3,ref,clim);axis image;colormap gray
title('Reference')

subplot(2,4,5)
orig_adj1=histmatch(orig,ref,'method','partial');
imagesc(x*1e3,z*1e3,orig_adj1,clim);axis image;colormap gray
title('Partial matched')

subplot(2,4,6)
orig_adj2=histmatch(orig,ref,'method','full');
imagesc(x*1e3,z*1e3,orig_adj2,clim);axis image;colormap gray
title('Full matched')

subplot(2,4,7)
orig_adj3=histmatch(orig,ref,'method','point');
imagesc(x*1e3,z*1e3,orig_adj3,clim);axis image;colormap gray
title('Point matched')

subplot(2,4,8)
orig_adj4=histmatch(orig,ref,'method','adaptive');
imagesc(x*1e3,z*1e3,orig_adj4,clim);axis image;colormap gray
title('Adaptive matched')

subplot(2,2,2)
h=[];
data={orig,ref,orig_adj1,orig_adj2,orig_adj3,orig_adj4};
for i=1:2
    h(end+1)=histogram(data{i},'Normalization','probability','EdgeColor','none','FaceAlpha',.5,'FaceColor',colors(i,:));
    hold all
end
for i=3:length(data)
    h(end+1)=histogram(data{i},'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colors(i,:));
end
hold off
ylabel('Probability')
xlabel('Pixel value')
legend(labels,'Location','NorthWest')

%% Demo 2 - Compounding
figure(2);clf
set(gcf,'Name','Spatial compounding')

clim=[-60 0];

% Identify a speckle region to match
[x1,z1]=meshgrid(x,z);

subplot(2,4,1)
orig=db(env_comp);
imagesc(x*1e3,z*1e3,orig,clim);axis image;colormap gray
title('Original (compounded)')
xlabel('Lateral (mm)')
ylabel('Axial (mm)')

subplot(2,4,2)
ref=db(env);
imagesc(x*1e3,z*1e3,ref,clim);axis image;colormap gray
title('Reference')

subplot(2,4,5)
orig_adj1=histmatch(orig,ref,'method','partial');
imagesc(x*1e3,z*1e3,orig_adj1,clim);axis image;colormap gray
title('Partial matched')

subplot(2,4,6)
orig_adj2=histmatch(orig,ref,'method','full');
imagesc(x*1e3,z*1e3,orig_adj2,clim);axis image;colormap gray
title('Full matched')

subplot(2,4,7)
orig_adj3=histmatch(orig,ref,'method','point');
imagesc(x*1e3,z*1e3,orig_adj3,clim);axis image;colormap gray
title('Point matched')

subplot(2,4,8)
orig_adj4=histmatch(orig,ref,'method','adaptive');
imagesc(x*1e3,z*1e3,orig_adj4,clim);axis image;colormap gray
title('Adaptive matched')

subplot(2,2,2)
h=[];
data={orig,ref,orig_adj1,orig_adj2,orig_adj3,orig_adj4};
for i=1:2
    h(end+1)=histogram(data{i},'Normalization','probability','EdgeColor','none','FaceAlpha',.5,'FaceColor',colors(i,:));
    hold all
end
for i=3:length(data)
    h(end+1)=histogram(data{i},'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colors(i,:));
end
hold off
ylabel('Probability')
xlabel('Pixel value')
legend(labels,'Location','NorthWest')

%% Demo 3 - SLSC to B-mode
figure(3);clf
set(gcf,'Name','SLSC to B-mode')

clim=[-60 0];

subplot(2,4,1)
orig=slsc;
imagesc(x*1e3,z*1e3,orig,[0 1]);axis image;colormap gray
title('Original (SLSC)')
xlabel('Lateral (mm)')
ylabel('Axial (mm)')

subplot(2,4,2)
ref=db(env);
imagesc(x*1e3,z*1e3,ref,clim);axis image;colormap gray
title('Reference (B-mode)')

subplot(2,4,5)
orig_adj1=histmatch(orig,ref,'method','partial');
imagesc(x*1e3,z*1e3,orig_adj1,clim);axis image;colormap gray
title('Partial matched')

subplot(2,4,6)
orig_adj2=histmatch(orig,ref,'method','full');
imagesc(x*1e3,z*1e3,orig_adj2,clim);axis image;colormap gray
title('Full matched')

subplot(2,4,7)
orig_adj3=histmatch(orig,ref,'method','point');
imagesc(x*1e3,z*1e3,orig_adj3,clim);axis image;colormap gray
title('Point matched')

subplot(2,4,8)
orig_adj4=histmatch(orig,ref,'method','adaptive');
imagesc(x*1e3,z*1e3,orig_adj4,clim);axis image;colormap gray
title('Adaptive matched')

subplot(2,2,2)
h=[];
data={orig,ref,orig_adj1,orig_adj2,orig_adj3,orig_adj4};
for i=1:2
    h(end+1)=histogram(data{i},'Normalization','probability','EdgeColor','none','FaceAlpha',.5,'FaceColor',colors(i,:));
    hold all
end
for i=3:length(data)
    h(end+1)=histogram(data{i},'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colors(i,:));
end
hold off
ylabel('Probability')
xlabel('Pixel value')
legend(labels,'Location','NorthWest')

%% Demo 4 - B-mode to SLSC
figure(4);clf
set(gcf,'Name','B-mode to SLSC')

clim=[0 1];

subplot(2,4,1)
orig=db(env);
imagesc(x*1e3,z*1e3,orig,[-60 0]);axis image;colormap gray
title('Original (B-mode)')
xlabel('Lateral (mm)')
ylabel('Axial (mm)')

subplot(2,4,2)
ref=slsc;
imagesc(x*1e3,z*1e3,ref,clim);axis image;colormap gray
title('Reference (SLSC)')

subplot(2,4,5)
orig_adj1=histmatch(orig,ref,'method','partial');
imagesc(x*1e3,z*1e3,orig_adj1,clim);axis image;colormap gray
title('Partial matched')

subplot(2,4,6)
orig_adj2=histmatch(orig,ref,'method','full');
imagesc(x*1e3,z*1e3,orig_adj2,clim);axis image;colormap gray
title('Full matched')

subplot(2,4,7)
orig_adj3=histmatch(orig,ref,'method','point');
imagesc(x*1e3,z*1e3,orig_adj3,clim);axis image;colormap gray
title('Point matched')

subplot(2,4,8)
orig_adj4=histmatch(orig,ref,'method','adaptive');
imagesc(x*1e3,z*1e3,orig_adj4,clim);axis image;colormap gray
title('Adaptive matched')

subplot(2,2,2)
h=[];
data={orig,ref,orig_adj1,orig_adj2,orig_adj3,orig_adj4};
for i=1:2
    h(end+1)=histogram(data{i},'Normalization','probability','EdgeColor','none','FaceAlpha',.5,'FaceColor',colors(i,:));
    hold all
end
for i=3:length(data)
    h(end+1)=histogram(data{i},'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colors(i,:));
end
hold off
ylabel('Probability')
xlabel('Pixel value')
legend(labels,'Location','NorthWest')

%% Demo 5 - analytic CDF
figure(5);clf
set(gcf,'Name','Analytic CDF')

clim=[-60 0];

% Identify a speckle region to match
[x1,z1]=meshgrid(x,z);
zi=30e-3;
dz=10e-3;
roi=abs(z1-zi)<dz;

subplot(241)
ref=db(env);
imagesc(x*1e3,z*1e3,ref,clim);axis image;colormap gray
hold all
plot(1e3*[x(1) x(end) x(end) x(1) x(1)]*.95,1e3*(zi+[dz dz -dz -dz dz]),'w--','LineWidth',2)
hold off
title('Reference')
xlabel('Lateral (mm)')
ylabel('Axial (mm)')

subplot(4,4,9)
histogram(ref(roi),'EdgeColor','none','Normalization','probability')
ylabel('Probability')
xlabel('Pixel value')

cdf_param=.3;
cdf=[];
cdf.x=linspace(-60,0,256);
cdf.pdf=raylpdf(10.^(cdf.x/20),cdf_param);
cdf.cdf=cumsum(cdf.pdf)/sum(cdf.pdf);

subplot(242)
orig_adj1=histmatch(ref,cdf,'method','analytic','roi',roi);
imagesc(x*1e3,z*1e3,orig_adj1,clim);axis image;colormap gray
title(sprintf('Rayleigh: %.1f',cdf_param))

subplot(4,4,10)
h=histogram(orig_adj1(roi),'EdgeColor','none','Normalization','probability');
hold all
plot(cdf.x,cdf.pdf/sum(cdf.pdf)*h.BinWidth/diff(cdf.x(1:2)),'r-','LineWidth',2)
hold off

cdf_param=[.4,.15];
cdf.pdf=pdf('norm',10.^(cdf.x/20),cdf_param(1),cdf_param(2));
cdf.cdf=cumsum(cdf.pdf)/sum(cdf.pdf);

subplot(243)
orig_adj2=histmatch(ref,cdf,'method','analytic','roi',roi);
imagesc(x*1e3,z*1e3,orig_adj2,clim);axis image;colormap gray
title(sprintf('Normal: %.1f, %.2f', cdf_param(1), cdf_param(2)))

subplot(4,4,11)
h=histogram(orig_adj2(roi),'EdgeColor','none','Normalization','probability');
hold all
plot(cdf.x,cdf.pdf/sum(cdf.pdf)*h.BinWidth/diff(cdf.x(1:2)),'r-','LineWidth',2)
hold off

cdf_param=[-2,1];
cdf.pdf=pdf('logn',10.^(cdf.x/20),cdf_param(1),cdf_param(2));
cdf.cdf=cumsum(cdf.pdf)/sum(cdf.pdf);

subplot(244)
orig_adj3=histmatch(ref,cdf,'method','analytic','roi',roi);
imagesc(x*1e3,z*1e3,orig_adj3,clim);axis image;colormap gray
title(sprintf('Lognormal: %d, %d',cdf_param(1), cdf_param(2)))

subplot(4,4,12)
h=histogram(orig_adj3(roi),'EdgeColor','none','Normalization','probability');
hold all
plot(cdf.x,cdf.pdf/sum(cdf.pdf)*h.BinWidth/diff(cdf.x(1:2)),'r-','LineWidth',2)
hold off