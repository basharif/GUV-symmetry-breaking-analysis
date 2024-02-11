    clc; clearvars; close all;
rng('default')
% see data from figure 2A and 2B 
%https://elifesciences.org/articles/01433#s4
nframes0 = 8;
nframes1 = 80;
npoints = 360; % for each degree
nallfr = nframes0 + nframes1;

tsize = 14;

truememb0 = 1*ones(npoints,nframes0);
truememb1 = 1*ones(npoints,nframes1);

% before input, nts time points
memb0 = truememb0+ 0.1*randn(npoints,nframes0);
% after input, nts time points
memb1 = truememb1+ 0.1*randn(npoints,nframes1);


% before input, nts time points
simActA0 = (truememb0 + 0.1*randn(npoints,nframes0));%./memb0;
% after input, nts time points
simActA1 = (2 + truememb1 + 0.1*randn(npoints,nframes1));%./memb1;

x = [0:npoints-1]';
y = 150*normpdf(x,npoints/2,40)+1;

% before input, nts time points
simactin0 = (truememb0 + 0.1*randn(npoints,nframes0));%./memb0;
%take ratio to acta which is like membrane marker
trueactin1 = repmat(y,1,nframes1);
simactin1 = (trueactin1 + 0.1*randn(npoints,nframes1));%./memb1;

simmembfull = [memb0,memb1];
simActAfull = [simActA0,simActA1];
simactinfull = [simactin0,simactin1];

figure;
imagesc(simmembfull);
colormap("parula");colorbar;
title('Membrane kymograph','FontSize',tsize);
xline(nframes0,'LineWidth',2)
f1 = gcf; 
f1.Position(3:4) = [500,300];


figure;
imagesc(simActAfull);
colormap("parula");colorbar;
xline(nframes0,'LineWidth',2)
title('ActA kymograph','FontSize',tsize);
f1 = gcf; 
f1.Position(3:4) = [500,300];

figure;
imagesc(simactinfull);
colormap("parula");colorbar;
xline(nframes0,'LineWidth',2)
title('actin kymograph','FontSize',tsize);
f1 = gcf; 
f1.Position(3:4) = [500,300];


figure;
plot(x,memb0(:,1),'b'); hold on
plot(x,memb1(:,1),'r')
legend('pre-Rapa','post-Rapa','AutoUpdate','off');
plot(x,memb0(:,2:end),'b'); 
plot(x,memb1(:,2:end),'r')
ylim([0 5])
hold off
title('membrane aligned','FontSize',tsize);

figure;
plot(x,simActA0(:,1),'b'); hold on
plot(x,simActA1(:,1),'r')
legend('pre-Rapa','post-Rapa','AutoUpdate','off');
plot(x,simActA0(:,2:end),'b'); 
plot(x,simActA1(:,2:end),'r')
ylim([0 5])
hold off
title('ActA aligned','FontSize',tsize);

figure;
plot(x,simactin0(:,1),'b'); hold on
plot(x,simactin1(:,1),'r')
legend('pre-Rapa','post-Rapa','AutoUpdate','off');
plot(x,simactin0(:,2:end),'b'); 
plot(x,simactin1(:,2:end),'r')
ylim([0 5])
hold off
title('actin aligned','FontSize',tsize);

figure;
binscatter(simActA0(:),simactin0(:),250);
aacorr0 = corr(simActA0(:),simactin0(:));
text(3,4,strcat('correlation: ',num2str(aacorr0)))
colormap(gca,'winter')
axis([0 5 0 5])
title('before input, ActA-actin correlations','FontSize',tsize);

figure;
binscatter(simActA1(:),simactin1(:),250);
aacorr1 = corr(simActA1(:),simactin1(:));
text(3,4,strcat('correlation: ',num2str(aacorr1)))
colormap(gca,'spring')
axis([0 5 0 5])
title('after input, ActA-actin correlations','FontSize',tsize);

aacorrfull = nan(nallfr,1);


for i = 1:nallfr

    aacorrfull(i) = corr(simActAfull(:,i),simactinfull(:,i));

end
figure;
plot(aacorrfull,'LineWidth',2); hold on
ylim([-1 1])
xline(nframes0,'LineWidth',2)
yline(0)
hold off
title('single frame ActA-actin correlations','FontSize',2)
xlabel('frames','FontSize',tsize)
ylabel('pearson corr coeff.','FontSize',tsize)

fh = findall(0,'Type','Figure');
txt_obj = findall(fh,'Type','text');
set(txt_obj,'FontSize',tsize);