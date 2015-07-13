

% HydrolightData.mat has the following elements
% Rrs:  LUT in remote-sensing reflectance (Rrs) units [1/sr]. Each column
%       is a different Hydrolight runs and the row represents wavelength 
%       (120 wavelengths).
%
% c:    Cell that contains the Hydrolight inputs used to generate Rrs. The
%       different cell are:
%       c{1}:   IOP inputs. Either from LONG pond (input140929ONTOS) or 
%                 the ONTario lake (input140929LONGS)
%       c{2}:   Chlorophyll-a concentration
%   	c{3}:   Total Suspended Solids (Sediments)
%       c{4}:   Absorption coefficient at 440nm for CDOM
%      	c{5}:   backscatter fraction for the discrete phase function (dpf)
%
% wavelength:   wavelength in um (120 wavelengths).

% load the spectra text file 
load HydrolightData.mat;

%% plot all curves
figure
fs = 15;
set(gcf,'color','white')
set(gca,'fontsize',fs)
plot(wavelength,Rrs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('R_{rs} [1/sr]','fontsize',fs)
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Isolate Suspended Solids (Sediments) spectra

% specify fixed values of other 4 parameters
iop = 'input140929LONGS';
chlA = 0.01;
cdom = 0.0954;
backscatter = 'FFbb005.dpf';

% find indices of spectral curves based on above parameters
indices_tss = find((strcmp(c{1},iop))&(c{2}==chlA)...
                    &(c{4}==cdom)&(strcmp(c{5},backscatter)));

% create empty array to hold spectra values
tss_spectra_all = zeros(size(wavelength,1),numel(indices_tss));

% plot subset of spectral curves 
figure; hold on
for i=1:numel(indices_tss)
    all_tss = plot(wavelength,Rrs(:,indices_tss(i)));
    tss_spectra_all(:,i) = Rrs(:,indices_tss(i));
end
fs = 15; set(gcf,'color','white'); set(gca,'fontsize',fs);
xlabel('wavelength [\mum]','fontsize',fs);
ylabel('R_{rs} [1/sr]','fontsize',fs);
if iop=='input140929LONGS'; water_body = ' Long Pond';
    else water_body = ' Ontario Lake';
end
title({'Hydrolight-Derived Total Suspended Sediment Spectral Curves';
    strcat('for ',water_body)})


% compute average spectral curve
tss_spectra_avg = mean(tss_spectra_all,2);
avg_tss = plot(wavelength, tss_spectra_avg,'r');
legend([all_tss, avg_tss],{'All','Average'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Isolate Chlorophyll-A spectra

% specify fixed values of other 4 parameters
iop = 'input140929LONGS'; %iop = 'input140929ONTOS';
tss= 0.01;
cdom = 0.0954;
backscatter = 'FFbb005.dpf';

% find indices of spectral curves based on above parameters
indices_chlA = find((strcmp(c{1},iop))&(c{3}==tss)...
                     &(c{4}==cdom)&(strcmp(c{5},backscatter)));

% create empty array to hold spectra values
chlA_spectra_all = zeros(size(wavelength,1),numel(indices_chlA));

% plot subset of spectral curves 
figure; hold on
for i=1:numel(indices_chlA)
    all_chlA = plot(wavelength,Rrs(:,indices_chlA(i)));
    chlA_spectra_all(:,i) = Rrs(:,indices_chlA(i));
end
fs = 15; set(gcf,'color','white'); set(gca,'fontsize',fs);
xlabel('wavelength [\mum]','fontsize',fs);
ylabel('R_{rs} [1/sr]','fontsize',fs);
if iop=='input140929LONGS'; water_body = ' Long Pond';
    else water_body = ' Ontario Lake';
end
title({'Hydrolight-Derived Chlorophyll-A Spectral Curves';
    strcat('for ',water_body)})


% compute average spectral curve
chlA_spectra_avg = mean(chlA_spectra_all,2);
avg_chlA = plot(wavelength, chlA_spectra_avg,'r');
legend([all_chlA, avg_chlA],{'All','Average'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Isolate CDOM spectra

% specify fixed values of other 4 parameters
iop = 'input140929LONGS'; %iop = 'input140929ONTOS';
tss= 0.01;
chlA = 0.01;
backscatter = 'FFbb005.dpf';

% find indices of spectral curves based on above parameters
indices_cdom = find((strcmp(c{1},iop))&(c{2}==chlA)...
                     &(c{3}==tss)&(strcmp(c{5},backscatter)));

% create empty array to hold spectra values
cdom_spectra_all = zeros(size(wavelength,1),numel(indices_cdom));

% plot subset of spectral curves 
figure; hold on
for i=1:numel(indices_cdom)
    all_cdom = plot(wavelength,Rrs(:,indices_cdom(i)));
    cdom_spectra_all(:,i) = Rrs(:,indices_cdom(i));
end
fs = 15; set(gcf,'color','white'); set(gca,'fontsize',fs);
xlabel('wavelength [\mum]','fontsize',fs);
ylabel('R_{rs} [1/sr]','fontsize',fs);
if iop=='input140929LONGS'; water_body = ' Long Pond';
    else water_body = ' Ontario Lake';
end
title({'Hydrolight-Derived CDOM Spectral Curves';
    strcat('for ',water_body)})


% compute average spectral curve
cdom_spectra_avg = mean(cdom_spectra_all,2);
avg_cdom = plot(wavelength, cdom_spectra_avg,'r');
legend([all_cdom, avg_cdom],{'All','Average'});



%%  write average spectra to file
% 
% fid_tss = fopen('tss_spectra_avg.txt', 'w');
% fprintf(fid_tss,'%6d\n', tss_spectra_avg)
% fclose(fid_tss)
% 
% 
% fid_chlA = fopen('chlA_spectra_avg.txt', 'w');
% fprintf(fid_chlA,'%6d\n', chlA_spectra_avg)
% fclose(fid_chlA)
% 
% 
% fid_cdom = fopen('cdom_spectra_avg.txt', 'w');
% fprintf(fid_cdom,'%6d\n', cdom_spectra_avg)
% fclose(fid_cdom)
% 
% fid_wavelength = fopen('wavelength.txt', 'w');
% fprintf(fid_wavelength,'%4f\n', wavelength)
% fclose(fid_wavelength)

%% write out curve series to file

mkdir('tss_spectra')

for r = 1:size(tss_spectra_all,2)
    fid_tss_tmp = fopen(['tss_spectra/tss_spectra_',num2str(r),'.txt'], 'w');
    fprintf(fid_tss_tmp,'%6d\n', tss_spectra_all(:,r))
    fclose(fid_tss_tmp)
end     


mkdir('cdom_spectra')

for r = 1:size(cdom_spectra_all,2)
    fid_cdom_tmp = fopen(['cdom_spectra/cdom_spectra_',num2str(r),'.txt'], 'w');
    fprintf(fid_cdom_tmp,'%6d\n', cdom_spectra_all(:,r))
    fclose(fid_cdom_tmp)
end     

mkdir('chlA_spectra')

for r = 1:size(chlA_spectra_all,2)
    fid_chlA_tmp = fopen(['chlA_spectra/chlA_spectra_',num2str(r),'.txt'], 'w');
    fprintf(fid_chlA_tmp,'%6d\n', chlA_spectra_all(:,r))
    fclose(fid_chlA_tmp)
end     

