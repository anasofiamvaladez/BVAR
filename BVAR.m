%% Financial Shocks and the Macroeconomy

clear all; close all; format compact; clc;

addpath('\\bmfiles\macros\DGIE\DAM\Proyectos Especiales\Proyectos\COPOM DAM\2020\1. Febrero\GAM\_BVAR Fs\_data')

 

%% I.- Importar datos

dataMaster = datastore('datos_ene21_j.csv');
data.TextscanFormats{1,1} = '%{uuuu/MM/dd}D';

% preview(data)

%Select variables
dataMaster.SelectedVariableNames = {'fecha','dsp500','dtcn','ymxgap','R_nom','ygap','fgap','ifb_gap','tpp_gr_nem','ff_av','y10y2','y10y2_mx','dipc'};
t = readall(dataMaster);
t.fecha = datetime(t.fecha,'InputFormat', 'uuuu-MM-dd');

sdate = datetime('2004-01-01');
fdate = datetime('2019-09-01'); 

%selecting smaple
ind = t.fecha < sdate| ...
      t.fecha > fdate; 

t(ind,:) = [];

head(t)
%% Todas las variables (Gráfica)
fecha = t.fecha;
Data = [t.ygap t.ff_av t.y10y2 t.dsp500 t.ymxgap t.ifb_gap t.R_nom t.fgap t.tpp_gr_nem t.dipc t.dtcn t.y10y2_mx];

Names_en = {'US: Output gap','US: Shadow policy rate','US: term premium','US: SP500'...
    'Producto, brecha','Inversión, brecha','R','F, brecha','Tasa de créditos nuevos (empr. grandes)','IPC','dtcn','MX: Term premium'};
%% Correlaciones dinámicas
%[fig2, corrXboot, corrXbootC] = corr_BS_CI(t.Financiamiento_b, t.CS_2009, sdate, fdate);

%% II.- Plot data

%English
do_it = 1;


if do_it 
fullscreen = get(0,'ScreenSize');
%TAMAÑO PARA POWERPOINT
 F1= figure('Position',[0 0 fullscreen(3) fullscreen(4)]); 
 set(gcf,'color','w')

 for i = 1:size(Data,2)
    subplot(3,4,i)
    plot(fecha, Data(:,i));
    title(Names_en{i});
    end


save_graph = 1;
fn = ['AllDataM1_en'];

if save_graph
img = getframe(gcf);
imwrite(img.cdata, [fn, '.png']);
end
end

 
%% 1.- Define the model in SUR representation
close all
p = 6;
constant = 1;
[Y,X,~]=SUR(Data,p,constant); % Y es Ts*n; X es TS*np+1
[Ts, n] = size(Y);   % number of observations in sample, number of variables/equations
m = size(X, 2); % number of regressors per equation (n*p+1 if there is constant)

%% 2.- Priors and starting values
%0) Define is SS priors
ss_priors =1;

%a)Compute the Minnesota priors for coeficients
    %Specify the parameters of the Minnesota prior
    lambda1 = 0.1; %controls the s.d. of prior on own lags. Lower => Lower variance of the prior (tighter)
    lambda2 = 0.5; %controls the s.d. of prior on lags of other variables. Lower => Lower variance of the prior (tighter)
    lambda3 = 1; %controls the the degree to which coefficients on lags higher than 1, are likely to be zero. Higher => Lower variance of the prior (tighter)
    lambda4 = 10^5;
    lambda0 = 0.5; %in case of SS priors
    
 %Minnesota prior withou exogenous variables:  
    var_exo = [1 2 3 4];
    [B0, H0, S0] = minnesota_prior_EV(Data, p, constant, ss_priors, lambda1, lambda2, lambda3, lambda4, var_exo);
    
%b)Define priors for parameters (S and alpha) of the VAR covariance matrix
%(inverse Wishart distribution)
alpha = n+1; %degrees of freedom

%c) Define priors for long term mean (if its the case)
mu0 = [0 3 mean(t.y10y2) mean(t.dsp500) 0 0 5.5 0 mean(t.tpp_gr_nem) mean(t.dipc) mean(t.dtcn) mean(t.y10y2_mx) ]'; %
Hmu0 = diag(diag(S0).*(lambda0*ones(1,n))');  % Variance

%d)Define starting values, iterations and burn
Sigma = eye(n);
reps = 10000;
burn = reps - 1000;

%save Model_2020 t Data Y X constant Names_en ss_priors m n Ts  B0 H0 S0 p fecha

%% 3.- Estimate Model by Gibbs Sampling
warning off

%For model with SS priors
[Beta, Omega, mu] = BayesianVAR_SSpriors(Y, X, p, B0, H0, S0, mu0, Hmu0, alpha, Sigma, reps, burn);


%Plot mean english
do_it = 1;

if do_it 
fullscreen = get(0,'ScreenSize');
%TAMAÑO PARA POWERPOINT
 F1= figure('Position',[0 0 fullscreen(3)*4 fullscreen(4)]); 
 set(gcf,'color','w')    
%A) mean
 for i = 1:size(Data,2)
    subplot(3,4,i)
    histogram(mu(i,:), 50);
    title(Names_en{i});
 end

save_graph = 1;
fn = ['SSM1_EBPen'];

if save_graph
img = getframe(gcf);
imwrite(img.cdata, [fn, '.png']);
end
end

save Gibbs_M1 Beta Omega mu 
  
%% 4.- Identification and IRF by Cholesky decomposition
close all
hor_irf = 40;
choque = -1;
clear IRF_chol IRF_plot_chol

[IRF_chol, IRF_plot_chol] = Bayesian_IRF(Beta, Omega, hor_irf, constant, ss_priors, choque, p, n, m);

%save IRFchol_M1 IRF_chol IRF_plot_chol hor_irf

%% Gráficas IRF chol: EBP SHOCK
%English
close all
clear Yfi YFf
YFf = IRF_plot_chol;  %%CAMBIAR AL MÉTODO DE IDENTIFICACION QUE QUEREMOS EXPLORAR
nameshock = {'xxx'};
fullscreen = get(0,'ScreenSize');
x = linspace(0,hor_irf-1,hor_irf); 

I = []; %variables for cumsum 

%TAMAÑO PARA POWERPOINT
 F1= figure('Position',[0 0 fullscreen(3)*4 fullscreen(4)]); 
 set(gcf,'color','w');
%TAMAÑO PARA PANTALLA 
% F1= figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
% set(gcf,'color','w');

j=1;

%Variables a graficar
pos_graph = [8 5 6];

%Posición del choque
r = 8;


f_f   = prctile(squeeze(YFf(8,:,8,:))',[50])';


for i = pos_graph
                
    
    %Yfi     = squeeze(YFf(i,:,r,:)); %muestra la distribucion de la IRF de la variable i al choque r    

    Yfi     = squeeze(YFf(i,:,r,:))/f_f(1)*(-1);
    % acumula IRF para variables seleccionadas
        if ismember(i,I); Yfi = cumsum(Yfi); end
    %----------
    
        pYfi    = percentile(Yfi',[16 84])'; 
        medYfi  = percentile(Yfi',50)';

        subplot(3,1,j)    
                
        fill( [x fliplr(x)],  [pYfi(:,1)' fliplr(pYfi(:,2)')],rgb('Lavender'), 'edgecolor', 'none'); % shaded area between 16-84 percentiles       
        hold on
        
        plot(x,medYfi,'color',rgb('Purple'),'LineWidth', 1.5); % mediana
        hold on
        
        plot(x,zeros(1,hor_irf),'.','color','k','LineWidth', .2) % eje x 
        %alpha(.8);


        xlim([0 40]);        
        set(gca,'XTick',0:10:40)
    
     % límites en las gráficas
        if  j == 1 
            ylim([-1 1]); 
            set(gca,'YTick',-1:1:1  )
           
        elseif  j == 2 
           ylim([-.2 .2]); 
           set(gca,'YTick',-.2:.1:.2)
  
        elseif j == 3
           ylim([-0.3 0.3]); 
           set(gca,'YTick',-0.3:.1:0.3)
        
                          
        else   
           ylim([-0.4 0.4]); 
           set(gca,'YTick',-0.4:.1:0.4)
        
            end               
        hold off
        
         set(gca,'fontsize',16,'fontname','calibri');
         t = title(Names_en{i});
         set(t,'fontname','Calibri','fontsize',16','fontweight','normal');%,'interpreter','latex')
         ylabel({'%'},'fontsize',16,'fontname','calibri','FontAngle','italic')
         xlabel({'meses'},'fontsize',16,'fontname','calibri','FontAngle','italic')

         if j == 1; ylabel({''},'fontsize',14,'fontname','calibri','FontAngle','italic'); 
         end

        
         
         j=j+1;

end

save_graph = 1;
fn = ['IRFcholM1_EBPen'];

if save_graph
img = getframe(gcf);
imwrite(img.cdata, [fn, '.png']);
end




