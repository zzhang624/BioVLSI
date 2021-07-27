clear
clc
close all

%change namelist, Calc(,dcd), number&for_loop

%Constants
global k T dt eps_inf q nframe
k=1.38064852E-23;
eps_inf=8.854187E-12;
T=300;
q=1.602176634E-19;
%Dalton_factor=1/3.336E-30;
dt=1.0E-12;
nframe=10000;

%namelist={'TIP3P1.5X', 'Eps1.5X',  'TIP3P2X', 'Eps2X'};
%namelist={'Eps3X'};
%for name_num=1:length(namelist)
%    name=namelist{name_num};
%    Calc(name,name);
%end

namelist={'TIP3P-50mg:mL-ALa',  'TIP3P-100mg:mL-ALa', 'TIP3P-150mg:mL-ALa',  'Eps-50mg:mL-Ala', 'Eps-100mg:mL-Ala','Eps-150mg:mL-Ala'}; %two rep names
for name_num=1:length(namelist)
    name=namelist{name_num};
    rep_name={'-rep1', '-rep2'};
    for rep=1:length(rep_name)
    Calc(name,strcat(name,rep_name{rep}));
    end
end


function Calc(psf_name,dcd_name)

global k T eps_inf q 

%psf=strcat('dcd/',psf_name,'.psf');
%dcd=strcat('dcd/',dcd_name,'.dcd');
%dt=20E-15;
%dcd='output/TIP3P3XNoAla/TIP3P3XNoAla_07.dcd';

%Reading in files
B=readpsf(strcat('dcd/',psf_name,'.psf'));

M_water = [];
M_nonwater = [];
n_dcd = {'-1','-2','-3','-4'};
%timetoshow = {'-10ns', '-20ns','-30ns', '-40ns'};
%i_timetoshow = 0;

for number=1:length(n_dcd)
    
%C=readdcd_large(strcat('dcd/',dcd_name, n_dcd{number}, '.dcd'),B.atom_id,seg/10*nframe+1,seg/10*nframe+1000);
%C=readdcd_large(strcat('dcd/',dcd_name, '.dcd'),B.atom_id,seg/10*nframe+1,seg/10*nframe+1000);
[C, box]=readdcd(strcat('dcd/',dcd_name, n_dcd{number}, '.dcd'),B.atom_id);
water_index=(all(B.residue_name=='TIP3'|B.residue_name=='TIP4',2));
nonwater_index=(water_index==0);


%Unit conversions
charges=B.charge*q;
C.x=C.x*1E-10;
C.y=C.y*1E-10;
C.z=C.z*1E-10;
box=box(1,:)*1E-10;
V=box(1)*box(2)*box(3);
%x=C(:,1:3:end);
%y=C(:,2:3:end);
%z=C(:,3:3:end);
%clear C;



[M_water_seg,~]=sys_dipole(charges(water_index),C.x(:,water_index),C.y(:,water_index),C.z(:,water_index),V);
[M_nonwater_seg,~]=sys_dipole(charges(nonwater_index),C.x(:,nonwater_index),C.y(:,nonwater_index),C.z(:,nonwater_index),V);
clear C;

M_water = cat(1,M_water,M_water_seg);
M_nonwater = cat(1,M_nonwater,M_nonwater_seg);

end 

Phi_WW=correlation(M_water,M_water);
Phi_XW=correlation(M_nonwater,M_water);
Phi_XX=correlation(M_nonwater,M_nonwater);


[~,WW_params]=fit_LM(Phi_WW,1);
[~,XW_params]=fit_LM(Phi_XW,1);
[~,XX_params]=fit_LM(Phi_XX,1);


freq_range=[1E6 1E13];
Num_freq=5000;
omega=2*pi*(10.^(linspace(log10(freq_range(1)),log10(freq_range(2)),Num_freq)))';

Chi_Fit_WW=eps_func(omega,WW_params)/(3*k*T*V);
Chi_Fit_XW=eps_func(omega,XW_params)/(3*k*T*V);
Chi_Fit_XX=eps_func(omega,XX_params)/(3*k*T*V);
Chi_Fit=Chi_Fit_WW+2*Chi_Fit_XW+Chi_Fit_XX;


eps_Fit_WW=1+(4*pi*Chi_Fit_WW)/eps_inf;
eps_Fit_XW=1+(4*pi*Chi_Fit_XW)/eps_inf;
eps_Fit_XX=1+(4*pi*Chi_Fit_XX)/eps_inf;
eps_Fit=1+(4*pi*Chi_Fit)/eps_inf;

save(strcat('eps/', dcd_name,'.mat'), 'eps_Fit', 'eps_Fit_WW', 'eps_Fit_XW' , 'eps_Fit_XX', 'omega', '-v7.3')
%save(strcat('eps/', dcd_name, timetoshow{i_timetoshow} ,'.mat'), 'eps_Fit_WW', 'omega', '-v7.3')
end
%delete(strcat('dcd/',dcd_name, n_dcd{number}, '.dcd'));

%save(strcat('dipole/', dcd_name ,'_100ns.mat'), 'M_water','-v7.3')


 


function[Phi]=correlation(M_i,M_j)
Phi=zeros(length(M_i),1);
%M_i_mean=mean(M_i,1);
%M_j_mean=mean(M_j,1);

for tt=1:length(M_i)
    tau_max=(length(M_i)-tt);
    %tt
    for tau=0:tau_max
        
        %M_tilda_i=M_i(1+tau,:);
        %M_tilda_j=M_j(tt+tau,:);
        %Phi(tt)=Phi(tt)+(M_tilda_i*M_tilda_j.');
        
        M_tilda_i=M_i(1+tau,:);
        M_tilda_j=M_j(tt+tau,:);
        Phi(tt)=Phi(tt)+(M_tilda_i*M_tilda_j.');
        
    end


    Phi(tt)=Phi(tt)/(tau_max + 1);
    %Phi(tt)=Phi(tt)/length(M_i);
    %Phi(tt)=Phi(tt)/denominator;
end
end

function [M,Chi_0]=sys_dipole(charges,x,y,z,V)
global k T
M=zeros(size(x,1),3);

avg_charge=sum(charges)/length(charges);
for ii=1:size(x,1)
    M(ii,1)=double((x(ii,:)*(charges-avg_charge)));
    M(ii,2)=double((y(ii,:)*(charges-avg_charge)));
    M(ii,3)=double((z(ii,:)*(charges-avg_charge)));
end
%Normalize!!!!!
M=M/sqrt(4*pi*1);
%M_mean=mean(M,1);
%M_tilda_sq=sum(mean((M-M_mean).^2));
M1=mean(sum((M.^2),2));
M2=sum(mean(M,1).^2);
Chi_0=1/(3*k*T*V)*(M1-M2);
end

function[Fit,x_out]=fit_LM(Phi,NumExp)
global dt
Normalization=Phi(1);
x_guess=[1;1E-3];
x_out=LevMarq(Phi/Normalization,dt,NumExp,x_guess,true);
x_out(1:2:end)=x_out(1:2:end)*Normalization;

Fit=zeros(length(Phi),1);
for ii=1:NumExp
    Fit=Fit+x_out(2*ii-1)*exp(-x_out(2*ii)*dt*[0:(length(Phi)-1)]).';
end

end

function[Chi_Fit]=eps_func(omega,x_out)
NumExp=length(x_out)/2;
Chi_prime_Fit=zeros(length(omega),1);
Chi_double_prime_Fit=zeros(length(omega),1);

for ii=1:NumExp
    Amp=x_out(2*ii-1);
    tau_exp=1/x_out(2*ii);
    Chi_prime_Fit=Chi_prime_Fit+omega.*((Amp*omega*tau_exp.^2)./(1+(tau_exp*omega).^2));
    Chi_double_prime_Fit=Chi_double_prime_Fit+omega.*((Amp*tau_exp)./(1+(tau_exp*omega).^2));
end
Chi_prime_Fit=max(Chi_prime_Fit)-Chi_prime_Fit;
Chi_Fit=(Chi_prime_Fit+1i.*Chi_double_prime_Fit);
end
