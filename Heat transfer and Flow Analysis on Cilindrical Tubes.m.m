clear all; close all; clc;
format long
% Escoamento sobre Banco de Tubos cilíndricos
% Cálculo para validação da simulação

% Informações da geometria e condições do escoamento
re = 50;
d = 0.55;
n = 20;
nt = 5;
st = 2.5;
l = 1;
tin = 140;
ts = 20;

% Propriedades do ar a T = 120ºC
rho = 0.8977;
cp = 1011;
k = 0.03235;
ni = 2.522*10^-5;
pr = 0.7073;
%Pr para temperatura de superfície
prs = 0.7309;
%Fator de correção do número de Nusselt para Nl = 5
f = 0.93;

%Cálculo das incógnitas desejadas 
vmax = (re * ni)/d
v = vmax * ((st-d)/st)
m = rho * v * nt * st* l
as = n * pi * d * l
nu = f * 0.9 * re^0.4 * pr^0.36 * (pr/prs)^0.25
h = ( nu * k)/ d
ntu = (h * as)/(m * cp)
te = ts - (ts -tin)*exp(-ntu)
Q = m * cp * (tin - te)

%%%%%%% Codigos de refino %%%%%%
clc,clear all,clear global,format short,close all



 [pde_fig,ax]=pdeinit;
 pdetool('appl_cb',1);
 set(ax,'DataAspectRatio',[1 0.77929197994987476 1]);
 set(ax,'PlotBoxAspectRatio',[35 7.69929648241206 1]);
 set(ax,'XLimMode','auto');
 set(ax,'YLimMode','auto');
 set(ax,'XTickMode','auto');
 set(ax,'YTickMode','auto');
 
%  Geometry description:
%  retangulo

%  circulo
esp_y=1.25;
a=.55;
pdeellip(esp_y*0,-esp_y*2,a,a,0,'E1');
pdeellip(esp_y*2,-esp_y*2,a,a,0,'E2');
pdeellip(esp_y*4,-esp_y*2,a,a,0,'E3');
pdeellip(esp_y*6,-esp_y*2,a,a,0,'E4');
pdeellip(esp_y*8,-esp_y*2,a,a,0,'E5');

% 
pdeellip(esp_y*0,0,a,a,0,'E6');
pdeellip(esp_y*2,0,a,a,0,'E7');
pdeellip(esp_y*4,0,a,a,0,'E8');
pdeellip(esp_y*6,0,a,a,0,'E9');
pdeellip(esp_y*8,0,a,a,0,'E10');


pdeellip(esp_y*0,+esp_y*2,a,a,0,'E11');
pdeellip(esp_y*2,+esp_y*2,a,a,0,'E12');
pdeellip(esp_y*4,+esp_y*2,a,a,0,'E13');
pdeellip(esp_y*6,+esp_y*2,a,a,0,'E14');
pdeellip(esp_y*8,+esp_y*2,a,a,0,'E15');

pdeellip(esp_y*0,+esp_y*4,a,a,0,'E16');
pdeellip(esp_y*2,+esp_y*4,a,a,0,'E17');
pdeellip(esp_y*4,+esp_y*4,a,a,0,'E18');
pdeellip(esp_y*6,+esp_y*4,a,a,0,'E19');
pdeellip(esp_y*8,+esp_y*4,a,a,0,'E20');

pderect([-esp_y*2.0 esp_y*17.5 -esp_y*3 esp_y*5],'R1');

%  combinação de geometrias
 set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1-(E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+E16+E17+E18+E19+E20)')
%  
% %  Mesh generation:
 setappdata(pde_fig,'Hgrad',1.3); %1.1 para a malha mais refinada
 setappdata(pde_fig,'refinemethod','regular');
 setappdata(pde_fig,'jiggle',char('on','mean',''));
 setappdata(pde_fig,'MesherVersion','preR2013a');
 pdetool('initmesh')
 pdetool('refine')%refino grossa
 pdetool('refine')%refino comum
 pdetool('refine')%refino fina
 
 %%%%%%% CALCULO DO H GLOBAL
 
 clc,clear all,close all

load RESULTADOS_SIMULACAO_RE_N
load trocador_pdetool_alinhado_REFINO_X
format long e

VS.Temp(:,1)=MTemp;
% A descrição desta parte deve constar detalhadamente no trabalho

%Cálculo dos gradiente associados a todos os nós do domínio
[dT_dxi,dT_dyi] = gradiente(VS.Temp(:,1),P,T);
%Cálculo das áreas associadas aos nós do contorno de interesse:
[Aglobal,Aglobalx,Aglobaly] = area_contorno(P,E,5);% O contorno "5" é o de interesse.
%cálculo da derivada normal ao cilindro
dT_dn=(dT_dxi.*Aglobalx+dT_dyi.*Aglobaly)./Aglobal;
% Lembrando que PP.N_cilindro=>nós no cilindro
% cálculo de H local
T_inf=140;
H_L=-PF.k.mix(PP.N_cilindro).*dT_dn(PP.N_cilindro)./(VS.Temp(PP.N_cilindro,1)-T_inf);
%calculo de H_G, coeficiente convectivo global
Acilindro=Aglobal(PP.N_cilindro);
H_G=(1./sum(Acilindro)).*sum(H_L.*Acilindro)

% Re=20
% Pr=0.7155;
 D=1.1;
% PF.k.mix=0.02953;   %condutividade térmica
% 
% Nu=0.911*Re^0.385*Pr^(1/3)
% 
 Nu=H_G/PF.k.mix*D;
%  Nu = mean(Nu)

[Aglobal_01] = area_contorno(P,E,1);%%%%saída
Ts=20;
Te=sum(MTemp.*Aglobal_01)/sum(Aglobal_01);%%%%saída
Ti=140;
DTlm=((Ts-Te)-(Ts-Ti))/(log((Ts-Te)/(Ts-Ti)))
[As] = area_contorno(P,E,5);%%%%entrada

Q=sum(H_G.*As.*DTlm)
 
 
 
 %%%%%% GRAFICOS DO NUSSELT
 clc,clear,close all
Re=linspace(.1,70,1200);
Pr=0.7155;
D=.55;

PF.k.mix=0.02953;   %condutividade térmica

Nu=0.989*Re.^0.330*Pr.^(1/3);
Nu(Re>4)=0.911*(Re(Re>4)).^0.385*Pr.^(1/3);
Nu(Re>40)=0.683*(Re(Re>40)).^0.466*Pr.^(1/3);

Nu2=0.3+((0.62*(Re.^0.5).*Pr.^(1/3))./((1+(0.4/Pr)^(2/3)).^(1/4))).*((1+(Re/282000).^(5/8))).^(4/5);

Nu3=(0.43+0.5*Re.^0.5)*Pr^0.38;% Holman (1981):
% Nu1=0.989*Re.^0.330*Pr.^(1/3);
% Nu2=0.911*Re.^0.385*Pr.^(1/3);
% Nu3=0.683*Re.^0.466*Pr.^(1/3);
% 
% h=Nu*PF.k.mix/D;
% Re=40,Nu=3.3840, malha1, 500ite
% Re=20,Nu=2.6354, malha1, 500ite
% Re=20,Nu=2.5143, malha1, 1000ite ou 
Re1 = 20; Nu_1=2.4290;
Re2 = 30; Nu_2=2.8567;
Re3 = 60; Nu_3= 3.7695;
Re4_1 =50; Nu4_1 =2.5056;
Re4_2 =50;Nu4_2=2.9535;
Re4_3 =50;Nu4_3=3.4521;
% malha pequena Nu=2.5022  
% malha media   Nu=2.5143
% malha grande  Nu=2.5292
% malha m grande  Nu=2.5810
figure(1)
plot(Re,Nu,'-b',[Re1 Re2 Re3 Re4_1 Re4_2 Re4_3],[Nu_1 Nu_2 Nu_3 Nu4_1 Nu4_2 Nu4_3],'pk',Re,Nu2,'-m','LineWidth',1.0,'MarkerSize',10),hold on
% plot(Re,Nu1,'-b',Re,Nu2,'-r',Re,Nu3,'-m','LineWidth',1.0)
% legend('Nu','Nu1','Nu2','Nu3',4)
xlabel('Re'),ylabel('Nu')
 