clc
clear all
close all

u = msgbox( {'Supervisor: Dr. Prasanta Sahoo';'Developed by: M. Barooni, T. Vishwakarma, D. Wilkinson';'';'';'The purpose of this program is to provide the hydrostatic properties of a vessel using the table of offsets.';'';'';'Instructions for use:';'';'1. Prepare the table of offsets in an excel file as directed by the program user manual.';'2. Run the script and select the table of offsets excel file.';'';'';'NOTE-The program should receive the table of offsets in metric units (meters).'},'Vessel Hydrostatics','modal'); 
uiwait(u);

A=uigetfile;
A=readmatrix(A);

[r,c]=size(A);
a=A(2:r,2:c);

num_waterlines = 49;  % set the number of rows interpolated rows you want
num_stations = 49;   % set the number of columns interpolated rows you want
[x, y] = meshgrid(1:size(a,2), 1:size(a,1));
[xq, yq] = meshgrid(linspace(1, size(a, 2), num_stations), linspace(1, size(a, 1), num_waterlines));
ToF = interp2(x, y, a, xq, yq,'spline');
ToF3=ToF.*ToF.*ToF;

waterlines= linspace(A(2,1),A(r,1),49);
w_span=abs(waterlines(1)-waterlines(49));
w_spacing=w_span/48;
waterlines=0:w_spacing:w_span;

stations= linspace(A(1,2),A(1,c),49);
s_span=abs(stations(1)-stations(49));
s_spacing=s_span/48;
stations=0:s_spacing:s_span;

SMS= [1 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 1];
S_levers= [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48];
SMW= [1 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 2 4 1];
W_levers= [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48] ;

T_max=waterlines(1);
hs=stations(49)/48;
LWL=stations(49);

for i=1:49;
    tof=ToF((50-i):49,1:49);
    B_max = max(tof,[],'all');
    B(i)=B_max*2;
    T(i)=waterlines(i);
    hw=T(i)/48;
    
    %assessing for waterline length
    k=find(tof(1,1:49));
    k=size(k,2);
    if k==0
        LWL(i)=0;
    elseif (k+1)<=(num_stations-1)    
        LWL(i)=(k+1)*hs;
    else
        LWL(i)=(num_stations-1)*hs;
    end  
         
    %Calculate Waterplane Area,LCF(forward of aft perp),IT,IL,Cw
    for s = 1:49
            WFA(s) = tof(1,s)*SMS(s);
            SFM(s) = tof(1,s)*SMS(s)*S_levers(s);
            FIT(s) = ToF3((50-i),s)*SMS(s);
            FI0(s) = tof(1,s)*SMS(s)*S_levers(s)*S_levers(s);
    end
    WA(i)= 2*(1/3)*hs*sum(WFA);
    LCF(i)= hs*sum(SFM)/sum(WFA);
    IT(i)=2*(1/9)*hs*sum(FIT);
    I0(i)=2*(1/3)*(hs^3)*sum(FI0);
    IL(i)=I0(i)-WA(i)*LCF(i)^2;
    Cw(i)=WA(i)/(LWL(i)*B(i));
    TPC(i)=WA(i)*1.025/100;    
    
    if size(tof ,1)==1 
        tof=tof;
        VOL(i)=0;
        LCB(i)=0;
        VCB(i)=0;
        BMT(i)=0;
        BML(i)=0;
        Cb(i)=0;
        Cp(i)=0;
        Cm(i)=0;
        Cw(i)=0;
        LCF(i)=0;
               
    else 
    %resize tof    
    [x, y] = meshgrid(1:size(tof,2), 1:size(tof,1));
    [xq, yq] = meshgrid(linspace(1, size(tof, 2), num_stations), linspace(1, size(tof, 1), num_waterlines));
    tof = interp2(x, y, tof, xq, yq);
    tof3=tof.*tof.*tof;
     
    %calculate station areas, volumetric displacement, LCB(forward aft
    %perp)
    for s=1:49
        for w=1:49
            SFA(w)=tof(w,s)*SMW(w);
        end
        SA(s)=2*(1/3)*sum(SFA)*hw;
        FV(s)=SA(s)*SMS(s);
        FVM(s)=FV(s)*S_levers(s);        
    end
    VOL(i)=(1/3)*hs*sum(FV);
    LCB(i)=hs*sum(FVM)/sum(FV);
    
    %calculate waterplane areas and VCB
    for w=1:49
        for s=1:49
            WFA(s)=tof((50-w),s)*SMS(s);  %(14-w) bias ref to baseline
        end
        WA(w)=2*(1/3)*hs*sum(WFA);
        FWM(w)=WA(w)*W_levers(w)*SMW(w);
    end
    Moment=(1/3)*(hw^2)*sum(FWM);
    VCB(i)=Moment/VOL(i);
    
    BMT(i)=IT(i)/VOL(i);
    BML(i)=IL(i)/VOL(i);
    MCTC(i)=1.025*IL(i)/(100*LWL(i));
    Cb(i)=VOL(i)/(LWL(i)*B(i)*T(i));
    Cm(i)=SA(12)/(tof(1,12)*2*T(i));
    Cp(i)=Cb(i)/Cm(i);
    Disp(i)=VOL(i)*1.025;
      
    end
        
end

figure
plot(T,Cb,T,Cm,T,Cp,T,Cw);
title('Form Coefficients')
xlabel('Draft (m)')
ylabel('Coefficient Value')
legend('Cb','Cm','Cp','Cw')
movegui('northwest')

figure
plot(T,LCF,T,LCB);
title('LCF & LCB')
xlabel('Draft (m)')
ylabel('Position (m) fwd station 0')
legend('LCF','LCB')
movegui('west')

figure
plot(T,VCB);
title('VCB')
xlabel('Draft (m)')
ylabel('Position (m) above baseline')
movegui('southwest')

figure
plot(T,Disp);
title('Displacement')
xlabel('Draft (m)')
ylabel('Tons')
movegui('northeast')

figure
plot(T,BMT);
title('BMT')
xlabel('Draft (m)')
ylabel('BMT (m)')
movegui('east')

figure
plot(T,BML);
title('BML')
xlabel('Draft (m)')
ylabel('BML (m)')
movegui('southeast')

figure
plot(T,TPC);
title('Tons per cm immersion')
xlabel('Draft (m)')
ylabel('Tons')
movegui('north')

figure
plot(T,MCTC);
title('Moment to change trim by 1 cm')
xlabel('Draft (m)')
ylabel('Tons-m')
movegui('south')

T11=0:w_span/10:w_span;
BMT11 = interp1(T,BMT,T11,'spline');
BML11 = interp1(T,BML,T11,'spline');
Cb11 = interp1(T,Cb,T11,'spline');
Cm11 = interp1(T,Cm,T11,'spline');
Cw11 = interp1(T,Cw,T11,'spline');
Cp11 = interp1(T,Cp,T11,'spline');
Disp11 = interp1(T,Disp,T11,'spline');
Vol11 = interp1(T,VOL,T11,'spline');
LCB11 = interp1(T,LCB,T11,'spline');
LCF11 = interp1(T,LCF,T11,'spline');
MCTC11 = interp1(T,MCTC,T11,'spline');
WA11 = interp1(T,WA,T11,'spline');
VCB11 = interp1(T,VCB,T11,'spline');
TPC11 = interp1(T,TPC,T11,'spline');

Draft= round(T11.',3);
BMT=round(BMT11.',3);
BML=round(BML11.',3);
Cb=round(Cb11.',3);
Cm=round(Cm11.',3);
Cw=round(Cw11.',3);
Cp=round(Cp11.',3);
Displacement=round(Disp11.',1);
VOL=round(Vol11.',1);
LCB=round(LCB11.',3);
LCF=round(LCF11.',3);
MCTC=round(MCTC11.',3);
WaterplaneArea=round(WA11.',3);
VCB=round(VCB11.',3);
TPC=round(TPC11.',3);
KMT=VCB+BMT;
KML=VCB+BML;

Out=table(Draft,VOL,Displacement,WaterplaneArea,Cb,Cm,Cw,Cp,LCB,LCF,MCTC,TPC,VCB,BMT,KMT,BML,KML);
C={'Draft (m)','Volume in SW (m^3)', 'Disp in SW (t)','Waterplane Area (m^2)','CB','CM','CW','CP','LCB (m fwd transom)','LCF (m fwd transom)','MCTC (tons-m/cm)','TPC (tons/cm)','KB (m)','BMT (m)','KMT (m)','BML (m)','KML (m)'};
writetable(Out,'GradsHydrostatics.xlsx','Range','A2','WriteVariableNames',false);
writecell(C,'GradsHydrostatics.xlsx');