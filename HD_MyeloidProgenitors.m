n=3000;m=1000;
for i=1:n
    a0(i)=0.01;mua(i)=0.1;
    a1(i)=rand;a2(i)=rand;
    a3(i)=rand;a4(i)=rand;
    a5(i)=rand;a6(i)=rand;

    b0(i)=0.01;mub(i)=0.1;
    b1(i)=rand;b2(i)=rand;b3(i)=rand;b4(i)=rand;
    b5(i)=rand;b6(i)=rand;b7(i)=rand;

    c0(i)=0.01;muc(i)=0.1;
    c1(i)=rand;c2(i)=rand;c3(i)=rand;

    d0(i)=0.01;mud(i)=0.1;
    d1(i)=rand;d2(i)=rand;d3(i)=rand;

    e0(i)=0.01;mue(i)=0.1;
    e1(i)=rand;e2(i)=rand;e3(i)=rand;
    
    f0(i)=0.01;muf(i)=0.1;
    f1(i)=rand;f2(i)=rand;f3(i)=rand;
    
    g0(i)=0.01;mug(i)=0.1;
    g1(i)=rand;g2(i)=rand;g3(i)=rand;
    g4(i)=rand;g5(i)=rand;g6(i)=rand; 

    h0(i)=0.01;muh(i)=0.1;
    h1(i)=rand;h2(i)=rand;h3(i)=rand;

    i0(i)=0.01;mui(i)=0.1;
    i1(i)=rand;i2(i)=rand;i3(i)=rand;i4(i)=rand;i5(i)=rand;
    
    j0(i)=0.01;muj(i)=0.1;
    j1(i)=rand;j2(i)=rand;j3(i)=rand;

    for j=1:m
        %Initial value
        iniGATA2(i,j)=10*rand;
        iniGATA1(i,j)=10*rand;
        iniEKLF(i,j)=10*rand;
        iniFli1(i,j)=10*rand;
        iniSCL(i,j)=10*rand;
        iniCEBPa(i,j)=10*rand;
        iniPU1(i,j)=10*rand;
        inicJun(i,j)=10*rand;
        iniEgrNab(i,j)=10*rand;
        iniGfi1(i,j)=10*rand;

        %GATA2=x(1);GATA1=x(2);EKLF=x(3);Fli1=x(4);SCL=x(5);CEBPa=x(6);PU1=x(7);cJun=x(8);EgrNab=x(9);Gfi1=x(10);
        MyeloidProg = @(t,x) [(a0(i)+a1(i)*x(1))/(1+a1(i)*x(1)+a2(i)*x(2)+a3(i)*x(7)+ ...H1
            a4(i)*x(1)*x(2)+a5(i)*x(1)*x(7)+a6(i)*x(2)*x(7))-mua(i)*x(1);
            (b0(i)+b1(i)*x(1)+b2(i)*x(2)+b3(i)*x(4))/(1+b1(i)*x(1)+b2(i)*x(2)+b3(i)*x(4)+b4(i)*x(7)+ ...H2
            (b5(i)*x(1)+b6(i)*x(2)+b7(i)*x(4))*x(7))-mub(i)*x(2);
            (c0(i)+c1(i)*x(2))/(1+c1(i)*x(2)+c2(i)*x(4)+c3(i)*x(2)*x(4))-muc(i)*x(3);
            (d0(i)+d1(i)*x(2))/(1+d1(i)*x(2)+d2(i)*x(3)+d3(i)*x(2)*x(3))-mud(i)*x(4);
            (e0(i)+e1(i)*x(2))/(1+e1(i)*x(2)+e2(i)*x(7)+e3(i)*x(2)*x(7))-mue(i)*x(5);
            (f0(i)+f1(i)*x(6))/(1+f1(i)*x(6)+f2(i)*x(2)*x(5)+f3(i)*x(2)*x(5)*x(6))-muf(i)*x(6);
            (g0(i)+g1(i)*x(6)+g2(i)*x(7)+g3(i)*x(6)*x(7))/(1+g1(i)*x(6)+g2(i)*x(7)+g3(i)*x(6)*x(7)+g4(i)*x(1)+g5(i)*x(2)+g6(i)*x(10))-mug(i)*x(7);
            (h0(i)+h1(i)*x(7))/(1+h1(i)*x(7)+h2(i)*x(10)+h3(i)*x(7)*x(10))-muh(i)*x(8);
            (i0(i)+i1(i)*x(7)*x(8))/(1+i1(i)*x(7)*x(8)+i2(i)*x(10)+i3(i)*x(7)*x(10)+i4(i)*x(8)*x(10)+i5(i)*x(7)*x(8)*x(10))-mui(i)*x(9);
            (j0(i)+j1(i)*x(6))/(1+j1(i)*x(6)+j2(i)*x(9)+j3(i)*x(6)*x(9))-muj(i)*x(10)];
        inival=[iniGATA2(i,j),iniGATA1(i,j),iniEKLF(i,j),iniFli1(i,j),iniSCL(i,j),iniCEBPa(i,j),iniPU1(i,j),inicJun(i,j),iniEgrNab(i,j),iniGfi1(i,j)];
        tspan=[0 8000];

        %Solve the ODE
        [t,x]=ode45(MyeloidProg,tspan,inival);
        GATA2=x(:,1);GATA1=x(:,2);EKLF=x(:,3);Fli1=x(:,4);SCL=x(:,5);
        CEBPa=x(:,6);PU1=x(:,7);cJun=x(:,8);EgrNab=x(:,9);Gfi1=x(:,10);

        %Filter out oscillatory states
        GATA2_EndPart=GATA2(t>=6000&t<=8000);
        PKS=findpeaks(GATA2_EndPart);TF=islocalmin(GATA2_EndPart);VLS=GATA2_EndPart(TF);
        Amplitude=median(PKS)-median(VLS);
        %L_PKS=length(PKS)
        if Amplitude>=0.0001
            ssGATA2(i,j)=0;ssGATA1(i,j)=0;ssEKLF(i,j)=0;ssFli1(i,j)=0;ssSCL(i,j)=0;
            ssCEBPa(i,j)=0;ssPU1(i,j)=0;sscJun(i,j)=0;ssEgrNab(i,j)=0;ssGfi1(i,j)=0;
        else
            ssGATA2(i,j)=GATA2(end);ssGATA1(i,j)=GATA1(end);ssEKLF(i,j)=EKLF(end);ssFli1(i,j)=Fli1(end);ssSCL(i,j)=SCL(end);
            ssCEBPa(i,j)=CEBPa(end);ssPU1(i,j)=PU1(end);sscJun(i,j)=cJun(end);ssEgrNab(i,j)=EgrNab(end);ssGfi1(i,j)=Gfi1(end);
        end
    end
end


%ssForParamt: n sets of parameters correspond to n tables of number
paramt=zeros(10,m,n);ssForParamt=zeros(10,m,n);
p=[];P=zeros(n,m);
for k=1:n
    paramt(:,:,k)=[ssGATA2(k,:);ssGATA1(k,:);ssEKLF(k,:);ssFli1(k,:);ssSCL(k,:);
        ssCEBPa(k,:);ssPU1(k,:);sscJun(k,:);ssEgrNab(k,:);ssGfi1(k,:);];
    ssForParamt(:,:,k)=AproxMerge(paramt(:,:,k),0.01);   
    Aprxm_ssGATA2(k,:)=ssForParamt(1,:,k);
    Aprxm_ssGATA1(k,:)=ssForParamt(2,:,k);
    Aprxm_ssEKLF(k,:)=ssForParamt(3,:,k);
    Aprxm_ssFli1(k,:)=ssForParamt(4,:,k);
    Aprxm_ssSCL(k,:)=ssForParamt(5,:,k);
    Aprxm_ssCEBPa(k,:)=ssForParamt(6,:,k);
    Aprxm_ssPU1(k,:)=ssForParamt(7,:,k);
    Aprxm_sscJun(k,:)=ssForParamt(8,:,k);
    Aprxm_ssEgrNab(k,:)=ssForParamt(9,:,k);
    Aprxm_ssGfi1(k,:)=ssForParamt(10,:,k);  
    Logic_SS(k,:)=logical(ssForParamt(1,:,k));
end

data_ssGATA2=Aprxm_ssGATA2(Aprxm_ssGATA2~=0);length(data_ssGATA2);
data_ssGATA1=Aprxm_ssGATA1(Aprxm_ssGATA1~=0);length(data_ssGATA1);
data_ssEKLF=Aprxm_ssEKLF(Aprxm_ssEKLF~=0);length(data_ssEKLF);
data_ssFli1=Aprxm_ssFli1(Aprxm_ssFli1~=0);length(data_ssFli1);
data_ssSCL=Aprxm_ssSCL(Aprxm_ssSCL~=0);length(data_ssSCL);
data_ssCEBPa=Aprxm_ssCEBPa(Aprxm_ssCEBPa~=0);length(data_ssCEBPa);
data_ssPU1=Aprxm_ssPU1(Aprxm_ssPU1~=0);length(data_ssPU1);
data_sscJun=Aprxm_sscJun(Aprxm_sscJun~=0);length(data_sscJun);
data_ssEgrNab=Aprxm_ssEgrNab(Aprxm_ssEgrNab~=0);length(data_ssEgrNab);
data_ssGfi1=Aprxm_ssGfi1(Aprxm_ssGfi1~=0);length(data_ssGfi1);


%Arrange these data in columns
data=[data_ssGATA2,data_ssGATA1,data_ssEKLF,data_ssFli1, ...
    data_ssSCL,data_ssCEBPa,data_ssPU1,data_sscJun,data_ssEgrNab,data_ssGfi1];
P=[a1',a2',a3',a4',a5',a6',...
    b1',b2',b3',b4',b5',b6',b7',...
    c1',c2',c3',...
    d1',d2',d3',...
    e1',e2',e3',...
    f1',f2',f3',...
    g1',g2',g3',g4',g5',g6',...
    h1',h2',h3',...
    i1',i2',i3',i4',i5',...
    j1',j2',j3'];

%Append the data results
%writematrix(data,'HD_MyeloidProgenitors_data2.txt','WriteMode','append')
%writematrix(P,'HD_MyeloidProgenitors_P2.txt','WriteMode','append')
%writematrix(Logic_SS,'HD_MyeloidProgenitors_Logic_SS2.txt','WriteMode','append')