
function faraday
clear all
%close all
eps1=1;
eps3=1;

epsa1=5.6;
epsa2=5.6;
gamaa=0.1;
epsa3=5.6;

epsb1=2.1;
epsb2=2.1;
gamab=0;
epsb3=2.1;

Na=4;

epsa=[epsa1 0 -1i*gamaa;0 epsa2 0;1i*gamaa 0 epsa3 ];

epsb=[epsb1 0 -1i*gamab;0 epsb2 0;1i*gamab 0 epsb3 ];

%epsa=epsb;

color='-or';

%memory

unt=1;

% a=1.17*unt;
% R=.504*unt;

a=1*unt;
R=a/4;

a1=a;
a2=a;
d1=0*a;
d2=a2/2-R+0*a;

t1=cputime;

theta=0;

nGx=4;
nGy=120;

wn1=.3;

wn2=.9;


ndiv=18*(wn2-wn1)*10;

dwn=(wn2-wn1)/ndiv;

plotWave=0;

cf=(plotWave==0)

for p=1:1*cf*ndiv+1
    
    p
    
    Fn(p)=wn1+dwn*(p-1);
    k1=2*pi*Fn(p)/a;
    
    
    [Ts Rs, Fs]=tranRef(epsa,epsb,eps1,eps3,a1,a2,R,d1,d2,Na,nGx,nGy,k1,p,theta,plotWave,color);
           
      if(abs(Ts)>1) 
        Ts=1;
      end
      
    T(p)=real(Fs);
    T2(p)=real(Ts);
    
    uu=Ts+Rs;
    
    uu
end

t2=cputime;

comptation_time=t2-t1;
comptation_time;
        
        T'
        T2'
        
        if(length(T)>1)
                figure(1)
             plot(Fn,T,color);
 
             axis([wn1,wn2,-90,90]);
             hold on
            
                figure(2)
               plot(Fn,T2,'-or');
                 axis([wn1,wn2,min(T2),max(T2)]);
                %axis([wn1,wn2,0,1]);
                 hold on
            
        end

end



function [Ts Rs Fr ]=tranRef(epsa,epsb,eps1,eps3,a1,a2,R,d1,d2,Na,nGx,nGy,k1,p,theta,plotWave,colorAng)

twoMat=1;

d=d1+d2;

triang=0;

global MM;

global Kapa;
global bE;

w2c2=(k1)^2;

L=(Na-1)*a2+2*(R+d);

Lx1=2*nGx+1;
Lx2=3*nGx+2;

Ly1=2*nGy+1;


numG=(Lx1)*nGy;

dd=Lx1;

E0=[0 0 1]';


kx=k1*sin(theta);
k1y=k1*cos(theta);
k3=sqrt(eps3*w2c2);
bx=2*pi/a1;

by=pi/L;

pph=2/pi;


nk=size(epsb,2);

if(p==1)
    disp('Computing the Fourier series ..');
    Kapa=zeros(4*nGx+1,4*nGy+1,4);
    if(triang)
        FillKapaTriang(nGx,nGy,epsa,epsb,L,R,Na,a1,a2);  %triangular
        %lattivce not implemented properly.
    else
       %FillKapaCylinderAntiSym(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1);
        FillKapaAntiSymNum(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1);
        %FillKapa1D(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1);
    end
    
    if(2>11)
        
        nd=50;
        
        x=linspace(-a1/2,a1/2,nd);
        y=linspace(-L,L,nd);
        
        [xx yy]=ndgrid(x,y);
        
        zz=zeros(size(xx));
        
        for ix=1:nd
            x1=x(ix);
            for iy=1:nd
                y1=y(iy);
                tt=0;
                for n=-nGx:nGx
                    for m=-nGy:nGy
                        Gn=bx*n;
                        Gm=by*m;
                        
                        TT=Kapa(n+2*nGx+1,m+2*nGy+1,4);
                        
                        tt=tt+TT*exp(1i*(Gn*x1+Gm*y1));
                    end
                end
                
                zz(ix,iy)=real(tt);
            end
        end
        
        
        X=abs(Kapa(:,:,1));
        
        figure(4);
        surf(xx,yy,zz);
        %  surf(X);
        axis equal;
        set(gca,'DataAspectRatio',[1 1 .05]);
        return;
    end
    disp('Computing matrix, step 1...');
    
    dimx=nk*(numG+2*(Lx1));
   
    MM=  single(zeros(dimx,dimx));
    

        
    bE=zeros(nk*(numG+2*(Lx1)),1);
    
    countG=0;
    
    
    for Gx=-nGx:nGx
        for Gy=1:nGy
            
            countG=countG+1;
            countG1=0;
            
            for Gx1=-nGx:nGx
                
                dGxp=Lx1+Gx-Gx1;
                Gx1p=nGx+1+Gx1;
                Gx1pp=Lx2+Gx1;


                
                kxn1=kx+Gx1*bx;
               % kxn1=0;
                
                kxn2=kxn1^2;
             %    kxn2=1;
                
                tempT=zeros(nk,nk);
                tempR=zeros(nk,nk);
                tempE=zeros(nk,nk);
                
                for Gy1=1:nGy
                    kym=Gy1*by;
                    dGyp=Ly1+(Gy-Gy1);
                    sGyp=Ly1+Gy+Gy1;
                    countG1=countG1+1;
                         
       
                    
                    Kapas(1:4)=Kapa(dGxp,sGyp,1:4);
      
                                        
                    Kapad(1:4)=Kapa(dGxp,dGyp,1:4);
                    
                    KapTenss=[Kapas(1) 0 1i*Kapas(4); 0 Kapas(2) 0 ;
                        -1i*Kapas(4) 0  Kapas(3)];

                    
                    KapTensd=[Kapad(1) 0 1i*Kapad(4); 0 Kapad(2) 0 ;
                        -1i*Kapad(4) 0  Kapad(3)];

                    
                    Akg=zeros(3,3);
                    
                    
                    Akg(1,1)=kym^2;
                    Akg(1,2)=1i*kxn1*kym;
                    Akg(2,1)= Akg(1,2);
                    Akg(2,2)=kxn2;
                    
                    Akg(3,3)=Akg(1,1)+ Akg(2,2);

                    AkdKap=(KapTenss*Akg-KapTensd*conj(Akg));
                    
                    Akg1=zeros(3,3);
                    
                    if(mod(Gy1,2)==0)
                        Akg1(1,2)=-1i*kxn1*pph/Gy1;
                    else
                        Akg1(1,2)=1i*kxn1*pph/Gy1;
                    end
                    
                    Akg1(2,1)=Akg1(1,2);
                    
                    if(mod(Gy1,2)==0)
                        Akg1(2,2)=kxn2*pph/Gy1;
                    else
                        Akg1(2,2)=-kxn2*pph/Gy1;
                    end
                    
                    
                    Akg1(3,3)=Akg1(2,2);
                    
                    
                    dKapTens1=-1*(KapTenss*Akg1-KapTensd*conj(Akg1));
                    
              
                    tempT=tempT+dKapTens1;

                    
                   Akg2=zeros(3,3);
                    
                    Akg2(1,2)=-1i*kxn1*pph/Gy1;
         
                    Akg2(2,1)=Akg2(1,2);
                    Akg2(2,2)=kxn2*pph/Gy1;
                   Akg2(3,3)=Akg2(2,2);


                   dKapTens2=-1*(KapTenss*Akg2-KapTensd*conj(Akg2));
           
                    
                   tempR=tempR-dKapTens2;
                    
                    if(Gx1==0)
                        
                       tempE=tempE-dKapTens2;
                    end
                    
                    r1=(countG-1)*nk;
                    c1=(countG1-1)*nk;
                    for j=1:nk
                        for k=1:nk
                            MM(r1+j,c1+k)= AkdKap(j,k);
                            
                        end
                    end
                    
                    
                end
                
                
                
                if(Gx1==0)
                    v=tempE*E0;
                    for k=1:nk
                        bE((countG-1)*nk+k)=v(k);
                    end
                    
                end
                
                
                for j=1:nk
                    for k=1:nk
                        
                        MM((countG-1)*nk+j,(numG+Gx1p-1)*nk+k)=tempT(j,k);
                        MM((countG-1)*nk+j,(numG+dd+Gx1p-1)*nk+k)=tempR(j,k);
                        
                    end
                end
                
                
                
                if(Gx==Gx1)
                    for k=1:nk
                        if(k~=2)
                            
                            MM((Gx1p+numG-1)*nk+k,((Gx1p-1)*nGy+Gy-1)*nk+k)=pi*Gy;
                        end
                    end
                    
                    if(mod(Gy,2)==0)
                        for k=1:nk
                            if(k~=2)
                                MM((Gx1pp+numG-1)*nk+k,((Gx1p-1)*nGy+Gy-1)*nk+k)=pi*Gy;
                            end
                        end
                    else
                        for k=1:nk
                            if(k~=2)
                                MM((Gx1pp+numG-1)*nk+k,((Gx1p-1)*nGy+Gy-1)*nk+k)=-pi*Gy;
                            end
                        end
                    end
                    
                end
                
            end
            
            
        end
        
        
    end
    
end


bN=zeros(nk*(numG+2*(Lx1)),1);

disp('Computing matrix, step 2...');

if(twoMat)
    NN=single(zeros(nk*(numG+2*(Lx1))));
    
    
    
    for Gx=1:size(NN,1)
        for Gy=1:size(NN,2)
            NN(Gx,Gy)=single(1e-10+1i*1e-10);
        end
    end
    
end

countG=0;

for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    Gxpp=Lx2+Gx;
    
   
    for Gy=1:nGy
        
        countG=countG+1;
        for k=1:nk
            if(twoMat)
                NN((countG-1)*nk+k,(countG-1)*nk+k)= w2c2;
            else
                MM((countG-1)*nk+k,(countG-1)*nk+k)= MM((countG-1)*nk+k,(countG-1)*nk+k)+w2c2;
            end
        end
        
        if(mod(Gy,2)==0)
            for k=1:nk
                if(twoMat)
                    NN((countG-1)*nk+k,(numG+Gxp-1)*nk+k)=-pph/Gy*w2c2;
                else
                    MM((countG-1)*nk+k,(numG+Gxp-1)*nk+k)=MM((countG-1)*nk+k,(numG+Gxp-1)*nk+k)-pph/Gy*w2c2;
                end
            end
        else
            for k=1:nk
                if(twoMat)
                    NN((countG-1)*nk+k,(numG+Gxp-1)*nk+k)=pph/Gy*w2c2;
                else
                    MM((countG-1)*nk+k,(numG+Gxp-1)*nk+k)= MM((countG-1)*nk+k,(numG+Gxp-1)*nk+k)+pph/Gy*w2c2;
                end
            end
        end
        
        for k=1:nk
            if(twoMat)
                NN((countG-1)*nk+k,(numG+Gxpp-1)*nk+k)=pph/Gy*w2c2;
            else
                MM((countG-1)*nk+k,(numG+Gxpp-1)*nk+k)=MM((countG-1)*nk+k,(numG+Gxpp-1)*nk+k)+pph/Gy*w2c2;
            end
            
            
            
            
            if(Gx==0)
                
                bN((countG-1)*nk+k)=-pph/Gy*w2c2*E0(k);
                
            end
        end
        
    end
end


real(MM);

for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    Gxpp=Lx2+Gx;
    
    kxn1=kx+Gx*bx;
    kxn2=kxn1^2;
    
    
    if(k1>=abs(kxn1))
        krny=-sqrt(k1^2-kxn2);
    else
        krny=-1i*sqrt(kxn2-k1^2);
    end
    
    
    if(k3>=abs(kxn1))
        ktny=sqrt(k3^2-kxn2);
    else
        ktny=1i*sqrt(kxn2-k3^2);
    end
    

    
    for k=1:nk
        
        if(twoMat)
            if(k~=2)
                NN((numG+Gxp-1)*nk+k,(numG+Gxp-1)*nk+k)=1;
                
                
                NN((numG+Gxp-1)*nk+k,(numG+Gxpp-1)*nk+k)=-(1i*L*krny+1);
                
                
                NN((numG+Gxpp-1)*nk+k,(numG+Gxp-1)*nk+k)=-(1i*L*ktny-1);
                
                NN((numG+Gxpp-1)*nk+k,(numG+Gxpp-1)*nk+k)=-1;
            else
                
                NN((numG+Gxp-1)*nk+k,(numG+Gxp-1)*nk+k)=1;
            
                NN((numG+Gxpp-1)*nk+k,(numG+Gxpp-1)*nk+k)=1;
            end
        else
            if(k~=2)
                MM((numG+Gxp-1)*nk+k,(numG+Gxp-1)*nk+k)=1;
                
                
                MM((numG+Gxp-1)*nk+k,(numG+Gxpp-1)*nk+k)=-(1i*L*krny+1);
                
                
                NN((numG+Gxpp-1)*nk+k,(numG+Gxp-1)*nk+k)=-(1i*L*ktny-1);
                
                MM((numG+Gxpp-1)*nk+k,(numG+Gxpp-1)*nk+k)=-1;
            else
                
                MM((numG+Gxp-1)*nk+k,(numG+Gxp-1)*nk+k)=1;
            
                MM((numG+Gxpp-1)*nk+k,(numG+Gxpp-1)*nk+k)=1;
            end
            
        end
        
        if(Gx==0)
                      
            if(k~=2)
            bN((numG+Gxp-1)*nk+k)= (1i*L*k1y+1)*E0(k);
            bN((numG+Gxpp-1)*nk+k)=E0(k);
            end
        end
    end
end


bN=bN+bE;

disp('solving matrix...');

if(twoMat)
    NN=NN+MM;   
    x=linsolve(NN,bN);
else

    x=linsolve(MM,bN);
end
% figure(12);
% surf(imag(NN(1:end,1:end)));
% plot(real(x),'g');
% hold on

Anm=zeros(2*nGx+1,nGy,3);
Tn=zeros(2*nGx+1,3);
Tn2=zeros(2*nGx+1,1);
Rn=zeros(2*nGx+1,3);
Rn2=zeros(2*nGx+1,1);

ix=0;
for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    for Gy=1:nGy
        for k=1:nk
            ix=ix+1;
            Anm(Gxp,Gy,k)=x(ix);
        end
    end
end


kp=numG*nk;
kpp=(numG+2*nGx+1)*nk;
for Gx=-nGx:nGx
    k=nGx+Gx;
    for j=1:nk
        Tn(k+1,j)=x(kp+k*nk+j);
        Tn2(k+1)= Tn2(k+1)+abs(Tn(k+1,j))^2;
        
        Rn(k+1,j)=x(kpp+k*nk+j);
        Rn2(k+1)= Rn2(k+1)+abs(Rn(k+1,j))^2;
    end
end



Ts=0;
Rs=0;

for Gx=-nGx:nGx
    k=nGx+1+Gx;
    kxn1=kx+Gx*bx;
    kxn2=kxn1^2;
    
    if(k1>=abs(kxn1))
        krny=-sqrt(k1^2-kxn2);
    else
        krny=-1i*sqrt(kxn2-k1^2);
    end
    
    
    if(k3>=abs(kxn1))
        ktny=sqrt(k3^2-kxn2);
    else
        ktny=1i*sqrt(kxn2-k3^2);
    end
    
    
  %  Ts=Ts+ktny/k3*sqrt(eps3/eps1)*Tn2(k)/cos(theta);
        Ts=Ts+ktny/k3*sqrt(eps3/eps1)*abs(Tn(k,3))^2/cos(theta);

    
    Rs=Rs+abs(krny)/k1*Rn2(k)/cos(theta);
end

Fr=0;


nL=max(int32(L*k1/2/pi)*17*2,35);



yy=linspace(0,L,nL);

Nx=1;

if(Nx==1)
    xx=zeros(1,1);
else
    xx=linspace(-a1/2,a1/2,Nx);
end

ff=zeros(Nx,nL,nk);
si=zeros(Nx,nL,nk);

for ix=1:Nx
    
    x=xx(ix);
    
    for k=1:nL
        y=yy(k);
        
        for Gx=-nGx:nGx
            
            kxn=kx+Gx*bx;
            
            Gxp=nGx+1+Gx;
            del=(Gx==0);
            for Gy=1:nGy
                for j=1:nk
                    si(ix,k,j)=si(ix,k,j)+Anm(Gxp,Gy,j)*sin(Gy*by*(y))*exp(1i*kxn*x);
                end
            end
            for j=1:nk
                ff(ix,k,j)=ff(ix,k,j)+(y*Tn(Gxp,j)+(L-y)*(Rn(Gxp,j)+del*E0(j)));
            end
        end
        
    end
end

E1=ff/L+si;

E2=real(E1);

E2(:,:,1);

%writeMeshAndField(Nx,nL,1,E2,2,Na);



%surf(x2,y2,E2(:,:,2));
plotWave2=false;

if(plotWave2)
    figure(5)
    
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-2 2 0 L -2 2]);
    az = 40;
    el = 30;
    %  az=90;
    %  el=0;
    view(az, el);
    
    
    
    hold on
    
    for ix=1:Nx
        x=xx(ix);
        Vx=E2(ix,1,1);
        Vy=E2(ix,1,2);
        Vz=E2(ix,1,3);
        
        
        for k=1:nL
            
            y=yy(k);
            Vxp=Vx;
            Vx=E2(ix,k,1);
            
            Vyp=Vy;
            Vy=E2(ix,k,2);
            
            Vzp=Vz;
            Vz=E2(ix,k,3);
            
            
            color = 'r';
            
            
            if(k>1)
                h1(k)=plot3([Vxp+x Vx+x],[Vyp+yy(k-1) Vy+y],[Vzp Vz],'LineWidth',2,'Color','b');
                h2(k)= plot3([-Vxp+x -Vx+x],[yy(k-1) y],[-Vzp -Vz],'LineWidth',2,'Color','b');
                h3= plot3([0 0],[0 L],[0 0],'LineWidth',2,'Color','k');
                
                
            end
            
            
            arrow1(k) = arrow3d([x x+Vx],[y y+Vy],[0 Vz],.92,.01,.02,color);
            arrow1(k) = arrow3d([x x-Vx],[y y-Vy],[0 -Vz],.92,.01,.02,color);
            
            
            
        end
        
    end;
end



tn0=E2(1,1,1)/E2(1,1,3);

ang0=atan(tn0);
nan=0;

tans1=zeros(nL,1);
for k=1:nL
    
    Vx=E2(1,k,1);
    Vz=E2(1,k,3);
    
    magn=sqrt(Vx^2+Vz^2);
    
    if( magn>.1)
        tans1(k)=Vx/Vz;
        nan=nan+1;
    else
        tans1(k)=1e10;
    end
end

tans=zeros(nan,1);
y=zeros(nan,1);
j=0;
for k=1:nL
    if(tans1(k)~=1e10)
        j=j+1;
        tans(j)=tans1(k);
        y(j)=yy(k);
    end
end;

angs=zeros(nan,1);
angs(1)=atan(tans(1));

for k=2:nan
    
    tnp=tans(k-1);
    tn=tans(k);
    angs(k)=atan(tn)*180/pi;
    %  angs(k)=E2(1,k,3);
    % angs(k)=angs(k-1)+atan((tn-tnp)/(1+tn*tnp))*180/pi;
    
    
end



Fr=(angs(nan)-ang0);

%Eout=sqrt(E2(1,nL,1)^2+E2(1,nL,2)^2+E2(1,nL,3)^2);

if(plotWave)
    
    figure(7)
    
    plot(y,angs,colorAng);
    hold on
    
    angsHomog=-0.5*y*k1*imag(epsb(1,3))/sqrt(epsb(1,1))*180/pi;
    plot(y,angsHomog,'-k');
    
end

end




function FillKapa1D(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;

invepsa1=inv(epsa);

invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];

invepsb1=inv(epsb);

invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

nk=length(invepsa);

nGx=0;
rx=a2/2;
ry=R;
ngrid=30;
Lh=Na*a1-a1/2;

invep=zeros(ngrid,Na*2*ngrid,3);
%The following loop carries out the definition of the unit cell.
for n=0:2*Na-1
    nx=1;
    for countX=-a1/2:a1/ngrid:a1/2
        ny=1;
        for countY=-a2/2:a2/ngrid:a2/2
            %Following condition allows to define the circle with of radius r
            inside=0;
            
            inside=abs(countX/rx)<1 && abs(countY/ry)<1;
            
            
            if(inside)
                
                ip=invepsa;
                
            else
                
                ip=invepsb;
            end
            
            
            
            for k=1:nk
                inveps(nx,ny+ngrid*n,k)=ip(k);
            end
            
            xSet(nx)=countX;
            ySet(ny+ngrid*n)=countY+a2*n-Lh;
            ny=ny+1;
        end
        nx=nx+1;
    end
end

for n=1:length(ySet)/2
    inveps(:,n,nk)=-inveps(:,n,nk);
end


MtNt=(length(xSet)-1)*(length(ySet)-1);
%The next loop computes the Fourier expansion coefficients
bx=2*pi/a1;

by=pi/L;


for dGx=-2*nGx:2*nGx
    
    for dGy=-2*nGy:2*nGy
        
        
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        for nx=1:length(xSet)-1
            for ny=1:length(ySet)-1
                x=xSet(nx);
                y=ySet(ny);
                tt=dGx*bx*x+dGy*by*y;
                
                for k=1:nk
                    Kapa(dGxp,dGyp,k)=Kapa(dGxp,dGyp,k)+inveps(nx,ny,k)*exp(-1i*(tt));
                end
            end
        end
        
        
    end
end

Kapa=  Kapa/MtNt;


end


function FillKapaCylinderAntiSym(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;


invepsa1=inv(epsa);

invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];

invepsb1=inv(epsb);

invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

nk=size(invepsb,2);

ff=Na*pi*R*R/(a1*L);
d11=d1/L;

for dGx=-2*nGx:2*nGx
    Gn=2*pi*dGx/a1;
    for dGy=-2*nGy:2*nGy
        Gm=dGy*pi/L;
        GnmR=sqrt(Gn*Gn+Gm*Gm)*R;
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        if(dGx==0&& dGy==0)
            for k=1:nk
                if(k~=4)
                    Kapa(dGxp,dGyp,k)=1*(ff*invepsa(k)+(1-ff)*invepsb(k))+2*d11*(invepsa(k)-invepsb(k));
                end
                
            end
            
        elseif(dGy==0)
            tt=2*ff*besselj(1,GnmR)/GnmR;
            
            for k=1:nk
                if(k~=4)
                    Kapa(dGxp,dGyp,k)=tt*(invepsa(k)-invepsb(k));
                end
            end
            
        else
            tt=a2*dGy*pi/(2*L);
            
            
            for k=1:nk
                if(k==4)
                    factm=-1i*sin(dGy*pi/2)*sin(Na*tt)/sin(tt)/Na;
                    vv=2*ff*besselj(1,GnmR)/GnmR;
                    
                    Kapa(dGxp,dGyp,k)=1*factm*vv*(invepsa(k)-invepsb(k));
                    if(dGx==0)
                        vvL=2*sin(dGy*pi*.5*(L-2*d1)/L)/(dGy*pi);
                        factm=-1i*sin(dGy*pi*.5);
                        
                        Kapa(dGxp,dGyp,k)=Kapa(dGxp,dGyp,k)+factm*vvL*invepsb(k);
                    end
                else
                    factm=cos(dGy*pi/2)*sin(Na*tt)/sin(tt)/Na;
                    vv=2*ff*factm*besselj(1,GnmR)/GnmR;
                    
                    Kapa(dGxp,dGyp,k)=1*vv*(invepsa(k)-invepsb(k));
                end
                
                
                if(dGx==0 && d11>0 && k~=4)
                    aa= pi*dGy*d11;
                    
                    zz=1*sin(aa)/aa*d11+2*cos(pi*dGy*(L-d1/2)/L)*sin(aa/2)/aa*d11;
                    Kapa(dGxp,dGyp,k)= Kapa(dGxp,dGyp,k)+zz*(invepsa(k)-invepsb(k));
                end
                
            end
            
        end
        
        
    end
end
end



function FillKapaCylinderSym(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;

invepsa1=inv(epsa);

invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];

invepsb1=inv(epsb);

invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];


nk=length(invepsa);
ff=Na*pi*R*R/(a1*L);
d11=d1/L;

for dGx=-2*nGx:2*nGx
    Gn=2*pi*dGx/a1;
    for dGy=-2*nGy:2*nGy
        Gm=dGy*pi/L;
        GnmR=sqrt(Gn*Gn+Gm*Gm)*R;
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        if(dGx==0&& dGy==0)
            for k=1:nk
                
                Kapa(dGxp,dGyp,k)=1*(ff*invepsa(k)+(1-ff)*invepsb(k))+2*d11*(invepsa(k)-invepsb(k));
                
                
            end
            
        elseif(dGy==0)
            tt=2*ff*besselj(1,GnmR)/GnmR;
            
            for k=1:nk
                
                Kapa(dGxp,dGyp,k)=tt*(invepsa(k)-invepsb(k));
                
            end
            
        else
            tt=a2*dGy*pi/(2*L);
            
            
            for k=1:nk
                
                factm=cos(dGy*pi/2)*sin(Na*tt)/sin(tt)/Na;
                vv=2*ff*factm*besselj(1,GnmR)/GnmR;
                
                Kapa(dGxp,dGyp,k)=1*vv*(invepsa(k)-invepsb(k));
                
                
                
                if(dGx==0 && d11>0)
                    aa= pi*dGy*d11;
                    
                    zz=sin(aa)/aa*d11+2*cos(pi*dGy*(L-d1/2)/L)*sin(aa/2)/aa*d11;
                    Kapa(dGxp,dGyp,k)= Kapa(dGxp,dGyp,k)+zz*(invepsa(k)-invepsb(k));
                end
                
            end
            
        end
        
        
    end
end
end
function FillKapaAntiSymNum(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)

invepsa1=inv(epsa);

invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];

invepsb1=inv(epsb);

invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

nk=length(invepsa);

global Kapa;
elementShape=1;
%0 for ellipse, 1 for rectangle
rx=1*R;
ry=2.5*R;

ngrid=30;
Lh=Na*a1-a1/2;
d11=d1/L;
invep=zeros(ngrid,Na*2*ngrid,nk);
%The following loop carries out the definition of the unit cell.
zz=zeros(ngrid+1,Na*2*ngrid+1);

for n=0:2*Na-1
    nx=1;
    for countX=-a1/2:a1/ngrid:a1/2
        ny=1;
        for countY=-a2/2:a2/ngrid:a2/2
            %Following condition allows to define the circle with of radius r
            
            
            if(elementShape==0)
                inside=(countX/rx)^2+(countY/ry)^2<1;
            else
                inside=abs(countX/rx)<1 && abs(countY/ry)<1;
            end
            
            if(inside)
                
                ip=invepsa;
                
                zz(nx,n*ngrid+ny)=1;
            else
                
                ip=invepsb;
            end
            
            
            
            for k=1:nk
                inveps(nx,ny+ngrid*n,k)=ip(k);
            end
            
            xSet(nx)=countX;
            ySet(ny+ngrid*n)=countY+a2*n-Lh;
            ny=ny+1;
        end
        nx=nx+1;
    end
end

for n=1:length(ySet)/2
    inveps(:,n,nk)=-inveps(:,n,nk);
end


%=========== The following lines plots the lattice shape. starts from here

plot =1; % if plot=0. it does nothing. If plot =1 it plots the graph and 

if(plot==1)
[xx yy]=ndgrid(xSet,ySet);

figure(4);
        surf(xx,yy,zz);
        %  surf(X);
        axis equal;
        set(gca,'DataAspectRatio',[1 1 1]);

 error('Shape plotted. Program ends stops without calculation. For calculation, set plot==0 ( in line 1128).');
end
 %=========== End of plot 

MtNt=(length(xSet)-1)*(length(ySet)-1);
%The next loop computes the Fourier expansion coefficients
bx=2*pi/a1;

by=pi/L;


for dGx=-2*nGx:2*nGx
    
    for dGy=-2*nGy:2*nGy
        
        
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        for nx=1:length(xSet)-1
            for ny=1:length(ySet)-1
                x=xSet(nx);
                y=ySet(ny);
                tt=dGx*bx*x+dGy*by*y;
                
                for k=1:nk
                    Kapa(dGxp,dGyp,k)=Kapa(dGxp,dGyp,k)+inveps(nx,ny,k)*exp(-1i*(tt));
                end
            end
        end
        
        
    end
end

Kapa=  Kapa/MtNt;


end

function [h]=arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
 %The function plotting 3-dimensional arrow
%
% h=arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
% The inputs are:
%       x,y,z =  vectors of the starting point and the ending point of the
%           arrow, e.g.:  x=[x_start, x_end]; y=[y_start, y_end];z=[z_start,z_end];
%       head_frac = fraction of the arrow length where the head should  start
%       radii = radius of the arrow
%       radii2 = radius of the arrow head (defult = radii*2)
%       colr =   color of the arrow, can be string of the color name, or RGB vector  (default='blue')
%
% The output is the handle of the surfaceplot graphics object.
% The settings of the plot can changed using: set(h, 'PropertyName', PropertyValue)
%
% example #1:
%        arrow3d([0 0],[0 0],[0 6],.5,3,4,[1 0 .5]);
% example #2:
%        arrow3d([2 0],[5 0],[0 -6],.2,3,5,'r');
% example #3:
%        h = arrow3d([1 0],[0 1],[-2 3],.8,3);
%        set(h,'facecolor',[1 0 0])
%
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% July 2010 (C)

if nargin==5
    radii2=radii*2;
    colr='blue';
elseif nargin==6
    colr='blue';
end
if size(x,1)==2
    x=x';
    y=y';
    z=z';
end

x(3)=x(2);
x(2)=x(1)+head_frac*(x(3)-x(1));
y(3)=y(2);
y(2)=y(1)+head_frac*(y(3)-y(1));
z(3)=z(2);
z(2)=z(1)+head_frac*(z(3)-z(1));
r=[x(1:2)',y(1:2)',z(1:2)'];

N=10;
dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;

normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;

X1=[];Y1=[];Z1=[];
j=1;
for theta=([0:N])*2*pi./(N);
    R1=Pc+radii*cos(theta).*(P1-Pc) + radii*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(2:3,j)=R1(:,1);
    Y1(2:3,j)=R1(:,2);
    Z1(2:3,j)=R1(:,3);
    j=j+1;
end

r=[x(2:3)',y(2:3)',z(2:3)'];

dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;

normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;

j=1;
for theta=([0:N])*2*pi./(N);
    R1=Pc+radii2*cos(theta).*(P1-Pc) + radii2*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(4:5,j)=R1(:,1);
    Y1(4:5,j)=R1(:,2);
    Z1(4:5,j)=R1(:,3);
    j=j+1;
end

X1(1,:)=X1(1,:)*0 + x(1);
Y1(1,:)=Y1(1,:)*0 + y(1);
Z1(1,:)=Z1(1,:)*0 + z(1);
X1(5,:)=X1(5,:)*0 + x(3);
Y1(5,:)=Y1(5,:)*0 + y(3);
Z1(5,:)=Z1(5,:)*0 + z(3);

h=surf(X1,Y1,Z1,'facecolor',colr,'edgecolor','none');

light('Color','r')

lighting phong
end



function writeMeshAndField(nx,ny,nz, E,mode,Na)

factor=1;


Ne=nx*ny*nz;

Nn=Ne*8;
ie=0;

elVert=zeros(Ne,8);

for kz=1:nz
    for ky=1:ny
        for kx=1:nx
            ie=ie+1;
            for j=1:8
                elVert(ie,j)=(ie-1)*8+j;
            end
        end
    end
end




% 	Nn=(ny+1)*(nx+1)*(nz+1);
% 		Nx=nx+1;
% 		Nxy=(nx+1)*(ny+1);
% 		Ne=nx*ny*nz;
%
%
%     frx=1;
%     frxy=1;
%     ie=0;
%
%      elVert=zeros(Ne,8);

%     while(frxy*Nxy<Nn)
%    ie=ie+1;
%         elVert(ie,1)=frxy*Nxy+frx;
%           elVert(ie,2)=frxy*Nxy+1+frx;
%             elVert(ie,3)=frxy*Nxy+Nx+1+frx;
%               elVert(ie,4)=frxy*Nxy+Nx+frx;
%               elVert(ie,5)=(frxy-1)*Nxy+frx;
%                 elVert(ie,6)=(frxy-1)*Nxy+1+frx;
%                   elVert(ie,7)=(frxy-1)*Nxy+Nx+1+frx;
%                     elVert(ie,8)=(frxy-1)*Nxy+Nx+frx;
%
%
%
%
%         if (mod(frx+1,Nx)==0 &&mod(frx+Nx+1,Nxy)>0 )
%             frx=frx+2;
%         elseif (mod(frx+Nx+1,Nxy)==0 )
%
%             frxy=frxy+1;
%             frx=1;
%
%         else frx=frx+1;
%         end
%     end

if(mode>1)
    
    fid = fopen('C:\\Users\\Hassan\\\Documents\MATLAB\\myMesh.txt','wt');  % Note the 'wt' for writing in text mode
    
    
    
    fprintf(fid,'hexahedron\n');
    fprintf(fid,'//numberOfNodes\n');
    fprintf(fid,'%d\n',Nn);
    fprintf(fid,'//numberOfElements\n');
    fprintf(fid,'%d\n',Ne);
    fprintf(fid,'//numberOfRegs\n');
    fprintf(fid,'%d\n',1);
    fprintf(fid,'//factor\n');
    fprintf(fid,'%f\n',factor);
    
    for ie=1:Ne
        
        for k=1:7
            fprintf(fid,'%d,', elVert(ie,k));
        end
        fprintf(fid,'%d\n', elVert(ie,8));
    end
    
    
    
    
    y=linspace(0,Na,ny+1);
    
    %          x=linspace(0,.001,nx+1);
    %
    %          z=linspace(0,.001,nz+1);
    dx=.001;
    dy=y(2)-y(1);
    dz=.001;
    ddx=1.0/nx;
    D=(nx-1)*ddx/2;
    
    for kz=1:nz
        for ky=1:ny
            for kx=1:nx
                fprintf(fid,'%f,%f,%f\n',(kx-1)*ddx-D,y(ky),dz);
                fprintf(fid,'%f,%f,%f\n',dx+(kx-1)*ddx-D,y(ky),dz);
                fprintf(fid,'%f,%f,%f\n',dx+(kx-1)*ddx-D,y(ky)+dy,dz);
                fprintf(fid,'%f,%f,%f\n',(kx-1)*ddx-D,y(ky)+dy,dz);
                
                fprintf(fid,'%f,%f,%f\n',(kx-1)*ddx-D,y(ky),0);
                fprintf(fid,'%f,%f,%f\n',dx+(kx-1)*ddx-D,y(ky),0);
                fprintf(fid,'%f,%f,%f\n',dx+(kx-1)*ddx-D,y(ky)+dy,0);
                fprintf(fid,'%f,%f,%f\n',(kx-1)*ddx-D,y(ky)+dy,0);
            end
        end
    end
    
    fprintf(fid,'%d,%d,xxx',1,Ne);
    
    
    fclose(fid);
end



nL=size(E,2);
fid = fopen('C:\\Users\\Hassan\\\Documents\MATLAB\\myFields.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'displacement\n');
fprintf(fid,'3\n');
fprintf(fid,'%d\n',Nn);

ie=0;

for ky=1:ny
    for kx=1:nx
        
        ie=ie+1;
        fprintf(fid,'%d\t%f\t%f\t%f\n',elVert(ie,1),E(kx,ky,1),E(kx,ky,2),E(kx,ky,3));
        
        fprintf(fid,'%d\t%f\t%f\t%f\n',elVert(ie,5),-E(kx,ky,1),-E(kx,ky,2),-E(kx,ky,3));
    end
end
fclose(fid);


end



