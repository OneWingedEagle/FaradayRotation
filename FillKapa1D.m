
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