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