trs=dlmread('FQTS_0.5nm_extrapolation.dat',"\t");%-----------------------------------------------------------------------------------------------Read data from tranmission data file
d=size(trs,1);%------------------------------------------------------------------------------------------------------------------Declaration of size of transmission Matrix
l0=trs(1,1)-1;%------------------------------------------------------------------------------------------------------------------wavelength index declaration.
A=zeros(d,1);en=zeros(d,1);n=zeros(d,1); l=zeros(d,1); kl=zeros(d,1); dn=zeros(d,d);%--------------------------------------------Declaration of Matrices with zero initial values
% The Physical Quantities represented: en---Energy(eV) ; l----wavelength(nm) ; A----decadic Absorbance ; kl----Extinction coefficient ; n----refractive index
ds= 0.5;%------------------------------------------------------------------------------------------------------------------------Step size declaration
kinf=8.8e-06;
for i=1:d
l(i,1) = trs(i,1);
en(i,1)= 1239.9219/(l(i,1));
A(i,1) = -2.303*log(trs(i,2)/100);
kl(i,1)= (A(i,1)*l(i,1)*(1e-9))/(4*pi*1e-3);
end
a=l(1,1);b=l(d,1);
for i=1:d
    for j=1:d
        if l(i,1)==l0 + ds*j
            dn(i,j) = 0;
        else
            dn(i,j) = 1e9*(kl(j,1)/(l0+j*ds))*(l(i,1)^2/(l(i,1)^2-(l0+j*ds)^2));%----------------------------------------------------Declaration of Integrand Matrix
        end
    end
end
s=sum(dn,2);%--------------------------------------------------------------------------------------------------------------------Summing the all elements in each row

for i=1:d
n(i,1)=1+(2/pi)*((l(d,1)-l(d-1,1))*(1e-9)*(s(i,1)-dn(i,1)/2-dn(i,d)/2)) + (l(i,1)*(1e-9)*kinf)/(pi)*log(abs(b^2)/(b^2-l(i,1)^2));%------------Integrating with Trapezoidal rule
end

figure; plot(l,n,'b'); axis square; grid on; xlabel('wavelength (nm)'); ylabel('refractive index');legend('Fused Quartz 0.5 nm resolution');ylim([0 3]);
%figure; plot(en,n,'b'); axis square; grid on; xlabel('Energy (eV)'); ylabel('refractive index');legend('Fused Quartz 0.5 nm resolution');ylim([0 3]);
figure; plot(l,kl,'r'); axis square; grid on; xlabel('wavelength'); ylabel('extinction coefficient');legend('Fused Quartz 0.5 nm resolution');
%figure; plot(en,kl,'r'); axis square; grid on; xlabel('Energy (eV)'); ylabel('extinction coefficient');legend('Fused Quartz 0.5 nm resolution');
