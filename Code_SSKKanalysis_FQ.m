trs=dlmread('FQTS_0.5nm.dat');%----------------------------------------------------------------------------------read data from tranmission data file
d=size(trs,1);
l0=trs(1,1)-1;
A=zeros(d,1);en=zeros(d,1);n=zeros(d,1);l=zeros(d,1); kl=zeros(d,1); dn=zeros(d,d);s=zeros(d,1);%--------------Declaration of Matrices
% The Physical Quantities represented: en---Energy(eV) ; l----wavelength(nm) ; A----decadic Absorbance ; kl----Extinction coefficient ; n----refractive index

n1=1.4433;%known refractive index at wavelength of 800.99 nm
wl1=800.99;

ds=0.5;%-------------------------------------------------------------------------------------------------------Step size declaration

for i=1:d
l(i,1) = trs(i,1);
en(i,1)= 1239.9219/(l(i,1));
A(i,1) = -log(trs(i,2)/100);
kl(i,1)= (2.303*A(i,1)*l(i,1)*(1e-9))/(4*pi*1e-3);
end

for i=1:d
    for j=1:d
        if l0+ds*j==l(i,1)||l0+ds*j==wl1
            dn(i,j) = 0;
        else
            dn(i,j) = (kl(j,1)*(l0+j*ds))/(((l0+j*ds)^2-l(i,1)^2).*((l0+j*ds)^2-wl1^2));
        end
    end
end


s=sum(dn,2);

for i=1:d
n(i,1)= n1 +(2/pi)*(l(i,1)^2-wl1^2)*(l(d,1)-l(d-1,1))*((s(i,1)-dn(i,1)/2-dn(i,d)/2));
end


%figure; plot(l,n,'b'); axis square; grid on; xlabel('wavelength'); ylabel('refractive index');;ylim([0 3]);
figure; plot(en,n,'b'); axis square; grid on; xlabel('Energy (eV)'); ylabel('refractive index');;ylim([0 3]);
%figure; plot(l,kl,'r'); axis square; grid on; xlabel('wavelength'); ylabel('extinction coefficient');
figure; plot(en,kl,'r'); axis square; grid on; xlabel('Energy (eV)'); ylabel('extinction coefficient');

