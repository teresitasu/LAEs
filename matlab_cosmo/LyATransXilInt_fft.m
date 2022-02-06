% STILL TO WRITE
%
%************************************************************
%*  Yk = cdenLyATransXilInt_fft(fk,fkk,Pk_alpha_z,lorder)  *
%************************************************************
%************************************************************
% Computes integrand for correlation function.
%
% ARGUMENTS
%  fk                  Input wavenumber (fk in comoving h/ Mpc).
%  fkk                 Input wavenumber for FFT array (comoving h/ Mpc).
%  Pk_alpha_z     Flux power spectrum at specific redshift defined at lgfka.
%  lorder            Legendre order
%
% COMPATIBILITY: Matlab(?), Octave
%
% REQUIREMENTS:
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  15 05 18 Creation date.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Yk = LyATransXilInt_fft(fk,fkk,Pk_alpha_z,lorder);
Pk_alpha_z_arr = interp1(fk,Pk_alpha_z,fkk,'pchip','extrap');
n = length(fkk);
maskp = find(fkk > 0);
Yk = zeros(5,n);
lok = 0;
if(lorder==0)
    Yk(1,maskp) = Pk_alpha_z_arr(maskp);
    lok = 1;
end
if(lorder==2)
    Yk(1,maskp) = 3*Pk_alpha_z_arr(maskp)./ fkk(maskp).^2;
    Yk(2,maskp) = -Pk_alpha_z_arr(maskp);
    Yk(3,maskp) = -3*Pk_alpha_z_arr(maskp)./ fkk(maskp);
    lok = 1;
end
if(lorder==4)
    Yk(1,maskp) = 105*Pk_alpha_z_arr(maskp)./ fkk(maskp).^4;
    Yk(2,maskp) = -45*Pk_alpha_z_arr(maskp)./ fkk(maskp).^2;
    Yk(3,maskp) = Pk_alpha_z_arr(maskp);
    Yk(4,maskp) = -105*Pk_alpha_z_arr(maskp)./ fkk(maskp).^3;
    Yk(5,maskp) = 10*Pk_alpha_z_arr(maskp)./ fkk(maskp);
    lok = 1;
end
if(lok==0)
    disp('LyATransXilInt: lorder must be 0, 2 or 4');
    return;
end
for i=1:5
    Yk(i,maskp) = fkk(1,maskp).*Yk(i,maskp)/ (2*pi);
end
for j = n/2+2:n
    if(lorder==0)
        Yk(1,j) = -Yk(1,n+2-j);
        %if (j == n/2+2) | (j == n)
        %    fprintf("n %i %i %i %e\n",n,j,n+2-j,Yk(1,j));
        %end
    end
    if(lorder==2)
        Yk(1,j) = -Yk(1,n+2-j);
        Yk(2,j) = -Yk(2,n+2-j);
        Yk(3,j) = Yk(3,n+2-j);
    end
    if(lorder==4)
        Yk(1,j) = -Yk(1,n+2-j);
        Yk(2,j) = -Yk(2,n+2-j);
        Yk(3,j) = -Yk(3,n+2-j);
        Yk(4,j) = Yk(4,n+2-j);
        Yk(5,j) = Yk(5,n+2-j);
    end
end
