function x=topgeo(c1,c2,lar,lor,opz)

% Conversione da coordinate geografiche a coordinate topografiche (origine in S) e
% viceversa
% se opz=0 o è omesso:  c1,c2 = latitudine e longitudine dei punti da trasformare
%                       lar, lor = latitudine e longitudine dell'origine
% se opz =/ 0           c1,c2 = coordinate topografiche cartesiane (x, y) y=nord
%                       lar, lor = latitudine e longitudine dell'origine

Rt=6371;

if nargin < 5, opz=0; end

if opz == 0
   x=[Rt*tan(c2-lor).*cos(c1)./cos(c1-lar),Rt*tan(c1-lar)];
else
   x1=lar+atan(c2/Rt);
   y1=lor+atan(c1.*cos(x1-lar)./cos(x1)/Rt);
   x=[x1,y1];
end

return