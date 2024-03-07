% Script para representar los datos obtenidos
% josem-canom
clc
clear all
close all
rvt=load('rvt.txt');
nc=(length(rvt(1,:))-1)/6;
hold on
for i=1:nc
plot3(rvt(:,6*i-5),rvt(:,6*i-4),rvt(:,6*i-3))
end
% for i=1:nc
% text(rvt(1,6*i-5),rvt(1,6*i-4),rvt(1,6*i-3),sprintf('%d',i))
% end
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
legend('�rbita del Sol','�rbita de la Luna','�rbita de la Tierra','�rbita de J�piter','�rbita de Saturno')
title('Trayectorias')
axis(10*1.496e11*[-1 1 -1 1])
grid
hold off