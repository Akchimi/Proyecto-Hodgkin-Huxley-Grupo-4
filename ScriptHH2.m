% Proyecto: Simulación del módelo de Hodgkin Huxley 
%
% Valeria Bonilla Rojas
% Johnny Borbón Valverde
% Angello Crawford Clark
% Daniel Espinoza Castro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parametros de la membrana
%
% Conductancias maximas (ms/cm^2) (1=K, 2=Na, 3=L)
g(1)=36; 
g(2)=120; 
g(3)=0.3;
% Area (cm2) usando el radio (a) mencionado en el paper de H-H
A = (0.0238)^2*pi;
% Corriente maxima (mA)
Imax = 0.1;
% Potenciales de equilibrio (mV) E:(1=K, 2=Na, 3=L)
E(1)=-12; 
E(2)=115; 
E(3)=10.613;
% Capacitancia (uF/cm2)
C=1;
% Condiciones iniciales: V0 (mV), x:(1->n0, 2->m0, 3->h0)
I_ext=0.0;
V0=-10;
x=zeros(1,3); 
x(1) = 0.318; 
x(2) = 0.053; 
x(3)=0.6;
t_i=0;
% Paso para realizar la integración
dt=0.01;
% Se establece la duracion total de la simulación (y por tanto, el eje x del gráfico)
t_f = 200;
% Implementacion del metodo de Euler para integrar el sistema
V = V0;
for t=0:dt:t_f
% Estimulo inicia en t=0ms y termina en 50ms, cuando Imax es diferente de 0
if t==0; 
    I_ext=Imax/A; 
end
if t==100; 
    I_ext=0.0; 
end
% Calculo de las funciones alfa
alfa(1)=(10-V)/(100*(exp((10-V)/10)-1));
alfa(2)=(25-V)/(10*(exp((25-V)/10)-1));
alfa(3)=0.07*exp((-V)/20);
% Calculo de las funciones beta
beta(1)=0.125*exp((-V)/80);
beta(2)=4*exp((-V)/18);
beta(3)=1/(exp((30-V)/10)+1);
% Actualización del vector de parametros n, m y h
x = x.*(1-dt.*(alfa+beta)) + alfa.*dt;
% Calculo de las conductancias actualizadas
gnmh(1)=g(1)*x(1)^4;
gnmh(2)=g(2)*x(2)^3*x(3);
gnmh(3)=g(3);
% Calculo de la suma de corrientes total
I=gnmh.*(V-E);
% Actualización del potencial de membrana con el método de Euler
V=V+dt/C*(I_ext-sum(I));
% Guardado de las variables para la gráfica
t_i=t_i+1;
x_plot(t_i)=t;
y_plot(t_i)=V;
end
% Gráfica
figure(1);hold off;plot(x_plot,y_plot,'red');
xlabel('Tiempo (ms)'); ylabel('Potencial de membrana (mV)');
title('Potencial de la membrana vs Tiempo');

