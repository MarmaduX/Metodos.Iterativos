%El presente archivo contiene el algoritmo desarrollado en GNU Octave para el calculo de la geodesica en el paraboloide hiperbolico (de tipo f(u,v)=u*v)
%Si se busca ejecutar el script recordar que debe encontrarse en la misma carpeta donde se ha guardado
%Es interesante realizar varias ejecuciones con distintos valores de a, b, M y el vector y0
%Notar que cuanto mayor es M, menor sera el error_porcentual


%Calculo de los valores de u(t) y v(t) a traves del algoritmo de Euler implementado para la parametrizacion del paraboloide hiperbolico
a=0; %Inicio del intervalo
b=10; %Fin del intervalo
M=10; %Cantidad de iteraciones
y0=[0.50148,0.50148,0.70293,0.70293]; %Valores iniciales de u, v, p, q respectivamente. Los valores por defecto cumplen que la norma de la derivada de la geodesica en este vector es 1
h=(b-a)/M; %Tamano del paso
Y=zeros(4,M);
Y(1,1)=y0(1);
Y(2,1)=y0(2);
Y(3,1)=y0(3);
Y(4,1)=y0(4);
for k=1:M-1 
  Y(1,k+1)=Y(1,k)+h*Y(3,k);
  Y(2,k+1)=Y(2,k)+h*Y(4,k);
  Y(3,k+1)=Y(3,k)+h*((-2*Y(2,k)*Y(3,k)*Y(4,k))/(Y(1,k)^2+Y(2,k)^2+1));
  Y(4,k+1)=Y(4,k)+h*((-2*Y(1,k)*Y(3,k)*Y(4,k))/(Y(1,k)^2+Y(2,k)^2+1));
endfor
u_euler=Y(1,1:M) %Valores de u segun Euler
v_euler=Y(2,1:M) %Valores de v segun Euler


%Calculo de los valores de u(t) y v(t) utilizando el comando "lsode" propio de GNU Octave para la resolucion de EDO's de primer orden
t=[a:h:b];
x0=[y0(1),y0(3),y0(2),y0(4)]; 
fun=@(x,t)[x(2);((-2*x(3)*x(2)*x(4))/(x(1)^2+x(3)^2+1));x(4);((-2*x(1)*x(2)*x(4))/(x(1)^2+x(3)^2+1))]; %x(1) es u, x(2) es p, x(3) es v y x(4) es q. Por esta razon cambia el orden de x0 respecto de y0 a pesar de usarse los mismos valores
Z=lsode(fun,x0,t); %Columna 1 es u, 2 es p, 3 es v y 4 es q. Cada fila es el paso h
u_octave=Z(1:M,1)' %Valores de u segun Octave
v_octave=Z(1:M,3)' %Valores de v segun Octave


%Valores de la parametrizacion
g_euler=zeros(M,3);
for i=1:M
g_euler(i,1)=u_euler(i);
g_euler(i,2)=v_euler(i);
g_euler(i,3)=u_euler(i).*v_euler(i);
endfor
g_octave=zeros(M,3);
for i=1:M
g_octave(i,1)=u_octave(i);
g_octave(i,2)=v_octave(i);
g_octave(i,3)=u_octave(i).*v_octave(i);
endfor
g_euler %Valor de la parametrizacion en cada paso utilizando u y v obtenidos por Euler (g_euler(i,j) : valor de la coordenada j de la parametrizacion en el paso i)
g_octave %Valor de la parametrizacion en cada paso utilizando u y v obtenidos por Octave (g_octave(i,j) : valor de la coordenada j de la parametrizacion en el paso i)


%Calculo de errores
g_error=(g_euler-g_octave); 
error=zeros(1,M);
for i=1:M
error(1,i)=norm(g_error(i,1:3),2); %Es la norma de la diferencia de los vectores de la parametrizacion en cada paso utilizando los valores obtenidos por Euler y Octave 
endfor
norma_octave=zeros(1,M);
for i=1:M
norma_octave(1,i)=norm(g_octave(i,1:3),2); %Norma de cada vector de los valores de octave en cada paso
endfor
error_porcentual=(error./norma_octave) %Error porcentual. Es el mostrado en pantalla 
p=mean(error_porcentual)

%Grafica parabolide hiperbolico y su curva geodesica obtenida a traves de Euler
%L=[-10:1:10];
%K=L;
%[Mx,My]=meshgrid(L,K);
%Mz=Mx.*My;
%mesh(Mx,My,Mz)
%hold
plot3(u_euler,v_euler,u_euler.*v_euler)
hold
plot3(u_octave,v_octave,u_octave.*v_octave)
