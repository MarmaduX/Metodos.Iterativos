%Se realizan los calculos del PVI del Paraboloide Hiperbolico mediante distintos metodos y los valores obtenidos se contrastan con los valores obtenidos al utilizar el comando "lsode" propio de Octave.
%El usuario puede modificar cualquiera de los valores iniciales.


function principal
%Valores iniciales
a=0; %Inicio del intervalo
b=10; %Fin del intervalo
M=10; %Cantidad de iteraciones
y0=[0.50148,0.50148,0.70293,0.70293]; %Vector de valores iniciales
tol=0.001; %Tolerancia utiliada para controlar las iteraciones dentro de los metodos de euler hacia atras y el trapecio

%Calculo de los valores de u(t) y v(t) utilizando el comando "lsode" propio de GNU Octave para la resolucion de EDO's de primer orden
t=[a:((b-a)/M):b];
x0=[y0(1),y0(3),y0(2),y0(4)]; 
fun=@(x,t)[x(2);((-2*x(3)*x(2)*x(4))/(x(1)^2+x(3)^2+1));x(4);((-2*x(1)*x(2)*x(4))/(x(1)^2+x(3)^2+1))]; %x(1) es u, x(2) es p, x(3) es v y x(4) es q. Por esta razon cambia el orden de x0 respecto de y0 a pesar de usarse los mismos valores
Z=lsode(fun,x0,t); %Columna 1 es u, 2 es p, 3 es v y 4 es q. Cada fila es el paso h
u_octave=Z(1:M,1)'; %Valores de u segun Octave
v_octave=Z(1:M,3)'; %Valores de v segun Octave

%Calculo de los valores de u(t) y v(t) utilizando los metodos clasicos de resolucion de EDO's implementados
[u_euler_adelante,v_euler_adelante]=EulerAdelante(a,b,M,y0);
[u_euler_atras,v_euler_atras]=EulerAtras(a,b,M,y0,tol);
[u_punto_medio,v_punto_medio]=PM(a,b,M,y0);
[u_heun,v_heun]=Heun(a,b,M,y0);
[u_trapecio,v_trapecio]=Trapecio(a,b,M,y0,tol);
[u_rungekutta,v_rungekutta]=RK4(a,b,M,y0);

%Calculo de los valores de la parametrizacion segun cada metodos
g_octave=valor_parametrizacion(u_octave,v_octave,M);
g_euler_adelante=valor_parametrizacion(u_euler_adelante,v_euler_adelante,M);
g_euler_atras=valor_parametrizacion(u_euler_atras,v_euler_atras,M);
g_punto_medio=valor_parametrizacion(u_punto_medio,v_punto_medio,M);
g_heun=valor_parametrizacion(u_heun,v_heun,M);
g_trapecio=valor_parametrizacion(u_trapecio,v_trapecio,M);
g_rungekutta=valor_parametrizacion(u_rungekutta,v_rungekutta,M);

%Calculo de errores pocentuales
ge_euler_adelante=g_euler_adelante-g_octave;
ge_euler_atras=g_euler_atras-g_octave;
ge_punto_medio=g_punto_medio-g_octave;
ge_heun=g_heun-g_octave;
ge_trapecio=g_trapecio-g_octave;
ge_rungekutta=g_rungekutta-g_octave;
norma_octave=norma(g_octave,M);
norma_euler_adelante=norma(ge_euler_adelante,M);
norma_euler_atras=norma(ge_euler_atras,M);
norma_punto_medio=norma(ge_punto_medio,M);
norma_heun=norma(ge_heun,M);
norma_trapecio=norma(ge_trapecio,M);
norma_rungekutta=norma(ge_rungekutta,M);
%Valores mostrados en pantalla
error_porcentual_euler_adelante=mean((norma_euler_adelante./norma_octave))
error_porcentual_euler_atras=mean((norma_euler_atras./norma_octave))
error_porcentual_punto_medio=mean((norma_punto_medio./norma_octave))
error_porcentual_heun=mean((norma_heun./norma_octave))
error_porcentual_trapecio=mean((norma_trapecio./norma_octave))
error_porcentual_rungekutta=mean((norma_rungekutta./norma_octave))
endfunction

%Funcion que devuelve una matriz donde cada fila es el valor (vector de 3 entradas) de la parametrizacion en un instante de tiempo
function G = valor_parametrizacion(u,v,M)
  G=zeros(M,3);
  for i=1:M 
    G(i,1)=u(i);
    G(i,2)=v(i);
    G(i,3)=u(i).*v(i);
  endfor
endfunction
%Funcion que calcula la norma de cada vector (fila) de una matriz
function G=norma(H,M)
  G=zeros(1,M);
  for i=1:M
    G(1,i)=norm(H(i,1:3),2);  
  endfor
endfunction

%Implementacion de Euler Hacia Adelante
function [u,v]= EulerAdelante(a,b,M,y0)
	h=(b-a)/M; 
	Y=zeros(4,M);
	Y(1:4,1)= y0(1,1:4);
	for k=1:M-1
		Y(1:4,k+1)=Y(1:4,k)+h.*(f(Y(1,k),Y(2,k),Y(3,k),Y(4,k)));
	endfor
	u=Y(1,1:M);
	v=Y(2,1:M);
end
%Implementacion de Euler Hacia Atras utilizando aproximacion por MIG de punto fijo
function [u,v]= EulerAtras(a,b,M,y0,tol)
	h=(b-a)/M; 
	Y=zeros(4,M);
	Y(1:4,1)= y0(1,1:4);
	for k=1:M-1
    P=zeros(4,3);
		P(1:4,1)=Y(1:4,k)+h.*(f(Y(1,k),Y(2,k),Y(3,k),Y(4,k))); 
    P(1:4,2)=Y(1:4,k)+h.*(f(P(1,1),P(2,1),P(3,1),P(4,1)));
    while(abs(P(1,2)-P(1,1))>=tol)
      P(1:4,1)=P(1:4,2);
      P(1:4,3)=Y(1:4,k)+h.*(f(P(1,1),P(2,1),P(3,1),P(4,1)));
      P(1:4,2)=P(1:4,3);
    endwhile
    Y(1:4,k+1)=P(1:4,2);
	endfor
	u=Y(1,1:M);
	v=Y(2,1:M);
end
%Implemetnacion de Punto Medio
function [u,v]= PM(a,b,M,y0)
  h=(b-a)/M; 
	Y=zeros(4,M);
	Y(1:4,1)= y0(1,1:4);
  Y(1:4,2)= Y(1:4,1)+h.*(f(Y(1,1),Y(2,1),Y(3,1),Y(4,1)));
  for k=3:M
		Y(1:4,k)=Y(1:4,k-2)+2*h.*(f(Y(1,k-1),Y(2,k-1),Y(3,k-1),Y(4,k-1)));
	endfor
  u=Y(1,1:M);
	v=Y(2,1:M);
end 
%Implementacion de Heun
function [u,v]= Heun(a,b,M,y0)
  h=(b-a)/M; 
  Y=zeros(4,M);
  P=zeros(4,1);
  Y(1:4,1)= y0(1,1:4);
  for k=1:M-1
    P(1:4,1)=Y(1:4,k)+h.*(f(Y(1,k),Y(2,k),Y(3,k),Y(4,k)));
	  Y(1:4,k+1)=Y(1:4,k)+(h/2).*(f(Y(1,k),Y(2,k),Y(3,k),Y(4,k))+f(P(1,1),P(2,1),P(3,1),P(4,1)));
  endfor
  u=Y(1,1:M);
  v=Y(2,1:M);
end
%Implementacion de Trapecio utilizando aproximacion por MIG de punto fijo
function [u,v]= Trapecio(a,b,M,y0,tol)
  tol=0.001;
	h=(b-a)/M; 
	Y=zeros(4,M);
	Y(1:4,1)= y0(1,1:4);
	for k=1:M-1
    P=zeros(4,3);
		P(1:4,1)=Y(1:4,k)+h.*(f(Y(1,k),Y(2,k),Y(3,k),Y(4,k))); 
    P(1:4,2)=(h/2).*f(P(1,1),P(2,1),P(3,1),P(4,1))+(Y(1:4,k)+(h/2).*(f(Y(1,k),Y(2,k),Y(3,k),Y(4,k)))); %El ultimo termino esta "fijo"
    while(abs(P(1,2)-P(1,1))>=tol)
      P(1:4,1)=P(1:4,2);
      P(1:4,3)=(h/2).*f(P(1,1),P(2,1),P(3,1),P(4,1))+(Y(1:4,k)+(h/2).*(f(Y(1,k),Y(2,k),Y(3,k),Y(4,k)))); %El ultimo termino esta "fijo"
      P(1:4,2)=P(1:4,3);
    endwhile
    Y(1:4,k+1)=P(1:4,2);
	endfor
	u=Y(1,1:M);
	v=Y(2,1:M);
end
%Implementacion de Runge Kutta de cuarto orden
function [u,v]= RK4(a,b,M,y0)
  h=(b-a)/M; 
  Y=zeros(4,M);
  F1=zeros(4,1);
  F2=zeros(4,1);
  F3=zeros(4,1);
  F4=zeros(4,1);
  Y(1:4,1)= y0(1,1:4);
  for k=1:M-1
    F1(1:4,1)=f(Y(1,k),Y(2,k),Y(3,k),Y(4,k));
    F2(1:4,1)=f(Y(1,k)+(h*0.5)*F1(1,1),Y(2,k)+(h*0.5)*F1(2,1),Y(3,k)+(h*0.5)*F1(3,1),Y(4,k)+(h*0.5)*F1(4,1));
    F3(1:4,1)=f(Y(1,k)+(h*0.5)*F2(1,1),Y(2,k)+(h*0.5)*F2(2,1),Y(3,k)+(h*0.5)*F2(3,1),Y(4,k)+(h*0.5)*F2(4,1));
    F4(1:4,1)=f(Y(1,k)+h*F3(1,1),Y(2,k)+h*F3(2,1),Y(3,k)+h*F3(3,1),Y(4,k)+h*F3(4,1));
	  Y(1:4,k+1)=Y(1:4,k)+(h/6).*(F1(1:4,1)+2.*(F2(1:4,1)+F3(1:4,1))+F4(1:4,1));
  endfor
  u=Y(1,1:M);
  v=Y(2,1:M);
end 
%Definicion de la funcion f(t,y(t))
function Z = f(w,x,y,z)
  Z=zeros(4,1);
	Z(1,1)=y;
	Z(2,1)=z;
	Z(3,1)=(-2*x*y*z)/(w^2+x^2+1);
	Z(4,1)=(-2*w*y*z)/(w^2+x^2+1);
end
