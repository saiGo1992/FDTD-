clear all;
clc;
%pi/5取样时二阶精度和四阶精度的一阶导数
N_a=31;
x_a = linspace (0 ,6 * pi , N_a);
f_a = sin ( x_a ).* exp( -0.3* x_a );
f_derivative_a = cos( x_a ).* exp ( -0.3* x_a)-0.3* sin ( x_a ).* exp (-0.3* x_a );
dx_a = pi /5;
f_derivative_forward_a = zeros (1 ,N_a);
f_derivative_backward_a = zeros (1 ,N_a);
f_derivative_central_a = zeros (1 ,N_a);

f_derivative_forward_a(1:N_a -2) = (4*f_a(2:N_a-1)-f_a(3:N_a)-3*f_a(1:N_a-2))/(2*dx_a);
f_derivative_backward_a(3:N_a) = (3*f_a(3:N_a)-4*f_a(2:N_a-1)+f_a(1:N_a-2))/(2*dx_a);
f_derivative_central_a(3:N_a-2) = (-f_a(5:N_a)+8*f_a(4:N_a-1)-8*f_a(2:N_a-3)+f_a(1,N_a-4))/(12*dx_a);

figure(1);
plot( x_a , f_derivative_a , 'k' ,...
      x_a(1:N_a-2),f_derivative_forward_a(1:N_a-2), 'b--' , ...
      x_a(3:N_a),f_derivative_backward_a(3:N_a), 'r-.', ...
      x_a(3:N_a-2), f_derivative_central_a(3:N_a-2),': ms', ...
      'MarkerSize',4 ,'linewidth',1.5);
set(gca ,  'FontSize' ,12 , 'fontweight' , 'demi') ;
axis([0 6* pi -1 1]);
grid on ;
legend ('exact',' forward difference' ,...
        'backward difference' ,' central difference') ;
% xlabel ( ' $x$ ' , ' Interpreter ' , ' latex ' , ' FontSize ' ,16);
% ylabel ( '$f ', ' ( x )$ ' , ' Interpreter ' , ' latex ' , ' FontSize ' ,16);
%  text ( pi ,0.6 , ' $\ Delta x = \ pi /5$ ' , ' Interpreter ' , ...
%         'latex' , ' fontsize ' ,16 , ' BackgroundColor ' , 'w ' , ' EdgeColor ' , 'k ' ) ;

%pi/5取样时误差
% plot error for finite difference derivatives
% using pi /5 sampling period
error_forward_a =f_derivative_a - f_derivative_forward_a ;
error_backward_a=f_derivative_a - f_derivative_backward_a ;
error_central_a =f_derivative_a - f_derivative_central_a ;

figure(2);
plot(x_a(1:N_a-2), error_forward_a(1:N_a-2), 'b--' , ...
     x_a(3:N_a) , error_backward_a(3:N_a), 'r-.' , ...
     x_a(3:N_a-2), error_central_a (3:N_a-2), ': ms' , ...
     'MarkerSize' ,4 , 'linewidth' ,1.5);
set(gca ,'FontSize' ,12 , 'fontweight' , 'demi' ) ;
axis([0 6*pi -0.2 0.2 ]) ;
grid on ;
legend('forward difference' ,  'backward difference' , ...
       'central difference'); 
