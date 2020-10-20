%% 1- get our Data

%linear accelration x y z
%linear_acc = moh_Pos_sensorData(:,18:20);
linear_acc = moh_Pos_sensorData(:,18:20);
% GPS position x y z
gps_data= moh_Pos_xyz;                     

% get time step
dt= (moh_Pos_sensorData(2,1)-moh_Pos_sensorData(1,1)) * 1e-3 ;

%get intial position values from GPS data
xo=gps_data(1,1);
yo=gps_data(1,2);

% our time vector
time =moh_Pos_sensorData(:,1);

% get our gps data in two seprate vectors one for X and one for Y
X_gps=moh_Pos_xyz(:,1);
Y_gps=moh_Pos_xyz(:,2);

% get our linear acc data in two seprate vectors one for X and one for Y
X_acc=linear_acc(:,1);
Y_acc=linear_acc(:,2);




%% 2- Design our constant velocity model

%%%%%%%%%%%%%%%%%%%% for x data %%%%%%%%%%%%%%%%%%%%%%%

% Fit line to data using polyfit
c_x = polyfit(time,X_gps,1);
% Display evaluated equation y = m*x + b
disp(['Equation is x = ' num2str(c_x(1)) '*t + ' num2str(c_x(2))]);
x_model= c_x(1).*time + c_x(2);

%now let`s plot our fit data the slope will be our v_x
figure;plot(time,X_gps,time,x_model);
title('x with time');

%%%%%%%%%%%%%%%%%%%% for y data %%%%%%%%%%%%%%%%%%%%%%%

% Fit line to data using polyfit
c_y = polyfit(time,Y_gps,1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c_y(1)) '*t + ' num2str(c_y(2))]);

y_model= c_y(1).*time + c_y(2);
%now let`s plot our fit data the slope will be our v_x
figure;plot(time,Y_gps,time,y_model);
title('y with time');
%our velocity for model are
v_x = c_x(1);
v_y = c_y(1);


%transfer our X_acc and Y_acc into X_pos and Y_pos by integrating 2 times
X_acc=linear_acc(:,1);
Y_acc=linear_acc(:,2);


% create arrays to hold transformation value from accelration of linear acc to position
vel_x(1:length(X_acc))=zeros;
vel_y(1:length(X_acc))=zeros;
X_pos(1:length(X_acc))=zeros;
Y_pos(1:length(X_acc))=zeros;

%inital values for linear acc data in X
vel_x(1)=v_x + X_acc(1)*dt;
X_pos(1)=xo + vel_x(1)*dt;

%inital values for linear acc data in Y
vel_y(1)=v_y + Y_acc(1)*dt;
Y_pos(1)=yo + vel_y(1)*dt;

for i=2:length(X_acc)

    %inital values for linear acc data in X
    vel_x(i)=vel_x(i-1)+ X_acc(i-1)*dt;
    X_pos(i)=X_pos(i-1) + vel_x(i)*dt;

    %inital values for linear acc data in Y
    vel_y(i)=vel_y(i-1)+ Y_acc(i-1)*dt;
    Y_pos(i)=Y_pos(i-1) + vel_y(i)*dt;

end


%% 3- Design our matrcies Phi , Q , H , P , R
phi= [1 dt 0 0; 
      0 1 0 0;
      0 0 1 dt;
      0 0 0 1];

Q=   [1 0 0 0; 
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];
  
H=   [1 0 0 0; 
      1 0 0 0;
      0 0 1 0;
      0 0 1 0];

R=   [25 0 0 0; 
      0 250 0 0;
      0 0 25 0;
      0 0 0 250];
  
   
P=   [400 0 0 0; 
      0 400 0 0;
      0 0 400 0;
      0 0 0 400];
  
  
  
%% 4- now let's start our Kalman filter loop

%inital value for X vector
Xo=[xo ;v_x ;yo; v_y];
X(:,1)=Xo;

for i=1:length(Y_acc)
   Zk=[X_gps(i) ;X_pos(i); Y_gps(i); Y_pos(i)];
   
   %prediction step
   Xp= phi* X(:,i);
   Pp= phi * P * phi' + Q;
   
   %update step
   S= H * Pp * H' + R;
   K= Pp * H' * inv(S);
   Y= Zk - H * Xp;
   
   X(:,i+1)= Xp + K * Y;
   P= (eye(4) - K*H) * Pp * (eye(4) - K*H)' + K* R *K';
   
end

%% 5- now let's plot our graphs
figure
plot ( X_gps,Y_gps,'g')

hold on
plot(X_pos,Y_pos,'k')

plot( X(1,:),X(3,:),'b')
plot(x_model,y_model,'r');
legend('GPS','Linear Acc','Kalman Filter','Our Model')
hold off
