%% **************************** part 1 ************************* %%

% ************************** 1 - givens *************************** %

% A matrix
A=[0.7 0.5 0;
   -0.5 0.7 0;
   0 0 0.9];

% B vector
B=[1;1;1];

% C vector
C=[0 -1 1];


% N
N=0:19;

% y ref
y_ref=[1 0 0 4 4 1 0 2 0 1 3 4 4 2 1 2 4 3 2 2]';

% xo
xo=[0.1;0.2;0.3];

% ******************** 2 - Calculate Q *************************** %

% let's compute the Qk and put them in a vector
% Q_vector= [Qo Q1 Q2 .... Qn-1]
Q_vector=zeros(1,20);

% Qo = 0.5
Q_vector(1,1)=0.5;


for i=2:20
    Q_vector(1,i)= C* A^(i-2) * B;
end


% let's estaplish our Q matrix

Q=zeros(20,20);
for ii=1:20
    L=length(Q_vector)-(ii-1);
    vec= Q_vector(1,ii)*ones(1,L);
    
    Q_k= diag(vec,-(ii-1));
    Q= Q+Q_k;
end

% **************** 3 - calculate phi matrix & xk ******************* %

phi=zeros(20,3);
xk=zeros(3,20);

for i=1:20
   phi(i,:)= C * A^(i-1) ;
    xk(:,i)=  A^(i-1)*xo ;
end


% **************** 4 - calculate u_opt ******************* %

u_opt = (Q' * Q) \ Q' * (y_ref - phi*xo);

% here we limit our control output from -100 to 100
u_limit=u_opt;
u_limit(u_limit>100)=100;
u_limit(u_limit<-100)=-100;

% **************** 5 - calculate yk & yk_limit ******************* %

yk= phi * xo + Q * u_opt;
yk_limit= phi * xo + Q * u_limit;


% **************** 6 - plot ****************************** %

% plot yk and y_ref on same graph
figure;hold on;
plot(N,yk,'b');
plot(N,y_ref,'g');

% plot yk_limit and y_ref 
figure;hold on;
plot(N,yk_limit);
plot(N,y_ref,'r')

%plot the control outputs u_optimal and u_limit
figure;hold on;
plot(N,u_opt,'r');
plot(N,u_limit ,'b')






