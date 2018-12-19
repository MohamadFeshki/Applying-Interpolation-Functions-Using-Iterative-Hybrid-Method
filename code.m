%% mohamad feshki - 810195155

clc;
clear;

%% initialize:

SNR_1=zeros(10,5);
SNR_2=zeros(10,5);
SNR_3=zeros(10,5);


SNR_11=zeros(10,10);
SNR_22=zeros(10,10);
SNR_33=zeros(10,10);

% repeating loop for 10 random signals:

%% creating random signal and input

clc;
clear;



for c=1:10
    
Lopt=1; %94; %optimum landa based on article
N=2; %number of mixed harmonics
T=5; %sampling interval
Itr=5; % number of iteration
t=1:500; %sample of original signal
x1= rand(1,500);

X1=fft(x1);
X2=zeros(size(X1));
X2(1:15)=X1(1:15);
x2=abs(ifft(X2));

x=zeros(size(x1));
x(1:5:500)=x2(1:5:500);
x_org=x;

X3=fft(x2);
X4=zeros(size(X1));
X4(1:15)=X3(1:15);
x5=abs(ifft(X4));



s=1:N;



%% hybrid method:

for j=1:Itr
    
    
    if j==1
        x_k=x;
    else
        x_k=x_k_1;
    end
        
S1t=zeros(1,size(x,2)+T);
S2t=zeros(1,size(x,2)+T);


for n=1:(size(x,2)/T)
    rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
    rec(rec==0.5)=1;
    S1t= S1t+(rec.* (repmat ( x((n-1)*T+1),[ 1 size(x,2)+T ]))) ;
    S2t= S2t+(rec.* (repmat ( x_k((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

end

S1=S1t(1:size(x,2));
S2=S2t(1:size(x,2));


ss=repmat(s,[size(x,2) 1]);
x_m=(1+2*sum((cos(2*pi* ss' .*(t)/T)),1)).* S1;
X_m_f=fft(x_m);
X_m_f(60:500)=0;
x_m_f=abs(ifft(X_m_f));

x_mk=(1+2*sum((cos(2*pi* ss' .*(t)/T)),1)).* S2;
X_mk_f=fft(x_mk);
X_mk_f(60:500)=0;
x_mk_f=abs(ifft(X_mk_f));


X_k_f=fft(x_k);
X_k_f(60:500)=0;
x_k_f=abs(ifft(X_k_f));


x_k_1= (Lopt* x_m_f + x_k_f - x_mk_f);

SNR_1(c,j)=snr(x_k_1(1,50:450),x_k_1(1,50:450)-x5(1,50:450));

[Cor_1(c),Pval_1(c)]=corr(x_k_1(1,50:450)',x5(1,50:450)');


end

%% accelerated hybrid method:


A=abs((min(x(1:10:100))^2))/(sum((x(1:10:100)).^2));
B=140*(max(x)^2)/(sum((x).^2));

ro=(B-A)/(B+A);
K=2/(A+B);
L=0.8;

for j=2:Itr
    
  
    
 St=zeros(1,size(x,2)+T);

for n=1:(size(x,2)/T)
    rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
    rec(rec==0.5)=1;
    St= St+(rec.* (repmat ( x2((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

end
  
    S=St(1:size(x,2));
    Xt_f=fft(S);
    Xt_f(60:500)=0;
    xt_f=abs(ifft(Xt_f));
    

    
    x0=K*xt_f;
    if j==2
        x_n_1=x0;
        x_n_2=xt_f;
    else
        x_n_2=x_n_1;
        x_n_1=x_n;
    end
    
    St=zeros(1,size(x,2)+T);

for n=1:(size(x,2)/T)
    rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
    rec(rec==0.5)=1;
    St= St+(rec.* (repmat ( x_n_1((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

end

S=St(1:size(x,2));

X_n1_f=fft(S);
X_n1_f(60:500)=0;
x_n1_f=abs(ifft(X_n1_f));

x_n=( ( x0 + x_n1_f -K*x_n1_f -x_n_2 ) * L ) + x_n_2;

L=1/(1-((ro^2)/4)*L);

dif=mean(x_n)-mean(x5);
x5=x5+dif;

SNR_2(c,j)=snr(x_n(100:400),x5(100:400)-x_n(100:400));

[Cor_2(c),Pval_2(c)]=corr(x_n(1,50:450)',x5(1,50:450)');

end

%% interrative

L=0.8;

for j=1:Itr
        
     St=zeros(1,size(x,2)+T);

    for n=1:(size(x,2)/T)
        rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
        rec(rec==0.5)=1;
        St= St+(rec.* (repmat ( x2((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

    end
  
    S=St(1:size(x,2));
    Xt_f=fft(S);
    Xt_f(20:500)=0;
    x0=abs(ifft(Xt_f));
    

    
  if j==1
      xx_k=x0;
  else
      xx_k=xx_k_1;
  end   
      
        
     St=zeros(1,size(x,2)+T);

    for n=1:(size(x,2)/T)
        rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
        rec(rec==0.5)=1;
        St= St+(rec.* (repmat ( xx_k((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

    end
  
    S=St(1:size(x,2));
    Xt_f=fft(S);
    Xt_f(20:500)=0;
    xx_k_f=abs(ifft(Xt_f));
    
    
    xx_k_1=L*x0+0.2*xx_k-L*xx_k_f;
   
    dif=mean(xx_k_1)-mean(x5);
    x5=x5+dif;
   
SNR_3(c,j)=snr(xx_k_1(50:450),x5(50:450)-xx_k_1(50:450));
[Cor_3(c),Pval_3(c)]=corr(xx_k_1(1,50:450)',x5(1,50:450)');

    
end

ITr=10;

%% hybrid method:

for j=1:Itr
    
    noise= rand(1,500)/4;
   
    
    if j==1
        x_k=x_org+noise;
    else
        x_k=x_k_1+noise;
    end
        
S1t=zeros(1,size(x,2)+T);
S2t=zeros(1,size(x,2)+T);


for n=1:(size(x,2)/T)
    rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
    rec(rec==0.5)=1;
    S1t= S1t+(rec.* (repmat ( x((n-1)*T+1),[ 1 size(x,2)+T ]))) ;
    S2t= S2t+(rec.* (repmat ( x_k((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

end

S1=S1t(1:size(x,2));
S2=S2t(1:size(x,2));

 

ss=repmat(s,[size(x,2) 1]);
x_m=(1+2*sum((cos(2*pi* ss' .*(t)/T)),1)).* S1;
X_m_f=fft(x_m);
X_m_f(60:500)=0;
x_m_f=abs(ifft(X_m_f));

x_mk=(1+2*sum((cos(2*pi* ss' .*(t)/T)),1)).* S2;
X_mk_f=fft(x_mk);
X_mk_f(60:500)=0;
x_mk_f=abs(ifft(X_mk_f));


X_k_f=fft(x_k);
X_k_f(60:500)=0;
x_k_f=abs(ifft(X_k_f));


x_k_1= (Lopt* x_m_f + x_k_f - x_mk_f);

SNR_11(c,j)=snr(x_k_1(1,50:450),x_k_1(1,50:450)-x5(1,50:450));

[Cor_11(c),Pval_11(c)]=corr(x_k_1(1,50:450)',x5(1,50:450)');


end

%% accelerated hybrid method:


A=abs((min(x(1:10:100))^2))/(sum((x(1:10:100)).^2));
B=140*(max(x)^2)/(sum((x).^2));

ro=(B-A)/(B+A);
K=2/(A+B);
L=0.8;

for j=2:Itr
    
  noise= rand(1,500)/4;
  
  x2noise=x2+noise;
    
 St=zeros(1,size(x,2)+T);

for n=1:(size(x,2)/T)
    rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
    rec(rec==0.5)=1;
    St= St+(rec.* (repmat ( x2noise((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

end
  
    S=St(1:size(x,2));
    Xt_f=fft(S);
    Xt_f(60:500)=0;
    xt_f=abs(ifft(Xt_f));
    

    
    x0=K*xt_f;
    if j==2
        x_n_1=x0;
        x_n_2=xt_f;
    else
        x_n_2=x_n_1;
        x_n_1=x_n;
    end
    
    St=zeros(1,size(x,2)+T);

for n=1:(size(x,2)/T)
    rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
    rec(rec==0.5)=1;
    St= St+(rec.* (repmat ( x_n_1((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

end

S=St(1:size(x,2));

X_n1_f=fft(S);
X_n1_f(60:500)=0;
x_n1_f=abs(ifft(X_n1_f));

x_n=( ( x0 + x_n1_f -K*x_n1_f -x_n_2 ) * L ) + x_n_2;

L=1/(1-((ro^2)/4)*L);

dif=mean(x_n)-mean(x5);
x5=x5+dif;

SNR_22(c,j)=snr(x_n(100:400),x5(100:400)-x_n(100:400));
[Cor_22(c),Pval_22(c)]=corr(x_n(1,50:450)',x5(1,50:450)');

end

%% itrrative

noise= rand(1,500)/4;

L=0.8;

for j=1:Itr
        
     St=zeros(1,size(x,2)+T);

    for n=1:(size(x,2)/T)
        rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
        rec(rec==0.5)=1;
        St= St+(rec.* (repmat ( x2((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

    end
  
    S=St(1:size(x,2));
    Xt_f=fft(S);
    Xt_f(20:500)=0;
    x0=abs(ifft(Xt_f));
    

    
  if j==1
      xx_k=x0+noise;
  else
      xx_k=xx_k_1+noise;
  end   
      
        
     St=zeros(1,size(x,2)+T);

    for n=1:(size(x,2)/T)
        rec=rectangularPulse((n-1)*T+1+0.5-T/2,n*T+0.5-T/2,1:size(x,2)+T);
        rec(rec==0.5)=1;
        St= St+(rec.* (repmat ( xx_k((n-1)*T+1),[ 1 size(x,2)+T ]))) ;

    end
  
    S=St(1:size(x,2));
    Xt_f=fft(S);
    Xt_f(20:500)=0;
    xx_k_f=abs(ifft(Xt_f));
    
    
    xx_k_1=L*x0+0.2*xx_k-L*xx_k_f;
   
    dif=mean(xx_k_1)-mean(x5);
    x5=x5+dif;
   
SNR_33(c,j)=snr(xx_k_1(50:450),x5(50:450)-xx_k_1(50:450));
[Cor_33(c),Pval_33(c)]=corr(xx_k_1(1,50:450)',x5(1,50:450)');

    
end




end

M1=mean(SNR_1,1);
M2=mean(SNR_2,1);
M3=mean(SNR_3,1);
figure();
plot(M1); hold on; plot(M2); hold on; plot(M3)
xlabel('iteration')
ylabel('SNR (dB)')
legend('hybrid method', 'accelerated hybrid method', 'iterative method')


M11=mean(SNR_11,1);
M22=mean(SNR_22,1);
M33=mean(SNR_33,1);
figure();
plot(M11); hold on; plot(M22); hold on; plot(M33)
xlabel('iteration')
ylabel('SNR with noise (dB)')
legend('hybrid method', 'accelerated hybrid method', 'iterative method')



   
