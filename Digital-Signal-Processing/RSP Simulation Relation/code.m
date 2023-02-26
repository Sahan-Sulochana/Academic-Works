%%%%%%%%%%%%%%===Defining Variables===%%%%%%%%%%%%%%%%%%%%%%%%

L = 100000; %Binary sequence of length
D = zeros(1,L); %10^3 or 10^5 zeros sequence
p = randperm(L,L/2); %Selecting the placement indexes for ones(1-1000)
D(p) = ones(1,L/2); %Replaces the zeros to ones reguarding to above placement indexes
A = 1; %Amplitude

%%%%%%%%%%%%%%===Generating the Transmitting signals===%%%%%%%%%%%%%%%%%%%%%%%%

S = zeros(1,L); %10^3 or 10^5 zeros sequence
for i = 1:L
    if D(i) == 1 %if D = 0 amplitude set => A
        S(i) = A;
    else
        S(i) = -1*A; %if D = 1 amplitude set => -A
    end
end

%%%%%%%%%%%%%%===Generating the AWGN===%%%%%%%%%%%%%%%%%%%%%%%%
M = 0;
sigma = 1;
N = M + sigma*randn(1,L);

%%%%%%%============Question 6 interference part======%%%%%%%%%% 
%generating intterference
MI = 0;
sigmaI = 1;
I = MI + sigmaI*randn(1,L);
%%%%%%%%==============================================%%%%%%%%%%%

alpha = 10;
%%%%%%%%%%%%%%===Generating the Recived signals===%%%%%%%%%%%%%%%%%%%%%%%%
%R = S + N; %for Question 5
%R = S + N + I; %for Question 6
R = alpha*S + N; %for Question 7
%%%%%%%%%===========================================%%%%%%%%%%
figure;
stairs([1:L],R);
title("Received Signal");




%%%%%%%%%%%%%%===Decoded signals===%%%%%%%%%%%%%%%%%%%%%%%%

tau = 0;
Y = zeros(1,L);
for i = 1:L
    if R(i) > tau
        Y(i) = A;
    else
        Y(i) = -1*A;
    end
end

%%==plotting transmitted signal==%
figure;
subplot(2,1,1);
stairs([1:L],S);
title("Transmitted Signal");
xlim([0 L]);
ylim([-1*A-1 A+1]);

%%==plotting Decoded signal==%
subplot(2,1,2);
stairs([1:L],Y);
title("Decoded Signal");
xlim([0 L]);
ylim([-1*A-1 A+1]);


%%%%%%%%%%%%%%===Bin sequence===%%%%%%%%%%%%%%%%%%%%%%%%

bin0 = 100; %Number of bins
Rmin = min(R); %minimum value of list for width calculation
Rmax = max(R); %maximum value of list for width calculation
width = (Rmax-Rmin)/(bin0-1); %width calculation
bins = [Rmin-width/2:width:Rmax]; %Histogram value list

y = zeros(1,bin0);% Bins y values
for i = 1:L
    for j = 1:bin0
        if (R(i) >= bins(j)-width/2) && (R(i) < bins(j)+width/2)
            y(j) = y(j) + 1;
        end
    end
end

new = y/width;

%histogram plot
figure;
bar(bins,new);
title("Histogram of R");

%using the MATLAB buit in function hist()
figure;
hist(R,bin0);
title("Histogram of R (Using built-in function)");


%%%%%%%%%==========plotting the PDF of f_R|S(r|S=A)=========%%%%%%%%

list1 = []; %creating a list for assign values when S=A
index = 1;

for i = 1:L
    if S(i) == A
        list1(index) = R(i);
        index = index + 1;
    end
end

bin = 100;
Rmax1 = max(list1);%maximum value of list for width calculation
Rmin1 = min(list1);%minimum value of list for width calculation
width1 = (Rmax1-Rmin1)/(bin-1); %width calculation
bins1 = [Rmin1-width1/2:width1:Rmax1]; %Histogram value list
[y1,x1] = hist(list1,bins1); %plotting the histogram
y1 = y1/((index-1)*width1);
figure;
bar(x1,y1);
hold on;
plot(x1,y1,'r'); %plotting the pdf
title("f_{R|S}(r|S=A)");



%%%%%%%%%==========plotting the PDF of f_R|S(r|S=-A)=========%%%%%%%%

list2 = []; %creating a list for assign values when S=-A
index = 1;

for i = 1:L
    if S(i) == -1*A
        list2(index) = R(i);
        index = index + 1;
    end
end


Rmax2 = max(list2);%maximum value of list for width calculation
Rmin2 = min(list2);%minimum value of list for width calculation
width2 = (Rmax2-Rmin2)/(bin-1); %width calculation
bins2 = [Rmin2-width2/2:width2:Rmax2]; %Histogram value list
[y2,x2] = hist(list2,bins2); %plotting the histogram
y2 = y2/((index-1)*width2);
figure;
bar(x2,y2);
hold on;
plot(x2,y2,'r'); %plotting the PDF
title("f_{R|S}(r|S=-A)");



%%%%%%%%%==========plotting the PDF of f_R|S(r|S=-A)=========%%%%%%%%

%E[R|S=A] calculation
ERlSA = 0;
for i = 1:bin
    ERlSA = ERlSA + (x1(i)*y1(i)*width1); %mean function
end
fprintf("E[R|S=A] = %f\n",ERlSA);

%E[R|S=-A] calculation
ERlS_A = 0;
for i = 1:bin
    ERlS_A = ERlS_A + (x2(i)*y2(i)*width2); %mean function
end
fprintf("E[R|S=-A] = %f\n",ERlS_A);

%E[R] calculation
[y,x] = hist(R,bins);
y = y/(L*width);
ER = 0;
for i = 1:bin0
    ER = ER + (x(i)*y(i)*width); %mean function
end
fprintf("E[R] = %f\n",ER);

%plotting the PDF of f_R(r)
figure;
bar(x,y);
hold on;
plot(x,y,'r');
title("f_R(r)");

