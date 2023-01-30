%% Bifurcation for ELF3 network
clear all
[x,v,s,h,f] =elf3_bifur; 
a = x(9,:); %bifurcation parameter
b = x(2,:); 
c = a./1000;

%% Based on eigenvalues to judge stable vs. unstable states
ind = zeros(1,4);
snum = size(f);
num = snum(2);
j = 1;
n = 1;

for n = 1:1:(num-1)
    x11 = find(f(:,n) > 0);
    x12 = find(f(:,n+1) > 0);
    if isempty(x11) && ~isempty(x12)
        ind(j) = n + 1;
        j = j + 1;
    elseif ~isempty(x11) && isempty(x12)
        ind(j) = n + 1;
        j = j + 1;
    end
end

%%

amat=[c(ind(1)) c(ind(2)) c(ind(7)) c(ind(end))];
figure1 = figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
% % 
% % ax2 = subplot(2,2,2);
plot(c(1:ind(1)),b(1:ind(1)),'b');
hold on
plot(c(ind(1)+1:ind(2)),b(ind(1)+1:ind(2)),'r');
plot(c(ind(2)+1:ind(3)),b(ind(2)+1:ind(3)),'b'); 
plot(c(ind(3)+1:ind(end)),b(ind(3)+1:ind(end)),'r');
plot(c(ind(end)+1:end),b(ind(end)+1:end),'b');
% 
% 
xlim([0 200]);
xlabel('I ext (10^3 molecules)');
ylabel('ZEB mRNA');
sound(sin(1:3000));
% Bifurcation for core network

[x1,v1,s1,h1,f1] =core_bifur; 
a1 = x1(8,:); %bifurcation parameter
b1 = x1(2,:); 
c1 = a1./1000;

% Based on eigenvalues to judge stable vs. unstable states
ind1 = zeros(1,4);
snum1 = size(f1);
num1 = snum1(2);
j = 1;
n = 1;

for n = 1:1:(num1-1)
    x21 = find(f1(:,n) > 0);
    x22 = find(f1(:,n+1) > 0);
    if isempty(x21) && ~isempty(x22)
        ind(j) = n + 1;
        j = j + 1;
    elseif ~isempty(x21) && isempty(x22)
        ind(j) = n + 1;
        j = j + 1;
    end
end
 
%

amat1=[c1(ind(1)) c1(ind(2)) c1(ind(5)) c1(ind(end))];

plot(c1(1:ind(1)),b1(1:ind(1)),'g');
hold on
plot(c1(ind(1)+1:ind(2)),b1(ind(1)+1:ind(2)),'k');
plot(c1(ind(2)+1:ind(3)),b1(ind(2)+1:ind(3)),'g'); 
plot(c1(ind(3)+1:ind(end)),b1(ind(3)+1:ind(end)),'k');
plot(c1(ind(end)+1:end),b1(ind(end)+1:end),'g');
% 
sound(sin(1:3000));
