%% Generate initial DNA's
global DNAl
DNAl = 200;

global nGen
nGen = 12;

global DNAn
DNAn = 500;

DNAs = initDNAs;
disp('end of initDNAs');

Qdnas = zeros(DNAn,1);

%% Body
nEpoch = 50;
Vlist = zeros(100,1);
for epoch = 1:nEpoch
    disp(['=========== epoch:', num2str(epoch),' ===========']);
    %% Evolve
    for i=1:DNAn
%         if mod(i, 100) == 0
%             disp(['DNA num: ',num2str(i)]);
%         end
        p = (i-1) * 2 + 1;
        dna = DNAs(p:p+1,:);
        Qdnas(i) = evolve(dna, 0);
    end
    
%     disp(['end of evolution of epoch:', num2str(epoch)]);

    %% Update
    if epoch ~= nEpoch
        [Qdnas,Ps] = sort(Qdnas, 'descend');
        for i = 250:DNAn
            p = Ps(i);
            
            p = (p-1) * 2 + 1;
            DNAs(p:p+1,:) = genDNA(DNAs, Ps);
        end
    end
    
    Vlist(epoch) = Qdnas(1);
    %% Exit on convergence
%     if Vmx < Qdnas(1)
%         Vmx = Qdnas(1);
%     else
%         break;
%     end
    
    %% Visualization
    if mod(epoch, 5) == 0
        p = (Ps(1)-1) * 2 + 1;
        dna = DNAs(p:p+1,:);
        evolve(dna, 1);
    end

end

%% Result visualization
i = 1;
p = Ps(i);
p = (p-1) * 2 + 1;
dna = DNAs(p:p+1,:);
evolve(dna, 1);

%% Convergence
figure(101);
x=1:1:100;
plot(x, Vlist);
axis([1,99,0.89,1]);
xlabel('Epoch')
ylabel('Qvalue of the best DNA')
saveas(gcf,'convergence.png');


%% Color denotation
clrmp = zeros(13, 3);
clrmp = uint8(clrmp);
clrmp(1,:) = [135,206,250];
clrmp(2,:) = [95,158,160];  
clrmp(3,:) = [176,224,230]; 
clrmp(4,:) = [0,191,255]; 
clrmp(5,:) = [30,144,255]; 
clrmp(6,:) = [255,215,0];
clrmp(7,:) = [0,255,255]; 
clrmp(8,:) = [255,140,0];
clrmp(9,:) = [255,99,71]; 
clrmp(10,:) = [173,255,47];
clrmp(11,:) = [127,255,0]; 
clrmp(12,:) = [0,100,0];
clrmp(13,:) = [100,149,237];

img = zeros(256, 256, 3);
img = uint8(img);
for i=1:256
    for j=1:256
        img(i,j,:) = [255,255,255];
    end
end
for i = 1:13
    xs = 15*(i-1) + 30;
    for x = xs:xs+10
        for y = 50:70
            img(x,y,:) = clrmp(i,:);
        end
    end
end
figure(600);
imshow(img);
    
%% Utility functions
function DNAs = initDNAs
    global DNAn
    global DNAl
    DNAs = ones(DNAn * 2, DNAl);
    
    % obtain object category
    for i=1:DNAn
        p = 2 * (i-1) + 1;
        for j=1:DNAl
            obj = floor(rand() * 12) + 1;
            DNAs(p,j) = obj;
        end
        dna = DNAs(p,1:25);
        dna = sort(dna);
        DNAs(p,1:25) = dna;
    end
    
    % obtain pos's
    for i=1:DNAn
        p = i*2;
        for j=1:DNAl
            posRRate = 0.1;
            pos = 0;
            if rand() < posRRate
                pos = 1;
            end
            DNAs(p,j) = pos;
        end
    end    
    
    % Validation check for DNAs
    for i=1:DNAn
        p = 2 * (i-1) + 1;
        dna = DNAs(p:p+1, :);
        checkDNA(dna);
    end
end

function dna = genDNA(dnas, Ps)
    global DNAl
    sppt = floor(rand() * (DNAl-1)) + 1; % split point
    if sppt <= 1
        sppt = sppt + 1;
    elseif sppt > DNAl-2
        sppt = sppt - 1;
    end
    
    dp1 = floor(rand() * 250)+1;
    dp2 = floor(rand() * 250)+1;
    dp1 = Ps(dp1);
    dp1 = 2*(dp1-1)+1;
    dp2 = Ps(dp2);
    dp2 = 2*(dp2-1)+1;
    d1 = dnas(dp1:dp1+1, :);
    d2 = dnas(dp2:dp2+1, :);
    dna = [d1(:,1:sppt), d2(:,sppt+1:DNAl)];
    dna = mutate(dna);
    
    checkDNA(dna);
end

function dna = mutate(dnaIn)
    global DNAl
    mutationRate = 0.1;
    dna = dnaIn;
    for i = 1:DNAl
        if rand() < mutationRate
            sppt = i;
            obj = floor(rand() * 12) + 1;
            dna(1, sppt) = obj;
        end
    end
end


function bans = checkDNA(dna)
    bans = 1;
    for i = dna
        if i(1) == 0
            bans = 0;
            break;
        end
    end
    if bans == 0
        disp('ERROR: DNA invalid !!!');
    end
end


