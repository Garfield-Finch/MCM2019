function Qvalue = evolve(dna, mode)
%     disp('in evolve()');
%% init
    m = zeros(231, 92);
    bitm = zeros(1, 12);
    
    global sR
    sR = zeros(12, 3);
    sR(1,:) = [45,45,25];  %A
    sR(2,:) = [30,30,22];
    sR(3,:) = [60,50,30];
    sR(4,:) = [25,20,25];
    sR(5,:) = [25,20,27];
    sR(6,:) = [40,40,25];
    sR(7,:) = [32,32,17];  %G
    sR(8,:) = [8,10,14];  %type1
    sR(9,:) = [24,20,20];  %type2
    sR(10,:) = [14,7,5];  %med1
    sR(11,:) = [5,8,5];
    sR(12,:) = [12,7,4];
    
    global numout
    numout = 0;
    
    xylist = [1;1];

%% body
    dnast = 0;
    for i = dna

        dnast = dnast + 1;
%         disp(['============== geneNum: ', num2str(dnast),' ================']);
        
%         pos = decidePos(x, y, i);
        boxtype = i(1);
        pos = i(2);
        
        [xylist, x,y] = decidexy(m, xylist, boxtype, pos);
        if x == -1
            continue;
        end
%         disp('out decidexy()');
        
        sizeBox = sR(boxtype,:);
        bl = sizeBox(1);
        bw = sizeBox(2);
        bh = sizeBox(3);

        [m, bitm, xylist] = setM(m,bitm,xylist, boxtype, x, y, pos, bl, bw, bh);
%         disp('out setM()');
    end    
    
    V = calV(bitm);
    Qvalue = getQ(bitm, V);
    
    %% Visualization
    if mode ~= 0
        disp(['Qvalue: ',num2str(Qvalue)]);
%         disp('bitm');
        bitm
        visM(m, V, bitm);
    end
end


%% Utility functions
function [m, bitm, xylist] = setM(m,bitm,xylist, key, x, y, pos, dl, dw, dh)
%     disp('in setM');
    if pos ~= 0
        pos = 1;
    end
    if pos == 0
        m(x:x+dl-1,y:y+dw-1) = key;
        
        m(x:x+dl-1,y) = 13;
        m(x:x+dl-1,y+dw-1) = 13;
        m(x,y:y+dw-1) = 13;
        m(x+dl-1,y:y+dw-1) = 13;

        newco = [x,x+dl;y+dw,y];
    else
        m(x:x+dw-1,y:y+dl-1) = key;
        
        m(x:x+dw-1,y) = 13;
        m(x:x+dw-1,y+dl-1) = 13;
        m(x,y:y+dl-1) = 13;
        m(x+dw-1,y:y+dl-1) = 13;
        
        newco = [x,x+dw;y+dl,y];
    end
    
    bitm(1, key) = bitm(1, key) + floor(94/dh);

    xylist = [xylist,newco];
    xylist = xylist';
    xylist = sortrows(xylist, 1);
    xylist = xylist';
end

function [xylist, x, y] = decidexy(m, xylist, boxtype, pos)
    x = -1;y = -1;
    global sR
    sizeBox = sR(boxtype,:);
    if pos == 0
        bl = sizeBox(1);
        bw = sizeBox(2);
    else
        bw = sizeBox(1);
        bl = sizeBox(2);
    end
    
%     disp('In decidexy()');
    
    for i = 1:size(xylist,2)
        xi = xylist(1, i);
        yi = xylist(2, i);
        bans = checkBox(m, xi, yi, bl, bw);
        if bans == 1
            x = xi;
            y = yi;
            xylist(:,i) = [];
            break;
        end
    end
end

function bans = checkBox(m, x, y, l, w)
    if x+l > 231 || y+w > 92
        bans = 0;
        return;
    end
    boxm = m(x:x+l,y:y+w);
    sans = sum(sum(boxm));
    if sans == 0
        bans = 1;
    else
        bans = 0;
    end
end

function Vans = calV(bitm)
    dntr = 231 * 92;
    global sR
    vsep = zeros(12, 1);
    for i=1:12
        vsep(i) = sR(i,1) * sR(i,2);
    end
    
    pbitm = bitm;
    for i=1:12
        dh = sR(i,3);
        pbitm(i) = bitm(i) / floor(94/dh);
    end
    
    vsum = pbitm * vsep;
    Vans = (dntr - vsum) / dntr;
end

function visM(m, V, bitm)
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
    
    m = m';
    
    mx = size(m,1);
    my = size(m,2);
    img = zeros(mx,my,3);
    img = uint8(img);
    
    % calculate img
    for i=1:mx
        for j=1:my
            mv = m(i,j);
            if mv ~= 0
                img(i,j,:) = clrmp(mv,:);
            end
        end
    end
    
    disp(['V: ',num2str(V)]);

    % convert bitm to str
    bitms = num2str(bitm(1));
    for i=2:12
        bitms = [bitms,',',num2str(bitm(i))];
    end
    
    close;
    figure(1);
    
    title(['bitm:',bitms]);
    xlabel(['V:',num2str(V*100),'%']);
    
    imshow(img);
    
    global numout
    numout = numout + 1;
    imwrite(img,strcat(num2str(numout),'.png'));
end

function Qans = getQ(bitm, V)
    Qans = 0;
    if bitm(3) == 0 || bitm(9) == 0
        return;
    end
    Qans = QV(V) * QM(bitm) * QD(bitm);
end

%% Seperate Q's
function qans = QV(V)
%     qm = V;
%     Kv = 1; 
%     qans = 1/qm * Kv;
    qans = 1 - V;
end

function qans = QM(bitm)
    med1 = bitm(1, 10);
    med2 = bitm(1, 11);
    med3 = bitm(1, 12);
    ndMed1 = 7;
    ndMed2 = 2;
    ndMed3 = 4;
    qm = min([floor(med1/ndMed1), floor(med2/ndMed2), floor(med3/ndMed3)]);
    Km = 1/10;
    qans = 1 - exp(- Km * qm);
end

function qans = QD(bitm)
    nA = bitm(1);
    nB = bitm(2);
    nC = bitm(3);
    nD = bitm(4);
    nE = bitm(5);
    nF = 0;
    nG = bitm(7);
    Kd = 0.5;
    qans = 1 - exp(-Kd*nB*0.155*21);
    qans = qans + 1 - exp(-Kd*(nB+nC)*0.156*21);
%     qans = qans + 1 - exp(-Kd*(nB+nC+nF)*0.211*21);
    qans = qans + 1 - exp(-Kd*(nB+nC+nF+nA)*0.857*4);
    qans = qans + 1 - exp(-Kd*(nB+nC+nF+nA+nD)*1.667*4);
    qans = qans + 1 - exp(-Kd*(nB+nC+nF+nA+nD+nG)*1.882*4);
    qans = qans + 1 - exp(-Kd*(nB+nC+nF+nA+nD+nG+nE)*2*4);
    
    qans =qans / 6;
end

