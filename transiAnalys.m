%% He-energy matr
%read data
thin = char(8201);
addpath heliumbös
heL=importdata('heliumL.txt');
heJ=importdata('heliumJ.txt');
heS=importdata('heliumS.txt');
heLS=importdata('heliumLS.txt');
heConf=importdata('heliumConfShort.txt');
heElevel=importdata('heliumenergy.txt');
antalData=length(heElevel);
heEnergyDelta=ones(antalData,antalData);%rad minus
for i=1:antalData
    for j=1:antalData
        heEnergyDelta(i,j)=heElevel(i)-heElevel(j);
        if heEnergyDelta(i,j)<=0 %We want a "triangular" matrix,
            heEnergyDelta(i,j)=1e38; %Something large
        end
    end
end
toppar=[];
%% Fitting
%Runs Oskars fit peaks script

kvant_deluxe
switch model 
    case 9    
        toppartemp=fitted_peaks(:,2);
    case 11  
        toppartemp=fitted_peaks(:,1);

end
 %% probably not really needed
   toppar=horzcat(toppar, toppartemp);

%% Match dipol

konvfaktor=10000000 ; %don't ask
topparcm=1./toppar.*konvfaktor;
topparkoll=1./topparcm.*konvfaktor;
antalToppar=length(toppar);
inteDipol=true; %while loop flag
k=1; %counter for "false" matches

TransmissionForbiddenIndex={};
TransmissionForbiddenString={};
TransmissionAllowedString={};
TransmissionAllowedIndex={};
falsecounter=0*[1:antalToppar];
energylevels=[];
energylevelsindex=[];
TransmissionAllowedString{1}=['#','  ','Measurement',' ','Theory','  ', 'confH',' ', 'confL', ' ','LS',thin, 'L',' LS',thin, 'H', ' ', '   #false ', 'diff2false','  diff2theory     '];
for i=1:antalToppar
    k=1;
    E_First=0;
    inteDipol=true;
    %we need some temps, because that's cool
    forbiddenStringTemp={};
    indexTemp={};
    heEnergyDeltaTemp=heEnergyDelta; %
    while inteDipol
        energyDiff=heEnergyDeltaTemp-topparcm(i);
        %WE want 2 index(indeii?), I for for upper energy state, J for the
        %lower, the min function is a bit tricky with matrii in older
        %matlab
        [energyMin1 energyIndexradtemp]=(min(abs(energyDiff)));
        [energyMin2 energyIndexcolontemp]=min(min(abs(energyDiff)));
        energyI=energyIndexradtemp(energyIndexcolontemp);
        energyJ=energyIndexcolontemp;
        if energyJ==1
        end
        nmteory=1/heEnergyDelta(energyI,energyJ).*konvfaktor; %"theoretical" wavelength
        
        %Rules for dipol transission
        deltaS=heS(energyI)-heS(energyJ);
        deltaJ=heJ(energyI)-heJ(energyJ);
        %this took way to many tries, but it's just an if statement broken
        %up in a few steps, matlab did some wierd things with one of the
        %statments in an earlier version
        %"If it works, it works" <- applies to this wole .m file
        deltaJstupid=sum([deltaJ==0 deltaJ==1 deltaJ==-1]);
        JnotZEROstupid=((abs(heJ(energyI))+abs(heJ(energyJ)))~=0);
        deltaSstupid=(deltaS==0);
        if (JnotZEROstupid && deltaJstupid && deltaSstupid)
            TransmissionAllowedIndex{i}=[energyI energyJ]; %saves index
            if E_First~=0.0000
                deltaEFirstPicked=toppar(i)-E_First; %if we've failed the above if statement and have a false transition
            else
                deltaEFirstPicked=0;
            end
            deltaEAllowedTeory=(toppar(i)-nmteory);
            %builds a cell of a bunch of strings, if statements for nicer
            %layout
            space1='  '; 
            space2='        ';
            if (k-1)==0
                space=space2;
            else
                space=space1;
            end
            if (k-1)>9
                sp=' ';
            else
                sp='  ';
            end
            if i<10
                TransmissionAllowedString{i+1}=[num2str(i),'  ',num2str(toppar(i)),thin,'nm',' ',num2str(nmteory),thin,'nm','  ', heConf{energyI},' ', heConf{energyJ}, ' ', heLS{energyI},num2str(heJ(energyI)), ' ', heLS{energyJ}, num2str(heJ(energyJ)),'      ',num2str(k-1),sp,num2str(deltaEFirstPicked),thin,'nm',space,num2str(deltaEAllowedTeory),thin,'nm'];
            else
                TransmissionAllowedString{i}=[num2str(i),' ',num2str(toppar(i)),thin,'nm',' ',num2str(nmteory),thin,'nm','  ', heConf{energyI},' ', heConf{energyJ}, ' ', heLS{energyI},num2str(heJ(energyI)), ' ', heLS{energyJ}, num2str(heJ(energyJ)),'      ',num2str(k-1),sp,num2str(deltaEFirstPicked),thin,'nm',space ,num2str(deltaEAllowedTeory),thin,'nm'];
            end
            inteDipol=0;
        else %saves "false" transitions in a separate cellarray
            indexTemp{k}=[energyI energyJ];
            if k==1
                E_First=nmteory;
            end
            %just saves away which rule we broke
            Srule=(deltaS==0);
            Jrule=(deltaJ==1 || deltaJ==-1 || deltaJ==0);
            Jrule0=(heJ(energyI~=0) || heJ(energyJ)~=0);
            if k <10
                forbiddenStringTemp{k}=[num2str(k),'   ',num2str(toppar(i)),thin,'nm',' ',num2str(nmteory),thin,'nm',' ', heConf{energyI},'  ', heConf{energyJ}, ' ', heLS{energyI},num2str(heJ(energyI)), ' ', heLS{energyJ}, num2str(heJ(energyJ)),' ',num2str(Srule),' ',num2str(Jrule),' ',num2str(Jrule0)];
            else
                forbiddenStringTemp{k}=[num2str(k),'  ',num2str(toppar(i)),thin,'nm',' ',num2str(nmteory),thin,'nm',' ', heConf{energyI},'  ', heConf{energyJ}, ' ', heLS{energyI},num2str(heJ(energyI)), ' ', heLS{energyJ}, num2str(heJ(energyJ)),' ',num2str(Srule),' ',num2str(Jrule),' ',num2str(Jrule0)];
            end
            TransmissionForbiddenIndex{k}=[energyI energyJ];
            heEnergyDeltaTemp(energyI,energyJ)=1e38; %we dont want to check the same transition over and over
            k=k+1;
            if k>antalData
                inteDipol=0;
                disp('DET SKET SIG') 
            end
        end
    end
    TransmissionForbiddenIndex{i}=indexTemp;
    TransmissionForbiddenString{i}=forbiddenStringTemp;
    falsecounter(i)=k-1;
    energylevels=[energylevels  heElevel(energyI) heElevel(energyJ)];
    energylevelsindex=[energylevelsindex energyI energyJ];
end
falsePrintInfo(falsecounter); %mostly for debugging
%printCellCust(TransmissionForbiddenString);
printCellCust(TransmissionAllowedString) %gives nice printout


%% Plot diagram
% a bunch of this is constants and numbers that are very "godtyckliga" and
% hardcoded for nicer looking diagram
figure
h=4.135667662e-15;
cspeed=299792458; 
mellanlambdameter=1e-3./energylevels;
energylevelsplot=energylevels;

hold on
deltax=0.8; 
deltaxarrow=0;
singlett=0; % for if statement
for i=2:2:length(energylevelsplot)

   if ((heS(energylevelsindex(i))==singlett) && (heS(energylevelsindex(i-1))==singlett))
%     if i==length(energylevelsplot)
%         pause
%     end
    deltax=0.8;
   else
       deltax=0.8+2.4+0.6; %change this if you get more spectral lines of singlet states
   end
    if abs(heL(energylevelsindex(i-1))-heL(energylevelsindex(i)))==0
        deltaxarrow=i*0.0009; %godtyckligt
    else
        deltaxarrow=0;
    end
    
    
    
    mylineLow=line([0 0], [energylevelsplot(i) energylevelsplot(i)]);  %create a line objext length 0
    xlow=get(mylineLow,'XData'); %these bunch of 0 will be manipulated in some wierd way
    ylow=get(mylineLow,'YData'); %well, same for these
    xlinelow=[floor(deltax+heL(energylevelsindex(i))) deltax+heL(energylevelsindex(i))]; %lower energylvl
    xpointlow=sum(xlinelow)/2; %transition "arrow" endpoint
    
    mylinehigh=line([0 0], [energylevelsplot(i-1) energylevelsplot(i-1)]);%line pbject for the upper energy level
    xhigh=get(mylinehigh,'XData');%fetch data
    yhigh=get(mylinehigh,'YData');
    xlinehigh=[floor(deltax+heL(energylevelsindex(i-1))) deltax+heL(energylevelsindex(i-1))]; %upper energy level
    xpointhigh=sum(xlinehigh)/2;%transition arrow start
    
    %draws a bunch of lines
    line(xlinehigh,yhigh,'color','black') ; 
    line(xlinelow,ylow,'color','black'); 
    line([xpointhigh+deltaxarrow xpointlow],[energylevelsplot(i-1) energylevelsplot(i)], 'linestyle','--'); %"arrow"
    
    %the next 2 lines dnt work after the singlett flag was implemented
    LineStringInPlot=num2str(i/2); %markers
   % text((xpointhigh+xpointlow)/2+0.2, energylevelsev(i)+ (energylevelsev(i-1)-energylevelsev(i))/2,LineStringInPlot);
    
   % pause
    
end

%what follows is just to make the axis and stuff look proper
spdf={texlabel('^(1)S');texlabel('^(1)P');texlabel('^(1)D'); ...
    texlabel('^(3)S');texlabel('^(3)P');texlabel('^(3)D');texlabel('^(3)F');texlabel('^(3)G')};
set(gca,'xtick',[0.4:7.4],'xticklabel',spdf)
midlineX=(2.4/2+3.4/2);
xmidline=[midlineX midlineX];



ymidline=[1.55e5 2e5];
plot(xmidline, ymidline,'k','linewidth',1.25)
%ylabelstring=('Energy [' texlabel('cm^(-1)') ')')1009.146162695160
ylabelstring=strcat('Energi  [',texlabel('cm^(-1)'),']');
ylabel(ylabelstring,'fontsize',12)
xlabel('LS-termer ','fontsize',12)
%d är på 2.4

