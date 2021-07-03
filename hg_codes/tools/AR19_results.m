%% thickness plot

% Res = readtable('G:\HAW_tip\data\Sizing_Results.xlsx','Sheet','AR19_FWT');

% Eta=0.85
x1=[0.027340383	0.028603848	0.029725752	0.029974759	0.031034137	0.032735202	0.039609329	0.043529875	0.046493235	0.04608511	0.045938318	0.045959811	0.045988237	0.041956391	0.041430907	0.040835445	0.039913865	0.038766877	0.033458223	0.036063169	0.037005598	0.03764346	0.034429109	0.028950006	0.016087067	0.009141008	0.009203225	0.009273596	0.009372603	0.00950832	0.009661973	0.009631135	0.009746061	0.009932049	0.00963159	0.009315086	0.00898109	0.00861313	0.008270483	0.007891068	0.007472047	0.007165947	0.006975031	0.0067466	0.006363138	0.005839796	0.005179279	0.004282995	0.003211139	0.001331628	0.000117651	0.00012046	0.000123778	0.000128853	0.000136692	0.000146527	0.000142149	0.000148416	0.000159236	0.000135214	0.000116392	0.000101625	8.95E-05	8.08E-05	7.31E-05	6.55E-05	6.02E-05	5.69E-05	5.34E-05	4.73E-05	3.99E-05	3.14E-05	2.18E-05	1.16E-05	1.84E-06
];

Y1=[2	2.758168763	3.516337526	4.274506289	5.032675052	5.790843815	6.549012578	7.307181341	8.065350104	8.879679516	9.694008929	10.50833834	11.32266775	12.13699717	12.95132658	13.76565599	14.5799854	15.39431481	16.20864423	17.02297364	17.83730305	18.65163246	19.46596187	20.28029129	21.0946207	21.43158459	21.76854849	22.10551238	22.44247628	22.77944017	23.11640407	23.45336796	23.79033186	24.12729575	24.46425965
];

Disp1=[-0.002772512	-0.000361963	0.007838251	0.022330484	0.043739537	0.072838848	0.110581921	0.157932249	0.21599734	0.293564661	0.385739256	0.492061201	0.612012237	0.74502119	0.890559837	1.048048064	1.216694698	1.395507845	1.583326153	1.778940906	1.981027635	2.188109567	2.39878092	2.611615968	2.817451885	2.973944539	3.12234841	3.270657232	3.418852542	3.566936181	3.714922818	3.86283594	4.010702444	4.158546122	4.306383683
];


% Eta=0.8

x2=[0.027271727	0.028170607	0.029030873	0.028870509	0.030299687	0.031496228	0.037214804	0.0402045	0.042006046	0.045425073	0.045066543	0.044433855	0.044061764	0.043657371	0.042976798	0.04090401	0.037839473	0.036353178	0.034373859	0.034904096	0.035869532	0.035481756	0.032934135	0.026823544	0.01575832	0.009067587	0.009076604	0.009115128	0.00916492	0.009227016	0.009302398	0.009159568	0.009178833	0.009206502	0.009262611	0.008962744	0.008651363	0.008316654	0.00795962	0.007575908	0.007172649	0.006811135	0.006602081	0.006327993	0.005931433	0.00542962	0.004795852	0.004005556	0.002994765	0.001739658	0.00011314	0.000113382	0.000114783	0.000116651	0.000119396	0.000123002	0.00011471	0.000114644	0.000114304	0.000115537	0.000102167	9.15E-05	8.25E-05	7.49E-05	6.78E-05	6.07E-05	5.47E-05	5.13E-05	4.71E-05	4.14E-05	3.46E-05	2.69E-05	1.85E-05	9.75E-06	3.55E-06
];

Y2=[2	2.758168763	3.516337526	4.274506289	5.032675052	5.790843815	6.549012578	7.307181341	8.065350104	8.809478705	9.553607306	10.29773591	11.04186451	11.78599311	12.53012171	13.27425031	14.01837891	14.76250751	15.50663611	16.25076471	16.99489331	17.73902191	18.48315051	19.22727912	19.97140772	20.42069291	20.8699781	21.3192633	21.76854849	22.21783368	22.66711887	23.11640407	23.56568926	24.01497445	24.46425965
];

Disp2=[-0.002775043	-0.000410158	0.007657646	0.021953467	0.043118093	0.071925086	0.109340207	0.156369206	0.214203058	0.284412842	0.366935981	0.461345343	0.567172399	0.683903712	0.810980024	0.947817144	1.093704893	1.247775943	1.409114856	1.576742918	1.749619962	1.926682672	2.106864275	2.289115479	2.465393195	2.661636811	2.850651071	3.039325681	3.227596409	3.415465025	3.602980019	3.790219934	3.977277005	4.16424602	4.351198481
];



% Eta=0.75

x3=[0.026391931	0.027028704	0.028364051	0.029047692	0.030500032	0.031145441	0.039680499	0.042544402	0.044152565	0.04826052	0.046406844	0.045722449	0.045627934	0.044862853	0.043556921	0.041773117	0.039681596	0.038048044	0.035831663	0.033968792	0.033726868	0.032390219	0.029145185	0.02342323	0.019152407	0.008938093	0.008980212	0.009029353	0.009086682	0.009176732	0.009288384	0.009213429	0.009276953	0.009355275	0.009462643	0.0091693	0.008850428	0.008507324	0.008148999	0.007771577	0.007394087	0.00696701	0.006505875	0.006135004	0.005723503	0.005187593	0.004547727	0.003756063	0.002756568	0.002075519	0.00010777	0.000109278	0.000111098	0.000113313	0.000117207	0.000122294	0.000117327	0.00011976	0.000122454	0.000125624	0.000110453	9.77E-05	8.71E-05	7.86E-05	7.13E-05	6.45E-05	5.73E-05	5.00E-05	4.44E-05	3.87E-05	3.18E-05	2.44E-05	1.65E-05	8.67E-06	4.60E-06
];

Y3=[2	2.673927789	3.347855579	4.021783368	4.695711157	5.369638947	6.043566736	6.717494526	7.391422315	8.065350104	8.784206413	9.503062722	10.22191903	10.94077534	11.65963165	12.37848796	13.09734426	13.81620057	14.53505688	15.25391319	15.9727695	16.69162581	17.41048212	18.12933843	18.84819473	19.40980123	19.97140772	20.53301421	21.0946207	21.65622719	22.21783368	22.77944017	23.34104666	23.90265315	24.46425965
];

Disp3=[-0.002761397	-0.000991654	0.005114411	0.015875878	0.031660212	0.052914022	0.080168804	0.113904401	0.154695228	0.203364147	0.266220594	0.339346614	0.422374742	0.51486478	0.616339098	0.726295493	0.844215184	0.969520316	1.101536682	1.239465085	1.382417877	1.529442511	1.67958997	1.831926141	1.979828575	2.208548217	2.430952915	2.652284724	2.872458212	3.091601487	3.309934583	3.527704425	3.745139373	3.962418806	4.17965875
];


% Eta=0.7

x4=[0.026303128	0.026698557	0.027761331	0.0281162	0.029918642	0.03027628	0.038068198	0.040086413	0.040972039	0.043097354	0.045482365	0.044731783	0.04385713	0.045680034	0.044760645	0.04372247	0.040786127	0.03942214	0.03784525	0.035936987	0.033021589	0.031101034	0.027548656	0.021815492	0.0190545	0.008825121	0.008833359	0.008844332	0.008860649	0.008876047	0.008937498	0.008791227	0.008753988	0.008745263	0.008720004	0.008733003	0.008428895	0.008115727	0.007752967	0.0073929	0.007002846	0.006606562	0.006196882	0.005749686	0.005267912	0.004802838	0.004224674	0.003555158	0.002643589	0.002276389	0.000103202	0.000103268	0.000103402	0.000103715	0.000103946	0.000105987	9.96E-05	9.78E-05	9.70E-05	9.50E-05	9.45E-05	8.56E-05	7.82E-05	7.10E-05	6.46E-05	5.79E-05	5.16E-05	4.54E-05	3.92E-05	3.29E-05	2.73E-05	2.11E-05	1.44E-05	7.64E-06	5.55E-06
];

Y4=[2	2.673927789	3.347855579	4.021783368	4.695711157	5.369638947	6.043566736	6.717494526	7.391422315	8.065350104	8.709325547	9.353300991	9.997276434	10.64125188	11.28522732	11.92920276	12.57317821	13.21715365	13.86112909	14.50510454	15.14907998	15.79305542	16.43703087	17.08100631	17.72498175	18.39890954	19.07283733	19.74676512	20.42069291	21.0946207	21.76854849	22.44247628	23.11640407	23.79033186	24.46425965
];

Disp4=[-0.002775533	-0.001062239	0.004882421	0.015377841	0.030785877	0.051535808	0.078125319	0.111071069	0.150972053	0.198649701	0.253393643	0.316249859	0.386865067	0.464821478	0.549690205	0.641008091	0.738312032	0.841095577	0.948836748	1.06098471	1.176932908	1.296015432	1.417534727	1.540798577	1.660339404	1.927406146	2.188794624	2.448563652	2.706402135	2.962297461	3.216472368	3.469289952	3.721195712	3.972638715	4.223950452
];


% Eta= 0.65

x5=[0.026974484	0.027004133	0.026407057	0.027396789	0.028970195	0.031420536	0.031307892	0.039883561	0.04331768	0.044247968	0.047795639	0.046652932	0.047430507	0.046165846	0.041411996	0.04332418	0.042467943	0.039947531	0.041392436	0.038570762	0.03555692	0.032640138	0.028832423	0.024926993	0.021774466	0.008622199	0.008646185	0.00867164	0.008713096	0.008738759	0.008783291	0.008880871	0.008733315	0.008733182	0.00875697	0.008785893	0.00849097	0.008159432	0.007832171	0.007484035	0.007110287	0.006700303	0.006298552	0.005832914	0.00535928	0.004841409	0.00425523	0.003585636	0.002816592	0.002527369	9.57E-05	9.63E-05	9.69E-05	9.80E-05	9.86E-05	9.98E-05	0.000102223	9.66E-05	9.62E-05	9.62E-05	9.58E-05	8.69E-05	7.88E-05	7.25E-05	6.63E-05	5.98E-05	5.31E-05	4.68E-05	3.97E-05	3.34E-05	2.72E-05	2.13E-05	1.51E-05	8.46E-06	7.07E-06
];

Y5=[2	2.60653501	3.213070021	3.819605031	4.426140042	5.032675052	5.639210063	6.245745073	6.852280083	7.458815094	8.065350104	8.675094295	9.284838485	9.894582675	10.50432687	11.11407106	11.72381525	12.33355944	12.94330363	13.55304782	14.16279201	14.7725362	15.38228039	15.99202458	16.60176877	17.38801786	18.17426694	18.96051603	19.74676512	20.53301421	21.3192633	22.10551238	22.89176147	23.67801056	24.46425965
];

Disp5=[-0.0027818	-0.001506023	0.00300761	0.010965348	0.022585114	0.038121129	0.057864423	0.082178169	0.111341567	0.145650828	0.185586348	0.232823667	0.286334609	0.345804578	0.410927085	0.481403209	0.556844854	0.636839142	0.720956041	0.808692766	0.899534154	0.992998009	1.088580708	1.185735761	1.280108588	1.576294809	1.867546092	2.156604994	2.443068256	2.726889771	3.008333927	3.287857268	3.56604605	3.843512166	4.120757945
];



% eTA=0.6

x6=[0.0280464	0.028003793	0.027415338	0.027310235	0.027530586	0.028346271	0.028093968	0.029109149	0.039468121	0.042631008	0.043789303	0.046898291	0.04630414	0.04668308	0.045149366	0.043588885	0.041963957	0.040287201	0.038818888	0.039254899	0.036100515	0.033167253	0.030305365	0.027096041	0.02005283	0.008409679	0.008409547	0.008412361	0.008448422	0.008459648	0.008506125	0.00853697	0.008620353	0.008395352	0.008353848	0.008345607	0.008329362	0.008022131	0.007684721	0.007345015	0.006986021	0.006601387	0.006186542	0.005756827	0.005291025	0.004796036	0.004235451	0.003577262	0.00283683	0.002691354	8.89E-05	8.87E-05	8.87E-05	8.94E-05	8.95E-05	9.06E-05	9.11E-05	9.33E-05	8.62E-05	8.46E-05	8.42E-05	8.27E-05	7.60E-05	6.98E-05	6.38E-05	5.77E-05	5.16E-05	4.52E-05	3.90E-05	3.26E-05	2.68E-05	2.10E-05	1.49E-05	8.65E-06	8.57E-06
];

Y6=[2	2.551395464	3.102790928	3.654186392	4.205581856	4.75697732	5.308372784	5.859768248	6.411163712	6.962559176	7.51395464	8.065350104	8.635596695	9.205843286	9.776089877	10.34633647	10.91658306	11.48682965	12.05707624	12.62732283	13.19756942	13.76781601	14.33806261	14.9083092	15.47855579	16.37712617	17.27569656	18.17426694	19.07283733	19.97140772	20.8699781	21.76854849	22.66711887	23.56568926	24.46425965
];

Disp6=[-0.002788299	-0.001826516	0.00165509	0.007789471	0.016725916	0.028639263	0.04372769	0.0622207	0.084374135	0.110345105	0.140268809	0.174447697	0.215390644	0.261216692	0.311652023	0.366414489	0.425213733	0.487735922	0.553632525	0.62252015	0.693982072	0.767602564	0.842977856	0.919652056	0.994210016	1.320208844	1.641295067	1.95886029	2.272236022	2.581620572	2.887716705	3.191331333	3.493181127	3.793968496	4.094268064
];


% thickness plot 
figure 
plot(Y1(1:25)-2,x1(1:25),'bs-','MarkerFaceColor','b') % 0.85

% hold on 
% plot(Y2(1:25),x2(1:25),'ks-','MarkerFaceColor','k') % 0.8
hold on 
plot(Y3(1:25)-2,x3(1:25),'ks-','MarkerFaceColor','k') % 0.75
% hold on 
% plot(Y4(1:25),x4(1:25),'gs-','MarkerFaceColor','g')  % 0.7
hold on 
plot(Y5(1:25)-2,x5(1:25),'rs-','MarkerFaceColor','r')  % 0.65
% hold on 
% plot(Y6(1:25),x6(1:25),'rs-','MarkerFaceColor','r') % 0.60
xlabel('Spanwise Distance (m)','Interpreter','latex','FontSize',12)
ylabel('Spar thickness (m)','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')
legend('Fold length = 15 \%','Fold length = 25 \%','Fold length = 35 \%','Interpreter','latex','FontSize',12)


figure 
plot(Y1(1:25)-2,x1(26:50),'bs-','MarkerFaceColor','b') % 0.85

% hold on 
% plot(Y2(1:25),x2(1:25),'ks-','MarkerFaceColor','k') % 0.8
hold on 
plot(Y3(1:25)-2,x3(26:50),'ks-','MarkerFaceColor','k') % 0.75
% hold on 
% plot(Y4(1:25),x4(1:25),'gs-','MarkerFaceColor','g')  % 0.7
hold on 
plot(Y5(1:25)-2,x5(26:50),'rs-','MarkerFaceColor','r')  % 0.65
% hold on 
% plot(Y6(1:25),x6(1:25),'rs-','MarkerFaceColor','r') % 0.60
xlabel('Spanwise Distance (m)','Interpreter','latex','FontSize',12)
ylabel('Skin thickness (m)','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')
legend('Fold length = 15 \%','Fold length = 25 \%','Fold length = 35 \%','Interpreter','latex','FontSize',12)


figure 
plot(Y1(1:25)-2,x1(51:75),'bs-','MarkerFaceColor','b') % 0.85

% hold on 
% plot(Y2(1:25),x2(1:25),'ks-','MarkerFaceColor','k') % 0.8
hold on 
plot(Y3(1:25)-2,x3(51:75),'ks-','MarkerFaceColor','k') % 0.75
% hold on 
% plot(Y4(1:25),x4(1:25),'gs-','MarkerFaceColor','g')  % 0.7
hold on 
plot(Y5(1:25)-2,x5(51:75),'rs-','MarkerFaceColor','r')  % 0.65
% hold on 
% plot(Y6(1:25),x6(1:25),'rs-','MarkerFaceColor','r') % 0.60
xlabel('Spanwise Distance (m)','Interpreter','latex','FontSize',12)
ylabel('Stringer area (m$^2$)','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')
legend('Fold length = 15 \%','Fold length = 25 \%','Fold length = 35 \%','Interpreter','latex','FontSize',12)

% Displacement

figure 
plot(Y1(1:end)-2,Disp1,'bs-','MarkerFaceColor','b') % 0.85

% hold on 
% plot(Y2(1:25),x2(1:25),'ks-','MarkerFaceColor','k') % 0.8
hold on 
plot(Y3(1:end)-2,Disp3,'ks-','MarkerFaceColor','k') % 0.75
% hold on 
% plot(Y4(1:25),x4(1:25),'gs-','MarkerFaceColor','g')  % 0.7
hold on 
plot(Y5(1:end)-2,Disp5,'rs-','MarkerFaceColor','r')  % 0.65
% hold on 
% plot(Y6(1:25),x6(1:25),'rs-','MarkerFaceColor','r') % 0.60
xlabel('Spanwise Distance (m)','Interpreter','latex','FontSize',12)
ylabel('Vertical displacement (m)','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')
legend('Fold length = 15 \%','Fold length = 25 \%','Fold length = 35 \%','Interpreter','latex','FontSize',12)

% Coast angle

fold_length=[0.15,0.2,0.25,0.3,0.35,0.4];

coast_angle1=atan((Disp1(end)-Disp1(end-10))/(Y1(end)-Y1(end-10)))*57;

coast_angle2=atan((Disp2(end)-Disp2(end-10))/(Y2(end)-Y2(end-10)))*57;

coast_angle3=atan((Disp3(end)-Disp3(end-10))/(Y3(end)-Y3(end-10)))*57;

coast_angle4=atan((Disp4(end)-Disp4(end-10))/(Y4(end)-Y4(end-10)))*57;

coast_angle5=atan((Disp5(end)-Disp5(end-10))/(Y5(end)-Y5(end-10)))*57;

coast_angle6=atan((Disp6(end)-Disp6(end-10))/(Y6(end)-Y6(end-10)))*57;

coast_angles=[coast_angle1,coast_angle2,coast_angle3,coast_angle4,coast_angle5,coast_angle6];

figure
plot(fold_length,coast_angles,'r-s','MarkerFaceColor','r','LineWidth',1)
xlabel('Fold length $(\%)$','Interpreter','latex','FontSize',12)
ylabel('Coast angle ($^{\circ}$)','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')


%% lift var

Lift_Res = readtable('G:\HAW_tip\data\Sizing_Results.xlsx','Sheet','AR19_FWT_L');

figure
plot(Lift_Res.Y7,Lift_Res.L7,'--','Color','k','LineWidth',1);
hold on 
plot(Lift_Res.Y1,Lift_Res.L1,'bs-','MarkerFaceColor','b');
hold on 
plot(Lift_Res.Y2,Lift_Res.L2,'ks-','MarkerFaceColor','k');
hold on 
plot(Lift_Res.Y3,Lift_Res.L3,'ms-','MarkerFaceColor','m');
hold on 
plot(Lift_Res.Y4,Lift_Res.L4,'ks-','MarkerFaceColor','g');
hold on 
plot(Lift_Res.Y5,Lift_Res.L5,'rs-','MarkerFaceColor','r');
hold on 
plot(Lift_Res.Y6,Lift_Res.L6,'ks-','MarkerFaceColor','y');

xlabel('Spanwise distance (m)','Interpreter','latex', 'FontSize',12);
ylabel('Lift per unit span (N/m)','Interpreter','latex', 'FontSize',12);
set(gcf,'Color','w')
legend('No FWT','Fold length = 0.15','Fold length = 0.20','Fold length = 0.25','Fold length = 0.30','Fold length = 0.35','Fold length = 0.40','Interpreter','latex','FontSize',12)


%% weight
Eta=[0,15,20,25,30,35,40];

Wing_Weight=[4315.9	4066.12	3910.3	3852.8	3705.8	3639.19	3506.878];

AOA=[0.1506	0.152	0.153	0.155	0.156	0.159];

figure 
plot(Eta,Wing_Weight,'bs-','MarkerFaceColor','b')
xlabel('Fold length ($\%$)','Interpreter','latex','FontSize',12)
ylabel('Wing weight (kg)','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')

figure 
plot(Eta(2:end),AOA,'bs-','MarkerFaceColor','b')
xlabel('Fold length ($\%$)','Interpreter','latex','FontSize',12)
ylabel('Root AoA (Rad)','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')


%% Lift dis

    fold_angle  = -10;   %[deg],
    flare_angle = 25;   %[deg],
    fold_eta=0.6;
    hinge_stiffness=1e-4;
    
%% Wing configurations for starboard wing
  
    Aspect_ratio=19; % Aspect ratio = 10.172 for A321 model
    
    Total_area=126;         % include two wing surface areas + floor size on the fuselage
    Fuselage_width=4;       % dimeter of the fuselage
  
    Wing_span = sqrt(Aspect_ratio*Total_area);
    BeamLoc = 0.4;          % choose a spar location: 0 --> 1
    Semi_span=(Wing_span-Fuselage_width)/2; % length of one wing: 16m for A321 model
    
    Root_chord =  Total_area/(1.064*Semi_span + 4);
    LE_sweep=27;            % deg
    
    Wing_area = (Total_area - Fuselage_width*Root_chord)/2;
    
    Mid_chord=0.63685*Root_chord;
    Tip_chord=0.2248*Root_chord;
    
    X0=Root_chord; 
    X1=0.27*Semi_span*tan(27*pi/180) + Mid_chord;
    X2=Semi_span*tan(27*pi/180) + Tip_chord;
    
    tan_TE_sweep1=(X1-X0)/(0.27*Semi_span);
    tan_TE_sweep2=(X2-X1)/(0.73*Semi_span);
    
    TE_sweep1=atan(tan_TE_sweep1)*180/pi; % deg
    TE_sweep2=atan(tan_TE_sweep2)*180/pi; % deg
      
 
    Taper_ratio=Tip_chord/Root_chord;
    
    Mean_cord_coefficient=(2/3)*(1+Taper_ratio+Taper_ratio^2)/(1+Taper_ratio);
    
    
    %% obtain wingbox geometric properties 

    Wingbox = awi.model.LiftingSurface;
    Wingbox.Origin=[20,2,0];

    %Use the Leading/Trailing edge sweep to define the planform
    Wingbox.ActiveSet = 'sSet';

    % Num of element
    Wingbox.NumBeamElem = 23;

    %Wing dimensions
    Wingbox.SpanVector  = 'Y';
    Wingbox.Span        = Semi_span;   %34.1/2;
    Wingbox.LESweep     = [LE_sweep, LE_sweep];
    Wingbox.LESweep_eta = [0, 1];
    Wingbox.TESweep     = [TE_sweep1, TE_sweep2, TE_sweep2];
    Wingbox.TESweep_eta = [0, 0.27, 1];
    Wingbox.RootChord   = Root_chord;
    
    build(Wingbox)
    
    
    FWT = insertWingFold(Wingbox, 'FlareAngle', flare_angle, 'FoldAngle', fold_angle,'EtaFold',fold_eta);
    FWT.HingeStiffness = [1e14 1e14 1e14 1e14 hinge_stiffness 1e14];
    
    


res_aeroF = mni.result.f06(strcat('D:\MATLAB_workspace\ALENA-master_v1\ALENA-master\hg_codes\Sizing_analysis\Result\AR19_FWT_eta60','/A321_36000ft_1g.f06')).read_aeroF;

% % index for each lifting surfaces: 0.85
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:660; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=911:1010; %tail wing _right

% % index for each lifting surfaces: 0.80
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:640; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=871:990; %tail wing _right

% % index for each lifting surfaces: 0.75
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:620; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=831:970; %tail wing _right

% % index for each lifting surfaces: 0.70
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:610; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=811:970; %tail wing _right

% % index for each lifting surfaces: 0.65
% 
% idx_w85_bef =411:470; % wing_right_bef_kink
% idx_w85_aft =471:590; % wing_right_aft_kink
% idx_c85=1:50; %conn_right
% idx_fwt85=771:940; %tail wing _right

% index for each lifting surfaces: 0.60

idx_w85_bef =411:470; % wing_right_bef_kink
idx_w85_aft =471:580; % wing_right_aft_kink
idx_c85=1:50; %conn_right
idx_fwt85=751:930; %tail wing _right

% right side of AC
lift_wing85_bef=res_aeroF.aeroFz(idx_w85_bef)';
lift_wing85_aft=res_aeroF.aeroFz(idx_w85_aft)';
lift_conn85=res_aeroF.aeroFz(idx_c85)';
lift_fwt85=res_aeroF.aeroFz(idx_fwt85)';

% calculate panel width for each segment: conn + wing_bef_kink +
% wing_aft_kink + FWT

R_chord=Wingbox.Chord(1);
semi_span=Wingbox.Span+FWT.Span;
conn=2;
kink_pos=Wingbox.YData(2);

width_conn=2/5;
width_wing1=kink_pos/6;

wing_sec2_num=ceil((Wingbox.YData(3)-Wingbox.YData(2))/(2.5*Wingbox.Chord(2)/10));
width_wing2=(Wingbox.YData(3)-Wingbox.YData(2))/wing_sec2_num;


fwt_num=ceil(FWT.Span/(2.5*Wingbox.Chord(3)/10));
width_fwt=FWT.Span/fwt_num;


% normalise lift by panel width to obtain lift per unit span

lift_wing=[lift_conn85;lift_wing85_bef;lift_wing85_aft;lift_fwt85];

lift_wing_matrix=reshape(lift_wing, 10, numel(lift_wing)/10);


lift_wing_matrix(:,1:5)=lift_wing_matrix(:,1:5)/width_conn;
lift_wing_matrix(:,6:11)=lift_wing_matrix(:,6:11)/width_wing1;
lift_wing_matrix(:,12:12+wing_sec2_num-1)=lift_wing_matrix(:,12:12+wing_sec2_num-1)/width_wing2;
lift_wing_matrix(:,12+wing_sec2_num:12+wing_sec2_num+fwt_num-1)=lift_wing_matrix(:,12+wing_sec2_num:12+wing_sec2_num+fwt_num-1)/width_fwt;

% lift per unit span
lift_wing_var=sum(lift_wing_matrix);

% find corresponding Y positions
Y_conn=0.5*width_conn:width_conn*0.999:2;

Y_wing1=2+0.5*width_wing1:width_wing1*0.999:2+Wingbox.YData(2);

Y_wing2=2+Wingbox.YData(2)+ 0.5*width_wing2:width_wing2*0.999:2+Wingbox.YData(3);

Y_fwt=2+Wingbox.YData(3)+0.5*width_fwt:width_fwt*0.999:2+Wingbox.Span+FWT.Span;

Y=[Y_conn,Y_wing1,Y_wing2,Y_fwt];



figure 
Y=Y';
lift_wing_var=lift_wing_var';
plot(Y,lift_wing_var,'s-')

%% polar
K=[0.035838281, 0.056235176, 0.057975428, 0.0599742, 0.065769787, 0.067821045,0.07413421];
cl=0:0.05:0.8;
cd_ref=K(1)*cl.^2 + 0.0174;

cd_fwt1=K(2)*cl.^2 + 0.0174;
cd_fwt2=K(3)*cl.^2 + 0.0174;
cd_fwt3=K(4)*cl.^2 + 0.0174;
cd_fwt4=K(5)*cl.^2 + 0.0174;
cd_fwt5=K(6)*cl.^2 + 0.0174;
cd_fwt6=K(7)*cl.^2 + 0.0174;

figure 

plot(cd_ref,cl,'k--','LineWidth',1)
hold on 
plot(cd_fwt1,cl,'b-s','LineWidth',0.8,'MarkerFaceColor','b')
hold on 
plot(cd_fwt2,cl,'k-s','LineWidth',0.8,'MarkerFaceColor','k')
hold on 
plot(cd_fwt3,cl,'m-s','LineWidth',0.8,'MarkerFaceColor','m')
hold on 
plot(cd_fwt4,cl,'k-s','LineWidth',0.8,'MarkerFaceColor','g')
hold on 
plot(cd_fwt5,cl,'r-s','LineWidth',0.8,'MarkerFaceColor','r')
hold on 
plot(cd_fwt6,cl,'k-s','LineWidth',0.8,'MarkerFaceColor','y')

xlabel('Drag coefficient $C_D$','Interpreter','latex','FontSize',12)
ylabel('Lift coefficient $C_L$','Interpreter','latex','FontSize',12)
set(gcf,'Color','w')
legend('No FWT','Fold length 15 $\%$','Fold length 20 $\%$','Fold length 25 $\%$','Fold length 30 $\%$','Fold length 35 $\%$','Fold length 40 $\%$','Interpreter','latex','FontSize',12)



% range estimation
pload=13000;
fmass=20000;

OEW=[51636.8,51137.24,50825.6,50710.6,50416.6,50283.38,50018.756];

TOW=OEW++pload+fmass;
Cl=TOW*9.81/(9700*154.1);

k=[0.076069,0.09684,0.09882, 0.1, 0.1069, 0.1091, 0.1156];
Cd=k.*Cl.^2 + 0.0174;

% range 

V=230*3.6;
Isp=1/0.565;
CLCD=Cl./Cd;

range=V.*(CLCD).*Isp.*log(TOW./(TOW-fmass));
Eta=[0, 0.15, 0.2, 0.25, 0.3, 0.35 ,0.4];

figure 

plot(Eta,range,'bs')

figure 
plot(Eta,CLCD,'bs')

figure 
plot(Eta,Cl,'bs')

figure 
plot(Eta,Cd,'bs')
























