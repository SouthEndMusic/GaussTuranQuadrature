# This file was automatically generated with rule_table_generation.jl
# (n, s) => (; X, W)
rules = Dict{Tuple{Int, Int}, @NamedTuple{X::Vector{Float64}, W::Matrix{Float64}}}()

rules[(1, 0)] = (; X = [0.5], W = [1.0;;])

rules[(1, 1)] = (; X = [0.5], W = [1.0 0.0 0.041666666666666664])

rules[(1, 2)] = (; X = [0.5],
    W = [1.0 0.0 0.041666666666666664 2.5679065925163143e-34 0.0005208333333333333])

rules[(1, 3)] = (; X = [0.5],
    W = [1.0 0.0 0.041666666666666664 2.5679065925163143e-34 0.0005208333333333333 2.6080301330243817e-35 3.1001984126984127e-6])

rules[(1, 4)] = (; X = [0.5],
    W = [1.0 0.0 0.041666666666666664 2.5679065925163143e-34 0.0005208333333333333 2.6080301330243817e-35 3.1001984126984127e-6 1.463882298224022e-36 1.076457782186949e-8])

rules[(2, 0)] = (; X = [0.2113248654051871, 0.7886751345948129], W = [0.5; 0.5;;])

rules[(2, 1)] = (; X = [0.1853944358250453, 0.8146055641749547],
    W = [0.5 0.024072942084497444 0.0036626496067172754;
         0.5 -0.024072942084497444 0.0036626496067172754])

rules[(2, 2)] = (; X = [0.1749861232620853, 0.8250138767379147],
    W = [0.5 0.03221761516310786 0.004896000318693137 0.00016374761828368364 6.92710342429444e-6;
         0.5 -0.03221761516310786 0.004896000318693137 -0.00016374761828368364 6.92710342429444e-6])

rules[(2, 3)] = (; X = [0.16924268152283997, 0.8307573184771601],
    W = [0.5 0.036465277184756055 0.005544389700939244 0.00025931431031748803 1.3480776074448803e-5 3.039463358775756e-7 5.876818996305826e-9;
         0.5 -0.036465277184756055 0.005544389700939244 -0.00025931431031748803 1.3480776074448803e-5 -3.039463358775756e-7 5.876818996305826e-9])

rules[(2, 4)] = (; X = [0.16555809445557323, 0.8344419055444268],
    W = [0.5 0.03910972351103327 0.00595041674363215 0.0003214045378443053 1.8321084404710193e-5 5.893337341791584e-7 1.6477030801634157e-8 2.5631763129986986e-10 2.817402174975488e-12;
         0.5 -0.03910972351103327 0.00595041674363215 -0.0003214045378443053 1.8321084404710193e-5 -5.893337341791584e-7 1.6477030801634157e-8 -2.5631763129986986e-10 2.817402174975488e-12])

rules[(3, 0)] = (; X = [0.11270166537925831, 0.5, 0.8872983346207417],
    W = [0.2777777777777778; 0.4444444444444444; 0.2777777777777778;;])

rules[(3, 1)] = (; X = [0.0927804072111184, 0.5, 0.9072195927888816],
    W = [0.2666582029608889 0.007791166649283887 0.0005134350911574401;
         0.46668359407822224 -1.1928679612878346e-21 0.002765885622271981;
         0.2666582029608889 -0.007791166649283887 0.0005134350911574401])

rules[(3, 2)] = (; X = [0.08538363096536741, 0.5, 0.9146163690346326],
    W = [0.26209446987372836 0.010065994810188755 0.0006982108470301828 1.2694196378915732e-5 2.2896517025727294e-7;
         0.47581106025254327 -1.6483955960009884e-21 0.003561493233522819 -7.006319150171902e-24 4.577443345063206e-6;
         0.26209446987372836 -0.010065994810188755 0.0006982108470301828 -1.2694196378915732e-5 2.2896517025727294e-7])

rules[(3, 3)] = (; X = [0.08145057382071315, 0.4999999999999997, 0.9185494261792867],
    W = [0.25957893901548573 0.011159029564292005 0.0007903103903145268 2.014235145117148e-5 4.839703228130958e-7 5.480256798303416e-9 4.4783524350585986e-11;
         0.4808421219690281 2.839912633057916e-17 0.00395327795585452 1.7132937809901752e-19 7.156927341541521e-6 1.724200518853433e-22 3.4309540649408188e-9;
         0.2595789390154862 -0.011159029564292038 0.000790310390314531 -2.014235145117162e-5 4.839703228131002e-7 -5.4802567983034775e-9 4.478352435058662e-11])

rules[(3, 4)] = (; X = [0.07898410779691255, 0.5000000000011252, 0.9210158922037434],
    W = [0.25797305116428926 0.011805263632966422 0.000845447121288389 2.4875907998806522e-5 6.765311769338742e-7 1.1065893231399942e-8 1.4289225389355437e-10 1.0600887936249058e-12 4.896985625968374e-15;
         0.48405389767316537 -1.1291705057005254e-13 0.004189326001808226 -7.896865294136057e-16 8.790394961329108e-6 -1.197263704637439e-18 6.347179267661228e-9 -4.717585919988653e-22 1.4591354805982023e-12;
         0.2579730511625454 -0.011805263632833213 0.0008454471212718469 -2.4875907998174097e-5 6.765311769112156e-7 -1.1065893230950558e-8 1.4289225388650317e-10 -1.0600887935637568e-12 4.8969856256339735e-15])

rules[(4, 0)] = (;
    X = [0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262],
    W = [0.17392742256872692; 0.32607257743127305; 0.3260725774312731; 0.17392742256872692;;])

rules[(4, 1)] = (;
    X = [0.05516159319722051, 0.320570736480558, 0.679429263519442, 0.9448384068027795],
    W = [0.16234344880949175 0.003004382093618984 0.00011246781956587063;
         0.33765655119050825 0.002547668023679286 0.0010166945613728222;
         0.33765655119050825 -0.002547668023679286 0.0010166945613728222;
         0.16234344880949175 -0.003004382093618984 0.00011246781956587063])

rules[(4, 2)] = (;
    X = [0.05008539366202137, 0.3170378226205417, 0.6829621773794601, 0.949914606337979],
    W = [0.15780210299761893 0.0037947981804902927 0.00015174706652647656 1.6712733952043003e-6 1.695760992815835e-8;
         0.34219789700238207 0.00338902338414494 0.001310020479917178 7.050156884234488e-6 8.238823596806355e-7;
         0.3421978970023811 -0.003389023384144983 0.0013100204799171673 -7.050156884234537e-6 8.238823596806235e-7;
         0.15780210299761785 -0.0037947981804902454 0.0001517470665264734 -1.6712733952042547e-6 1.695760992815774e-8])

rules[(4, 3)] = (;
    X = [0.04743808333375109, 0.3151668533577375, 0.6848331411116344, 0.9525619155030799],
    W = [0.15535210800252497 0.004154246081403507 0.00017041321735599316 2.627383892939806e-6 3.641320920743837e-8 2.426801057640091e-10 1.109738023343413e-12;
         0.34464788877091745 0.003813802350927886 0.0014517271548593765 1.134146913228312e-5 1.3159507620457354e-6 5.459161013126903e-9 3.001240286560581e-10;
         0.34464789182592864 -0.003813802200069096 0.0014517271924342479 -1.1341468886712821e-5 1.3159508198568138e-6 -5.459160999292307e-9 3.001240484416567e-10;
         0.15535211140062893 -0.004154246247962014 0.00017041322842726287 -2.6273841184074764e-6 3.641321325326741e-8 -2.426801384235476e-10 1.1097382041138087e-12])

rules[(4, 4)] = (;
    X = [0.04705443223503832, 0.3149998747601893, 0.6795827071732133, 0.9525168207619494],
    W = [0.15658128378790317 0.004462824183290979 0.00019046777300360236 3.4321403342771256e-6 5.573196077361811e-8 5.498095053264564e-10 4.205242038349158e-12 1.8553301040252444e-14 4.973827900952276e-17;
         0.3394880139896599 0.003492969003149812 0.0014420771063855062 1.1546387579010252e-5 1.462494865125023e-6 8.080732443170292e-9 4.943482751098804e-10 1.422298278758843e-12 5.055483635834853e-14;
         0.34533931133533124 -0.003408697707584355 0.0015201213060068491 -1.1741702033401443e-5 1.6037012113108796e-6 -8.610707167114485e-9 5.674160648626245e-10 -1.6000125051170713e-12 6.12558627000696e-14;
         0.15859139088710572 -0.004601212314057005 0.0001982217246646686 -3.624354381567177e-6 5.942912389715918e-8 -5.93134901957444e-10 4.573497105526363e-12 -2.035583341793736e-14 5.478263877690098e-17])

rules[(5, 0)] = (;
    X = [0.046910077030668004, 0.23076534494715845,
        0.5, 0.7692346550528415, 0.953089922969332],
    W = [0.11846344252809454; 0.23931433524968324; 0.28444444444444444; 0.23931433524968324; 0.11846344252809454;;])

rules[(5, 1)] = (;
    X = [0.03644106519505467, 0.21956629041918166,
        0.5, 0.7804337095808184, 0.9635589348049454],
    W = [0.10847688585191251 0.0013652858453592012 3.3076279639323965e-5;
         0.2439303314354623 0.0018767226274142023 0.0003778176406532087;
         0.2951855654252504 -3.523213392045612e-22 0.000669563074180224;
         0.2439303314354623 -0.0018767226274142023 0.0003778176406532087;
         0.10847688585191251 -0.0013652858453592012 3.3076279639323965e-5])

rules[(5, 2)] = (;
    X = [0.03440521217828337, 0.22452684217570137,
        0.5152621008258988, 0.797061551804066, 0.970048702452686],
    W = [0.10951046770130263 0.0018540826678916126 5.062613028683343e-5 3.8752411344887924e-7 2.659474629212364e-9;
         0.2538113095899912 0.002510907278064052 0.0005335304351891977 2.82921804644202e-6 1.7994804976841152e-7;
         0.30166608913789916 -0.00027854473667954113 0.0008841510776668243 -4.452587792835811e-7 4.259185302493216e-7;
         0.23818391136917175 -0.0025720620039383587 0.0004418352713474075 -2.5439892153134014e-6 1.2982923007343678e-7;
         0.09682822220163526 -0.0014819980883617814 3.4969376853553514e-5 -2.376504181728702e-7 1.387378797113466e-9])

rules[(5, 3)] = (;
    X = [0.02825500088663886, 0.20825318224052464,
        0.5023775664386808, 0.7950674944896856, 0.9724894207146948],
    W = [0.09723717049085707 0.0017409454849219394 4.2611276731710764e-5 4.2167576115494173e-7 3.468210216596606e-9 1.4185546689084742e-11 3.615349117398181e-14;
         0.24992427475830026 0.003076919771771085 0.0005582609076252262 4.774174060943447e-6 2.6514753901144954e-7 1.1818530859755733e-9 3.030846942246323e-11;
         0.3100575921679976 -5.141376880008349e-5 0.001040837928090127 -1.2303270150613418e-7 7.478076022771542e-7 -4.74360396742336e-11 1.3741495261535453e-10;
         0.247526066343117 -0.003077532790427243 0.0005427108371142863 -4.682329787703729e-6 2.528339666897305e-7 -1.1356345357431237e-9 2.82668045081824e-11;
         0.09525489623972808 -0.0016841674353872493 4.0184829798744325e-5 -3.90974800099482e-7 3.1340684161487328e-9 -1.2543092168775094e-11 3.0908437149820515e-14])

rules[(5, 4)] = (;
    X = [0.17121384607564052, 0.34140071709618397,
        0.5077958243722474, 0.6742235254466961, 0.8378587502118069],
    W = [2.7280464760882015e9 1.7017870728331098e8 4.818022058382893e6 80954.49801758672 884.1371700459172 6.436825046436324 0.030562300985927354 8.671895755548605e-5 1.1291164334608655e-7;
         3.2450468936650293e12 2.3417351712959567e11 7.740191209077458e9 1.5536444212495595e8 2.0217220873073249e6 18103.727686479953 96.55527791664251 0.318243871599901 -0.00041514890244388205;
         -3.656128882247107e12 2.5046489214073236e11 -1.713706068346678e10 2.729296873197766e8 -1.0656948580568925e7 63900.70179731948 -1884.9534768679384 3.8428248287647206 -0.09571554941112798;
         4.0394279773399884e11 -1.2355956712795467e10 -3.6312112194645196e8 2.5392182547216002e7 -723394.5772569334 10860.726943122076 -120.43949806708451 0.7197015948878546 -0.0034354538223852053;
         4.411144372990494e9 -2.5772040826849833e8 6.82648070814943e6 -107186.20455514971 1092.4625124916706 -7.411157254185841 0.032731469136768776 -8.621164507332034e-5 1.0394151176357649e-7])

rules[(6, 0)] = (;
    X = [0.03376524289842399, 0.16939530676686773, 0.38069040695840156,
        0.6193095930415985, 0.8306046932331322, 0.966234757101576],
    W = [0.08566224618958518; 0.1803807865240693; 0.23395696728634552; 0.23395696728634552; 0.1803807865240693; 0.08566224618958518;;])

rules[(6, 1)] = (;
    X = [0.026207026018493763, 0.16019721281018487, 0.37705621611733403,
        0.6236932001352493, 0.840621409958848, 0.9740211008509931],
    W = [0.07834778258529725 0.000716695678358457 1.2392973655668773e-5;
         0.1820072388861997 0.001174479525138785 0.0001557756259355574;
         0.23998104050295524 0.0005562401435952383 0.0003566369994885021;
         0.2401743479161354 -0.0005532919508455062 0.0003576630897364555;
         0.18171686279589314 -0.0011855122785907797 0.0001550209967998752;
         0.0777727273135193 -0.0007077394802436498 1.2102695076433981e-5])

rules[(6, 2)] = (;
    X = [0.025301971455268742, 0.16782024385157138, 0.3965221292912325,
        0.6374763507503035, 0.8428507561805891, 0.976056264797032],
    W = [0.08098587994753706 0.0010218472322094066 2.0456489435468772e-5 1.1588337810366034e-7 5.800175881531553e-10;
         0.19402979789841426 0.0016417648314118955 0.00023834296104686972 1.0745217966744247e-6 4.6429585290263225e-8;
         0.24563729998490813 0.0002680005835387971 0.00047673843429148286 2.823998585162851e-7 1.5176021595366036e-7;
         0.22495565046260488 -0.0008118194180335687 0.00036032951249129257 -7.003493382143163e-7 9.118442805352283e-8;
         0.17796648574010404 -0.0012230243615208424 0.00018410525780149985 -6.792943457481126e-7 3.0662855531693405e-8;
         0.07642488596643161 -0.0009067074638113536 1.7206094628854568e-5 -9.199282478701585e-8 4.3748155671223256e-10])

rules[(6, 3)] = (;
    X = [0.14083435526219712, 0.28350317678383097, 0.4282414523705482,
        0.5714917379182543, 0.7149799257476332, 0.8585924145860061],
    W = [1.391800959228562e7 701798.1981979774 15357.04885400442 187.01800865965333 1.339805648996253 0.005369247479566608 9.442357285896151e-6;
         5.042908445130297e10 2.869251106635573e9 7.241475019945417e7 1.0491796467197814e6 9161.633181948922 47.11309985047666 0.10010095846178606;
         -1.5414792543998848e12 -1.0057171503307417e11 -3.6004227703289905e9 -6.2033593510675095e7 -950131.0924769781 -6109.068590502996 -49.001437224098844;
         1.6690595997076882e12 -9.378317355601671e10 4.254870206767387e9 -6.098547986677406e7 1.26890696863112e6 -6408.568896169139 79.02264395630728;
         -1.7808459177087167e11 1.0734859299082018e10 -2.9245306478985786e8 4.632689719280404e6 -46439.455143275205 280.1558647651238 -0.9446623900954695;
         6.1244003172976926e7 -2.82145072064911e6 56411.90025473733 -627.5291174205253 4.103828238897965 -0.01499449220262126 2.3991576857028858e-5])

rules[(6, 4)] = (;
    X = [0.14402608943378706, 0.2878555758048066, 0.4327398720252977,
        0.5749248882170597, 0.7147645050236421, 0.8565018783685162],
    W = [-1.3665691918087051e10 -7.274402092068454e8 -1.753879518455488e7 -250302.13351283697 -2313.946282896039 -14.197410374360544 -0.05648985503997976 -0.0001333344369174792 -1.4299647840899878e-7;
         3.376753703404257e13 2.0022635141386694e12 5.459075782209691e10 8.898417775740446e8 9.628694153646288e6 69601.52545028849 347.7587062669465 1.008782147624832 0.0022075887963500197;
         -5.608555034321565e13 7.577077073601453e11 -2.1728995435832617e11 1.895625899705514e8 -1.1651310969419076e8 -105545.38000855535 -17892.769481035746 -14.202772071990921 -0.7777724899450001;
         2.4957645152521746e13 -7.843192186222891e11 3.418549293487724e10 5.6578648661452144e7 1.3278329783074211e7 16164.823271620116 3111.8146812094574 -4.16610546518103 0.1950874895791952;
         -2.635202222856013e12 1.0465885168512354e11 -1.5470635356993377e9 2.5294475655838125e6 127837.95123276798 -2034.6791799673151 -2.9377362978082595 0.06756913147008639 -0.0011306054754154076;
         9.236071426429087e9 -4.681861817100971e8 1.0743361234244796e7 -145880.94225151965 1283.2626469354482 -7.496372372041548 0.028434468612872868 -6.412646883891076e-5 6.596173517799792e-8])

rules[(7, 0)] = (;
    X = [0.025446043828620736, 0.12923440720030277, 0.2970774243113014,
        0.5, 0.7029225756886985, 0.8707655927996972, 0.9745539561713793],
    W = [0.06474248308443485; 0.13985269574463832; 0.19091502525255946; 0.2089795918367347; 0.19091502525255946; 0.13985269574463832; 0.06474248308443485;;])

rules[(7, 1)] = (;
    X = [0.015841152733811873, 0.11018496792070888, 0.28372940900939736,
        0.5011161161885772, 0.7165343111657518, 0.886854367406236, 0.9824129229743549],
    W = [0.05033559609066881 0.0003257599148959591 3.109864541625732e-6;
         0.1374264503735753 0.0008580810997640591 6.671928831773264e-5;
         0.20161948269853125 0.0006875068492781328 0.000211095494328221;
         0.22299328789151043 -3.763847998971467e-5 0.00028542726248230026;
         0.1982917070796926 -0.000700067599375663 0.00020045945109222656;
         0.13574531458639416 -0.0007939584849990602 6.407947172670347e-5;
         0.053588161279627455 -0.0003451322238580692 3.8793720419135074e-6])

rules[(7, 2)] = (;
    X = [0.0293481556098766, 0.0968669579818669, 0.25792339160695416, 0.5019003817731398,
        0.7427169526566737, 0.9040711135206501, 0.9731764858089322],
    W = [0.7465653085791529 0.018882765926061168 0.00022288434549367234 1.216491032742204e-6 3.388120920028985e-9;
         -0.6147529348703281 0.026323806823299056 -0.0005513434956560656 4.044232522506001e-6 -5.4295883895935644e-8;
         0.24556537030772269 -0.0003016046535451338 0.000426502765725777 1.3160802095986334e-7 1.1273635997731327e-7;
         0.24841973996660094 -6.41443610540175e-5 0.0004897491772918167 -6.747111599591425e-8 1.567031265345041e-7;
         0.22987528768776824 -0.0004950761522660313 0.00037169182758489916 -5.746185741425365e-7 9.403862478344469e-8;
         -0.21987668659693063 -0.012452246988201243 -0.0002471003350377893 -1.9627606856351604e-6 -2.746191710348599e-8;
         0.364203914926014 -0.008656382551474434 0.0001063205238552509 -5.762488024539545e-7 1.6927752368916777e-9])

rules[(7, 3)] = (;
    X = [0.12636402238444303, 0.25171043137344595, 0.3773231405279579,
        0.49850829789154505, 0.6272272918015391, 0.7518677629887717, 0.8793596222357084],
    W = [-3.2876144288899636e8 -1.2748182133736044e7 -213402.01313130456 -1974.6272915886725 -10.654118449229363 -0.031780065522240386 -4.092309633883123e-5;
         1.5848949225587512e12 8.304327015733109e10 1.9597752547624094e9 2.684337497509509e7 230869.1832854966 1193.4286577167427 3.334638172566183;
         -3.058247545205526e12 -7.52790184311041e10 -1.2594526277974642e10 -1.1973071168317963e8 -4.408795192978429e6 -16262.247524311668 -283.74096598615137;
         7.163438400003379e11 -6.227570552656937e10 -7.068494364216626e9 4.9959377888273135e7 -4.01283893724519e6 11155.005294606775 -345.3496710701244;
         7.325006889251498e11 -2.488824643362975e10 -1.9139162526711583e8 1.1921028368587838e7 -370112.8475874759 2826.908722944568 -27.550778303505886;
         2.4803289843491722e10 -7.577239349718779e8 3.5257984221828883e6 161820.16770523018 -3555.180932218658 30.016603902421615 -0.1370355516808071;
         3.3565321684413545e7 -1.382179004301502e6 24672.27704459596 -244.69076483223188 1.4243865022138433 -0.004624268478586563 6.5610419652819985e-6])

rules[(7, 4)] = (;
    X = [0.12482432404215345, 0.24959904268238084, 0.375823696745915, 0.5014929541890132,
        0.6263125220454192, 0.7523369792928704, 0.8797588506762075],
    W = [-3.4073251113093616e10 -5.073988918443246e8 1.0770194848073358e7 386795.9756670328 5046.504620185429 36.76494768368679 0.15967854094689438 0.00039157701992818356 4.226497006914564e-7;
         -1.033021435533527e14 -1.760315378181535e12 7.291723475045097e10 3.2881770530664105e9 6.1584126771609336e7 682067.2945164114 4942.238742824353 21.820038262614514 0.05658933234979985;
         1.1708905336287764e14 -9.96947500658793e12 6.192681111623231e11 -1.0856028054642136e10 3.8609572441611266e8 -1.840907755315186e6 63646.74467267885 -58.418147370145135 2.82478800181263;
         -1.0604347933458299e13 1.4753632072077551e12 6.04722879510074e10 5.541506434247907e9 5.9990093369520515e7 2.135410262047812e6 16363.211332957197 146.7984950819364 1.1349339506829421;
         -4.874149891569442e12 -8.554404919150209e11 5.393295376605481e9 2.6694217179182476e8 -3.3327822700147565e7 468821.4295569618 -7659.119370156147 39.02900549326592 -0.280135941058475;
         1.723713620322928e12 -3.378581783204494e10 -5.133560836621591e8 2.521007515600202e7 -403982.25281950284 3483.3176521834193 -19.07159518915512 0.055211028970167136 -0.00011959544771427846;
         1.947646293970534e9 -9.076215285347638e7 1.913359900283023e6 -23847.623647081353 192.3468773491461 -1.028946134868428 0.003568861044372065 -7.348036644405632e-6 6.8888420289048395e-9])

rules[(8, 0)] = (;
    X = [0.019855071751231884, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751,
        0.591717321247825, 0.7627662049581645, 0.8983332387068134, 0.9801449282487681],
    W = [0.05061426814518813; 0.11119051722668724; 0.15685332293894363; 0.181341891689181; 0.181341891689181; 0.15685332293894363; 0.11119051722668724; 0.05061426814518813;;])

rules[(8, 1)] = (;
    X = [0.01775618848344525, 0.1058052049316385, 0.2088129972887829, 0.37729659989303743,
        0.5730495913548818, 0.7555859575867768, 0.8996054226840804, 0.9838665220696268],
    W = [0.05295279675744153 0.00032667178575170486 3.843876645473482e-6;
         0.10613681568019244 5.017981590708254e-5 3.540510132642286e-5;
         0.12551271732749114 0.001245422517758422 5.2174797911006e-5;
         0.19028203081276843 0.0003237178518762323 0.00017790920946568709;
         0.1932625343021204 -0.00020926155406674295 0.0001849041203581214;
         0.16647570812537876 -0.0004908511936645379 0.00011806352825294597;
         0.11686665654067875 -0.0005379770066001594 4.09571140737481e-5;
         0.04851074045392854 -0.0002771146334314465 2.920746457034375e-6])

rules[(8, 2)] = (;
    X = [0.04937554101723252, 0.24057716769644547, 0.3828674535343338, 0.4910950198059525,
        0.5993771731295043, 0.6804879049398018, 0.748979362239717, 0.9075119225789154],
    W = [0.6143366602585483 0.03114232281221078 0.0007762850481528883 8.25303240470379e-6 4.348618887536117e-8;
         520430.8754962269 21825.45305133542 365.92178286268256 2.9259873489476447 0.009491698459205097;
         -3.555213326998642e9 -1.4703078904600456e8 -2.54271187753711e6 -22211.692113317466 -90.92997465679207;
         1.0133333755923878e11 4.485602997450006e9 1.0701692056703666e8 1.0068344656357617e6 8930.154897074784;
         -1.303506913426187e12 -4.475164092232371e10 -1.1126463657013495e9 -7.524303722052718e6 -89547.38871943261;
         1.3702114895940303e12 -4.001853586529572e10 9.584116328310881e8 -5.904956445177471e6 55655.0221865558;
         -1.6448323178034824e11 4.102815153773301e9 -4.162880233363787e7 205528.5324857647 -422.98367759656907;
         10949.774787200895 -398.94856323058985 5.7125752272561 -0.03828631657309529 0.00010235600311380764])

rules[(8, 3)] = (;
    X = [0.11621079188053773, 0.22584643331590482, 0.3376350183964973, 0.45295971096873067,
        0.5640576927876456, 0.6764019404751812, 0.7865994663541034, 0.8974911545759012],
    W = [-2.797486581520437e9 -9.569530073693766e7 -1.4144891368938382e6 -11570.662742386019 -55.27640137820975 -0.14630533932645243 -0.00016769709021474444;
         9.550000059671652e12 4.490093643359025e11 9.520862786830544e9 1.1718981347006989e8 908507.6179925709 4230.294962454733 10.80149592943708;
         -1.1055320762840953e13 2.1068476378728076e11 -3.932882120407685e10 -1.3967872454649884e8 -1.2213534450872963e7 -26479.649305752944 -703.8575964114632;
         -4.2027899856606494e12 -2.5505633327787534e11 -3.523335608544714e10 4.8607811679960206e7 -1.4385422101493107e7 16299.274431056374 -1035.2433057956723;
         5.782620704393232e12 -1.4220210129293463e11 -2.6687683689276057e8 6.414705720958276e7 -1.0885611699350683e6 9644.058806889128 -39.52406523625063;
         -7.779579040001065e10 3.7019887250582323e9 -1.5070371075708017e8 1.835232219616193e6 -43896.25750189282 203.15789971011978 -2.983848037403842;
         6.060939257678886e9 -8.960852140320647e7 -2.90694839029288e6 99513.7932228104 -1387.9010449387292 9.428351257267726 -0.03717554985975623;
         2.232216156983294e7 -800386.2957367868 12434.985416424723 -107.28084760339475 0.542902435983477 -0.001531075808575346 1.8854204433083025e-6])

rules[(8, 4)] = (;
    X = [
        0.11110357925609274, 0.22222085226237678, 0.33339816848895326, 0.44452474306613327,
        0.5556421302811534, 0.666780124373161, 0.7778536246824319, 0.8889839153549745],
    W = [1.1194564761980622e14 4.466936499300066e12 8.07109413600653e10 8.631063538932526e8 5.979164208879455e6 27497.389524333226 82.04726200725761 0.14534710312237267 0.00011714497818428859;
         -1.471343830376712e14 7.313788418616557e12 -1.077274782302173e12 -3.3914532838180153e10 -1.0139913346952852e9 -1.1249752479428759e7 -120127.1137539061 -526.5317134817392 -2.493030105342485;
         -9.089916663922666e13 4.503141108293013e12 2.6549533195893115e12 3.3593546259251793e10 3.827044593646768e9 1.3684312130769828e7 787681.023101893 1061.1951780964375 35.791799509697206;
         3.2906296333356275e14 -5.510329115010372e12 4.116589038947732e11 3.446256495762621e10 -5.679523777871045e8 9.61544102154318e6 -5614.114283860519 208.68308704208795 6.537731426484479;
         -2.589110439814556e14 -4.8091367638678e10 -2.5643063461816653e11 -1.8760177859691326e10 3.6859283830079064e7 -6.8006664110514065e6 7636.069715388714 -435.25950132855667 -0.47191846056963693;
         5.59906029192761e13 1.2520219763825545e11 -3.025564156037095e10 -1.7467832122076157e8 1.9702813745738268e7 -421619.3878612942 3858.577335863672 -24.27595964444314 0.051105646381587044;
         -5.5207984880617e10 8.03293906894343e9 7.854817155925852e7 -1.2398794280636638e7 323019.1727417581 -4307.876434514333 35.134066625662626 -0.16643679676464881 0.00044664303116121726;
         5.877705900135345e8 -4.430550029251748e7 1.1709743889515381e6 -16388.126306450573 139.46288185005747 -0.755982884466262 0.002582101085119266 -5.1235482654814715e-6 4.550516686184704e-9])

rules[(9, 0)] = (;
    X = [0.02057769382697001, 0.10421050026091837, 0.23853872138584103,
        0.39867313485382416, 0.5502486502919993, 0.6733169132381079,
        0.7989251363951934, 0.9120438771519325, 0.9827209048161779],
    W = [0.0523169606217991; 0.11239191190353619; 0.15211443151375162; 0.16239343232497178; 0.13520866801078665; 0.12040612699572774; 0.12589493419753972; 0.09527240775009659; 0.04400112668179061;;])

rules[(9, 1)] = (;
    X = [0.00930147028627406, 0.06811890283176837, 0.17612292089036016,
        0.31270352155272274, 0.47741025333844345, 0.6525201084494693,
        0.8095265718352896, 0.9260033329388501, 0.9887197918987861],
    W = [0.03052991198961334 0.00012605387247418466 6.82164784495857e-7;
         0.08627276316936303 0.00033793862803564126 1.658133581663282e-5;
         0.12400636675454435 0.00022932335016979168 4.8806851110263796e-5;
         0.1507068305084014 0.00037002602495473785 8.659116652344526e-5;
         0.17410345323921453 0.0001414949054526844 0.00013553026691511924;
         0.16975406349592514 -0.00024176857849049396 0.00012536845723645884;
         0.13935053723686308 -0.0004384731832518654 6.924309905265752e-5;
         0.09065437873988604 -0.00037878610030061134 1.8994149159237965e-5;
         0.03462169486618913 -0.00014559623826917775 1.037595388485697e-6])

rules[(9, 2)] = (;
    X = [0.09956243574752437, 0.19926400129243677, 0.299402326028075,
        0.3995778698951973, 0.49930365733990517, 0.5999291973070748,
        0.6999736138390834, 0.7998889361558479, 0.9007132452681615],
    W = [164498.07649853936 4819.323040438995 55.94372174288774 0.3060916865052064 0.0006701237425029998;
         -2.216448366663717e8 -7.996323768903481e6 -121567.59968096482 -930.6276772519794 -3.5304742401462703;
         4.228129522746808e10 1.6528630992603614e9 2.8434707070224836e7 238041.4228401048 1168.9539883845646;
         -9.135516386265774e11 -3.718567137526496e10 -8.635254013782058e8 -7.403737572493175e6 -65112.89890034294;
         2.1996932517766445e12 2.2114880462105637e10 2.791469052567297e9 5.2727511164224055e6 326099.8270806614;
         -1.5239043819434163e12 5.301593094739925e10 -1.684360834257054e9 1.1696535438343162e7 -162565.73194373018;
         1.9762778680494382e11 -8.112868115096329e9 1.5243971618661737e8 -1.355196889980091e6 7882.968515242498;
         -1.9253303010707018e9 6.81326230016e7 -998288.1233914937 7282.6522773359875 -24.51947942161411;
         497401.59788144147 -13754.92931670274 149.66353326979828 -0.7614975857371541 0.0015354324363939817])

rules[(9, 3)] = (;
    X = [0.09998548426763884, 0.19998401002357818, 0.2999906455608416,
        0.39998877783126263, 0.49998254437481937, 0.599992793709389,
        0.699987853134525, 0.7999849607639696, 0.8999920255965876],
    W = [1.6121907547109125e10 5.0808685416103643e8 6.915798830848382e6 52063.5567670159 228.73037208741704 0.5562481268508295 0.0005852335873980867;
         -8.049830516705811e12 -5.02912096839358e11 -1.2988154135387037e10 -1.8307566808688408e8 -1.5916769624068274e6 -8001.983546224669 -22.9048460192468;
         -1.497804643178219e14 -1.2800737207676658e13 -3.931556667894387e11 -6.9143208555190525e9 -7.919833984230818e7 -489025.22486122686 -2528.024301877015;
         2.319314524076298e14 -1.7701748470135261e12 -1.1220778261604536e11 1.1113570520813346e9 -1.1876242034723832e8 151765.28217950618 -9332.803227464807;
         -4.412610516652438e13 8.209964060225307e12 -9.37075471774712e10 4.19638799193567e9 927805.5438154216 256370.98454873636 1642.5285349717076;
         -3.039492516653858e13 6.692890087015717e11 1.5924400831638896e9 -2.6216534028650296e8 4.4112004280984895e6 -32807.67832444433 183.10166599503043;
         3.722297246619063e11 -2.136596338013129e10 1.3947057019011012e8 -1.2091710188198732e6 -118064.8786181477 626.2856571478371 -10.057036361313099;
         3.1418173552911346e10 -3.469321535282782e8 -1.489519740277935e7 421448.0646246944 -5113.049424605397 30.83023348314233 -0.10566287192092223;
         1.0295419995481846e8 -3.295127146360694e6 45642.803203972864 -350.594178261146 1.5770312759510892 -0.0039450893688020355 4.297543731043648e-6])

rules[(9, 4)] = (;
    X = [0.10001505836566953, 0.20092476962479977, 0.30207437709459567,
        0.4024370996621896, 0.5014809977923621, 0.603716659247225,
        0.7019272263540461, 0.8006580663409686, 0.9010418696453107],
    W = [1.6863737481481047e13 7.11338774612948e11 1.3561726590028292e10 1.5280768886303285e8 1.1141676352565242e6 5388.623746733286 16.899065100262483 0.03144973883004197 2.6618738699452154e-5;
         4.971380720928628e14 7.533369798890017e12 -2.7420357945678607e11 -7.725850474905943e9 -1.154677657650777e8 -745591.6938169638 -3126.5864127418586 0.5128555005525919 0.037415356524784296;
         -1.9613749167005362e15 -2.7666422336976582e13 -3.75036574266685e12 -1.9578477193555036e9 -1.1625961541984072e9 -83084.57914316743 -84916.27393410527 56.12840208484968 -0.75702364708872;
         2.1647373287491262e15 -1.0156954277341955e13 4.019531959246067e12 -1.4249911107800377e10 3.4897933339758265e8 -1.6666122241867155e7 -158332.8638063077 -1573.4939091098422 -11.348159709605444;
         -8.53656093799889e14 -2.446546524331762e12 2.5700431082718295e11 -1.3396865867997412e10 -1.409629007615303e8 1.0346019054875933e7 -117641.84304146402 1268.203453381615 -6.055778234549697;
         1.653459367957821e14 6.15852886649317e11 1.9116557692510933e10 1.037456505914787e10 4.421353682096979e7 1.4823623660373264e6 32631.51114043555 -27.39057632389983 2.161728283113508;
         -2.9344306340855574e13 -2.2835066651290253e11 -9.954101136947005e8 3.2401727231718993e8 -9.636000022459786e6 52834.76784600208 -349.6143779381547 -2.165846163162667 0.0013283537227849903;
         2.918230279555939e11 9.746339213213495e9 -2.5851438106206033e8 -3.0094679869655073e6 168847.7762615988 -2481.88693249556 20.63491279825864 -0.09315557986866038 0.0002461766612589245;
         -1.5813059260764732e9 3.483789146483143e7 -214236.1477839988 -1255.5129271365277 27.878147235253838 -0.1950579634982722 0.0007272239090677588 -1.4650670342279101e-6 1.270430879030874e-9])

rules[(10, 0)] = (;
    X = [0.014464270328525392, 0.06640087702397494, 0.1292033826376839,
        0.22816428211317738, 0.36880957304785944, 0.5269185259062504,
        0.6827877170410014, 0.8197084862297327, 0.9239103164231074, 0.9852659471402094],
    W = [0.03616354866452189; 0.06042993997002827; 0.07304210866278368; 0.12408958010825734; 0.15298863484203368; 0.1599933489128748; 0.14895506531136282; 0.12257865369707771; 0.08412392815065327; 0.037635191680406514;;])

rules[(10, 1)] = (;
    X = [0.013047950602590792, 0.07881161808375635, 0.14627587920307927,
        0.28218121854951245, 0.4431874968153275, 0.5961335161162916,
        0.7367561462448611, 0.8570099256460135, 0.9292861318981677, 0.9882423989833297],
    W = [0.03914446418773684 0.00018007822558481303 1.5409779615107856e-6;
         0.07270946993067252 -0.00012331039277158518 1.4052214053278504e-5;
         0.09597686667377826 0.0009423625636887865 2.872349447971619e-5;
         0.15379362857441217 0.00027491362546793526 9.35665463152031e-5;
         0.16101239083519195 -0.00012031294828115535 0.00010750293555155595;
         0.1442395361688392 -0.00013896273301540544 7.441038297471046e-5;
         0.1363061344524078 -0.00016382301961389724 6.573879478828958e-5;
         0.08984951966912685 -0.0006014894263637936 1.8929169502435608e-5;
         0.07177737609415148 -3.183694716677994e-5 1.0865067516364077e-5;
         0.035190613413682945 -0.00014508824944462273 1.1233704583387542e-6])

rules[(10, 2)] = (;
    X = [0.09090957651347137, 0.18181867942612825, 0.27272770131076596,
        0.363636803029375, 0.4545462430742926, 0.5454551099597403, 0.6363641450489766,
        0.7272726518937742, 0.818181986037891, 0.9090913714589799],
    W = [3.6916016185144023e6 89851.88666002316 855.8775914283002 3.7851616881370735 0.006566678843223363;
         -1.6531212535064499e10 -5.303071819996167e8 -7.0290213808954805e6 -46297.725033305425 -139.5258351708447;
         1.5666020749681443e12 5.922272242006088e10 1.0303030597795079e9 8.437159763086379e6 45717.02253090053;
         -6.0115182144826875e12 -1.3847916335176334e11 -6.724367511443964e9 -3.1383201510916036e7 -668085.1850138103;
         -1.7457510297067512e12 -4.709623677207471e11 -2.464385667627966e9 -9.36876416550545e7 -304933.0705534066;
         5.5555918102619e12 -3.044118524468017e11 5.206553899763597e9 -6.0786892072928965e7 405511.26993898925;
         7.363727009420698e11 -2.9682075328678143e10 9.107267503807712e8 -7.4265224202151345e6 85556.46082815964;
         -8.397198056942418e10 2.894506325915727e9 -3.7653660014580496e7 267992.3014554676 -359.2614475382738;
         -7.98307646088812e8 2.7996101631849535e7 -408522.13074075367 2965.89041550966 -10.212603965899419;
         467167.2835838014 -11848.768623177804 118.24204507597211 -0.5516664056405248 0.0010195277648492847])

rules[(10, 3)] = (;
    X = [0.09090916404062453, 0.181818249704436, 0.2727275159592682,
        0.3636368181024138, 0.454545824987883, 0.5454549147796255, 0.6363638086824235,
        0.7272726449086702, 0.8181816465485583, 0.909091283506985],
    W = [1.038077587406127e11 3.1351658089277215e9 4.102825471621039e7 298003.60055446764 1268.0146346320764 2.999199194328756 0.003083253575376709;
         7.95999681624534e13 2.6174707799028906e12 3.5709734063700005e10 2.5197812090119812e8 626411.155896279 -1490.222399652487 -21.153561291896303;
         1.0223877847196336e14 3.451609082110673e12 -3.636619304762432e11 -3.6606093630908175e9 -1.0145221152223048e8 -356391.02335845603 -4000.7104238029433;
         -6.100538042066915e14 -1.393064123080292e13 -1.6117600867312092e12 -5.728102653111093e9 -3.5761474241604596e8 -379232.2570090698 -14855.1072099858;
         5.0636942182950656e14 -1.3123959994820417e12 6.17221665356132e10 8.287867929497612e9 -4.26492281220514e7 783098.1253739216 -1883.6978412560493;
         -1.125105743303373e14 -1.1728698542483596e12 -3.232705863497267e10 -2.4740227712518787e9 7.643697260141823e6 -231275.81420637664 417.02846741983615;
         3.350810476401109e13 -8.588329805178354e11 -5.729839905480764e9 2.464021244792647e8 -6.078163273762397e6 34873.84631634885 -261.98149713001374;
         6.84755036854316e11 7.339030042455748e9 -1.2595098251082096e9 2.1136845206444424e7 -321647.4911352887 1569.3074691370598 -11.102506619349798;
         5.94239419775791e10 -1.0153932614707407e9 -7.156033013868795e6 337272.72174045356 -4123.4156254492555 23.440406955958842 -0.07348543217486538;
         1.1857152287505206e8 -3.4264426222352996e6 42821.9193830874 -296.5327911872476 1.2014102623583294 -0.0027042440361247635 2.6474488564510413e-6])

rules[(10, 4)] = (;
    X = [0.09090909102362564, 0.18181818203062813, 0.27272727295938853,
        0.3636363637515141, 0.45454545472306856, 0.5454545455380557,
        0.6363636364574026, 0.727272727219863, 0.8181818182397191, 0.9090909091061855],
    W = [-8.136496029130247e13 -1.8929500015912988e12 -1.4624777641382994e10 1.0607478509623539e7 992716.5717243744 7909.671705599158 31.45707477690163 0.06662561241550312 6.05401914349988e-5;
         6.773894720953282e14 -9.173252015900924e12 -7.6697366511308e11 -2.406646084048198e10 -9.363234807912178e7 -49354.134288219786 24516.728289547253 120.95570491889764 0.8138455981381884;
         -5.650741674249291e15 -7.909850050223064e13 -3.8580759449213403e12 8.595814164817575e10 3.5674086500763717e9 6.565880483978261e7 985349.0704652774 4776.650559722429 36.25144766783728;
         6.383855207899024e15 -1.1124687533291942e14 2.2404271395847516e13 5.95362423079036e10 1.5208651553538748e10 -2.0332733005193148e7 2.7522597957510604e6 -2812.2208125776606 118.21019808761874;
         -1.360211562250754e14 4.914037007083708e12 4.64407073665379e12 5.359145769285212e10 -3.3984174799119437e8 7.9153631060486855e6 -217907.68580066922 347.56054582763585 -16.010562922371875;
         -1.8508872633821152e15 -2.8232405070160805e13 -2.8572040768926746e11 -4.425219267338729e10 1.934233629494634e8 -3.5773381154633965e6 19400.50484212615 -55.91410751361578 0.6932790287886547;
         6.605229151922256e14 -1.7530647534241209e13 1.7989304201998288e11 4.736960498511042e9 3.080929404592648e7 721851.6648828905 17949.18846003322 -11.352531829388504 0.9957261841717919;
         -3.048807880446933e12 -7.406349026271543e11 4.2549616707343636e9 2.598674328276396e8 -4.4024241919612e6 15247.415135392052 308.98191781372293 -2.9680223429286205 0.017914830118083747;
         2.909684685503806e11 -9.688061105968662e9 4.5360066096525274e7 856987.4395072573 -16654.760743520128 24.40817818644479 0.3734957882485467 -0.005043919994235192 9.651763800202984e-6;
         5.2983731040061655e9 -1.807168042741506e8 2.7803807217464866e6 -25215.295727240242 147.51681522070194 -0.5704538203136557 0.0014250933756861336 -2.1048692745400624e-6 1.4091761698826774e-9])
