$setglobal type MINLP
$setglobal write 0

options
limrow = 0,
limcol = 0,
reslim = 2000,
optcr  = 0.000001,
SYSOUT = ON,
MINLP = BARON;

*Reading data from excel file for the input parameters in the case study
Parameter pm, pgrd, pte, pmc, puae, psfe, psdm, pflt1, pcnf, pdry1, pahy, pnt,pflt2, pcrys, pflt3, pnf, pchrm, pdry2, plbr, CON_mat, CMSTR, CMUNT;
*loading GDX from Excel
$GDXIN Data_IsoflavoneExtraction_v6_5-16-2021_addingstage4.gdx
$LOAD pm=P1 pgrd=P2 pte=P3 pmc=P4 puae=P4U psfe=P4S psdm=P5 pflt1=P6 pcnf=P7 pdry1=P8 pahy=P9 pnt=P10 pflt2=P11 pcrys=P12 pflt3=P13 pnf=P14 pchrm=P15 pdry2=P16 plbr=P19 CON_mat=Con_mat CMSTR=CMSTR CMUNT=CMUNT
$GDXIN

*Declaration of the sets for units, streams and components and subsets for the case study
Sets
         i               unit operations /split1,split2,GRD,TE,MC,UAE,SFE,mix1,CNF,FLT1,SDM,split3,mix2,split4,DRY1,AHY,NT,FLT2,split5,CRYS,FLT3,NF,CHRM,mix3,DRY2/
         i1(i)           unit ops with cost/GRD,TE,MC,UAE,SFE,CNF,FLT1,SDM,DRY1,AHY,NT,FLT2,CRYS,FLT3,NF,CHRM,DRY2/
         i2(i1)          units with Concentration Factor /SDM,FLT1,CNF,FLT2,FLT3,NF/
         i3(i)           unit ops with no cost /split1,split2,mix1,split3,mix2,split4,split5,mix3/
         i4(i1)          filter and membrane units /FLT1,FLT2,FLT3,NF/
         ibv(i)          units with binary variables/GRD,TE,MC,UAE,SFE,CNF,FLT1,SDM,DRY1,split4,split5,CRYS,NF,CHRM,mix3,DRY2/
         iv(i)           units not in extraction stage/split1,mix1,TE,MC,UAE,SFE,split2,CNF,FLT1,SDM,split3,mix2,DRY1,AHY,NT,FLT2,split5,CRYS,FLT3,NF,CHRM,mix3,DRY2/
         Nstg            number of stages/s1,s2,s3,s4/
         istg1(i1)       all units in stage-1/GRD/
         istg2(i1)       all units in stage-2/TE,MC,UAE,SFE,CNF,FLT1,SDM,DRY1/
         istg3(i1)       all units in stage-3/AHY,NT,FLT2/
         istg4(i1)       all units in stage-3/CRYS,FLT3,NF,CHRM,DRY2/

         j               streams/1,2,3,4,5,6,4U,5U,6U,4S,5S,6S,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43/

         jsplit1(j)              Splitter 1 streams/1,2/
         jsplit1F(j)             Splitter 1 feed stream/1/
         jsplit1E(j)             Splitter 1 effluent stream/2/

         jGRD(j)                 Grinding streams/2,3/
         jGRDF(j)                Grinding feed stream/2/
         jGRDE(j)                Grinding effluent stream/3/

         jsplit2(j)              Splitter 2 streams/3,4,7/
         jsplit2F(j)             Splitter 2 feed stream/3/
         jsplit2E(j)             Splitter 2 effluent stream/4,7/

         jTE(j)                  Turbo-Extraction streams/4,5,6/
         jTEF(j)                 Turbo-Extraction feed stream /4/
         jTES(j)                 Turbo-Extraction solvent stream in /5/
         jTEE(j)                 Turbo-Extraction output/6/

         jMC(j)                  Maceration streams/7,8,9/
         jMCF(j)                 Maceration feed stream /7/
         jMCS(j)                 Maceration solvent stream in /8/
         jMCE(j)                 Maceration output/9/

         jUAE(j)                 Ultrasound-Assisted Extraction streams/4U,5U,6U/
         jUAEF(j)                Ultrasound-Assisted Extraction stream /4U/
         jUAES(j)                Ultrasound-Assisted Extraction stream in /5U/
         jUAEE(j)                Ultrasound-Assisted Extraction output/6U/

         jSFE(j)                 Supercritical Fluid Extraction streams/4S,5S,6S/
         jSFEF(j)                Supercritical Fluid Extraction stream /4S/
         jSFES(j)                Supercritical Fluid Extraction stream in /5S/
         jSFEE(j)                Supercritical Fluid Extraction output/6S/

         jmix1(j)                Mixing Point 1 streams/6,9,6U,6S,10/
         jmix1F(j)               Mixing Point 1 feed stream/6,6U,6S,9/
         jmix1E(j)               Mixing Point 1 effluent stream/10/

         jsplit3(j)              Splitter 3 streams/10,11,14,17/
         jsplit3F(j)             Splitter 3 feed stream/10/
         jsplit3E(j)             Splitter 3 effluent stream/11,14,17/

         jSDM(j)                 Sedimentation streams/11,12,13/
         jSDMF(j)                Sedimentation feed stream/11/
         jSDMS(j)                Sedimentation solid phase leaving/12/
         jSDME(j)                Sedimentation liquid phase leaving/13/

         jFLT1(j)                Filtration 1 streams/14,15,16/
         jFLT1F(j)               Filtration 1 feed stream/14/
         jFLT1S(j)               Filtration 1 cake/15/
         jFLT1E(j)               Filtration 1 filtrate stream/16/

         jCNF(j)                 Centrifugation streams/17,18,19/
         jCNFF(j)                Centrifugation feed/17/
         jCNFS(j)                Centrifugation solid phase leaving/18/
         jCNFE(j)                Centrifugation liquid phase leaving/19/

         jmix2(j)                Mixing Point 2 streams/13,16,19,20/
         jmix2F(j)               Mixing Point 2 feed stream/13,16,19/
         jmix2E(j)               Mixing Point 2 effluent stream/20/

         jsplit4(j)              Splitter 4 streams/20,21/
         jsplit4F(j)             Splitter 4 feed stream/20/
         jsplit4E(j)             Splitter 4 effluent stream/21/

         jDRY1(j)                Drying 1 streams/21,22,23/
         jDRY1F(j)               Drying 1 feed stream /21/
         jDRY1S(j)               Drying 1 removed solvent /22/
         jDRY1E(j)               Drying 1 effluent stream (dried)/23/

         jAHY(j)                 Acid Hydrolysis streams/23,24,25/
         jAHYF(j)                Acid Hydrolysis feed streams/23,24/
         jAHYE(j)                Acid Hydrolysis effluent stream/25/

         jNT(j)                 Acid Hydrolysis streams/25,26,27/
         jNTF(j)                Acid Hydrolysis feed streams/25/
         jNTB(j)                Acid Hydrolysis base streams/26/
         jNTE(j)                Acid Hydrolysis effluent stream/27/

         jFLT2(j)                Filtration 2 streams/27,28,29/
         jFLT2F(j)               Filtration 2 feed stream/27/
         jFLT2S(j)               Filtration 2 cake/28/
         jFLT2E(j)               Filtration 2 filtrate stream/29/

         jsplit5(j)              Splitter 4 streams/29,30,38/
         jsplit5F(j)             Splitter 4 feed stream/29/
         jsplit5E(j)             Splitter 4 effluent stream/30,38/

         jCRYS(j)              Crystallization streams/30,31,32/
         jCRYSF(j)             Crystallization feed stream/30/
         jCRYSA(j)             Crystallization antisolvent stream/31/
         jCRYSE(j)             Crystallization effluent stream/32/

         jFLT3(j)                Filtration 3 streams/32,33,34/
         jFLT3F(j)               Filtration 3 feed stream/32/
         jFLT3E(j)               Filtration 3 filtrate stream/33/
         jFLT3S(j)               Filtration 3 product stream/34/

         jNF(j)                Nanofiltraiton streams/35,36,37/
         jNFF(j)               Nanofiltration feed stream/35/
         jNFE(j)               Nanofiltraiton filtrate stream/36/
         jNFR(j)               Nanofiltration rejection stream/37/

         jCHRM(j)               Chromatography streams/38,39,40/
         jCHRMF(j)               Chromatography feed stream /38/
         jCHRMR(j)               Chromatography removed solvent /39/
         jCHRME(j)               Chromatography effluent stream /40/

         jmix3(j)                Mixing Point 1 streams/34,40,41/
         jmix3F(j)               Mixing Point 1 feed stream/34,40/
         jmix3E(j)               Mixing Point 1 effluent stream/41/

         jDRY2(j)               Drying 2 streams/41,42,43/
         jDRYF(j)               Drying 2 feed stream /41/
         jDRYS(j)               Drying 2 removed solvent /42/
         jDRYE(j)               Drying 2 effluent stream (dried)/43/


         k               components/largesoy,soy,i-glucoside,i-aglycone,ethanol,water,CO2,HCl,Base,impurity/
         k1(k)           components in stream flows/largesoy,soy,i-glucoside,i-aglycone,ethanol,water,CO2,Base,HCl,impurity/
         k2(k)           non-soy components/i-glucoside,i-aglycone,ethanol,water,CO2,HCl,Base,impurity/
         k3(k)           solvent components/ethanol,water/
         k4(k)           SFE solvent components /ethanol,CO2/
         kP(k)           pre processing componenets/largesoy,soy,i-glucoside/
         kB(k)           initial basic components/largesoy,i-glucoside/
         kc(k)           extraction componenets/soy,i-glucoside,ethanol,water,largesoy/
         kw(k)           initial basic components/ethanol,water/
         kadd(k)         externally added components/ethanol,water,HCl,Base/
         kL4(k1)         light components in drying at stage 4 /ethanol,water/
         kSol(k1)        added solvents/ethanol,water,CO2,i-glucoside,impurity/
         kSld(k1)        solid components in/largesoy,soy/
         kNSDRY(k1)     non-solvent components in DRY/largesoy,soy,i-glucoside,i-aglycone/
         kNS(k1)         non-solvent components/largesoy,soy,i-glucoside,i-aglycone,HCl,Base,impurity/
         kFLT1(k)       relevant components in Filtration 1 /largesoy,soy,i-glucoside,ethanol,water,impurity/
         kNTSalt(k)      Components in Filtration 2 after Neutralization Reaction /i-aglycone,i-glucoside,ethanol,water/
         jI,jIin,jIout
         kI,kJ;

Alias(i,ii,iii); Alias(j,jj,jjj);Alias(k,kk,kkk);Alias(k1,kk1);Alias(kB,kkB);
Alias (kSol,kkSol);
Alias (kP,kPP,kPPP);
Alias (kB,kBB,kBBB);
Alias (kc,kcc,kccc);
Alias (kw,kww,kwww);
Alias (kFLT1,kFLT11,kFLT111);
Alias (Nstg,Nstg1);

Parameter CON(j,i),CMP1(k,i),CMP2(k,j),CMP22(j,k);

CON(j,i)=CON_mat(j,i);
CMP1(k,i)=CMUNT(k,i);
CMP22(j,k)=CMSTR(j,k);

*all stream for technology i, for when corresponding parameter is in CON(j,i) is not equal to 0.
*So whenever these streams have a value of either +1 or -1, that should be assigned to a particular technology i.

jI(j,i)$[CON(j,i)<>0]=Yes;
jIin(j,i)$[CON(j,i)=+1]=Yes;

*element -1 is out.
jIout(j,i)$[CON(j,i)=-1]=Yes;

kI(i,k)$[CMP1(k,i)<>0]=Yes;
CMP2(k,j)=CMP22(j,k);
kJ(k,j)$[CMP2(k,j)<>0]=Yes;

Parameters
*parameters for the material
Qf,g,pi,C_lbr,Rep_time,Den(k),MW(k),Cp(k),Hvap(k3),Hsub(k3),Cp_Air,C_soy,Csell_soy,C_air,C_stm,C_cwt,C_rfg,C_elec,C_water,C_ethanol,C_co2,C_acid,C_base,C_isoflavone,Tcw_in,Tcw_out,Trfg_out,Trfg_in,Tfrz,Tamb,LH_stm,Frac(k);

Qf=pm('Qf','Value');
g=pm('g','Value');
pi=pm('pi','Value');
C_lbr=pm('C_lbr','Value');
Rep_time=pm('Rep_time','Value');

Den('largesoy')=pm('Den_largesoy','Value');
Den('soy')=pm('Den_soy','Value');
Den('i-glucoside')=pm('Den_i-glucoside','Value');
Den('i-aglycone')=pm('Den_i-aglycone','Value');
Den('ethanol')=pm('Den_ethanol','Value');
Den('water')=pm('Den_water','Value');
Den('CO2')=pm('Den_CO2','Value');
Den('HCl')=pm('Den_HCl','Value');
Den('impurity')=pm('Den_impurity','Value');
Den('Base')=pm('Den_Base','Value');

*MW('largesoy')=pm('MW_largesoy','Value');
*MW('soy')=pm('MW_soy','Value');
MW('i-glucoside')=pm('MW_i-glucoside','Value');
MW('i-aglycone')=pm('MW_i-aglycone','Value');
MW('ethanol')=pm('MW_ethanol','Value');
MW('HCl')=pm('MW_HCl','Value');
MW('Base')=pm('MW_Base','Value');

Cp('largesoy')=pm('Cp_largesoy','Value');
Cp('soy')=pm('Cp_soy','Value');
Cp('i-glucoside')=pm('Cp_i-glucoside','Value');
Cp('i-aglycone')=pm('Cp_i-aglycone','Value');
Cp('ethanol')=pm('Cp_ethanol','Value');
Cp('HCl')=pm('Cp_HCl','Value');
Cp_Air=pm('Cp_Air','Value');

Hvap('ethanol')=pm('DeltaH_vap_ethanol','Value');
Hvap('water')=pm('DeltaH_vap_stm','Value');
Hsub('ethanol')=pm('DeltaH_sub_ethanol','Value');
Hsub('water')=pm('DeltaH_sub_water','Value');

C_soy=pm('C_soy','Value');
Csell_soy=pm('Csell_soy','Value');
C_air=pm('C_air','Value');
C_stm=pm('C_stm','Value');
C_cwt=pm('C_cwt','Value');
C_water=pm('C_water','Value');
C_ethanol=pm('C_ethanol','Value');
C_co2=pm('C_co2','Value');
C_rfg=pm('C_rfg','Value');
C_elec=pm('C_elec','Value');
C_water=pm('C_water','Value');
C_acid=pm('C_acid','Value');
C_base=pm('C_base','Value');
C_isoflavone=pm('C_isoflavone','Value');

Tcw_in=pm('Tcw_in','Value');
Tcw_out=pm('Tcw_out','Value');
Trfg_out=pm('Trfg_out','Value');
Trfg_in=pm('Trfg_in','Value');
Tfrz=pm('Tfrz','Value');
Tamb=pm('Tamb','Value');
LH_stm=pm('LH_stm','Value');

Frac('largesoy')=pm('largesoy_frac','Value');
Frac('soy')=0;
Frac('i-glucoside')=pm('i-glucoside_frac','Value');
Frac('i-aglycone')=0;
Frac('ethanol')=0;
Frac('water')=0;
Frac('HCl')=0;
Frac('impurity')=pm('impurity_frac','Value');

Parameters Q0(i1),Nlabr(i1),C0(i1),nc(i1),Wsp(i1);

Q0(i1)=plbr(i1,'Q0');
Nlabr(i1)=plbr(i1,'Nlabr');
C0(i1)=plbr(i1,'C0');
nc(i1)=plbr(i1,'nc');
Wsp(i1)=plbr(i1,'Wsp');

Parameters
*Parameters for stage-I (Pre-Processing)
*Grinder (Particle Size Reduction)
BWI,Y_grd, Frc_k,Tca_in_grd,Tca_out_grd,d1,d2,Theta_grd,
*Parameters for stage-II (Extraction)
*parameters for the Turbo Extraction (TE)
Tcw_in_TE,Tcw_out_TE,moist_soy_TE,P_agit_TE,Po_impeller_TE,D_soymeal,Delta_t_TE,
*parameters for the Maceration (MC)
moist_soy_MC,Delta_t_MC,
*Parameters for Ultrasound-assisted extraction (UAE)
Eff_UAE,tR_UAE,
*Parameters for Supercritical Fluid Extraction (SFE)
e_ratio_SFE,
*Parameters for solid separation following the extraction
*Sedimentation (SDM)
tT_SDM_TE,Eff_SDM,Pur_SDM,DpS_SDM,muL_SDM,
*Filtration (FLT,1)
Zeta_FLT1,CPM_FLT1,RFLT1(k),
*Centrifugation (CNF)
Eff_CNF,v_gas_CNF,Tcw_in_CNF,Tcw_out_CNF,
*Acid Hydrolysis (AHY)
Y_AHY,T_amb_AHY,T_op_AHY,Lambda_stm,phi_acid_AHY,Frc_k_AHY,
*Neutralization
theta_NT,
*Filtration 2 (FLT2)
Zeta_FLT2,CPM_FLT2, RFLT2(k),
*Thin-Film Evaporation 2 (TFE2)
Ts_TFE2,Tp_TFE2,Tf_TFE2,U_TFE2,DeltaH_vap_stm_TFE2,DeltaH_vap_ethanol_TFE2,
*Crystallization (CRYS)
phi_antisolv_CRYS,
*Filtration 3 (FLT3)
Zeta_FLT3,CPM_FLT3, RFLT3(k),
*Nanofiltration (NF)
Zeta_NF,CPM_NF,RNF(k),
*Chromatography (CHRM)
theta_CHRM,kc_CHRM,Width_CHRM,HETP_CHRM,Ratio_LD_CHRM,

*Drying 1 (DRY1)
dry_eff_DRY1,cp_rfg_DRY1,U_DRY1,
*Drying 2 (DRY2)
dry_eff_DRY2,cp_rfg_DRY2,U_DRY2

*Cost Per Membrane
*CPM(i4)
;
*Grinding
BWI=pgrd('BWI','Value');
Y_grd=pgrd('Y_grd','Value');
Frc_k=pgrd('Frc_k','Value');
Tca_in_grd=pgrd('Tca_in_grd','Value');
Tca_out_grd=pgrd('Tca_out_grd','Value');
d1=pgrd('d1','Value');
d2=pgrd('d2','Value');
Theta_grd=0.5;

*Turbo-Extraction
Tcw_in_TE=pte('Tcw_in_TE','Value');
Tcw_out_TE=pte('Tcw_out_TE','Value');
moist_soy_TE=pte('moist_soy_TE','Value');
P_agit_TE=pte('P_agit_TE','Value');
Po_impeller_TE=pte('Po_impeller_TE','Value');
D_soymeal=pte('D_soymeal','Value');
Delta_t_TE=pte('Delta_t_TE','Value');

*Maceration
moist_soy_MC=pmc('moist_soy_MC','Value');
Delta_t_MC=pmc('Delta_t_MC','Value');

*Sedimentation
tT_SDM_TE=psdm('tT_SDM','Value');
Eff_SDM=psdm('Eff_SDM','Value');
Pur_SDM=psdm('Pur_SDM','Value');
DpS_SDM = psdm('DpS_SDM','Value');
muL_SDM = psdm('muL_SDM','Value');


*Ultrasound-Assisted Extraction
Eff_UAE=puae('Eff_UAE','Value');
tR_UAE=puae('tR_UAE','Value');

*Supercritical Fluid Extraction
e_ratio_SFE=psfe('e_ratio_SFE','Value');

*Filtration,1
Zeta_FLT1=pflt1('Zeta_FLT1','Value');
CPM_FLT1=pflt1('CPM_FLT1','Value');
RFLT1('largesoy')=pflt1('RFLT1_largesoy','Value');
RFLT1('soy')=pflt1('RFLT1_soy','Value');
RFLT1('i-glucoside')=pflt1('RFLT1_i-glucoside','Value');
RFLT1('i-aglycone')=pflt1('RFLT1_i-aglycone','Value');
RFLT1('ethanol')=pflt1('RFLT1_ethanol','Value');
RFLT1('HCl')=pflt1('RFLT1_HCl','Value');
RFLT1('impurity')=pflt1('RFLT1_impurity','Value');

*Centrifugation
Eff_CNF=pcnf('Eff_CNF','Value');
Tcw_in_CNF=pcnf('Tcw_in_CNF','Value');
Tcw_out_CNF=pcnf('Tcw_out_CNF','Value');

*Drying 1
dry_eff_DRY1=pdry1('dry_eff_DRY1','Value');
cp_rfg_DRY1=pdry1('cp_rfg_DRY1','Value');
U_DRY1=pdry1('U_DRY1','Value');

*Acid Hydrolysis
Y_AHY=pahy('Y_AHY','Value');
T_amb_AHY=pahy('T_amb_AHY','Value');
T_op_AHY=pahy('T_op_AHY','Value');
Lambda_stm=pahy('Lambda_stm','Value');
phi_acid_AHY=pahy('phi_acid_AHY','Value');
Frc_k_AHY=pahy('Frc_k_AHY','Value');

*Neutralization
theta_NT=pnt('theta_NT','Value');

*Filtration,2
Zeta_FLT2=pflt2('Zeta_FLT2','Value');
CPM_FLT2=pflt2('CPM_FLT2','Value');
RFLT2('largesoy')=pflt2('RFLT2_largesoy','Value');
RFLT2('soy')=pflt2('RFLT2_soy','Value');
RFLT2('i-glucoside')=pflt2('RFLT2_i-glucoside','Value');
RFLT2('i-aglycone')=pflt2('RFLT2_i-aglycone','Value');
RFLT2('ethanol')=pflt2('RFLT2_ethanol','Value');
RFLT2('HCl')=pflt2('RFLT2_HCl','Value');
RFLT2('Base')=pflt2('RFLT2_Base','Value');
RFLT2('impurity')=pflt2('RFLT2_impurity','Value');

*Crystallization
phi_antisolv_CRYS=pcrys('phi_antisolv_CRYS','Value');


*Filtration,3
Zeta_FLT3=pflt3('Zeta_FLT3','Value');
CPM_FLT3=pflt3('CPM_FLT3','Value');
RFLT3('largesoy')=pflt3('RFLT3_largesoy','Value');
RFLT3('soy')=pflt3('RFLT3_soy','Value');
RFLT3('i-glucoside')=pflt3('RFLT3_i-glucoside','Value');
RFLT3('i-aglycone')=pflt3('RFLT3_i-aglycone','Value');
RFLT3('ethanol')=pflt3('RFLT3_ethanol','Value');
RFLT3('HCl')=pflt3('RFLT3_HCl','Value');
RFLT3('water')=pflt3('RFLT3_water','Value');
RFLT3('impurity')=pflt3('RFLT3_impurity','Value');

*Nanofiltration
Zeta_NF=pnf('Zeta_NF','Value');
CPM_NF=pnf('CPM_NF','Value');
RNF('largesoy')=pnf('RNF_largesoy','Value');
RNF('soy')=pnf('RNF_soy','Value');
RNF('i-glucoside')=pnf('RNF_i-glucoside','Value');
RNF('i-aglycone')=pnf('RNF_i-aglycone','Value');
RNF('ethanol')=pnf('RNF_ethanol','Value');
RNF('HCl')=pnf('RNF_HCl','Value');
RNF('water')=pnf('RNF_water','Value');
RNF('impurity')=pnf('RNF_impurity','Value');

*Chromatography
theta_CHRM=pchrm('theta_CHRM','Value');
kc_CHRM=pchrm('kc_CHRM','Value');
Width_CHRM=pchrm('Width_CHRM','Value');
HETP_CHRM=pchrm('HETP_CHRM','Value');
Ratio_LD_CHRM=pchrm('Ratio_LD_CHRM','Value');

*Drying 2
dry_eff_DRY2=pdry2('dry_eff_DRY2','Value');
cp_rfg_DRY2=pdry2('cp_rfg_DRY2','Value');
U_DRY2=pdry2('U_DRY2','Value');


Parameters M1(k),MIsoflavoneF;
*Big M values for logical constraints
M1('largesoy')=1000000000;
M1('soy')=1000000000;
M1('i-glucoside')=100000000;
M1('ethanol')=1000000000;
M1('water')=1000000000;
M1('Base')=1000000000;
M1('HCl')=1000000000;
M1('impurity')=1000000000;
MIsoflavoneF=1;

Parameter CRF,Tann,Clbhr,BMC_mult;
BMC_mult=5.4;
CRF=0.11;
*330 days x 24 hrs operating time annually
Tann = 7920;
Clbhr = 30;

Positive Variables
         M(j,k)        mass flowrate of component k in streams j
*         V_feed_soy    Volume of soy feed
         CF(i1)        concentration factor of unit i
         A(i1)         area of unit i
         V(i1)         volume of unit i
         Qc(i1)        standard for cost estimation of unit i
         Cc(i1)        purchase cost of unit i
         PW(i1)        power required in unit i
         Mcw(i1)       cooling water required
         Mstm(i1)      steam required for heating
         M_rfg_DRY1         mass of refrigerant of dryer 1
         M_rfg_DRY2         mass of refrigerant of dryer 2

         Nlbr(i1)      number of laborers for unit i
         Cpr(k)     cost of added components
         Cons(i4)      consumable cost for membrane units
         U(i1)         sigma factor or overall heat transfer coefficient of technolohy;

*single variables
Positive variables
         Feed_C          Cost of entering feed in the system;

Binary variables
         y(ibv)          binary variable corresponding to the unit ibv;

Variable
         Obj             Objective function;

*Upperbounds
M.up('1',k)=2000;
M.up('2',k)=2000;
M.up('3',k)=2000;
M.up('4',k)=2000;
M.up('5',k)=15000;
M.up('6',k)=20000;
M.up('7',k)=2000;
M.up('8',k)=15000;
M.up('9',k)=20000;

Equations
Initial_mass(k)  flow assignment
InitFeed         initial feed cost
CMB(iv,k1)       component mass balance in all units
Cost_units(i1)   cost of units
Nlabr_units(i1)  labors required in units;

Initial_mass(k)..  M('1',k) =e= Qf*Frac(k);
InitFeed.. Feed_C=e=M('1','largesoy')*C_soy;
CMB(iv,k1)$(kI(iv,k1)).. sum(j$jIin(j,iv),M(j,k1))=e=sum(j$jIout(j,iv),M(j,k1));
Cost_units(i1).. Cc(i1)/C0(i1)=e=((Qc(i1)/Q0(i1))**nc(i1));
Nlabr_units(i1).. Nlbr(i1)=e=Nlabr(i1)*(Qc(i1)/Q0(i1));


*Stage-I Pre-Processing
*Grinding
Positive variables
E_grd            Energy per mass (kJ*kg^-1)
PW_grd           Power required  (kW)
AnnualCost       Annualized Cost ($*yr^-1)
Mca_grd          Mass of cooling air (kg*h^-1)
theta            residence time (sec)
Nlbr_new         Calculated number of labor required based on capacity
Cpr_Air_GRD      cost of air ($*hr^-1);

Equations
*CMB_GRD1(k2)                     Component balance
Comp_Released_GRDeqn(jGRDE)      Component Released Equation
Residual_Comp_GRDeqn(jGRDE)      Residual large components after grinding
Costing_variable_GRDeqn          Costing variable of grinding operation (m^3)
Power_required_GRDeqn(j,k)       Power required
Cooling_required_GRDeqn          Cooling air required for removing heat generated through grinding process
Air_cost_GRDeqn                  Cost of air;


*CMB_GRD1(k2).. M('1',k2)=e=M('3',k2);
Comp_Released_GRDeqn(jGRDE).. M(jGRDE,'soy')=e=Y_grd*(M('2','largesoy'));
Residual_Comp_GRDeqn(jGRDE).. M(jGRDE,'largesoy')=e=M('2','largesoy')-M(jGRDE,'soy');
Costing_variable_GRDeqn.. Qc('GRD')=e=Theta_grd*sum(kP,((M('3',kP))/(Den(kP))));
Power_required_GRDeqn(j,k)..   PW('GRD')=e=Wsp('GRD')*sum(kk,M('2',kk));
Cooling_required_GRDeqn.. Mca_grd*Cp_Air*(Tca_out_grd-Tca_in_grd)=e=0.6*PW('GRD');
Air_cost_GRDeqn.. Cpr_Air_GRD=e=Mca_grd*C_Air;
M.fx('3','impurity')=20;
Equations
*Big-M constraints:
*Flowrates
logTE(j,k) logical eqn-TE
logMC(j,k) logical eqn-MC
logUAE(j,k) logical eqn-UAE
logSFE(j,k) logical eqn-SFE
*unit selection
select1 selection of units;

logTE(j,k)$(jI(j,'TE') and (CMP22(j,k)$(jI(j,'TE'))<>0)).. M(j,k)-M1(k)*y('TE')=l=0;
logMC(j,k)$(jI(j,'MC') and (CMP22(j,k)$(jI(j,'MC'))<>0)).. M(j,k)-M1(k)*y('MC')=l=0;
logUAE(j,k)$(jI(j,'UAE') and (CMP22(j,k)$(jI(j,'UAE'))<>0)).. M(j,k)-M1(k)*y('UAE')=l=0;
logSFE(j,k)$(jI(j,'SFE') and (CMP22(j,k)$(jI(j,'SFE'))<>0)).. M(j,k)-M1(k)*y('SFE')=l=0;
select1.. y('TE') + y('MC')+y('UAE')+y('SFE')=e=1;
*y.fx('MC')=0;
y.fx('SFE')=0;
y.fx('UAE')=0;

Positive Variables
*Design variables
         D_tank_TE - diameter of the mixing tank (m)
         D_agit_TE - diameter of agitator (m)
         c_agit_TE - off bottom impeller clearance (m)
         Z_liq_TE - height of liquid in mixing tank (m)
         V_liq_TE - Volume flow of liquid (m^3)
         H_tank_TE - Height of mixing tank (m)
         V_tank_TE  - Volume of tank (m)
         D_baffle_TE - Baffle width (m)
         QS_tank_TE - Costing variable for coned tank (m^3)
         QS_agit_TE - Costing variable for agitator (W) ;

Equations
*Mixing Design

Ethanol_Req_TE           Amount of solvent required
Water_Req_TE             Amount of water required for extraction
V_liq_TEeqn(j,kc)        Volume of Liquid
Z_liq_TEeqn              Liquid Height
D_tank_TEeqn             Tank Diameter

*Costing
QS_tank_TEeqn      Costing variable of cone tank excluding an agitator (m^3)
QS_agit_TEeqn     Costing variable of agitator (Watts) - used for power requirement
Qc_eqn_TE         Combined Costing variable of Turbo-Extraction;



Ethanol_Req_TE..    M('5','ethanol')=e=sum(k,(M('4',k)*8.774));
Water_Req_TE..      M('5','water')=e=sum(k,(M('4',k)*2.780));
V_liq_TEeqn(j,kc).. V_liq_TE=e=sum(kcc,((M('4',kcc)/Den(kcc)))+sum(kccc,(M('5',kccc)/Den(kccc))))*Delta_T_TE;
Z_liq_TEeqn.. V_liq_TE*4=e=(pi)*(D_tank_TE**2)*Z_liq_TE;
D_tank_TEeqn.. QS_tank_TE*4=e=pi*(D_tank_TE**2)*(H_tank_TE);
QS_tank_TEeqn.. QS_tank_TE=e=1.3* V_liq_TE;
QS_agit_TEeqn.. QS_agit_TE=e=P_agit_TE;
Qc_eqn_TE.. Qc('TE')=e=QS_tank_TE;


M.fx('5',kNS)=0;
M.fx('5','i-glucoside')=0;

Positive variables
QH_TE         Heat produced by agitator
QC_TE         Cooling duty
Mcw_TE        Cooling water required
Cost_Util     Utility cost  ;

Equations
*utility requirements
PowReq_TEeqn        Power required for agitator
QH_TE_eqn           Heat produced by agitator
QC_TE_eqn           cooling duty
Mcw_TE_eqn          Cooling water required;

*Assuming an isothermal process (If cooling is enough to cancel heat produced by agitator)
PowReq_TEeqn..      PW('TE')=e=P_agit_TE*V_tank_TE;
QH_TE_eqn..         QH_TE=e=PW('TE');
QC_TE_eqn..         QC_TE =e= QH_TE;
Mcw_TE_eqn..        Mcw_TE*Cp('water')*(Tcw_out_TE-Tcw_in_TE) =e= QC_TE;

Equations
*Design Constraints as recommended by Grenville et al. (2015)
Z_liq_to_D_tank_TEmax    Maximum liquid height to Tank diameter ratio
Z_liq_to_D_tank_TEmin    Minimum liquid height to Tank diameter ratio
D_agit_to_D_tank_TEmax   Maximum Impeller diameter to Tank diameter ratio
D_agit_to_D_tank_TEmin   Minimum Impeller diameter to Tank diameter ratio
c_agit_to_D_tank_TEmax   Maximum off-bottom impeller clearance to Tank diameter ratio
c_agit_to_D_tank_TEmin   Minimum off-bottom impeller clearance to Tank diameter ratio;

*Design Constraints as recommended by Grenville et al. (2015)
Z_liq_to_D_tank_TEmin.. Z_liq_TE=l=D_tank_TE;
Z_liq_to_D_tank_TEmax.. Z_liq_TE=g=0.8*D_tank_TE;
*Design Constraints as recommended by Grenville et al. (2015)
D_agit_to_D_tank_TEmax.. D_agit_TE=l=D_tank_TE*0.50;
D_agit_to_D_tank_TEmin.. D_agit_TE=g=D_tank_TE*0.33;
c_agit_to_D_tank_TEmax.. 0.33*D_tank_TE=g=c_agit_TE;
c_agit_to_D_tank_TEmin.. 0.17*D_tank_TE=l=c_agit_TE;

D_tank_TE.lo=0.01;

*$ontext
*Choice 2: Maceration (MC) **********************************************************************************
Positive variables
*Design variables
         D_tank_MC - diameter of the tank (m)
         Z_liq_MC - height of liquid in the tank (m)
         V_liq_MC - Volume of liquid (m^3)
         H_tank_MC - Height of mixing tank (m)
         V_tank_MC  - Volume of tank (m)
         QS_tank_MC - Costing variable for coned tank (m^3)

Equations
*Tank Design

Ethanol_Req_MC(k3) Amount of solvent required
Water_Req_MC(k3)  Amount of water required for extraction
V_liq_MCeqn(j,kc) Volume of Liquid

QS_tank_MCeqn  Costing variable of cone tank excluding an agitator (m^3)
Qc_eqn_MC      Combined Costing variable of Maceration;


Ethanol_Req_MC(k3)..    M('8','ethanol')=e=8.774*sum(k,M('7',k));
Water_Req_MC(k3)..      M('8','water')=e=2.780*sum(k,M('7',k));
V_liq_MCeqn(j,kc)..     V_liq_MC=e=sum(kcc,((M('7',kcc)/Den(kcc)))+sum(kccc,(M('8',kccc)/Den(kccc))))*Delta_T_MC;

QS_tank_MCeqn.. QS_tank_MC=e=1.3*V_liq_MC;
Qc_eqn_MC.. Qc('MC')=e=QS_tank_MC;


*$ontext
*Ultrasound-Assisted Extraction (UAE)
Positive variables
V_UAE    Calculated volume of ultrasound-assisted extraction unit;

Equations
Ethanol_Req_UAE(k3) Amount of solvent required
Water_Req_UAE(k3)  Amount of water required for extraction
Retfac_UAE efficiency of isoflavone extraction

Confac_UAE CF equation for UAE
*Design equations
Vol_UAE equation for volume
SCap_UAE standard capacity of UAE
Power_req_UAE    power required in ultrasound-assisted extraction;


Ethanol_Req_UAE(k3)..    M('5U','ethanol')=e=8.774*sum(k,M('4U',k));
Water_Req_UAE(k3)..      M('5U','water')=e=2.780*sum(k,M('4U',k));
Retfac_UAE.. M('6U','i-glucoside')=e=Eff_UAE*M('4U','i-glucoside');
Confac_UAE.. CF('UAE')*(M('6U','i-glucoside')/Den('i-glucoside'))=e=(M('4U','i-glucoside')/Den('i-glucoside'));
*Design equations
Vol_UAE.. V_UAE=e=(sum(k,(M('4U',k)/Den(k))))*tR_UAE;
SCap_UAE.. Qc('UAE')=e= V_UAE;
*Power required in UAE
*PW.fx('UAE')=3250;

Power_req_UAE.. PW('UAE')=e=Wsp('UAE')*V_UAE;

*Bounds on CF:
CF.lo('UAE')=1.1;   CF.up('UAE')=15;

*$offtext

*Supercritical Fluid Extraction
Variables
textract_SFE     Extraction time of unit i    ;

Equations

Misoflavone_SFE     Mass of isoflavone extracted
ExtractRatio_SFE    Extraction ratio
Capacity_SFE        Capacity scaling
PowReq_SFE          Power required;



Misoflavone_SFE.. M('6S','i-glucoside')=e=M('4S','i-glucoside');

M.fx('5S','largesoy')=0;
M.fx('5S','soy')=0;
M.fx('5S','i-glucoside')=0;

ExtractRatio_SFE.. M('4S','i-glucoside')=e=e_ratio_SFE*M('5S','CO2');
textract_SFE.lo=0.5; textract_SFE.up=1.6;

*Capacity_SFE..Qc('SFE')=e=sum(k,M('2',k)/rho(k));
Capacity_SFE..Qc('SFE')=e=sum(k,M('6S',k)/Den(k))*textract_SFE;


PowReq_SFE..PW('SFE')=e=Qc('SFE')*Wsp('SFE');

Equations
*Purchase cost of raw materials in Stage 2
Cpur_water_eqn        Purchase cost of water
Cpur_ethanol_eqn      Purchase cost of ethanol
Cpur_CO2_eqn          Purchase cost of CO2
;

Cpur_water_eqn.. Cpr('water') =e= C_water*(M('5','water')+M('8','water')+M('5U','water'));
Cpur_ethanol_eqn.. Cpr('ethanol') =e= C_ethanol*(M('5','ethanol')+M('8','ethanol')+M('5U','ethanol'));
Cpur_CO2_eqn.. Cpr('CO2') =e= C_CO2*(M('5S','CO2'));



Equations
logSDM(j,k) logical eqn-SDM
logFLT1(j,k) logical eqn-FLT1
logCNF(j,k) logical eqn-CNF

*unit selection
select2 selection of units;
logSDM(j,k)$(jI(j,'SDM') and (CMP22(j,k)$(jI(j,'SDM'))<>0)).. M(j,k)-M1(k)*y('SDM')=l=0;
logFLT1(j,k)$(jI(j,'FLT1') and (CMP22(j,k)$(jI(j,'FLT1'))<>0)).. M(j,k)-M1(k)*y('FLT1')=l=0;
logCNF(j,k)$(jI(j,'CNF') and (CMP22(j,k)$(jI(j,'CNF'))<>0)).. M(j,k)-M1(k)*y('CNF')=l=0;
select2.. y('SDM')+y('FLT1')+y('CNF')=e=1;

y.fx('FLT1')=1;
y.fx('SDM')=0;
*If you don't fix the technology, it gives error in mass balance. Check the splitter#3 balances
*Check excel sheet again!   FLT - $634573.59   (SDM is infeasible)

*y.fx('CNF')=1;   (CNF is infeasible)

Variables
CSolvent_Recov    Solvent recovered after extraction
Req_Solv_Makeup   Required solvent makeup
CFeed_Recov       The amount of soy recovered and sold
;

*$ontext

Equations
CSolvent_Recovered(j,k)        Solvent recovered after extraction
*Req_Solvent_Makeup             Required solvent make-up after 5% purge
CFeed_Recovered(j,k)           Feed recovered after solid separation and sold
;

*This number will get subtracted from the hourly demand at the end. Ultimately, only 5% make-up is required hourly
CSolvent_Recovered(j,k)..  CSolvent_Recov=e=M('22','ethanol')*0.95*(C_ethanol)+M('22','water')*0.95*(C_water);
*Req_Solvent_Makeup(j,k)..
CFeed_Recovered(j,k).. CFeed_Recov =e= (M('12','soy')+M('12','largesoy')+M('15','soy')+M('15','largesoy')+M('18','soy')+M('18','largesoy'))*(Csell_soy);
*$offtext
*Evaluated parameters
Parameters Ug_SDM, SOR_SDM;
Ug_SDM=(g*(DpS_SDM**2)*(Den('soy')-Den('ethanol')))/(18*muL_SDM);
SOR_SDM=Eff_SDM*Ug_SDM;

Equations
effeqn_SDM_S(jSDMS,kSld) efficiency of solid removal
effeqn_SDM_W(jSDME,kSol) solvent retention efficiency
Confac_SDM(jSDME) CF equation for sdm
*Design equations
Area_SDM(jSDMF) equation for area
SCap_SDM standard capacity of sedimentation tank;

effeqn_SDM_S(jSDMS,kSld).. Eff_SDM*M('11',kSld)=e=M(jSDMS,kSld);
effeqn_SDM_W(jSDME,kSol).. Eff_SDM*M('11',kSol)=e=M(jSDME,kSol);
Confac_SDM(jSDME).. CF('SDM')*(M('13','soy')/Den('soy'))=e=(M(jSDME,'soy')/Den('soy'));
*Design equations
Area_SDM(jSDMF).. (SOR_SDM*3600)*A('SDM')=e=(M(jSDMF,'soy')/Den('soy'))+(M(jSDMF,'water')/Den('water'));
SCap_SDM.. Qc('SDM')=e=A('SDM');
PW.fx('SDM')=0;

*Filtration 1
Equations
Retent_FLT1(jFLT1E,k)   retention factor in filter
Conc_FLT1(jFLT1E)        concentration factor in filter
Flux_Bal_FLT1(jFLT1E)        flux balance around filter
PW_req_FLT1         power required;

Retent_FLT1(jFLT1E,k)..  M('15',k) =e= M('14',k)*RFLT1(k);
Conc_FLT1(jFLT1E)..   sum(k,M(jFLT1E,k)/Den(k))*CF('FLT1') =e= sum(k,M('14',k)/Den(k));
*Conc_FLT1(jFLT1E)..       CF('FLT1')*(sum(kk,M(jFLT1E,k)/Den(k)))=e= (sum(k,((M('13',k)/Den(k)))));
CF.lo('FLT1')=2;  CF.up('FLT1')=30;


*Note, jFLT1E should really be jFLT1F. However, there is a large amount of base that showed up in the material balance, which messes with the Qc('FLT1'). Therefore, we used jFLT1E instead for now.
Flux_Bal_FLT1(jFLT1E)..  Zeta_FLT1*Qc('FLT1')*CF('FLT1')=e= (sum(k,(M(jFLT1E,k)/Den(k))))*(CF('FLT1')-1);
PW_req_FLT1..         PW('FLT1') =e= Wsp('FLT1')*Qc('FLT1');
*M.up('14','Base')=0;

*Choice 3: Centrifugation (CNF)**************************************************
Equations
Eff_solid_CNF(jCNFS,kSld)          solid removal efficiency
Eff_solv_CNF(jCNFE,kSol)       solvent retention efficiency
CF_eqn_CNF(jCNFE)             concentration factor equation
Sigma_eqn_CNF          sigma factor equation
PW_eqn_CNF             power equation
PW_diss_CNF            power dissapation due to heat;

Eff_solid_CNF(jCNFS,kSld)..       Eff_CNF*M('17',kSld) =e= M(jCNFS,kSld);
Eff_solv_CNF(jCNFE,kSol)..        Eff_CNF*M('17',kSol) =e= M(jCNFE,kSol);
CF_eqn_CNF(jCNFE)..          CF('CNF')*sum(k,M(jCNFE,k)/Den(k)) =e= sum(k,M('17',k)/Den(k));
CF.lo('CNF')=2;CF.up('CNF')=25;
Sigma_eqn_CNF(jCNFF)..        Qc('CNF')*U('CNF') =e= sum(k,M(jCNFF,k)/Den(k));
PW_eqn_CNF(jCNFF)..        PW('CNF') =e= Wsp('CNF')*sum(k,M(jCNFF,k)/Den(k));
PW_diss_CNF..       Mcw('CNF')*Cp('water')*(Tcw_out_CNF-Tcw_in_CNF) =e= 0.4*PW('CNF');
U.up('CNF')=0.01;

*$offtext

*$ontext
Equations
Mixing2bal(k)  Mixing point 2 balance
Split4bal(k)   Splitting point 4 balance;

Mixing2bal(k).. M('13',k)+M('16',k)+M('19',k)=e=M('20',k);
Split4bal(k).. M('20',k)=e=M('21',k);
*$offtext


*Drying 1 Step after soy removal (Freeze drying)
*$ontext
Variable
Q_sub_DRY1    enthalpy of sublimation
A_dry1    area of dryer;

Equations
SolRemoved_DRY1(k3)   solvent removed
EnergyBal_DRY1        energy balance
Rfg_req_DRY1          amount of refridgerant required
Dry_area_DRY1         area of the dryer
Capacity_DRY1         sizing of the dryer
Pow_req_DRY1          power required;


SolRemoved_DRY1(k3).. M('22',k3) =e= dry_eff_DRY1 * M('21',k3);
EnergyBal_DRY1..      Q_sub_DRY1 =e= sum(k3,M('22',k3) * (Cp(k3) * (Tfrz-Tamb) + Hsub(k3)));
Rfg_req_DRY1..        Q_sub_DRY1 =e= M_rfg_DRY1 * cp_rfg_DRY1 * (Trfg_out - Trfg_in);
Dry_area_DRY1..       Q_sub_DRY1 =e= U_DRY1 * A_dry1 * (Tamb-Trfg_out);
Capacity_DRY1..       Qc('DRY1') =e= sum(k3,M('22',k3));
Pow_req_DRY1..        PW('DRY1') =e= Wsp('DRY1') * A_dry1;
*$offtext

M.fx('22','i-glucoside')=0;
M.fx('22','impurity')=0;
*$ontext
Equations

Comp_Released_eqn_AHY        Component Released Equation
*Residual_Comp_eqn        Residual components after release
Acid_added_eqn_AHY(k)           Acid added
Acid_cost_eqn_AHY            Cost of acid
Costing_variable_eqn_AHY     Costing variable of acid hydrolysis (m^3*h^-1)
Power_required_eqn_AHY       Power required
Steam_required_eqn_AHY       Steam required for heating to hydrolysis temperature;

Comp_Released_eqn_AHY.. M('25','i-aglycone')=e=Y_AHY*Frc_k_AHY*(M('23','i-glucoside'));
*The remaining isoflavone glucoside not hydrolyzed
*Residual_Comp_eqn.. M('3','b')=e=M('1','b')-sum(k,(M('3',k)));

*****I will implement an equation later to handle the automatic calculation for acid added. For now, I already specified acid added in the parameter section
Acid_added_eqn_AHY(k).. M('24','HCl')=e=phi_acid_AHY * sum(kk,(M('23',kk)));
Acid_cost_eqn_AHY.. Cpr('HCl')=e= C_acid*M('24','HCl');
Costing_variable_eqn_AHY.. Qc('AHY') =e= sum(k,((M('25',k))/(Den(k))));


*Utility Requirement
Power_required_eqn_AHY.. PW('AHY') =e= Wsp('AHY')*Qc('AHY');
Steam_required_eqn_AHY.. Mstm('AHY')*Lambda_stm =e= ((sum(k,(M('23',k))*Cp(k))+M('24','i-glucoside')*Cp('i-glucoside'))*(T_op_AHY-T_amb_AHY));
*$offtext


*Neutralization (NT)
Positive Variables
NinHCl_NT      Mole HCl input
Nbase_req_NT        Mole of base required
*V_NT_in        Volumetric flow rate of Incoming flow
*V_NT_out       Volumetric flow rate of Effluent flow

;


Equations
*FeedVolFlow_NT        Volumetric flow rate of the incoming feed
Acid_mole_in_NT(jNTF,k)       Mole flow
*OutVolFlow_NT(jNTE)     Volumetric flow rate of output flow
MolBaseRequired_NT      Number of Moles of base required
MassBaseRequired_NT    Mass of base required
PurchaseCost_Base(jNTB)        Base required cost
Flowsummation_NT(j,k)     Total flow out of Neutralization
Costing_variable_eqn(jNTE,k)  Costing variable of crystallization (m^3*h^-1)
Power_required_eqn       Power required
;

*pH calculation to determine the required base
*FeedVolFlow_NT.. V_NT_in=e=sum(k,(M('25',k)/Den(k)));
Acid_mole_in_NT(jNTF,k).. NinHCl_NT=e=1000*M('24','HCl')/MW('HCl');
*OutVolFlow_NT(jNTE).. V_NT_out=e=sum(k,(M('27',k)/Den(k)));
MolBaseRequired_NT.. Nbase_req_NT =e= NinHCl_NT;
MassBaseRequired_NT.. M('26','Base')=e=Nbase_req_NT*MW('Base')/1000;
PurchaseCost_Base(jNTB)..      Cpr('Base')=e=C_Base*M('26','Base');
FLowsummation_NT(j,k).. M('27',k)=e=M('26',k)+M('25',k);
Costing_variable_eqn(jNTE,k).. Qc('NT')=e=theta_NT*sum(kk,((M(jNTE,kk)/(Den(kk)))));


*Utility Requirement
Power_required_eqn.. PW('NT')=e=Wsp('NT')*Qc('NT');
*$ontext

*M.fx('23','i-aglycone')=0;

M.fx('24','i-aglycone')=0;

M.fx('25','largesoy')=0;
*M.fx('25','soy')=0;

M.fx('26','largesoy')=0;
M.fx('26','i-glucoside')=0;
M.fx('26','soy')=0;
M.fx('26','CO2')=0;
M.fx('26','HCl')=0;

*After dryer 2 is turned on
M.fx('26','water')=0;
M.fx('24','ethanol')=0;

M.up('25','water')=50;


*M.fx('27','largesoy')=0;
*M.fx('27','i-glucoside')=0;
*M.fx('27','soy')=0;
*M.fx('27','CO2')=0;







****NOTE: HCl + NaOH yield NaCl + H2O. We did not introduce salt into the component to avoid complications. Therefore, the mass of HCl + NaOH = the mass of NaCl + H2O (conservation of mass still applies)
***The filtration unit (2) will treat HCl and NaOH as 'salt and water' component that get removed.
*Filtration 2
Equations
Retent_FLT2(jFLT2E,k)   retention factor in filter
Conc_FLT2(jFLT2E)        concentration factor in filter
Flux_Bal_FLT2(jFLT2E)        flux balance around filter
PW_req_FLT2         power required;

Retent_FLT2(jFLT2E,k)..  M('28',k) =e= M('27',k)*RFLT2(k);
Conc_FLT2(jFLT2E)..   sum(k,M(jFLT2E,k)/Den(k))*CF('FLT2') =e= sum(k,M('27',k)/Den(k));
*Conc_FLT1(jFLT1E)..       CF('FLT1')*(sum(kk,M(jFLT1E,k)/Den(k)))=e= (sum(k,((M('13',k)/Den(k)))));
CF.lo('FLT2')=2;  CF.up('FLT2')=30;


*Note, jFLT1E should really be jFLT1F. However, there is a large amount of base that showed up in the material balance, which messes with the Qc('FLT1'). Therefore, we used jFLT1E instead for now.
Flux_Bal_FLT2(jFLT2E)..  Zeta_FLT2*Qc('FLT2')*CF('FLT2')=e= (sum(k,(M(jFLT2E,k)/Den(k))))*(CF('FLT2')-1);
PW_req_FLT2..         PW('FLT2') =e= Wsp('FLT2')*Qc('FLT2');

*$offtext


*Stage 4
*$ontext
Equations
*Big-M constraints:
*Flowrates
logCRYS(j,k) logical eqn-CRYS
logNF(j,k) logical eqn-NF
logCHRM(j,k) logical eqn-CHRM

*unit selection
select3 selection of units;

logCRYS(j,k)$(jI(j,'CRYS') and (CMP22(j,k)$(jI(j,'CRYS'))<>0)).. M(j,k)-100000*y('CRYS')=l=0;
logNF(j,k)$(jI(j,'NF') and (CMP22(j,k)$(jI(j,'NF'))<>0)).. M(j,k)-100000*y('NF')=l=0;
logCHRM(j,k)$(jI(j,'CHRM') and (CMP22(j,k)$(jI(j,'CHRM'))<>0)).. M(j,k)-100000*y('CHRM')=l=0;

select3.. y('CRYS') + y('NF')+y('CHRM')=e=1;

*$offtext





*$ontext
*Crystallization (CRYS)
Positive Variables
M_antisolv       mass flow rate of antisolvent in kg h^-1
Cpur_antisolv    Purchase cost of antisolvent ($)
;


Equations

Antisolvent_Req_CRYS(k)   Antisolvent required for crystallization (kg*hr^-1)
Antisolv_Stream_CRYS(k)       Antisolvent stream
Antisolvent_Cost_CRYS    Antisolvent cost
Costing_variable_CRYS(k)  Costing variable of crystallization (m^3*h^-1)
Power_required_CRYS       Power required;

Antisolvent_Req_CRYS(k).. M_antisolv=e=phi_antisolv_CRYS*sum(kk,M('30','i-aglycone'));
Antisolv_Stream_CRYS(k).. M('31','water')=e=M_antisolv;
Antisolvent_Cost_CRYS..   Cpur_antisolv=e=C_water*M('31','water');
Costing_variable_CRYS(k).. Qc('CRYS')=e=sum(kk,((M('32',kk)/(Den(kk)))));

*Utility Requirement
Power_required_CRYS.. PW('CRYS')=e=Wsp('CRYS')*Qc('CRYS');

M.fx('32','CO2')=0;
*M.fx('32','Base')=0;


*Filtration 3
Equations
Retent_FLT3(jFLT3E,k)   retention factor in filter
Conc_FLT3(jFLT3E)        concentration factor in filter
Flux_Bal_FLT3(jFLT3E)        flux balance around filter
PW_req_FLT3         power required;

Retent_FLT3(jFLT3E,k)..  M('33',k) =e= M('32',k)*RFLT3(k);
Conc_FLT3(jFLT3E)..   sum(k,M(jFLT3E,k)/Den(k))*CF('FLT3') =e= sum(k,M('32',k)/Den(k));
*Conc_FLT1(jFLT1E)..       CF('FLT1')*(sum(kk,M(jFLT1E,k)/Den(k)))=e= (sum(k,((M('13',k)/Den(k)))));
CF.lo('FLT3')=2;  CF.up('FLT3')=30;


*Note, jFLT1E should really be jFLT1F. However, there is a large amount of base that showed up in the material balance, which messes with the Qc('FLT1'). Therefore, we used jFLT1E instead for now.
Flux_Bal_FLT3(jFLT3E)..  Zeta_FLT3*Qc('FLT3')*CF('FLT3')=e= (sum(k,(M(jFLT3E,k)/Den(k))))*(CF('FLT3')-1);
PW_req_FLT3..         PW('FLT3') =e= Wsp('FLT3')*Qc('FLT3');

*$offtext

M.fx('30',k)=0;
M.fx('31',k)=0;
M.fx('32',k)=0;
M.fx('33',k)=0;
M.fx('34',k)=0;
M.fx('38',k)=0;
M.fx('39',k)=0;
M.fx('40',k)=0;

*Chromatography
Positive Variables
TauR_CHRM      Retention time (h)
Np_CHRM  Number of plates
L_CHRM   Length of chromatography column (m)
D_CHRM   Diameter of chromatography column (m)
V_CHRM   Volume of chromatography column (m^3)
Ncol_CHRM        Number of chromatography column

Qc_CHRM          Costing variable (m^3*hr^-1)
PW_CHRM          Power required  (kW);


Equations
Retention_Time_eqn_CHRM   Retention time
Number_of_Plates_eqn_CHRM     Number of plates
Length_of_column_eqn_CHRM     Length of column
Diameter_of_column_eqn_CHRM  Diameter of column
Volume_of_column_eqn_CHRM    Volume of column
Component_removed_eqn_CHRM      Solvents and other components removed by chromatography
Products_not_retained_eqn_CHRM Products not retained by chromatography
Number_of_columns_eqn_CHRM(k)    Number of columns required
Total_volume_eqn_CHRM         Total volume (costing variable) for chromatography
Power_required_eqn_CHRM       Power required
;

Retention_Time_eqn_CHRM.. TauR_CHRM=e=(1+kc_CHRM)*theta_CHRM;
Number_of_Plates_eqn_CHRM..     Np_CHRM=e=16*(TauR_CHRM/Width_CHRM)**2;
Length_of_column_eqn_CHRM.. L_CHRM=e=HETP_CHRM*Np_CHRM;
Diameter_of_column_eqn_CHRM.. L_CHRM*Ratio_LD_CHRM=e=D_CHRM;
Volume_of_column_eqn_CHRM.. 4*V_CHRM=e=pi*(D_CHRM**2)*L_CHRM;
Component_removed_eqn_CHRM.. M('39','impurity')=e=kc_CHRM*M('38','impurity');
Products_not_retained_eqn_CHRM.. M('40','i-aglycone')+M('40','i-glucoside')=e=M('38','i-aglycone')+M('38','i-glucoside');
Number_of_columns_eqn_CHRM(k).. Ncol_CHRM*V_CHRM=e=theta_CHRM*sum(kk,((M('38',kk)/(Den(kk)))));
Total_volume_eqn_CHRM.. Qc('CHRM')=e=Ncol_CHRM*V_CHRM;

*Utility
Power_required_eqn_CHRM.. PW('CHRM')=e=Wsp('CHRM')*Qc('CHRM');


*$ontext
*Nanofiltration
Equations
Retent_NF(jNF,k)   retention factor in nanofiltration
Conc_NF(jNFE)        concentration factor in nanofiltration
Flux_Bal_NF(jNFE)        flux balance around nanofiltration
PW_req_NF        power required;

Retent_NF(jNF,k)..  M('36',k) =e= M('35',k)*RNF(k);
Conc_NF(jNFE)..   sum(k,M(jNFE,k)/Den(k))*CF('NF') =e= sum(k,M('35',k)/Den(k));
*Conc_FLT1(jFLT1E)..       CF('FLT1')*(sum(kk,M(jFLT1E,k)/Den(k)))=e= (sum(k,((M('13',k)/Den(k)))));
CF.lo('NF')=1.05;  CF.up('NF')=30;


*Note, jFLT1E should really be jFLT1F. However, there is a large amount of base that showed up in the material balance, which messes with the Qc('FLT1'). Therefore, we used jFLT1E instead for now.
Flux_Bal_NF(jNFE)..  Zeta_NF*Qc('NF')*CF('NF')=e= (sum(k,(M(jNFE,k)/Den(k))))*(CF('NF')-1);
PW_req_NF..         PW('NF') =e= Wsp('NF')*Qc('NF');

M.fx('35','soy')=0;
M.fx('35','largesoy')=0;
*M.fx('35','HCl')=0;
*$offtext
*$ontext
*Drying Step after soy removal (Freeze drying)

Variable
Q_sub_DRY2    enthalpy of sublimation
A_dry2    area of dryer;

Equations
SolRemoved_DRY2(k3)   solvent removed
EnergyBal_DRY2        energy balance
Rfg_req_DRY2          amount of refridgerant required
Dry_area_DRY2         area of the dryer
Capacity_DRY2         sizing of the dryer
Pow_req_DRY2          power required;


SolRemoved_DRY2(k3).. M('42',k3) =e= dry_eff_DRY2 * M('41',k3);
EnergyBal_DRY2..      Q_sub_DRY2 =e= sum(k3,M('42',k3) * (Cp(k3) * (Tfrz-Tamb) + Hsub(k3)));
Rfg_req_DRY2..        Q_sub_DRY2 =e= M_rfg_DRY2 * cp_rfg_DRY2 * (Trfg_out - Trfg_in);
Dry_area_DRY2..       Q_sub_DRY2 =e= U_DRY2 * A_dry2 * (Tamb-Trfg_out);
Capacity_DRY2..       Qc('DRY2') =e= sum(k3,M('42',k3));
Pow_req_DRY2..        PW('DRY2') =e= Wsp('DRY2') * A_dry2;
*$offtext


M.fx('42','i-glucoside')=0;
M.fx('42','i-aglycone')=0;
M.fx('42','impurity')=0;
*M.fx('26','water')=0;
*M.fx('26','ethanol')=0;
*M.fx('24','i-glucoside')=0;
*M.fx('24','water')=0;
*M.fx('24','ethanol')=0;
*M.fx('38',k)=0;

*$ontext
*Defining Upperbounds to reduce calculation time
M.up('4U',k)=2000;
M.up('5U',k)=15000;
M.up('6U',k)=20000;
M.up('4S',k)=2000;
M.up('5S',k)=15000;
M.up('6S',k)=20000;
M.up('10',k)=15000;
M.up('11',k)=15000;
M.up('12',k)=15000;
M.up('13',k)=15000;
M.up('14',k)=30000;
M.up('15',k)=15000;
M.up('16',k)=15000;
M.up('17',k)=15000;
M.up('18',k)=15000;
M.up('19',k)=15000;
$ontext
M.up('20',k)=150000;
M.up('21',k)=150000;
M.up('22',k)=150000;
M.up('23',k)=50000;
M.up('24',k)=50000;
M.up('25',k)=50000;
M.up('26',k)=50000;
M.up('27',k)=50000;
M.up('28',k)=50000;

M.up('29',k)=50000;
M.up('30',k)=50000;
M.up('31',k)=50000;
M.up('32',k)=50000;
M.up('33',k)=50000;
M.up('34',k)=50000;
M.up('35',k)=50000;
M.up('36',k)=50000;
M.up('37',k)=50000;
M.up('38',k)=50000;
M.up('39',k)=50000;
M.up('40',k)=50000;
M.up('41',k)=50000;
M.up('42',k)=50000;
M.up('43',k)=50000;

$offtext

Qc.up('GRD')=5000;
Qc.up('TE')=5000;
Qc.up('MC')=5000;
Qc.up('UAE')=5000;
Qc.up('SDM')=5000;
Qc.up('FLT1')=5000;
Qc.up('CNF')=5000;
Qc.up('DRY1')=50000;
Qc.up('AHY')=5000;
Qc.up('NT')=5000;
Qc.up('FLT2')=5000;
Qc.up('CRYS')=5000;
Qc.up('FLT3')=5000;
Qc.up('NF')=5000;
Qc.up('CHRM')=5000;
Qc.up('DRY2')=50000;


Mca_grd.up=100;
*Cpr_Air_GRD.up=100;
*D_tank_TE.up=10;
D_agit_TE.up=10;
c_agit_TE.up=10;
*Z_liq_TE.up=10;
*V_liq_TE.up=100;
*H_tank_TE.up=10;
*V_tank_TE.up=100;
*QS_tank_TE.up=100;
*V_liq_MC.up=1000;
*QS_tank_MC.up=1000;
*V_UAE.up=1000;
*Q_sub_DRY1.up=2000000;
*A_dry1.up=600;
*Cpur_Base.up=100;
*NinHCl_NT.up=5000;
*Nbase_req_NT.up=5000;
M_antisolv.up=1000;
Cpur_antisolv.up=100000;
TauR_CHRM.up=200;
Np_CHRM.up=10000;
L_CHRM.up=50;
D_CHRM.up=10;
V_CHRM.up=500;
Ncol_CHRM.up=100;


*Equations for stagewise cost analysis:
Positive variables
CCAC(Nstg) Annualized capital cost in stages
CCRM(Nstg) Raw material costs in stages
CCCS(Nstg) Consumable cost in stages
CCLB(Nstg) Labor costs in stages
CCUT(Nstg) Utility costs in stages
CCTC(Nstg) Total cost in stages
CCOT(Nstg) Other costs in stages
CFAC Annual feed costs;

Equation
ACC_eq_s1,ACC_eq_s2,ACC_eq_s3,ACC_eq_s4  Annualized capital cost eqns
RMC_eq_s1,RMC_eq_s2,RMC_eq_s3,RMC_eq_s4 Raw material costs in stages eqns
CSC_eq_s1,CSC_eq_s2,CSC_eq_s3,CSC_eq_s4 Consumables costs in stages eqns
LBC_eq_s1,LBC_eq_s2,LBC_eq_s3,LBC_eq_s4  Labor costs in stages eqns
UTC_eq_s1,UTC_eq_s2,UTC_eq_s3,UTC_eq_s4 Utility costs in stages eqns

TC_eq(Nstg) Total cost in stages eqns
OTC_eq(Nstg) Other costs in stages eqns
FeedC_eq Annual cost of entering feed;


ACC_eq_s1.. CCAC('s1')=e=1.66*CRF*BMC_mult*(sum(istg1,Cc(istg1)));
ACC_eq_s2.. CCAC('s2')=e=1.66*CRF*BMC_mult*(sum(istg2,Cc(istg2)));
ACC_eq_s3.. CCAC('s3')=e=1.66*CRF*BMC_mult*(sum(istg3,Cc(istg3)));
ACC_eq_s4.. CCAC('s4')=e=1.66*CRF*BMC_mult*(sum(istg4,Cc(istg4)));

*$ontext
RMC_eq_s1.. CCRM('s1')=e=(Tann*(10**(-6))*(Feed_C+Cpr_Air_GRD));
RMC_eq_s2.. CCRM('s2')=e=(Tann*(10**(-6))*(Cpr('ethanol')+Cpr('water')+Cpr('CO2')-CSolvent_Recov-CFeed_Recov));
RMC_eq_s3.. CCRM('s3')=e=0;
RMC_eq_s4.. CCRM('s4')=e=(Tann*(10**(-6))*(Cpur_antisolv));
*RMC_eq_s4.. CCRM('s4')=e=0;
*$offtext



CSC_eq_s1.. CCCS('s1')=e=0;
CSC_eq_s2.. CCCS('s2')=e=(Tann*10**(-6))*(sum(i4,CPM_FLT1*Qc('FLT1')))/Rep_time;
CSC_eq_s3.. CCCS('s3')=e=((Tann*10**(-6))*((sum(i4,CPM_FLT2*Qc('FLT2')))/Rep_time)+(Tann*10**(-6))*((Cpr('HCl')+Cpr('Base'))));
CSC_eq_s4.. CCCS('s4')*Rep_time*1e6=e=Tann*(sum(i4,CPM_FLT3*Qc('FLT3')))+Tann*(sum(i4,CPM_NF*Qc('NF')));

LBC_eq_s1.. CCLB('s1')=e=Tann*(10**(-6))*Clbhr*sum(istg1,Nlbr(istg1));
LBC_eq_s2.. CCLB('s2')=e=Tann*(10**(-6))*Clbhr*sum(istg2,Nlbr(istg2));
LBC_eq_s3.. CCLB('s3')=e=Tann*(10**(-6))*Clbhr*sum(istg3,Nlbr(istg3));
LBC_eq_s4.. CCLB('s4')=e=Tann*(10**(-6))*Clbhr*sum(istg4,Nlbr(istg4));

UTC_eq_s1.. CCUT('s1')=e=Tann*(10**(-6))*((C_elec*sum(istg1,PW(istg1)))+(C_air*sum(istg1,Mca_grd)));
UTC_eq_s2.. CCUT('s2')=e=Tann*(10**(-6))*((C_elec*sum(istg2,PW(istg2)))+(C_cwt*sum(istg2,Mcw(istg2))));
UTC_eq_s3.. CCUT('s3')=e=Tann*(10**(-6))*((C_elec*sum(istg3,PW(istg3)))+(C_cwt*sum(istg3,Mcw(istg3))));
UTC_eq_s4.. CCUT('s4')=e=Tann*(10**(-6))*((C_elec*sum(istg4,PW(istg4)))+(C_cwt*sum(istg4,Mcw(istg4))));

TC_eq(Nstg).. CCTC(Nstg)=e=CCAC(Nstg)+CCRM(Nstg)+CCCS(Nstg)+CCUT(Nstg)+2.78*CCLB(Nstg);
OTC_eq(Nstg).. CCOT(Nstg)=e=CCTC(Nstg)-(CCAC(Nstg)+CCRM(Nstg)+CCCS(Nstg)+CCLB(Nstg)+CCUT(Nstg));
FeedC_eq.. CFAC=e=Feed_C*(10**(-6));

M.lo('6','i-glucoside')=1;
M.up('15','HCl')=0;
M.up('14','i-aglycone')=0;

Positive variables
CCTAC - Annualized capital cost
CCTRM - raw material cost
CCTCS - Consumables cost
CCTUT - utility cost
CCTLB - labor cost
CCTOT - Other costs;

Variable
CCTPC - total process cost
Profit - profit
Revenue - revenue;

Equation
ACap_Cost annualized total capital cost
RM_Cost raw materials cost
CS_Cost Consumable cost eqn
Labor_Cost labor cost
Util_Cost utility cost
TPC_fun total process cost
OTH_cost total other costs
Ob_fun objective function
Profit_eqn(k) Amount of profit from the process
Revenue_eqn(k) Amount of revenue generated;

ACap_Cost.. CCTAC=e=sum(Nstg1,CCAC(Nstg1));
RM_Cost.. CCTRM=e=sum(Nstg1,CCRM(Nstg1));
CS_Cost.. CCTCS=e=sum(Nstg1,CCCS(Nstg1));
Labor_Cost.. CCTLB=e=sum(Nstg1,CCLB(Nstg1));
Util_Cost.. CCTUT=e=sum(Nstg1,CCUT(Nstg1));
OTH_cost.. CCTOT=e=sum(Nstg1,CCOT(Nstg1));
TPC_fun.. CCTPC=e=CCTAC+CCTRM+CCTCS+CCTLB+CCTUT+CCTOT+CFAC;
Ob_fun.. Obj=e=CCTPC;
Profit_eqn(k).. Profit=e=Revenue-CCTPC;
Revenue_eqn(k).. Revenue=e=Tann*(10**(-6))*C_isoflavone*M('43','i-aglycone');

Model IsoflavoneExtraction /all/;

Solve IsoflavoneExtraction minimizing Obj using MINLP;
EXECUTE_UNLOAD 'Base_case_Isoflavone_extraction_9-14-2020.gdx';
