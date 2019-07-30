#include "StdAfx.h"
const char * theNameVar_ParamNewRechPH[299] = {
"\n",
"<GenCpp  Difndef=\"Define_NotRechNewPH\">\n",
"\n",
"<Verbatim File=\".cpp\">\n",
"#include \"StdAfx.h\"\n",
"#include \"cParamNewRechPH.h\"\n",
"</Verbatim>\n",
"<Verbatim File=\".h.cpp\">\n",
"// NOMORE ... \n",
"</Verbatim>\n",
"<!--\n",
"    eTPR_Corner  = 2,\n",
"    eTPR_MaxLapl = 3,\n",
"    eTPR_MinLapl = 4,\n",
"    eTPR_NoLabel = 5\n",
"-->\n",
"   <enum Name=\"eTypePtRemark\">\n",
"      <eTPR_LaplMax>      </eTPR_LaplMax>\n",
"      <eTPR_LaplMin>      </eTPR_LaplMin>\n",
"\n",
"      <eTPR_GrayMax>      </eTPR_GrayMax>\n",
"      <eTPR_GrayMin>      </eTPR_GrayMin>\n",
"\n",
"      <eTPR_BifurqMax>    </eTPR_BifurqMax>\n",
"      <eTPR_BifurqMin>    </eTPR_BifurqMin>\n",
"\n",
"      <eTPR_NoLabel>  </eTPR_NoLabel>\n",
"      <!-- Apparement les sadl points sont peu fiables, en plus il\n",
"          semblent (???) couter cher en temps, donc on les met pas ,\n",
"          mais on laisse dans le code si on reprend + tard\n",
"      -->\n",
"\n",
"      <eTPR_GraySadl>      </eTPR_GraySadl>\n",
"      <eTPR_BifurqSadl>    </eTPR_BifurqSadl>\n",
"   </enum>\n",
"\n",
"   <!-- Sotcke les noms des differentes options de calcul d'invariants -->\n",
"   <enum Name=\"eTypeVecInvarR\">\n",
"         <eTVIR_Curve>     </eTVIR_Curve>\n",
"         <eTVIR_ACR0>      </eTVIR_ACR0>  <!-- Auto correl radiom basiques -->\n",
"         <eTVIR_ACGT>      </eTVIR_ACGT>  <!-- Auto correl gradient tangentiel-->\n",
"         <eTVIR_ACGR>      </eTVIR_ACGR>  <!-- Auto correl gradient radial-->\n",
"         <eTVIR_LogPol>      </eTVIR_LogPol>  <!-- NOT an invariant, BUT usefull for exporting -->\n",
"         <eTVIR_NoLabel>   </eTVIR_NoLabel>  <!-- Auto correl gradient radial-->\n",
"   </enum>\n",
"\n",
"<!--\n",
"         <CoeffRadiom   Nb=\"1\" Type=\"std::vector<double>\">           </CoeffRadiom>\n",
"         <CoeffRadiom2  Nb=\"1\" Type=\"std::vector<double>\">           </CoeffRadiom2>\n",
"         <CoeffGradRadial   Nb=\"1\" Type=\"std::vector<double>\">       </CoeffGradRadial>\n",
"         <CoeffGradTangent   Nb=\"1\" Type=\"std::vector<double>\">      </CoeffGradTangent>\n",
"         <CoeffGradTangentPiS2   Nb=\"1\" Type=\"std::vector<double>\">  </CoeffGradTangentPiS2>\n",
"         <CoeffGradTangentPi     Nb=\"1\" Type=\"std::vector<double>\">  </CoeffGradTangentPi>\n",
"         <CoeffGradCroise   Nb=\"1\" Type=\"std::vector<double>\">       </CoeffGradCroise>\n",
"         <CoeffGradCroise2  Nb=\"1\" Type=\"std::vector<double>\">       </CoeffGradCroise2>\n",
"         <CoeffDiffOpposePi   Nb=\"1\" Type=\"std::vector<double>\">     </CoeffDiffOpposePi>\n",
"         <CoeffDiffOppose2Pi  Nb=\"1\" Type=\"std::vector<double>\">     </CoeffDiffOppose2Pi>\n",
"         <CoeffDiffOpposePiS2   Nb=\"1\" Type=\"std::vector<double>\">   </CoeffDiffOpposePiS2>\n",
"         <CoeffDiffOppose2PiS2  Nb=\"1\" Type=\"std::vector<double>\">   </CoeffDiffOppose2PiS2>\n",
"-->\n",
"   <enum Name=\"eTypeInvRad\">\n",
"      <eTIR_Radiom>          </eTIR_Radiom>    <!-- 1 -->\n",
"      <eTIR_GradRad>         </eTIR_GradRad>\n",
"      <eTIR_GradCroise>      </eTIR_GradCroise>\n",
"      <eTIR_GradTan>         </eTIR_GradTan>\n",
"      <eTIR_GradTanPiS2>     </eTIR_GradTanPiS2>\n",
"      <eTIR_GradTanPi>       </eTIR_GradTanPi>   \n",
"      <eTIR_LaplRad>         </eTIR_LaplRad>\n",
"      <eTIR_LaplTan>         </eTIR_LaplTan>\n",
"      <eTIR_LaplCrois>       </eTIR_LaplCrois>\n",
"      <eTIR_DiffOpposePi>    </eTIR_DiffOpposePi>\n",
"      <eTIR_DiffOpposePiS2>   </eTIR_DiffOpposePiS2>  <!-- 11 -->\n",
"      <!--         BEGIN SQUARE -->\n",
"      <eTIR_Sq_Radiom>        </eTIR_Sq_Radiom>   \n",
"      <eTIR_Sq_GradRad>       </eTIR_Sq_GradRad>\n",
"      <eTIR_Sq_GradCroise>    </eTIR_Sq_GradCroise>\n",
"      <eTIR_Sq_GradTan>       </eTIR_Sq_GradTan>\n",
"      <eTIR_Sq_GradTanPiS2>   </eTIR_Sq_GradTanPiS2>    <!-- 16 -->\n",
"      <eTIR_Sq_GradTanPi>     </eTIR_Sq_GradTanPi>\n",
"      <eTIR_Sq_LaplRad>       </eTIR_Sq_LaplRad>\n",
"      <eTIR_Sq_LaplTan>       </eTIR_Sq_LaplTan>\n",
"      <eTIR_Sq_LaplCrois>     </eTIR_Sq_LaplCrois>\n",
"      <eTIR_Sq_DiffOpposePi>    </eTIR_Sq_DiffOpposePi>\n",
"      <eTIR_Sq_DiffOpposePiS2>  </eTIR_Sq_DiffOpposePiS2>  <!-- 22 -->\n",
"      <!--         BEGIN CUBE -->\n",
"      <eTIR_Cub_Radiom>        </eTIR_Cub_Radiom>   \n",
"      <eTIR_Cub_GradRad>       </eTIR_Cub_GradRad>\n",
"      <eTIR_Cub_GradCroise>    </eTIR_Cub_GradCroise>\n",
"      <eTIR_Cub_GradTan>       </eTIR_Cub_GradTan>\n",
"      <eTIR_Cub_GradTanPiS2>   </eTIR_Cub_GradTanPiS2>    <!-- 16 -->\n",
"      <eTIR_Cub_GradTanPi>     </eTIR_Cub_GradTanPi>\n",
"      <eTIR_Cub_LaplRad>       </eTIR_Cub_LaplRad>\n",
"      <eTIR_Cub_LaplTan>       </eTIR_Cub_LaplTan>\n",
"      <eTIR_Cub_LaplCrois>     </eTIR_Cub_LaplCrois>\n",
"      <eTIR_Cub_DiffOpposePi>  </eTIR_Cub_DiffOpposePi>\n",
"      <eTIR_Cub_DiffOpposePiS2></eTIR_Cub_DiffOpposePiS2>  <!-- 22 -->\n",
"\n",
"      <eTIR_NoLabel>     </eTIR_NoLabel>\n",
"   </enum>\n",
"\n",
"   <PtSc  Nb=\"*\" Class=\"true\" ToReference=\"true\">\n",
"       <Pt   Nb=\"1\" Type=\"Pt2dr\">     </Pt>\n",
"       <Scale   Nb=\"1\" Type=\"double\">  </Scale>\n",
"   </PtSc>\n",
"\n",
"   <OneInvRad  Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"         <!-- y :  NumVect , x : Rho -->\n",
"         <ImRad Nb=\"1\" Type=\"Im2D_INT1\"> </ImRad>\n",
"         <!-- Code Binaire : devrait etre 1D, mais flemme modifier xml_init et autres ...  -->\n",
"         <CodeBinaire  Nb=\"1\" Type=\"Im2D_U_INT2\"> </CodeBinaire>\n",
"   </OneInvRad>\n",
"\n",
"   <ProfilRad  Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"         <!-- y :  NumVect , x : Rho -->\n",
"         <ImProfil Nb=\"1\" Type=\"Im2D_INT1\"> </ImProfil>\n",
"         <!-- Code Binaire : devrait etre 1D, mais flemme modifier xml_init et autres ...  -->\n",
"   </ProfilRad>\n",
"\n",
"   <!-- Les invariants rotationnels obtenu par auto-correlation -->\n",
"   <RotInvarAutoCor Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"         <IR0 Nb=\"1\" Type=\"Im2D_INT1\"> </IR0>\n",
"         <IGT Nb=\"1\" Type=\"Im2D_INT1\"> </IGT>\n",
"         <IGR Nb=\"1\" Type=\"Im2D_INT1\"> </IGR>\n",
"   </RotInvarAutoCor>\n",
"\n",
"\n",
"   <OnePCarac   Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"         <Kind Nb=\"1\" Type=\"eTypePtRemark\">        </Kind>\n",
"         <Pt   Nb=\"1\" Type=\"Pt2dr\">     </Pt>\n",
"         <Pt0   Nb=\"1\" Type=\"Pt2dr\">     </Pt0>  <!-- Avant opt, pour test -->\n",
"         <NivScale   Nb=\"1\" Type=\"int\">  </NivScale>\n",
"         <Scale   Nb=\"1\" Type=\"double\">  </Scale>\n",
"         <ScaleStab   Nb=\"1\" Type=\"double\">  </ScaleStab> <!-- 4 MinMax, highest scale where they are visible -->\n",
"         <ScaleNature   Nb=\"1\" Type=\"double\">  </ScaleNature> <!-- May be != from scale when scale was forced -->\n",
"         <DirMS   Nb=\"1\" Type=\"Pt2dr\">     </DirMS>\n",
"         <DirAC   Nb=\"1\" Type=\"Pt2dr\">     </DirAC>\n",
"         <Contraste   Nb=\"1\" Type=\"double\">  </Contraste>\n",
"         <ContrasteRel   Nb=\"1\" Type=\"double\">  </ContrasteRel> <!-- Par raport au seuil, pour inspection-->\n",
"         <AutoCorrel   Nb=\"1\" Type=\"double\">  </AutoCorrel>\n",
"         <OK   Nb=\"1\" Type=\"bool\">  </OK>  <!-- Help 4 compute, should always be true -->\n",
"     <!--   codage binaire -->\n",
"    <!-- Pour visu -->\n",
"         <InvR Nb=\"1\" RefType=\"OneInvRad\">  </InvR>\n",
"         <MoyLP   Nb=\"1\" Type=\"double\">  </MoyLP> <!-- Moyenne de l'image log pol-->\n",
"         <ImLogPol Nb=\"1\" Type=\"Im2D_INT1\"> </ImLogPol>\n",
"         <VectRho   Nb=\"1\" Type=\"std::vector<double>\">  </VectRho>\n",
"         <ProfR Nb=\"1\" RefType=\"ProfilRad\">  </ProfR>\n",
"         <RIAC  Nb=\"1\" RefType=\"RotInvarAutoCor\">  </RIAC>\n",
"         <Id   Nb=\"1\" Type=\"int\">  </Id>  <!-- Tuning -->\n",
"         <HeapInd   Nb=\"1\" Type=\"int\">  </HeapInd>  <!-- Tuning -->\n",
"         <Prio   Nb=\"1\" Type=\"double\">  </Prio>  <!-- Tuning -->\n",
"   </OnePCarac>\n",
"\n",
"   <SetPCarac  Nb=\"1\" Class=\"true\">\n",
"      <OnePCarac  Nb=\"*\" RefType=\"OnePCarac\" Container=\"std::vector\"> </OnePCarac>\n",
"   </SetPCarac>\n",
"\n",
"   <SetRefPCarac Nb=\"1\" Class=\"true\">\n",
"         <SRPC_Truth Nb=\"*\" Container=\"std::vector\">\n",
"             <P1  Nb=\"1\" RefType=\"OnePCarac\"> </P1>\n",
"             <P2  Nb=\"1\" RefType=\"OnePCarac\"> </P2>\n",
"         </SRPC_Truth>\n",
"         <SRPC_Rand  Nb=\"*\" RefType=\"OnePCarac\" Container=\"std::vector\"> </SRPC_Rand>\n",
"   </SetRefPCarac>\n",
"\n",
"   <CBOneVect   Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"         <IndVec Nb=\"1\" Type=\"int\">     </IndVec> \n",
"         <CBOneBit Nb=\"*\" Container=\"std::vector\"> \n",
"            <Coeff  Nb=\"1\" Type=\"std::vector<double>\">  </Coeff> \n",
"            <IndInV Nb=\"1\" Type=\"std::vector<int>\">     </IndInV> \n",
"            <IndBit Nb=\"1\" Type=\"int\">     </IndBit> \n",
"         </CBOneBit>\n",
"   </CBOneVect>\n",
"\n",
"   <FullParamCB  Nb=\"1\" Class=\"true\">\n",
"        <CBOneVect  Nb=\"*\" RefType=\"CBOneVect\" Container=\"std::vector\"> </CBOneVect>\n",
"   </FullParamCB>\n",
"\n",
"   \n",
"   <CompCBOneBit Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"      <Coeff  Nb=\"1\" Type=\"std::vector<double>\">  </Coeff> \n",
"      <IndX   Nb=\"1\" Type=\"std::vector<int>\">   </IndX> \n",
"      <IndY   Nb=\"1\" Type=\"std::vector<int>\">   </IndY> \n",
"      <IndBit Nb=\"1\" Type=\"int\">     </IndBit> \n",
"   </CompCBOneBit>\n",
"\n",
"   <CompCB Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"       <BitThresh Nb=\"1\" Type=\"int\"> </BitThresh>\n",
"       <CompCBOneBit  Nb=\"*\" RefType=\"CompCBOneBit\" Container=\"std::vector\"> </CompCBOneBit>\n",
"   </CompCB>\n",
"\n",
"   \n",
"   <FitsOneBin Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"       <PrefName Nb=\"1\" Type=\"std::string\">     </PrefName> \n",
"       <PostName Nb=\"?\" Type=\"std::string\" Def=\"_Local.xml\">     </PostName> \n",
"       <CCB Nb=\"?\" RefType=\"CompCB\">            </CCB>\n",
"   </FitsOneBin>\n",
"\n",
"   <FitsOneLabel Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"       <KindOf Nb=\"1\" Type=\"eTypePtRemark\">     </KindOf> \n",
"       <BinIndexed Nb=\"1\" RefType=\"FitsOneBin\">  </BinIndexed>  <!-- Quick one, for indexation -->\n",
"       <BinDecisionShort Nb=\"1\" RefType=\"FitsOneBin\"> </BinDecisionShort>  <!-- Quick with lowe error -->\n",
"       <BinDecisionLong  Nb=\"1\" RefType=\"FitsOneBin\"> </BinDecisionLong>  <!-- Long with low error -->\n",
"   </FitsOneLabel>\n",
"\n",
"   <SeuilFitsParam  Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"          <SeuilCorrDR Nb=\"?\" Type=\"double\" Def=\"0.7\"> </SeuilCorrDR>  <!--  Seuil corr sur inv rad -->\n",
"          <SeuilInc    Nb=\"?\" Type=\"double\" Def=\"0.01\"> </SeuilInc>    <!-- Seuil Inco entre les div estim de shift -->\n",
"          <SeuilCorrLP Nb=\"?\" Type=\"double\" Def=\"0.93\"> </SeuilCorrLP> <!-- Seuil correl Log Polar-->\n",
"          <ExposantPdsDistGrad Nb=\"?\" Type=\"double\" Def=\"0.5\"> </ExposantPdsDistGrad> <!-- Seuil correl Log Polar-->\n",
"          <SeuilDistGrad Nb=\"?\" Type=\"double\" Def=\"0.5\"> </SeuilDistGrad> <!-- Seuil correl Log Polar-->\n",
"          <SeuilCorrelRatio12 Nb=\"?\" Type=\"double\" Def=\"0.6\"> </SeuilCorrelRatio12>\n",
"          <SeuilGradRatio12 Nb=\"?\" Type=\"double\" Def=\"0.6\">   </SeuilGradRatio12>\n",
"   </SeuilFitsParam>\n",
"   \n",
"   \n",
"   <FitsParam Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"       <DefInit  Nb=\"1\" RefType=\"FitsOneLabel\"> </DefInit>  <!-- 0 ou 1, sera complete autom -->\n",
"       <GenLabs  Nb=\"*\" RefType=\"FitsOneLabel\"> </GenLabs>  <!-- 0 ou 1, sera complete autom -->\n",
"       <SeuilOL Nb=\"1\" RefType=\"SeuilFitsParam\">  </SeuilOL>\n",
"       <SeuilGen Nb=\"1\" RefType=\"SeuilFitsParam\"> </SeuilGen>\n",
"   </FitsParam>\n",
"\n",
"<!--\n",
"mm3d Malt GeomImage \"Soldat-Hue-IMGP703[0-2].JPG\" AllRel Master=Soldat-Hue-IMGP7031.JPG ZoomF=4\n",
"mm3d PHom_RenameRef 6 Soldat-Hue-IMGP7031.JPG\n",
"-->\n",
"\n",
"   <XmlAimeParamApprentissage Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"       <AbsDir Nb=\"1\" Type=\"std::string\"> </AbsDir>  <!-- pour connaitre le nom absolu a passer aux sous prog-->\n",
"       <DefDoIt Nb=\"?\" Type=\"bool\" Def=\"true\"> </DefDoIt>\n",
"       <DefDoMatch Nb=\"?\" Type=\"bool\" Def=\"true\"> </DefDoMatch>\n",
"       <DefDoPtCar Nb=\"?\" Type=\"bool\" Def=\"true\"> </DefDoPtCar>\n",
"       <DefDoRef   Nb=\"?\" Type=\"bool\" Def=\"true\"> </DefDoRef>\n",
"       <DefDoApprComb   Nb=\"?\" Type=\"bool\" Def=\"true\"> </DefDoApprComb>\n",
"       <DefDoApprLocal1   Nb=\"?\" Type=\"bool\" Def=\"true\"> </DefDoApprLocal1>\n",
"       <DefDoApprLocal2   Nb=\"?\" Type=\"bool\" Def=\"true\"> </DefDoApprLocal2>\n",
"       <DefParamPtCar Nb=\"?\" Type=\"std::string\" Def=\"\"> </DefParamPtCar>\n",
"\n",
"       <XlmAimeOneDir Nb=\"*\">\n",
"            <DoIt Nb=\"?\" Type=\"bool\"> </DoIt>\n",
"            <DoMatch Nb=\"?\" Type=\"bool\"> </DoMatch>\n",
"            <DoPtCar Nb=\"?\" Type=\"bool\"> </DoPtCar>\n",
"            <DoRef   Nb=\"?\" Type=\"bool\"> </DoRef>\n",
"            <ZoomF Nb=\"?\" Type=\"int\" Def=\"4\"> </ZoomF>\n",
"            <NumMatch Nb=\"?\" Type=\"int\"> </NumMatch>\n",
"\n",
"            <Dir Nb=\"1\" Type=\"std::string\"> </Dir>\n",
"            <Ori Nb=\"1\" Type=\"std::string\"> </Ori>\n",
"\n",
"            <XAPA_OneMatch Nb=\"*\">\n",
"                <Master Nb=\"1\" Type=\"std::string\">        </Master>\n",
"                <!-- Pattern pour le match dense, plutot des B/H faibles -->\n",
"                <Pattern Nb=\"1\" Type=\"std::string\">       </Pattern> \n",
"                <!-- Pattern pour les appariements de points caracteristiques -->\n",
"                <PatternRef Nb=\"1\" Type=\"std::string\">    </PatternRef>\n",
"            </XAPA_OneMatch>\n",
"\n",
"            <!-- Ceux sur lesquel il faut lancer le calcul de pts car TestNewRechPH, union des pattern ...-->\n",
"            <XAPA_PtCar Nb=\"1\">\n",
"                <Pattern  Nb=\"1\" Type=\"std::string\" AccessorFils=\"false\"> </Pattern>\n",
"            </XAPA_PtCar>\n",
"       </XlmAimeOneDir>\n",
"\n",
"       <XlmAimeApprent Nb=\"1\">\n",
"             <NbExEt0 Nb=\"1\" Type=\"int\">       </NbExEt0>\n",
"             <NbExEt1 Nb=\"1\" Type=\"int\">       </NbExEt1>\n",
"             <TimeOut Nb=\"?\" Type=\"double\" Def=\"300.0\">       </TimeOut>\n",
"             <XlmAimeOneApprent  Nb=\"*\">\n",
"                  <PdsW Nb=\"1\" Type=\"double\">  </PdsW>\n",
"                  <NbBB Nb=\"1\" Type=\"int\">     </NbBB>\n",
"                  <BitM Nb=\"?\" Type=\"int\">     </BitM>\n",
"             </XlmAimeOneApprent>\n",
"       </XlmAimeApprent>\n",
"\n",
"   </XmlAimeParamApprentissage>\n",
"\n",
"   <Xml2007Pt Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"       <Pt   Nb=\"1\" Type=\"Pt2dr\">       </Pt>\n",
"       <NumOct   Nb=\"1\" Type=\"int\">     </NumOct>\n",
"       <NumIm    Nb=\"1\" Type=\"int\">     </NumIm>\n",
"       <ScaleInO   Nb=\"1\" Type=\"double\">   </ScaleInO>\n",
"       <ScaleAbs   Nb=\"1\" Type=\"double\">   </ScaleAbs>\n",
"   </Xml2007Pt>\n",
"\n",
"   <Xml2007FilePt  Nb=\"1\" Class=\"true\" ToReference=\"true\">\n",
"       <Pts  Nb=\"*\" RefType=\"Xml2007Pt\"> </Pts> \n",
"       <IsMin Nb=\"1\" Type=\"bool\">        </IsMin>\n",
"       <TypePt Nb=\"1\" Type=\"int\">        </TypePt>\n",
"       <NameTypePt Nb=\"1\" Type=\"std::string\">    </NameTypePt>\n",
"   </Xml2007FilePt>\n",
"\n",
"<Verbatim File=\".h.cpp\">\n",
"// };\n",
"</Verbatim>\n",
"\n",
"</GenCpp>\n",
"","//#_-=+{}@$##$##@"};
