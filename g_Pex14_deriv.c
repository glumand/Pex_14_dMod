#include <R.h>
 #include <math.h>
 void g_Pex14_deriv_k968ekox ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (p[1])*(x[11+i**k]) ;
y[1+i**l] = (p[1])*(x[12+i**k]) ;
y[2+i**l] = (p[1])*(x[13+i**k]) ;
y[3+i**l] = (p[1])*(x[14+i**k]) ;
y[4+i**l] = (p[1])*(x[15+i**k]) ;
y[5+i**l] = (p[1])*(x[16+i**k]) ;
y[6+i**l] = (p[1])*(x[17+i**k]) ;
y[7+i**l] = (p[1])*(x[18+i**k]) ;
y[8+i**l] = (p[1])*(x[19+i**k]) ;
y[9+i**l] = (p[1])*(x[22+i**k]) ;
y[10+i**l] = (p[1])*(x[23+i**k]) ;
y[11+i**l] = (p[1])*(x[24+i**k]) ;
y[12+i**l] = (p[1])*(x[25+i**k]) ;
y[13+i**l] = (p[1])*(x[26+i**k]) ;
y[14+i**l] = (p[1])*(x[27+i**k]) ;
y[15+i**l] = (p[1])*(x[28+i**k]) ;
y[16+i**l] = (p[1])*(x[29+i**k]) ;
y[17+i**l] = (p[1])*(x[30+i**k]) ;
y[18+i**l] = (p[1])*(x[33+i**k]) ;
y[19+i**l] = (p[1])*(x[34+i**k]) ;
y[20+i**l] = (p[1])*(x[35+i**k]) ;
y[21+i**l] = (p[1])*(x[36+i**k]) ;
y[22+i**l] = (p[1])*(x[37+i**k]) ;
y[23+i**l] = (p[1])*(x[38+i**k]) ;
y[24+i**l] = (p[1])*(x[39+i**k]) ;
y[25+i**l] = (p[1])*(x[40+i**k]) ;
y[26+i**l] = (p[1])*(x[41+i**k]) ;
y[27+i**l] = (p[1])*(x[44+i**k]) ;
y[28+i**l] = (p[1])*(x[45+i**k]) ;
y[29+i**l] = (p[1])*(x[46+i**k]) ;
y[30+i**l] = (p[1])*(x[47+i**k]) ;
y[31+i**l] = (p[1])*(x[48+i**k]) ;
y[32+i**l] = (p[1])*(x[49+i**k]) ;
y[33+i**l] = (p[1])*(x[50+i**k]) ;
y[34+i**l] = (p[1])*(x[51+i**k]) ;
y[35+i**l] = (p[1])*(x[52+i**k]) ;
y[36+i**l] = (p[1])*(x[55+i**k]) ;
y[37+i**l] = (p[1])*(x[56+i**k]) ;
y[38+i**l] = (p[1])*(x[57+i**k]) ;
y[39+i**l] = (p[1])*(x[58+i**k]) ;
y[40+i**l] = (p[1])*(x[59+i**k]) ;
y[41+i**l] = (p[1])*(x[60+i**k]) ;
y[42+i**l] = (p[1])*(x[61+i**k]) ;
y[43+i**l] = (p[1])*(x[62+i**k]) ;
y[44+i**l] = (p[1])*(x[63+i**k]) ;
y[45+i**l] = (p[1])*(x[66+i**k]) ;
y[46+i**l] = (p[1])*(x[67+i**k]) ;
y[47+i**l] = (p[1])*(x[68+i**k]) ;
y[48+i**l] = (p[1])*(x[69+i**k]) ;
y[49+i**l] = (p[1])*(x[70+i**k]) ;
y[50+i**l] = (p[1])*(x[71+i**k]) ;
y[51+i**l] = (p[1])*(x[72+i**k]) ;
y[52+i**l] = (p[1])*(x[73+i**k]) ;
y[53+i**l] = (p[1])*(x[74+i**k]) ;
y[54+i**l] = (p[1])*(x[77+i**k]) ;
y[55+i**l] = (p[1])*(x[78+i**k]) ;
y[56+i**l] = (p[1])*(x[79+i**k]) ;
y[57+i**l] = (p[1])*(x[80+i**k]) ;
y[58+i**l] = (p[1])*(x[81+i**k]) ;
y[59+i**l] = (p[1])*(x[82+i**k]) ;
y[60+i**l] = (p[1])*(x[83+i**k]) ;
y[61+i**l] = (p[1])*(x[84+i**k]) ;
y[62+i**l] = (p[1])*(x[85+i**k]) ;
y[63+i**l] = (p[1])*(x[88+i**k]) ;
y[64+i**l] = (p[1])*(x[89+i**k]) ;
y[65+i**l] = (p[1])*(x[90+i**k]) ;
y[66+i**l] = (p[1])*(x[91+i**k]) ;
y[67+i**l] = (p[1])*(x[92+i**k]) ;
y[68+i**l] = (p[1])*(x[93+i**k]) ;
y[69+i**l] = (p[1])*(x[94+i**k]) ;
y[70+i**l] = (p[1])*(x[95+i**k]) ;
y[71+i**l] = (p[1])*(x[96+i**k]) ;
y[72+i**l] = (p[1])*(x[99+i**k]) ;
y[73+i**l] = (p[1])*(x[100+i**k]) ;
y[74+i**l] = (p[1])*(x[101+i**k]) ;
y[75+i**l] = (p[1])*(x[102+i**k]) ;
y[76+i**l] = (p[1])*(x[103+i**k]) ;
y[77+i**l] = (p[1])*(x[104+i**k]) ;
y[78+i**l] = (p[1])*(x[105+i**k]) ;
y[79+i**l] = (p[1])*(x[106+i**k]) ;
y[80+i**l] = (p[1])*(x[107+i**k]) ;
y[81+i**l] = (p[1])*(x[110+i**k]) ;
y[82+i**l] = (p[1])*(x[111+i**k]) ;
y[83+i**l] = (p[1])*(x[112+i**k]) ;
y[84+i**l] = (p[1])*(x[113+i**k]) ;
y[85+i**l] = (p[1])*(x[114+i**k]) ;
y[86+i**l] = (p[1])*(x[115+i**k]) ;
y[87+i**l] = (p[1])*(x[116+i**k]) ;
y[88+i**l] = (p[1])*(x[117+i**k]) ;
y[89+i**l] = (p[1])*(x[118+i**k]) ;
y[90+i**l] = (p[1])*(x[121+i**k]) ;
y[91+i**l] = (p[1])*(x[122+i**k]) ;
y[92+i**l] = (p[1])*(x[123+i**k]) ;
y[93+i**l] = (p[1])*(x[124+i**k]) ;
y[94+i**l] = (p[1])*(x[125+i**k]) ;
y[95+i**l] = (p[1])*(x[126+i**k]) ;
y[96+i**l] = (p[1])*(x[127+i**k]) ;
y[97+i**l] = (p[1])*(x[128+i**k]) ;
y[98+i**l] = (p[1])*(x[129+i**k]) ;
y[99+i**l] = (p[1])*(x[132+i**k])+1.0 ;
y[100+i**l] = (p[1])*(x[133+i**k]) ;
y[101+i**l] = (p[1])*(x[134+i**k]) ;
y[102+i**l] = (p[1])*(x[135+i**k]) ;
y[103+i**l] = (p[1])*(x[136+i**k]) ;
y[104+i**l] = (p[1])*(x[137+i**k]) ;
y[105+i**l] = (p[1])*(x[138+i**k]) ;
y[106+i**l] = (p[1])*(x[139+i**k]) ;
y[107+i**l] = (p[1])*(x[140+i**k]) ;
y[108+i**l] = (p[1])*(x[143+i**k])+x[0+i**k] ;
y[109+i**l] = (p[1])*(x[144+i**k])+x[1+i**k] ;
y[110+i**l] = (p[1])*(x[145+i**k])+x[2+i**k] ;
y[111+i**l] = (p[1])*(x[146+i**k])+x[3+i**k] ;
y[112+i**l] = (p[1])*(x[147+i**k])+x[4+i**k] ;
y[113+i**l] = (p[1])*(x[148+i**k])+x[5+i**k] ;
y[114+i**l] = (p[1])*(x[149+i**k])+x[6+i**k] ;
y[115+i**l] = (p[1])*(x[150+i**k])+x[7+i**k] ;
y[116+i**l] = (p[1])*(x[151+i**k])+x[8+i**k] ;
y[117+i**l] = (p[1])*(x[154+i**k]) ;
y[118+i**l] = (p[1])*(x[155+i**k])+1.0 ;
y[119+i**l] = (p[1])*(x[156+i**k]) ;
y[120+i**l] = (p[1])*(x[157+i**k]) ;
y[121+i**l] = (p[1])*(x[158+i**k]) ;
y[122+i**l] = (p[1])*(x[159+i**k]) ;
y[123+i**l] = (p[1])*(x[160+i**k]) ;
y[124+i**l] = (p[1])*(x[161+i**k]) ;
y[125+i**l] = (p[1])*(x[162+i**k]) ;
y[126+i**l] = (p[1])*(x[165+i**k]) ;
y[127+i**l] = (p[1])*(x[166+i**k]) ;
y[128+i**l] = (p[1])*(x[167+i**k])+1.0 ;
y[129+i**l] = (p[1])*(x[168+i**k]) ;
y[130+i**l] = (p[1])*(x[169+i**k]) ;
y[131+i**l] = (p[1])*(x[170+i**k]) ;
y[132+i**l] = (p[1])*(x[171+i**k]) ;
y[133+i**l] = (p[1])*(x[172+i**k]) ;
y[134+i**l] = (p[1])*(x[173+i**k]) ;
y[135+i**l] = (p[1])*(x[176+i**k]) ;
y[136+i**l] = (p[1])*(x[177+i**k]) ;
y[137+i**l] = (p[1])*(x[178+i**k]) ;
y[138+i**l] = (p[1])*(x[179+i**k])+1.0 ;
y[139+i**l] = (p[1])*(x[180+i**k]) ;
y[140+i**l] = (p[1])*(x[181+i**k]) ;
y[141+i**l] = (p[1])*(x[182+i**k]) ;
y[142+i**l] = (p[1])*(x[183+i**k]) ;
y[143+i**l] = (p[1])*(x[184+i**k]) ;
y[144+i**l] = (p[1])*(x[187+i**k]) ;
y[145+i**l] = (p[1])*(x[188+i**k]) ;
y[146+i**l] = (p[1])*(x[189+i**k]) ;
y[147+i**l] = (p[1])*(x[190+i**k]) ;
y[148+i**l] = (p[1])*(x[191+i**k])+1.0 ;
y[149+i**l] = (p[1])*(x[192+i**k]) ;
y[150+i**l] = (p[1])*(x[193+i**k]) ;
y[151+i**l] = (p[1])*(x[194+i**k]) ;
y[152+i**l] = (p[1])*(x[195+i**k]) ;
y[153+i**l] = (p[1])*(x[198+i**k]) ;
y[154+i**l] = (p[1])*(x[199+i**k]) ;
y[155+i**l] = (p[1])*(x[200+i**k]) ;
y[156+i**l] = (p[1])*(x[201+i**k]) ;
y[157+i**l] = (p[1])*(x[202+i**k]) ;
y[158+i**l] = (p[1])*(x[203+i**k])+1.0 ;
y[159+i**l] = (p[1])*(x[204+i**k]) ;
y[160+i**l] = (p[1])*(x[205+i**k]) ;
y[161+i**l] = (p[1])*(x[206+i**k]) ;
y[162+i**l] = (p[1])*(x[209+i**k]) ;
y[163+i**l] = (p[1])*(x[210+i**k]) ;
y[164+i**l] = (p[1])*(x[211+i**k]) ;
y[165+i**l] = (p[1])*(x[212+i**k]) ;
y[166+i**l] = (p[1])*(x[213+i**k]) ;
y[167+i**l] = (p[1])*(x[214+i**k]) ;
y[168+i**l] = (p[1])*(x[215+i**k])+1.0 ;
y[169+i**l] = (p[1])*(x[216+i**k]) ;
y[170+i**l] = (p[1])*(x[217+i**k]) ;
y[171+i**l] = (p[1])*(x[220+i**k]) ;
y[172+i**l] = (p[1])*(x[221+i**k]) ;
y[173+i**l] = (p[1])*(x[222+i**k]) ;
y[174+i**l] = (p[1])*(x[223+i**k]) ;
y[175+i**l] = (p[1])*(x[224+i**k]) ;
y[176+i**l] = (p[1])*(x[225+i**k]) ;
y[177+i**l] = (p[1])*(x[226+i**k]) ;
y[178+i**l] = (p[1])*(x[227+i**k])+1.0 ;
y[179+i**l] = (p[1])*(x[228+i**k]) ;
y[180+i**l] = (p[1])*(x[231+i**k]) ;
y[181+i**l] = (p[1])*(x[232+i**k]) ;
y[182+i**l] = (p[1])*(x[233+i**k]) ;
y[183+i**l] = (p[1])*(x[234+i**k]) ;
y[184+i**l] = (p[1])*(x[235+i**k]) ;
y[185+i**l] = (p[1])*(x[236+i**k]) ;
y[186+i**l] = (p[1])*(x[237+i**k]) ;
y[187+i**l] = (p[1])*(x[238+i**k]) ;
y[188+i**l] = (p[1])*(x[239+i**k])+1.0 ;
y[189+i**l] = (p[1])*(x[242+i**k]) ;
y[190+i**l] = (p[1])*(x[243+i**k]) ;
y[191+i**l] = (p[1])*(x[244+i**k]) ;
y[192+i**l] = (p[1])*(x[245+i**k]) ;
y[193+i**l] = (p[1])*(x[246+i**k]) ;
y[194+i**l] = (p[1])*(x[247+i**k]) ;
y[195+i**l] = (p[1])*(x[248+i**k]) ;
y[196+i**l] = (p[1])*(x[249+i**k]) ;
y[197+i**l] = (p[1])*(x[250+i**k]) ;
y[198+i**l] = (p[1])*(x[253+i**k]) ;
y[199+i**l] = (p[1])*(x[254+i**k]) ;
y[200+i**l] = (p[1])*(x[255+i**k]) ;
y[201+i**l] = (p[1])*(x[256+i**k]) ;
y[202+i**l] = (p[1])*(x[257+i**k]) ;
y[203+i**l] = (p[1])*(x[258+i**k]) ;
y[204+i**l] = (p[1])*(x[259+i**k]) ;
y[205+i**l] = (p[1])*(x[260+i**k]) ;
y[206+i**l] = (p[1])*(x[261+i**k]) ;
y[207+i**l] = (p[1])*(x[264+i**k]) ;
y[208+i**l] = (p[1])*(x[265+i**k]) ;
y[209+i**l] = (p[1])*(x[266+i**k]) ;
y[210+i**l] = (p[1])*(x[267+i**k]) ;
y[211+i**l] = (p[1])*(x[268+i**k]) ;
y[212+i**l] = (p[1])*(x[269+i**k]) ;
y[213+i**l] = (p[1])*(x[270+i**k]) ;
y[214+i**l] = (p[1])*(x[271+i**k]) ;
y[215+i**l] = (p[1])*(x[272+i**k]) ;
y[216+i**l] = (p[1])*(x[275+i**k]) ;
y[217+i**l] = (p[1])*(x[276+i**k]) ;
y[218+i**l] = (p[1])*(x[277+i**k]) ;
y[219+i**l] = (p[1])*(x[278+i**k]) ;
y[220+i**l] = (p[1])*(x[279+i**k]) ;
y[221+i**l] = (p[1])*(x[280+i**k]) ;
y[222+i**l] = (p[1])*(x[281+i**k]) ;
y[223+i**l] = (p[1])*(x[282+i**k]) ;
y[224+i**l] = (p[1])*(x[283+i**k]) ;
y[225+i**l] = (p[1])*(x[286+i**k]) ;
y[226+i**l] = (p[1])*(x[287+i**k]) ;
y[227+i**l] = (p[1])*(x[288+i**k]) ;
y[228+i**l] = (p[1])*(x[289+i**k]) ;
y[229+i**l] = (p[1])*(x[290+i**k]) ;
y[230+i**l] = (p[1])*(x[291+i**k]) ;
y[231+i**l] = (p[1])*(x[292+i**k]) ;
y[232+i**l] = (p[1])*(x[293+i**k]) ;
y[233+i**l] = (p[1])*(x[294+i**k]) ;
y[234+i**l] = (p[1])*(x[297+i**k]) ;
y[235+i**l] = (p[1])*(x[298+i**k]) ;
y[236+i**l] = (p[1])*(x[299+i**k]) ;
y[237+i**l] = (p[1])*(x[300+i**k]) ;
y[238+i**l] = (p[1])*(x[301+i**k]) ;
y[239+i**l] = (p[1])*(x[302+i**k]) ;
y[240+i**l] = (p[1])*(x[303+i**k]) ;
y[241+i**l] = (p[1])*(x[304+i**k]) ;
y[242+i**l] = (p[1])*(x[305+i**k]) ;
y[243+i**l] = (p[1])*(x[308+i**k]) ;
y[244+i**l] = (p[1])*(x[309+i**k]) ;
y[245+i**l] = (p[1])*(x[310+i**k]) ;
y[246+i**l] = (p[1])*(x[311+i**k]) ;
y[247+i**l] = (p[1])*(x[312+i**k]) ;
y[248+i**l] = (p[1])*(x[313+i**k]) ;
y[249+i**l] = (p[1])*(x[314+i**k]) ;
y[250+i**l] = (p[1])*(x[315+i**k]) ;
y[251+i**l] = (p[1])*(x[316+i**k]) ;
y[252+i**l] = (p[1])*(x[319+i**k]) ;
y[253+i**l] = (p[1])*(x[320+i**k]) ;
y[254+i**l] = (p[1])*(x[321+i**k]) ;
y[255+i**l] = (p[1])*(x[322+i**k]) ;
y[256+i**l] = (p[1])*(x[323+i**k]) ;
y[257+i**l] = (p[1])*(x[324+i**k]) ;
y[258+i**l] = (p[1])*(x[325+i**k]) ;
y[259+i**l] = (p[1])*(x[326+i**k]) ;
y[260+i**l] = (p[1])*(x[327+i**k]) ;
y[261+i**l] = (p[1])*(x[330+i**k]) ;
y[262+i**l] = (p[1])*(x[331+i**k]) ;
y[263+i**l] = (p[1])*(x[332+i**k]) ;
y[264+i**l] = (p[1])*(x[333+i**k]) ;
y[265+i**l] = (p[1])*(x[334+i**k]) ;
y[266+i**l] = (p[1])*(x[335+i**k]) ;
y[267+i**l] = (p[1])*(x[336+i**k]) ;
y[268+i**l] = (p[1])*(x[337+i**k]) ;
y[269+i**l] = (p[1])*(x[338+i**k]) ;
y[270+i**l] = (p[1])*(x[341+i**k]) ;
y[271+i**l] = (p[1])*(x[342+i**k]) ;
y[272+i**l] = (p[1])*(x[343+i**k]) ;
y[273+i**l] = (p[1])*(x[344+i**k]) ;
y[274+i**l] = (p[1])*(x[345+i**k]) ;
y[275+i**l] = (p[1])*(x[346+i**k]) ;
y[276+i**l] = (p[1])*(x[347+i**k]) ;
y[277+i**l] = (p[1])*(x[348+i**k]) ;
y[278+i**l] = (p[1])*(x[349+i**k]) ;
y[279+i**l] = (p[1])*(x[352+i**k]) ;
y[280+i**l] = (p[1])*(x[353+i**k]) ;
y[281+i**l] = (p[1])*(x[354+i**k]) ;
y[282+i**l] = (p[1])*(x[355+i**k]) ;
y[283+i**l] = (p[1])*(x[356+i**k]) ;
y[284+i**l] = (p[1])*(x[357+i**k]) ;
y[285+i**l] = (p[1])*(x[358+i**k]) ;
y[286+i**l] = (p[1])*(x[359+i**k]) ;
y[287+i**l] = (p[1])*(x[360+i**k]) ;
y[288+i**l] = (p[1])*(x[363+i**k]) ;
y[289+i**l] = (p[1])*(x[364+i**k]) ;
y[290+i**l] = (p[1])*(x[365+i**k]) ;
y[291+i**l] = (p[1])*(x[366+i**k]) ;
y[292+i**l] = (p[1])*(x[367+i**k]) ;
y[293+i**l] = (p[1])*(x[368+i**k]) ;
y[294+i**l] = (p[1])*(x[369+i**k]) ;
y[295+i**l] = (p[1])*(x[370+i**k]) ;
y[296+i**l] = (p[1])*(x[371+i**k]) ;
y[297+i**l] = (p[1])*(x[374+i**k]) ;
y[298+i**l] = (p[1])*(x[375+i**k]) ;
y[299+i**l] = (p[1])*(x[376+i**k]) ;
y[300+i**l] = (p[1])*(x[377+i**k]) ;
y[301+i**l] = (p[1])*(x[378+i**k]) ;
y[302+i**l] = (p[1])*(x[379+i**k]) ;
y[303+i**l] = (p[1])*(x[380+i**k]) ;
y[304+i**l] = (p[1])*(x[381+i**k]) ;
y[305+i**l] = (p[1])*(x[382+i**k]) ;
y[306+i**l] = (p[1])*(x[385+i**k]) ;
y[307+i**l] = (p[1])*(x[386+i**k]) ;
y[308+i**l] = (p[1])*(x[387+i**k]) ;
y[309+i**l] = (p[1])*(x[388+i**k]) ;
y[310+i**l] = (p[1])*(x[389+i**k]) ;
y[311+i**l] = (p[1])*(x[390+i**k]) ;
y[312+i**l] = (p[1])*(x[391+i**k]) ;
y[313+i**l] = (p[1])*(x[392+i**k]) ;
y[314+i**l] = (p[1])*(x[393+i**k]) ;
y[315+i**l] = (p[1])*(x[396+i**k]) ;
y[316+i**l] = (p[1])*(x[397+i**k]) ;
y[317+i**l] = (p[1])*(x[398+i**k]) ;
y[318+i**l] = (p[1])*(x[399+i**k]) ;
y[319+i**l] = (p[1])*(x[400+i**k]) ;
y[320+i**l] = (p[1])*(x[401+i**k]) ;
y[321+i**l] = (p[1])*(x[402+i**k]) ;
y[322+i**l] = (p[1])*(x[403+i**k]) ;
y[323+i**l] = (p[1])*(x[404+i**k]) ;
y[324+i**l] = (p[1])*(x[407+i**k]) ;
y[325+i**l] = (p[1])*(x[408+i**k]) ;
y[326+i**l] = (p[1])*(x[409+i**k]) ;
y[327+i**l] = (p[1])*(x[410+i**k]) ;
y[328+i**l] = (p[1])*(x[411+i**k]) ;
y[329+i**l] = (p[1])*(x[412+i**k]) ;
y[330+i**l] = (p[1])*(x[413+i**k]) ;
y[331+i**l] = (p[1])*(x[414+i**k]) ;
y[332+i**l] = (p[1])*(x[415+i**k]) ;
y[333+i**l] = (p[1])*(x[418+i**k]) ;
y[334+i**l] = (p[1])*(x[419+i**k]) ;
y[335+i**l] = (p[1])*(x[420+i**k]) ;
y[336+i**l] = (p[1])*(x[421+i**k]) ;
y[337+i**l] = (p[1])*(x[422+i**k]) ;
y[338+i**l] = (p[1])*(x[423+i**k]) ;
y[339+i**l] = (p[1])*(x[424+i**k]) ;
y[340+i**l] = (p[1])*(x[425+i**k]) ;
y[341+i**l] = (p[1])*(x[426+i**k]) ;
y[342+i**l] = (p[1])*(x[429+i**k]) ;
y[343+i**l] = (p[1])*(x[430+i**k]) ;
y[344+i**l] = (p[1])*(x[431+i**k]) ;
y[345+i**l] = (p[1])*(x[432+i**k]) ;
y[346+i**l] = (p[1])*(x[433+i**k]) ;
y[347+i**l] = (p[1])*(x[434+i**k]) ;
y[348+i**l] = (p[1])*(x[435+i**k]) ;
y[349+i**l] = (p[1])*(x[436+i**k]) ;
y[350+i**l] = (p[1])*(x[437+i**k]) ; 
}
}