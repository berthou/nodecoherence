typedef struct {
	double x;
	double y;
}point;

point M1N_table[41] = {
	{0,0},
	{0.75,0.2558519645504121},
	{1.5,0.7615669471444275},
	{2.25,0.9995922605850198},
	{3.,0.7263309286461761},
	{3.75,0.2214407457767586},
	{4.5,0.001630292654197285},
	{5.25,0.291855314318987},
	{6.,0.7950972429526371},
	{6.75,0.9963343342149172},
	{7.5,0.6896189675463662},
	{8.25,0.1888460594267805},
	{9.,0.006510539200236254},
	{9.75,0.32921601109586485},
	{10.5,0.8267031592910089},
	{11.25,0.9898397269684284},
	{12.,0.6516704688070082},
	{12.75,0.1582804610113671},
	{13.5,0.014608914717737842},
	{14.25,0.3676904194030031},
	{15.,0.8561785885865996},
	{15.75,0.980150791287497},
	{16.5,0.6127329010630297},
	{17.25,0.12994327401278832},
	{18.,0.02587260831823429},
	{18.75,0.407027641059451},
	{19.5,0.8833313165359684},
	{20.25,0.9673307103747926},
	{21.,0.5730601828370937},
	{21.75,0.10401929006226401},
	{22.5,0.040228167533980785},
	{23.25,0.4469711513311938},
	{24.,0.9079842755674485},
	{24.75,0.9514630861652693},
	{25.5,0.5329110266935851},
	{26.25,0.08067756388219965},
	{27.,0.05758197731403136},
	{27.75,0.48726047177271764},
	{28.5,0.9299766995289873},
	{29.25,0.9326513941436805},
	{30.,0.4925472521298306},
};

point M2N1_table[41] = {
	{0.,0.},
	{0.75,0.13415556556308986},
	{1.5,0.4646313991661486},
	{2.25,0.8140868113613694},
	{3.,0.9949962483002229},
	{3.75,0.9102796786697807},
	{4.5,0.6053978997153903},
	{5.25,0.2439572613790798},
	{6.,0.01991485667481696},
	{6.75,0.053496827655461404},
	{7.5,0.32668234108248667},
	{8.25,0.6928739687261106},
	{9.,0.9555651309423386},
	{9.75,0.9737899019889971},
	{10.5,0.7377684639979968},
	{11.25,0.3741551749641241},
	{12.,0.0780730206337546},
	{12.75,0.008406276476204039},
	{13.5,0.2025396683450535},
	{14.25,0.5562968963169169},
	{15.,0.8798439564294102},
	{15.75,0.9995582933398672},
	{16.5,0.8511985287513576},
	{17.25,0.5143778151645948},
	{18.,0.16984164587796058},
	{18.75,0.0024757994818104534},
	{19.5,0.10209251509302733},
	{20.25,0.41523524542217},
	{21.,0.7738646301121337},
	{21.75,0.9855321574404183},
	{22.5,0.9366523200467594},
	{23.25,0.6534551268518676},
	{24.,0.2879104963315033},
	{24.75,0.03617781506970849},
	{25.5,0.033342443968038704},
	{26.25,0.28092590628139624},
	{27.,0.6460694043669164},
	{27.75,0.9328288082352075},
	{28.5,0.9873226378603295},
	{29.25,0.7803082911100546},
	{30.,0.4228742750562098},
};

point M2N10_table[41] = {
	{0,0},
	{0.75,0.1334015856765857},
	{1.5,0.4532583551296518},
	{2.25,0.774422402614494},
	{3.,0.9582162172585719},
	{3.75,0.9882440634055067},
	{4.5,0.8710735957379998},
	{5.25,0.6018562333045654},
	{6.,0.26447851638871916},
	{6.75,0.03543638775454617},
	{7.5,0.05140164068666749},
	{8.25,0.30378940827565215},
	{9.,0.6412641760543925},
	{9.75,0.890978136729808},
	{10.5,0.9859273022223775},
	{11.25,0.9289486697278156},
	{12.,0.7223477795006722},
	{12.75,0.4116347552211995},
	{13.5,0.13053709020468995},
	{14.25,0.036525241861020644},
	{15.,0.18840062319687162},
	{15.75,0.49379571684773416},
	{16.5,0.7852197914016029},
	{17.25,0.948177384385688},
	{18.,0.9526051951304715},
	{18.75,0.8065748028882014},
	{19.5,0.5447038925833655},
	{20.25,0.25680420301353274},
	{21.,0.08455377038024349},
	{21.75,0.12986700241078475},
	{22.5,0.3646208643062731},
	{23.25,0.6568782755293252},
	{24.,0.872859417942188},
	{24.75,0.9413237952935949},
	{25.5,0.8561486546475862},
	{26.25,0.6494505482887466},
	{27.,0.38441338352730614},
	{27.75,0.1734412753298811},
	{28.5,0.13047415125137501},
	{29.25,0.27808250202734025},
	{30.,0.5330002297522931},
};

point M3N1_table[41] = {
	{0.,0.0},
	{0.75,0.06867992295466127},
	{1.5,0.25585196455041237},
	{2.25,0.5100962792030017},
	{3.,0.7615669471444274},
	{3.75,0.9411800239763921},
	{4.5,0.9995922605850198},
	{5.25,0.9207566653307506},
	{6.,0.7263309286461767},
	{6.75,0.4697276289948972},
	{7.5,0.22144074577675882},
	{8.25,0.04967957503206141},
	{9.,0.0016302926541974655},
	{9.75,0.09049298269012045},
	{10.5,0.2918553143189865},
	{11.25,0.5503990898518292},
	{12.,0.7950972429526364},
	{12.75,0.9587263724129644},
	{13.5,0.9963343342149172},
	{14.25,0.8975894806823362},
	{15.,0.6896189675463669},
	{15.75,0.42955639008309215},
	{16.5,0.18884605942678098},
	{17.25,0.03361584343290294},
	{18.,0.006510539200236171},
	{18.75,0.11497650755422996},
	{19.5,0.32921601109586385},
	{20.25,0.5903732394368025},
	{21.,0.8267031592910085},
	{21.75,0.9732812879086107},
	{22.5,0.9898397269684287},
	{23.25,0.8718295471949504},
	{24.,0.6516704688070094},
	{24.75,0.38984452597041647},
	{25.5,0.15828046101136758},
	{26.25,0.02059348249168605},
	{27.,0.014608914717737675},
	{27.75,0.14197083630405166},
	{28.5,0.3676904194030005},
	{29.25,0.6297580497082136},
	{30.,0.8561785885865989},

};


point M3N5_table[41] = {
	{0.,0.0},
	{0.75,0.0682139624264106},
	{1.5,0.24872780841611797},
	{2.25,0.4789746944930624},
	{3.,0.6906510880862958},
	{3.75,0.8468341124937991},
	{4.5,0.9431782560103625},
	{5.25,0.9854533127702823},
	{6.,0.9781260635907736},
	{6.75,0.9176881205110966},
	{7.5,0.7936646228188851},
	{8.25,0.6068227966444502},
	{9.,0.3867412190927386},
	{9.75,0.19086027248191734},
	{10.5,0.07978769198446084},
	{11.25,0.08609892006494652},
	{12.,0.20102397805973},
	{12.75,0.37994884765692893},
	{13.5,0.5658557917593026},
	{14.25,0.7193111541367545},
	{15.,0.8275378659525284},
	{15.75,0.8900920248391855},
	{16.5,0.9066715499054412},
	{17.25,0.8754347776998737},
	{18.,0.7911142223935939},
	{18.75,0.6556413692084813},
	{19.5,0.4943484637728258},
	{20.25,0.34788777985207564},
	{21.,0.25342063941843374},
	{21.75,0.23459006744793717},
	{22.5,0.291621021649377},
	{23.25,0.3990677109515566},
	{24.,0.521436591242608},
	{24.75,0.6324943182352487},
	{25.5,0.7186519756381243},
	{26.25,0.7713858315840885},
	{27.,0.7850259416069214},
	{27.75,0.7561765793811756},
	{28.5,0.6849679488852258},
	{29.25,0.5844586510404053},
	{30.,0.4842216119680105},
};