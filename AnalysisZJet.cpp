#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"
#include <fstream>


using namespace std;
using namespace HepMC;
namespace Rivet {


  /// @brief MC validation analysis for Z + jets events
  class AnalysisZJet : public MC_JetAnalysis {
  public:

    /// Default constructor
    AnalysisZJet(string name = "AnalysisZJet")
    : MC_JetAnalysis(name, 4, "Jets")
	  {	  }


    /// @name Analysis methods
    //@{

    /// Initialize
    void init() 
    {
      i=1;
		  _dR=0.2;
      if (getOption("SCHEME") == "BARE")  _dR = 0.0;
		  _lepton=PID::ELECTRON;
      if (getOption("LMODE") == "MU")  _lepton = PID::MUON;

      FinalState fs;
      Cut cut = Cuts::abseta < 2.0 && Cuts::pT > 60*GeV; //etamax=3.5, originalmente
      //Encontrar partículas Z, na forma de 2 muoes
      ZFinder zfinder(fs, cut, PID::MUON, 65*GeV, 115*GeV, _dR, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES); 
      declare(zfinder, "ZFinder");
      FastJets jetpro(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.4);
      declare(jetpro, "Jets");

      book(_h_Z_jet1_deta ,"Z_jet1_deta", 50, -5, 5);
      book(_h_Z_jet1_dR ,"Z_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
      ZJetdata.open("ZJetdata10e3C0-5Njob0Pt45Jet60Z.dat");
      ZJetdata << "Angulo jprodr weight PtZ PtJet\n";
    }



    /// Do the analysis
    void analyze(const Event & e) {
      //cout<<"Evento "<<i<<endl;
      //++i;
      MSG_TRACE("AnalysisZJet: running ZFinder");
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) {
        //cout<<"Vetou"<<endl;
        vetoEvent;
      }

      const FourMomentum& zmom = zfinder.bosons()[0].momentum();
      //cout<<zmom.pT()<<endl;
      //cout<<"valor"<<i<<endl;
      MSG_TRACE("AnalysisZJet: have exactly one Z boson candidate");
      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(45*GeV); 
      //FourMomentum pmax(0,0,0,0);
      //cout<<"momento"<<jets[0].E()<<endl;
      if (jets.size() > 0) {
        MSG_TRACE("AnalysisZJet: have at least one valid jet");
        //preenche histogramas
        _h_Z_jet1_deta->fill(zmom.eta()-jets[0].eta()); 
        _h_Z_jet1_dR->fill(deltaR(zmom, jets[0].momentum()));
      }


    /*  if(jets.size() == 0){
        cout<<"Não se detetaram jets!"<<endl;
      }
      else{
      //FourMomentum pmax = jets[0].momentum();
      long unsigned int j=0;
      while(j<jets.size()){
        if(pmax.E()<jets[j].momentum().E()){
          pmax = jets[j].momentum();
        }
        ++j;
      }
      //cout<<"Momento de Jet mais energetico e Z  "<< pmax <<" "<< zmom <<endl;
      }
*/
      //cout<<"Peso do evento "<<e.weights()[0]<<endl;
       if(jets.size() == 0)
       {
        vetoEvent;
      }  //Garante que só aceito evento se houver jets
      
      else{
        //cout<<"Momento de Jet mais energetico e Z  "<< jets[0].momentum()<<" "<< zmom <<endl;
      }
      MC_JetAnalysis::analyze(e);

      //Definicao de variáveis uteis
      double angle = 0;
      double jprodr = 0;
      const HeavyIon *hi = e.genEvent()->heavy_ion();
      if(hi) {
        angle = hi->event_plane_angle();
        jprodr = hi->eccentricity();
        //cout << hi->event_plane_angle() << " " << hi->eccentricity() << endl;
      }

      double event_weight = e.weights()[0];
      double Zpt = zfinder.bosons()[0].pT();
      double jetpt = jets[0].E();
      
      //escrita em ficheiro
      ZJetdata << angle << " " << jprodr << " " << event_weight << " " << Zpt << " "<< jetpt <<"\n";
    }


    /// Finalize
    void finalize() {
      scale(_h_Z_jet1_deta, crossSection()/picobarn/sumOfWeights());
      scale(_h_Z_jet1_dR, crossSection()/picobarn/sumOfWeights());
      MC_JetAnalysis::finalize();
      ZJetdata.close(); //fechar o ficheiro
    }

    //@}


  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    //@{
    double _dR;
    PdgId _lepton;
    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Z_jet1_deta;
    Histo1DPtr _h_Z_jet1_dR;
    ofstream ZJetdata; //inicializar ficheiro
    //@}
    int i;

  };

  // The hooks for the plugin system
  DECLARE_RIVET_PLUGIN(AnalysisZJet);
}