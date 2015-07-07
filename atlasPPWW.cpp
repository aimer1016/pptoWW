#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <fastjet/ClusterSequence.hh>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <TStyle.h>
#include <TLatex.h>
//#include <TLorentzVector.h>
//#include <TLorentzRotation.h>


using namespace std;
using namespace fastjet;

const PseudoJet zero(0.,0.,0.,0.);


double deltaR(PseudoJet *obj1, PseudoJet *obj2)
{
  double deltaPhi=fabs(obj1->phi()-obj2->phi());
  if (deltaPhi > M_PI) deltaPhi=2*M_PI-deltaPhi;
  double deltaEta=fabs(obj1->eta()-obj2->eta());
  return sqrt(pow(deltaPhi,2)+pow(deltaEta,2));
}

void atlasppww(char rootFileLocation[]) // ATLAS-COM-CONF-2014-033
{


  Long64_t i,entryCount;
  Int_t j;
  unsigned int a,b,c;


  TFile *file = TFile::Open(rootFileLocation);
  TTree *tree = (TTree *)file->Get("dshepmc");
  /*
  cout << endl;
  ClusterSequence::print_banner();
  cout << endl;
  
  cout << endl << rootFileLocation
       << endl;
  */
  entryCount = tree->GetEntries();
  /*
  cout << endl << "entryCount = " << entryCount
       << endl;
  */  

  // Open file where to save the histograms
  TString output_name = "Output/";
  output_name += "hist_";
  output_name += rootFileLocation;
  TFile out (output_name,"RECREATE");

  // Initialize all the histograms
  TH1 *fptl1 = new TH1F("fptl1","p_{T}^{#it{l}_{1}}",23,25,140);
  TH1 *fptl2 = new TH1F("fptl2","p_{T}^{#it{l}_{2}}",16,20,100);
  TH1 *fptll = new TH1F("fptll","p_{T}^{#it{ll}}",20,20,120);
  TH1 *fptllem = new TH1F("fptllem","p_{T}^{#it{ll}}",24,0,120);
  TH1 *fdpll = new TH1F("fdpll","#Delta #phi (#it{ll})",20,0,3.2);
  TH1 *fmtllmet = new TH1F("fmtllmet","m_{T} (#it{ll} + E_{T}^{miss})",20,100,300);
  TH1 *fmtllmetem = new TH1F("fmtllmetem","m_{T} (#it{ll} + E_{T}^{miss})",23,70,300);
  TH1 *fptllmet = new TH1F("fptllmet","p_{T} (#it{ll} + E_{T}^{miss})",14,0,70);

  TH1 *fnu1 = new TH1F("fnu1","p_{T}^{#nu_{1}}",20,0,100);
  TH1 *fnu2 = new TH1F("fnu2","p_{T}^{#nu_{2}}",20,0,100);
  TH1 *fnu = new TH1F("fnu","p_{T} (E_{T}^{miss})",20,0,100);
  //TH1 *fne = new TH1F("fne","p_{T}^{#nu_{e}}",20,0,100);
  TH1 *fmll = new TH1F("fmll","mll",29,10,300);
  TH1 *fetmll = new TH1F("fetmll","EtmissRel",20,0,100);
  TH1 *fdr = new TH1F("fdr","DeltaR",20,0,1);


  Double_t maxWeight=0.;
  Int_t maxParticleCount=0;
  Double_t eventWeight=0.;
  Int_t eventParticleCount=0;
  Int_t eventProcessVertex=0;
  TBranch *eventWeightBranch=0;
  TBranch *eventParticleCountBranch=0;
  TBranch *eventProcessVertexBranch=0;
  
  tree->SetBranchAddress("eventWeight", &eventWeight, &eventWeightBranch);
  tree->SetBranchAddress("eventParticleCount", &eventParticleCount, &eventParticleCountBranch);
  tree->SetBranchAddress("eventProcessVertex", &eventProcessVertex, &eventProcessVertexBranch);
  
  for (i=0; i<entryCount; i++)
  {
    eventWeightBranch->GetEntry(i);
    eventParticleCountBranch->GetEntry(i);
    if (eventWeight > maxWeight) maxWeight=eventWeight;
    if (eventParticleCount > maxParticleCount) maxParticleCount=eventParticleCount;
  }
  
  Double_t minWeight=maxWeight;
  Int_t minParticleCount=maxParticleCount;
  
  for (i=0; i<entryCount; i++)
  {
    eventWeightBranch->GetEntry(i);
    eventParticleCountBranch->GetEntry(i);
    if (eventWeight < minWeight) minWeight=eventWeight;
    if (eventParticleCount < minParticleCount) minParticleCount=eventParticleCount;
  }
  /*
  cout << endl << "minParticleCount = " << minParticleCount
       << endl << "maxParticleCount = " << maxParticleCount
       << endl
       << endl << "minWeight = " << minWeight
       << endl << "maxWeight = " << maxWeight
       << endl;
  */
  vector<Int_t> *particleStatus=0;
  vector<Int_t> *particleIden=0;
  vector<Int_t> *particlePastVertex=0;
  vector<Int_t> *particleFutureVertex=0;
  vector<Double_t> *particlePx=0;
  vector<Double_t> *particlePy=0;
  vector<Double_t> *particlePz=0;
  vector<Double_t> *particleEnergy=0;
  
  TBranch *particleStatusBranch=0;
  TBranch *particleIdenBranch=0;
  TBranch *particlePastVertexBranch=0;
  TBranch *particleFutureVertexBranch=0;
  TBranch *particlePxBranch=0;
  TBranch *particlePyBranch=0;
  TBranch *particlePzBranch=0;
  TBranch *particleEnergyBranch=0;
  
  tree->SetBranchAddress("particleStatus", &particleStatus, &particleStatusBranch);
  tree->SetBranchAddress("particleIden", &particleIden, &particleIdenBranch);
  tree->SetBranchAddress("particlePastVertex", &particlePastVertex, &particlePastVertexBranch);
  tree->SetBranchAddress("particleFutureVertex", &particleFutureVertex, &particleFutureVertexBranch);
  tree->SetBranchAddress("particlePx", &particlePx, &particlePxBranch);
  tree->SetBranchAddress("particlePy", &particlePy, &particlePyBranch);
  tree->SetBranchAddress("particlePz", &particlePz, &particlePzBranch);
  tree->SetBranchAddress("particleEnergy", &particleEnergy, &particleEnergyBranch);

  double jetR=.4;
  JetDefinition jetDef(antikt_algorithm, jetR);
  // PseudoJet is a fastjet object
  // it has a fourmomentum and an integer called user_index
  // we will use this object for all particles, not just jets
  PseudoJet particle;
  vector<PseudoJet> pseudojets; // the actual pseudojets that go into the jetfinder
  vector<PseudoJet> jets;
  vector<PseudoJet> leptons;
  vector<PseudoJet> electrons;
  vector<PseudoJet> muons;
  vector<PseudoJet> photons;
  PseudoJet wboson;
  vector<PseudoJet> neutrinos;
  //PseudoJet miss;
  vector<PseudoJet> rLeptons;
  vector<PseudoJet> rElectrons;
  vector<PseudoJet> rMuons;
  vector<PseudoJet> rJets;
  // "r" stands for "reconstrcuted"
  //vector<PseudoJet> totalLeptonPT;
  
  Double_t eventWeightSum=0.;
  
  
  
  // each entry in an array is a different final state category (ee/mm/em or me)
  Double_t atlcut[3]={0.,0.,0.}; 

  // SM results including H->WW : central value of AWW from Table 5
  Double_t smeff[3]={0.025/0.291,0.044/0.471,0.116/0.511};

  // limits
  for (i=0; i<entryCount; 
       i++)
    {
      pseudojets.clear();
      jets.clear();
      leptons.clear();
      electrons.clear();
      muons.clear();
      photons.clear();
      
      neutrinos.clear();
      //miss=zero;
    
      //int leptonFromWCount=0;
    
      eventWeightBranch->GetEntry(i);
      eventParticleCountBranch->GetEntry(i);
      eventProcessVertexBranch->GetEntry(i);
      particleStatusBranch->GetEntry(i);
      particleIdenBranch->GetEntry(i);
      particlePastVertexBranch->GetEntry(i);
      particleFutureVertexBranch->GetEntry(i);
      particlePxBranch->GetEntry(i);
      particlePyBranch->GetEntry(i);
      particlePzBranch->GetEntry(i);
      particleEnergyBranch->GetEntry(i);
    
      eventWeightSum += eventWeight;

      
      for (j=0; j<eventParticleCount; j++)
	{
	  if (particleStatus->at(j) == 1)
	    {
	      particle.reset(particlePx->at(j), particlePy->at(j), particlePz->at(j), particleEnergy->at(j));
	      particle.set_user_index(particleIden->at(j)); // user_index is used to store the PID number
	      // find electrons
	      if (   abs(particleIden->at(j)) == 11
		     && ((fabs(particle.eta()) < 2.47 && fabs(particle.eta()) > 1.52) || fabs(particle.eta()) < 1.37) 
		     )
		{
		  leptons.push_back(particle);
		  electrons.push_back(particle);
		  pseudojets.push_back(particle);
		}
	      // find muons
	      else if (   abs(particleIden->at(j)) == 13
			  && fabs(particle.eta()) < 2.4)
		{
		  leptons.push_back(particle);
		  muons.push_back(particle);
		}
	      // find missing momentum
	      else if (  abs(particleIden->at(j)) == 12
			  || abs(particleIden->at(j)) == 14
			  || abs(particleIden->at(j)) == 16
			 //&& particle.pt() >= 20.
			 ) //miss+=particle;
		neutrinos.push_back(particle);
	      // find photons (QED radiation)
	      else if (  abs(particleIden->at(j)) == 22 )
		{
		  photons.push_back(particle);
		  pseudojets.push_back(particle);
		}
	      
	      else pseudojets.push_back(particle);
	    }

	} // end of event particle count
      
      
    
      // find jets
      // sort vectors by pt
      ClusterSequence cs(pseudojets, jetDef);
       
      jets = sorted_by_pt(cs.inclusive_jets());      
      
      leptons = sorted_by_pt(leptons);
            
      electrons = sorted_by_pt(electrons);
      muons = sorted_by_pt(muons);
      photons = sorted_by_pt(photons);

      neutrinos = sorted_by_pt(neutrinos);

      
      //lepton triggers
      if (leptons.size() < 2) continue;
      //single lepton trigger
      if (leptons[0].pt() < 24.) continue;  // trigger conditions are "or" : this is the weakest condition for a single lepton trigger
      //dilepton triggers
      if (electrons.size() == 2 && electrons[1].pt() < 12. ) continue;
      if (muons.size() == 2 && (muons[0].pt() < 18. || muons[1].pt() < 8.) ) continue;
      if (electrons.size() == 1 && muons.size() == 1 && (electrons[0].pt() < 12. || muons[0].pt() < 8.) ) continue;

      //single lepton triggers
      // if (leptons[0].pt() < 24.) continue;
      //if (electrons.size() >=1 && electrons[0].pt() < 60.) continue;
      //if (muons.size() >=1 && muons[0].pt() < 36.) continue;  // The key cut to reduce the VLL signal h -> e_4 \mu -> W \mu \nu_\mu -> e \mu \nu_\mu \nu_e
            
       
      // lepton isolation (very rough : instead of the detail track based isolation)
      for (a=0; a<leptons.size(); a++)
	{
	  for (b=a+1; b<leptons.size();)
            {
              if ( deltaR(&(*(leptons.begin()+a)),&(*(leptons.begin()+b))) < .3) leptons.erase(leptons.begin()+b);
              else b++;
            }
	}
      // important!!
      // a is not incremented if entry a is erased
      // since a will then correspond to what was the next entry
                
      
      // efficiencies loops.
      // the next three loops run over all combinations of leptons being kept or lost,
      // the eventWeight being multiplied by the correct efficiency in each case
      //unsigned int noltk; // number of leptons to keep
      //vector<unsigned int> leptonNumbers; // (to keep)
      //for (noltk=2; noltk<=leptons.size(); noltk++)
      //{
	  // keep the first ones on the first iteration of the central efficiecies loop (cel)
	  //leptonNumbers.clear();
	  //for (a=0; a<noltk; a++)
      //{
      //leptonNumbers.push_back(a);
      //}

      //bool fcel=false; // (finished cel)
      //while (!fcel)
      //{          
      Double_t weight=eventWeight;
      //rLeptons=leptons;
      rJets=jets;
      //rElectrons=electrons;
      //rMuons=muons;
      


      // lepton isolations (calorimetric isolation)
      // calorimetric isolation energy = sum of E_T of calorimeter cells in a cone of \Delta R < 0.3
      rLeptons.clear();
      rElectrons.clear();
      rMuons.clear();
      
      double ESum=0;
      double MSum1=0;
      double MSum2=0;
            
      for (a=0; a<leptons.size(); a++)
	{
	  for (b=0; b<rJets.size(); b++)
            {
              if (   abs(leptons[a].user_index()) == 11 
		     && deltaR(&(*(leptons.begin()+a)),&(*(rJets.begin()+b))) < .3 )
		{
		  //cout << "hello" << endl;
		  ESum+=rJets[b].Et();
		  //MSum1 = 100000.;
		  //MSum2 = 100000.;
		}
              else if (   abs(leptons[a].user_index()) == 13 
		     && deltaR(&(*(leptons.begin()+a)),&(*(rJets.begin()+b))) < .3 
		     && leptons[a].pt() > 25. ) 
		{
		  MSum1+=rJets[b].Et();
		  //ESum = 100000.;
		  //MSum2 = 100000.;
		}
              else if (   abs(leptons[a].user_index()) == 13 
		     && deltaR(&(*(leptons.begin()+a)),&(*(rJets.begin()+b))) < .3 
		     && leptons[a].pt() < 25. && leptons[a].pt() > 20. ) 
		{
		  MSum2+=rJets[b].Et();
		  //ESum = 100000.;
		  //MSum1 = 100000.;
		}
	      /*
	      else
		{
		  ESum = 100000.;
		  MSum1 = 100000.;
		  MSum2 = 100000.;
		}
	      */
            }//end of loop for jets

	  // energy of electron is subtracted from ESum : add electron E_T \approx p_T
	  if (abs(leptons[a].user_index()) == 11 && ESum < 1.28 * (leptons[a].pt())  )
	    {
	      rElectrons.push_back(leptons[a]);
	      rLeptons.push_back(leptons[a]);
	    }

	  //cout << "MSum1 = " << MSum1 << endl;
	  //cout << "MSum2 = " << MSum2 << endl;
	  	  
	  else if (abs(leptons[a].user_index()) == 13 && (MSum1 < 0.3 * (leptons[a].pt()) || MSum2 < 0.18 * (leptons[a].pt()) )  )
	    {
	      rMuons.push_back(leptons[a]);
	      rLeptons.push_back(leptons[a]);
	    }

	  //cout << "lepton pt = " << leptons[a].pt() << endl;

	  //cout << "ESum = " << ESum << endl;
	  //cout << "MSum1 = " << MSum1 << endl;
	  //cout << "MSum2 = " << MSum2 << endl;

	  
	}
            
      //cout << "number of electrons = " << rElectrons.size() << endl;
      //cout << "number of muons = " << rMuons.size() << endl;
      
      
      
      
      bool esFail=false;
      // fails event selection,
	      
      if (rLeptons.size() < 2) esFail=true;

      if (neutrinos.size() < 2) esFail=true;
          
      if (rJets.size() >= 1) // for jet veto : jets not included here survive after the jet veto
	{
	  if (rElectrons.size() >=1 )
	    {
	      for (a=0; a < rElectrons.size(); a++)
		{
		  for (b=0; b<rJets.size();)
		    {
		      if (deltaR(&(*(rElectrons.begin()+a)),&(*(rJets.begin()+b))) < .3  
			  ) 
			// electrons must not be vetoed.
			rJets.erase(rJets.begin()+b);
		      else b++;
		    }
		}
	    }
	  
	      //else vetoj++;
	}

      // remove low pt jets
      for (a=0; a<rJets.size();)
	{
	  if (rJets[a].pt() < 25. || fabs(rJets[a].eta()) > 4.5 ) rJets.erase(rJets.begin()+a);
	  else a++;
	}

      if (rJets.size() > 0) esFail=true;  // jet veto


      // extra requirement on the 3rd isolated lepton
      if (rLeptons.size() >= 3 && rLeptons[2].pt() > 7. ) esFail=true; 
      
      
      vector<PseudoJet> totalLeptons = rLeptons;	     
      if (!esFail)
	{
	  // Leptons originating (from the W only?) decay are corrected for QED radiation 
	  // by adding the four-momenta of photons with \Delta R = 0.1 of the decay leptons
	  	  
	  if (rLeptons.size() == 2 && photons.size() > 0)
	    {		  
	      for (a=0; a<rLeptons.size(); a++)
		{
		  for (b=0; b<photons.size(); b++)
		    {
		      if (  deltaR(&(*(rLeptons.begin()+a)),&(*(photons.begin()+b))) < .1  
			    && ( a > 0 && deltaR(&(*(rLeptons.begin()+a-1)),&(*(photons.begin()+b))) >= .1 ) 
			    )  
			{
			  totalLeptons[a]+= photons[b]; 	      
			}
		    }
		}
	    }
	  

	  
	  // define cut variables
	  PseudoJet miss = neutrinos[0] + neutrinos[1];
	  double mll = (totalLeptons[0] + totalLeptons[1]).m();
	  double mz = 91.1876;
	  double angle = min(fabs(totalLeptons[0].phi() - miss.phi()) , fabs(totalLeptons[1].phi()-miss.phi()) );
	  double ETmissRel = sin(min(angle, M_PI/2.0)) * miss.Et();


	  double mtllmet = sqrt(pow(totalLeptons[0].Et() + totalLeptons[1].Et() + miss.Et(), 2) - pow(fabs((totalLeptons[0]+totalLeptons[1]+miss).pt()), 2));


	  /*            
	  // set event category
	  unsigned int c;
	  if (zeroj) c=0;
	  else if (onej) c=1;
	  else c=2;
	  */
	  
	  // Multiply C_WW and cuts to each event category
		  
	  if (rElectrons.size() == 2 && rMuons.size() == 0 && totalLeptons[0].pt() > 25. && totalLeptons[1].pt() > 20.) 
	    {
	      //weight*=.291;
	      
	      if(mll > 15. && ETmissRel > 45. && miss.pt() > 45. && fabs(mll - mz) > 15.)
		{
		  atlcut[0]+=weight;

		  fptl1 -> Fill(totalLeptons[0].pt() );
		  fptl2 -> Fill(totalLeptons[1].pt() );
		  fptll -> Fill((totalLeptons[0] + totalLeptons[1]).pt() );
		  fdpll -> Fill(fabs(totalLeptons[0].phi() - totalLeptons[1].phi() ) );
		  fmtllmet -> Fill( mtllmet  );
		  fptllmet -> Fill( (totalLeptons[0]+totalLeptons[1]+miss).pt() );
		  
		  fnu1 -> Fill( neutrinos[0].pt() );
		  fnu2 -> Fill( neutrinos[1].pt() );
		  fnu -> Fill( miss.pt() );
		}
	    }
	  
	  else if (rElectrons.size() == 0 && rMuons.size() == 2 && totalLeptons[0].pt() > 25. && totalLeptons[1].pt() > 20.)
	    {
	      //weight*=.471;

	      if(mll > 15. && ETmissRel > 45. && miss.pt() > 45. && fabs(mll - mz) > 15.)
		{
		  atlcut[1]+=weight;

		  fptl1 -> Fill(totalLeptons[0].pt() );
		  fptl2 -> Fill(totalLeptons[1].pt() );
		  fptll -> Fill((totalLeptons[0] + totalLeptons[1]).pt() );
		  fdpll -> Fill(fabs(totalLeptons[0].phi() - totalLeptons[1].phi() ) );
		  fmtllmet -> Fill( mtllmet  );
		  fptllmet -> Fill( (totalLeptons[0]+totalLeptons[1]+miss).pt() );

		  fnu1 -> Fill( neutrinos[0].pt() );
		  fnu2 -> Fill( neutrinos[1].pt() );
		  fnu -> Fill( miss.pt() );

		  fmll -> Fill( mll );
		  fetmll -> Fill(ETmissRel);
		  fdr -> Fill( deltaR(&(*(totalLeptons.begin())),&(*(totalLeptons.begin()+1))) );


		}
	    }
	  
	  else if (rElectrons.size() == 1 && rMuons.size() == 1 && totalLeptons[0].pt() > 25. && totalLeptons[1].pt() > 20. ) 
	    {
	      
	      //weight*=.511; 
 
	      if(mll > 10. && ETmissRel > 15. && miss.pt() > 20.) 
		{
		  atlcut[2] += weight;

		  fptl1 -> Fill(totalLeptons[0].pt() );
		  fptl2 -> Fill(totalLeptons[1].pt() );
		  fptll -> Fill( (totalLeptons[0]+totalLeptons[1]).pt() );
		  fptllem -> Fill( (totalLeptons[0]+totalLeptons[1]).pt() );
		  fdpll -> Fill(fabs(totalLeptons[0].phi() - totalLeptons[1].phi() ) );
		  fmtllmet -> Fill( mtllmet  );
		  fmtllmetem -> Fill( mtllmet  );
		  fptllmet -> Fill( (totalLeptons[0]+totalLeptons[1]+miss).pt() );

		  fnu1 -> Fill( neutrinos[0].pt() );
		  fnu2 -> Fill( neutrinos[1].pt() );
		  fnu -> Fill( miss.pt() );

		  fmll -> Fill( mll );
		  fetmll -> Fill(ETmissRel);
		  fdr -> Fill( deltaR(&(*(totalLeptons.begin())),&(*(totalLeptons.begin()+1))) );


		}
	      
	    }
	  
	} // end analysis for selected events
            
	      
  
    }//end of limit
  
  


  ofstream myfile;    
  myfile.open("ppwwxi.dat");
  for(a=0; a<3; a++)
    {
      atlcut[a]/=eventWeightSum;
      myfile << atlcut[a] << "   ";
    }
  myfile.close();
  out.Write();

	  /*
    for (a=0; a<3; a++)
    {
      cout << "test = " << atlcut[a] << endl;
      //atlcut[a]/=eventWeightSum;
      //cout << "check = " << atlcut[a]/smeff[a] << endl;
    }
	  */
  
  
  

  return;
  
}// end of atlasppww



