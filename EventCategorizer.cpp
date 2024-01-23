/**
 *  @copyright Copyright 2019 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file EventCategorizer.cpp
 */

#include <JPetAnalysisTools/JPetAnalysisTools.h>
#include <JPetWriter/JPetWriter.h>
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetStatistics/JPetStatistics.h>
#include "EventCategorizer.h"
#include "../ModularDetectorAnalysis/EventCategorizerTools.h"
#include "../ModularDetectorAnalysis/HitFinderTools.h"
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <TMath.h>
#include <fstream>
#include <string>
#include <TVector3.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TCutG.h>
#include <cmath>

//#include "reconstructor.h"
#include "Hits/JPetBaseHit/JPetBaseHit.h"
#include "Hits/JPetPhysRecoHit/JPetPhysRecoHit.h"
#include "Hits/JPetMCRecoHit/JPetMCRecoHit.h"
#include "JPetRawMCHit/JPetRawMCHit.h"


using namespace jpet_options_tools;
using namespace std;

EventCategorizer::EventCategorizer(const char *name) : JPetUserTask(name) {}

bool EventCategorizer::init() {
  INFO("Event analysis started.");

 
 if (isOptionSet(fParams.getOptions(), kMC)) {
    fIsMC = getOptionAsBool(fParams.getOptions(), kMC);
 } 
     
  fOutputEvents = new JPetTimeWindow("JPetEvent");

  fReconstructor = Reconstructor::GetInstance();
  double scin_length = getParamBank().getScins().begin()->second->getLength();
  scin_length *= 0.1; // in cm
  fReconstructor->setBarrelLength( scin_length );
  fReconstructor->setChamberRadius( 10.0 );

  std::map<EvtType, std::string> EventTypeNames;
  EvtType event_type;

  initialiseHistograms();
  
  return true;
}


bool EventCategorizer::exec() {
  
  JPetTimeWindowMC *time_window_mc = nullptr;
  
  if (time_window_mc = dynamic_cast<JPetTimeWindowMC *const>(fEvent)) {
    fIsMC = true;
    EventTypeNames = EventTypeNames_MC;
  }
  else {
    fIsMC = false;
    EventTypeNames = EventTypeNames_data;
  }

   if (auto timeWindow = dynamic_cast<const JPetTimeWindow *const>(fEvent)) {
     
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      
      const auto &event = dynamic_cast<const JPetEvent &>(timeWindow->operator[](i));
     
      std::vector< const JPetPhysRecoHit*> hits, hits_tot;
      //std::vector< const JPetPhysRecoHit*> total_hits = event.getHits();
          
      getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(UNKNOWN)).c_str(), event.getHits().size());

      if ( event.getHits().size() == 1 ) {
      
	const JPetPhysRecoHit* Hit =  dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(0)); 
	if (fIsMC == false){

	  double tot_hits = Hit->getEnergy()/10000;
	
	  //getStatistics().fillHistogram("tot_data", tot_hits );
	  
	  if (tot_hits > 5.5){
	    pmtEmmTime = Hit->getTime()/1000 - Hit->getPos().Mag()/29.9792458;
	    prompt_time.push_back(pmtEmmTime);
	    
	  }
	}
      }
      
      EvtType event_type = UNKNOWN;

      getStatistics().fillHistogram( "hitMult_all", event.getHits().size());

      if(event.getHits().size() >= 3 ) {

	getStatistics().fillHistogram( "hitMult_3greater", event.getHits().size());

	hits_tot = tot( &event);
	getStatistics().fillHistogram("hitMult_tot", hits_tot.size() );
	
	hits = scatTest (hits_tot );
	getStatistics().fillHistogram( "hitMult_sTest", hits.size());
            	
	if( hits.size() == 3 ) {
	      
	  event_type = UNKNOWN;

	  getStatistics().fillHistogram( "HitMult", event_type);
	 	  	  
	  ReconHitsCalculation( hits, event_type = UNKNOWN, *time_window_mc);

	  //threeHits.clear();
	  hits.clear();
	  hits_tot.clear();

	}
		
      }

      lifeTime( prompt_time, threeHitTime);
      
    }
      
    prompt_time.clear();
    threeHitTime.clear();
    
   }


   else {
     return false;
   }

   return true;
}

	
	

void EventCategorizer::ReconHitsCalculation( vector<const JPetPhysRecoHit*> hits,  EvtType& event_type, const JPetTimeWindowMC& time_window_mc) {
  
  // Reconstructed hits -- Annihilation Point ( Trilateration Method ) (Alek's code)
  const JPetPhysRecoHit* hit1 = hits.at(0);
  const JPetPhysRecoHit* hit2 = hits.at(1);
  const JPetPhysRecoHit* hit3 = hits.at(2);
  
  fReconstructor->setGamma(0, hits.at(0));
  fReconstructor->setGamma(1, hits.at(1));
  fReconstructor->setGamma(2, hits.at(2));
  
  int error = 0;
  TVector3 sol[2];
  double t[2];
  error = fReconstructor->getSolution(sol[1], t[1], 1);
  getStatistics().fillHistogram("error_tri" , error);

  getStatistics().fillHistogram(("error_tri_"+ EventTypeNames.at(event_type)).c_str() , error);

  TVector3 decayPoint = sol[1];
  
  //TVector3 decayPoint;  
  //decayPoint = fReconstructor->RecoPos( hit1, hit2, hit3);
  
  std::vector<TVector3> momenta(3);	
  std::vector<pair<double, int>> k;
  std::vector<double> theta(3);
  std::vector<double> energies(3);
  TVector3 spin, normaltoPlane;
  double OperatorStudy[2];
  
  momenta.at(0) = hit1->getPos() - decayPoint;
  momenta.at(1) = hit2->getPos() - decayPoint;
  momenta.at(2) = hit3->getPos() - decayPoint;
  double radius = decayPoint.Perp();

  double transRad = pow( (pow(decayPoint.X(),2) + pow(decayPoint.Y(),2)),0.5);

  //p1_recs = hit1.getPos().Unit();
  //p1_uncer = p1_recs.Angle(p1_gen);
 
  std::vector<double> time;
  double vel =  29.9792458;
  time.push_back(momenta[0].Mag()/vel);
  time.push_back(momenta[1].Mag()/vel);
  time.push_back(momenta[2].Mag()/vel);
  meanTime_3Hits = (time[0]+time[1]+time[2])/3;
  
  Int_t combs[][2] = { {0,1}, {1,0}, {2,1}, {1,2}, {2,0}, {0,2} };
  TRandom3 rndm;
  int wc = rndm.Integer(3);

  theta[0] = momenta[0].Angle(momenta[1]);
  theta[1] = momenta[1].Angle(momenta[2]);
  theta[2] = momenta[0].Angle(momenta[2]);
  
  std::vector<double> theta_rndm;
  theta_rndm.push_back(theta[0]);
  theta_rndm.push_back(theta[1]);
  theta_rndm.push_back(theta[2]);

  std::random_shuffle (theta_rndm.begin(), theta_rndm.end());
  
  double mass_e = 510.99;
  energies[0] = -2*mass_e*(-cos(theta[2])+cos(theta[0])*cos(theta[1]))/((-1+cos(theta[0]))*(1+cos(theta[0])-cos(theta[2])-cos(theta[1])));
  energies[1] = -2*mass_e*(cos(theta[0])*cos(theta[2])-cos(theta[1]))/((-1+cos(theta[0]))*(1+cos(theta[0])-cos(theta[2])-cos(theta[1])));
  energies[2] = 2*mass_e*(1+cos(theta[0]))/(1+cos(theta[0])-cos(theta[2])-cos(theta[1]));
 
  std::vector<double> eng;       // deposited energy by gammas
  eng.push_back(hit1->getEnergy());
  eng.push_back(hit2->getEnergy());
  eng.push_back(hit3->getEnergy());
    
  k.push_back( {momenta[0].Mag(),0});
  k.push_back( {momenta[1].Mag(),1});
  k.push_back( {momenta[2].Mag(),2});
  std::sort( k.begin(), k.end());

  // using x and y  coordinates of hits
  std::vector<double> thetas;
  thetas.push_back(calculateAngle (hit1, hit2));
  thetas.push_back(calculateAngle (hit2, hit3));
  thetas.push_back(calculateAngle (hit3, hit1));
  std::sort (thetas.begin(), thetas.end());

  // using x y and z pos of hits
  std::vector<double> angle_3D;
  angle_3D.push_back(calculateAngle_3D (hit1, hit2));
  angle_3D.push_back(calculateAngle_3D (hit2, hit3));
  angle_3D.push_back(calculateAngle_3D (hit3, hit1));
 
  double sph_jacobian = pow(decayPoint.Mag(),2);

  std::random_shuffle (energies.begin(), energies.end());
  std::sort(theta.begin(), theta.end()); // from photons' momenta
  
  getStatistics().fillHistogram( "HitMult_beforeTri", event_type);
  getStatistics().fillHistogram(("scattTest_beforeTri_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( hits, event_type) );
  getStatistics().fillHistogram(("angle_dLOR_beforeTri_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( hits, event_type));
  getStatistics().fillHistogram(("E1E2_beforeTri_"+ EventTypeNames.at(event_type)).c_str(), energies.at(0), energies.at(1));
  getStatistics().fillHistogram(("theta12_beforeTri_"+ EventTypeNames.at(event_type)).c_str(), theta_rndm.at(0)*TMath::RadToDeg(), theta_rndm.at(1)*TMath::RadToDeg());
  getStatistics().fillHistogram(("hitDis_tDiff_min_beforeTri_"+ EventTypeNames.at(event_type)).c_str(), hitDistance_tDiff(hits).first, hitDistance_tDiff(hits).second );
  getStatistics().fillHistogram(("hitsDis_diff_beforeTri_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(hits).at(0), hitsDis_diff(hits).at(2) );
  getStatistics().fillHistogram(("sumDiff_angles_beforeTri_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());
  getStatistics().fillHistogram(("sumDiff_azmTheta_beforeTri_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));

  if ( error == 0){

    getStatistics().fillHistogram( "HitMult_afterTri", event_type);
    getStatistics().fillHistogram(("sumDiff_angles_afterTri_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());
    getStatistics().fillHistogram(("sumDiff_azmTheta_afterTri_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));
    getStatistics().fillHistogram(("angle_dLOR_afterTri_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( hits, event_type));      
    getStatistics().fillHistogram(("theta12_afterTri_"+ EventTypeNames.at(event_type)).c_str(), theta_rndm.at(0)*TMath::RadToDeg(), theta_rndm.at(1)*TMath::RadToDeg());    
    getStatistics().fillHistogram(("E1E2_afterTri_"+ EventTypeNames.at(event_type)).c_str(), energies.at(0), energies.at(1));
    getStatistics().fillHistogram(("hitDis_tDiff_min_afterTri_"+ EventTypeNames.at(event_type)).c_str(), hitDistance_tDiff(hits).first, hitDistance_tDiff(hits).second );
    getStatistics().fillHistogram(("hitsDis_diff_afterTri_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(hits).at(0), hitsDis_diff(hits).at(2) );
	
    if (energies[0] > 0 && energies[1] > 0 && energies[2] > 0){
	  
      getStatistics().fillHistogram( "HitMult_afterEng", event_type);
      //getStatistics().fillHistogram(("engDep_g1_"+ EventTypeNames.at(event_type)).c_str(), hit1.getEnergy()); 
         
      getStatistics().fillHistogram(("E1E2_afterEng_"+ EventTypeNames.at(event_type)).c_str(), energies.at(0), energies.at(1));
      getStatistics().fillHistogram(("theta12_afterEng_"+ EventTypeNames.at(event_type)).c_str(), theta_rndm.at(0)*TMath::RadToDeg(), theta_rndm.at(1)*TMath::RadToDeg());
      getStatistics().fillHistogram(("scattTest_afterEng_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( hits, event_type) );	
      getStatistics().fillHistogram(("sumDiff_angles_afterEng_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());
      getStatistics().fillHistogram(("sumDiff_azmTheta_afterEng_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));
      getStatistics().fillHistogram(("hitDis_tDiff_min_afterEng_"+ EventTypeNames.at(event_type)).c_str(), hitDistance_tDiff(hits).first, hitDistance_tDiff(hits).second );
      getStatistics().fillHistogram(("hitsDis_diff_afterEng_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(hits).at(0), hitsDis_diff(hits).at(2) );
      getStatistics().fillHistogram(("angle_dLOR_afterEng_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( hits, event_type));

      getStatistics().fillHistogram(("hit_xy_"+ EventTypeNames.at(event_type)).c_str(), hit1->getPos().X(), hit1->getPos().Y() );                                         
      getStatistics().fillHistogram(("hit_xz_"+ EventTypeNames.at(event_type)).c_str(), hit1->getPos().X(), hit1->getPos().Z() );                                          
      getStatistics().fillHistogram(("hit_x_"+ EventTypeNames.at(event_type)).c_str(), hit1->getPos().X() );                                                               
      getStatistics().fillHistogram(("hit_y_"+ EventTypeNames.at(event_type)).c_str(), hit1->getPos().Y() );                                                               
      getStatistics().fillHistogram(("hit_z_"+ EventTypeNames.at(event_type)).c_str(), hit1->getPos().Z() );                                                               
      getStatistics().fillHistogram(("hit_t_"+ EventTypeNames.at(event_type)).c_str(), hit1->getTime() );  

      if ( hit1->getEnergy()/1000000 < 4.5 && hit2->getEnergy()/1000000 < 4.5 && hit3->getEnergy()/1000000 < 4.5 ){

	getStatistics().fillHistogram( "HitMult_afterEngDep", event_type);
	getStatistics().fillHistogram(("scattTest_afterEngDep_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( hits, event_type) );
	getStatistics().fillHistogram(("angle_dLOR_afterEngDep_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( hits, event_type));
	getStatistics().fillHistogram(("E1E2_afterEngDep_"+ EventTypeNames.at(event_type)).c_str(), energies.at(0), energies.at(1));
	getStatistics().fillHistogram(("theta12_afterEngDep_"+ EventTypeNames.at(event_type)).c_str(), theta_rndm.at(0)*TMath::RadToDeg(), theta_rndm.at(1)*TMath::RadToDeg());
	getStatistics().fillHistogram(("sumDiff_angles_afterEngDep_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());
	getStatistics().fillHistogram(("sumDiff_azmTheta_afterEngDep_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));
	getStatistics().fillHistogram(("hitDis_tDiff_min_afterEngDep_"+ EventTypeNames.at(event_type)).c_str(), hitDistance_tDiff(hits).first, hitDistance_tDiff(hits).second );
	getStatistics().fillHistogram(("hitsDis_diff_afterEngDep_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(hits).at(0), hitsDis_diff(hits).at(2) );

	
   	if ( scatterTest( hits, event_type ) > 10 ){
	
	  getStatistics().fillHistogram( "HitMult_afterSTest", event_type);    
	  getStatistics().fillHistogram(("sumDiff_angles_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());
	  getStatistics().fillHistogram(("sumDiff_azmTheta_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));
	  getStatistics().fillHistogram(("angle_dLOR_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( hits, event_type));	         
	  getStatistics().fillHistogram(("theta12_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), theta_rndm.at(0)*TMath::RadToDeg(), theta_rndm.at(1)*TMath::RadToDeg());      
	  getStatistics().fillHistogram(("E1E2_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), energies.at(0), energies.at(1));
	  getStatistics().fillHistogram(("hitDis_tDiff_min_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), hitDistance_tDiff(hits).first, hitDistance_tDiff(hits).second );
	  getStatistics().fillHistogram(("hitsDis_diff_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(hits).at(0), hitsDis_diff(hits).at(2) );

	  
	  if ( (thetas[0]+ thetas[1]) > 181.5){

	    getStatistics().fillHistogram( "HitMult_after2DAngle", event_type);  
	    getStatistics().fillHistogram(("sumDiff_angles_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());
	    getStatistics().fillHistogram(("sumDiff_azmTheta_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));
	    getStatistics().fillHistogram(("angle_dLOR_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( hits, event_type));	         
	    getStatistics().fillHistogram(("theta12_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), theta_rndm.at(0)*TMath::RadToDeg(), theta_rndm.at(1)*TMath::RadToDeg());      
	    getStatistics().fillHistogram(("E1E2_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), energies.at(0), energies.at(1));
	    getStatistics().fillHistogram(("hitDis_tDiff_min_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), hitDistance_tDiff(hits).first, hitDistance_tDiff(hits).second );
	    getStatistics().fillHistogram(("hitsDis_diff_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(hits).at(0), hitsDis_diff(hits).at(2) );
	    getStatistics().fillHistogram(("scattTest_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( hits, event_type) );
	
      
	    
	    double angle_sum = thetas[0]+ thetas[1];
	    
	      getStatistics().fillHistogram(("distance_hits_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(hits).at(0) );
	      getStatistics().fillHistogram(("distance_hits_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(hits).at(1) );
	      getStatistics().fillHistogram(("distance_hits_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(hits).at(2) );
	    
	 	 
	      getStatistics().fillHistogram(("hitsDis_diff_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(hits).at(0), hitsDis_diff(hits).at(2) );
	    

	    
	      getStatistics().fillHistogram(("AnnhPointZComp_"+ EventTypeNames.at(event_type)).c_str(), decayPoint.Z());
	      getStatistics().fillHistogram(("AnnhPointXZComp_"+ EventTypeNames.at(event_type)).c_str(), decayPoint.X(), decayPoint.Z());
	      getStatistics().fillHistogram(("AnnhPointXYComp_"+ EventTypeNames.at(event_type)).c_str(), decayPoint.X(), decayPoint.Y());

	      //threeHitTime.push_back({decayTime_3hits (hits, event_type), event_type});

	      std::vector<int> z_position = {1, 11, 21, 22};
	      for (int n = 0; n < z_position.size(); n++){
		if ( sol[1].Z() < z_position.at(n) && sol[1].Z() > -z_position.at(n)){

		  getStatistics().fillHistogram( Form(("AnnhPointXY_Zto%d_"+ EventTypeNames.at(event_type)).c_str(), z_position.at(n)), sol[1].X(), sol[1].Y());
		       
		}
	      }
		 
	      double sph_jacobian = pow(decayPoint.Mag(),2);

	      normaltoPlane = momenta.at( k.at(2).second).Cross( momenta.at( k.at(1).second));
	      spin = decayPoint.Unit();

	      OperatorStudy[0]  = spin.Dot(momenta.at(k.at(2).second).Unit());
	      OperatorStudy[1]  = spin.Dot(normaltoPlane.Unit());

	      getStatistics().fillHistogram(("S_k1_"+ EventTypeNames.at(event_type)).c_str(), OperatorStudy[0] );
	      getStatistics().fillHistogram(("S_k1_k2_"+ EventTypeNames.at(event_type)).c_str(), OperatorStudy[1]);
	    }

	  }
	}
      }
     }		    
}



double EventCategorizer::scatterTest( vector<const JPetPhysRecoHit*> hits, EvtType event_type ) {

  double tDiff;
  double  c = 29.979246;   // cm/ns
  std::vector<double> deltaScatt;

  for(int i=0; i < hits.size(); ++i){
      for(int j=i+1; j < hits.size(); ++j){

	if (hits.at(i)->getTime() < hits.at(j)->getTime()) {
	  tDiff = (hits.at(j)->getTime() - hits.at(i)->getTime())/1000.0;
	} else {
	  tDiff = (hits.at(i)->getTime() - hits.at(j)->getTime())/1000.0;
	}
	
	float distance = sqrt(pow(hits.at(j)->getPosX()- hits.at(i)->getPosX(), 2) + pow(hits.at(j)->getPosY()- hits.at(i)->getPosY(), 2) + pow(hits.at(j)->getPosZ()- hits.at(i)->getPosZ(), 2));

	deltaScatt.push_back(  fabs(distance - tDiff*c ));
        
      }
  }
  
  std::sort(deltaScatt.begin(), deltaScatt.end());
  return deltaScatt.at(0);
  
}

	
std::vector<const JPetPhysRecoHit*> EventCategorizer::tot( const JPetEvent *event ){

  //const JPetEvent *temp;
  //temp = dynamic_cast<const JPetEvent *>(&event);
  //const JPetBaseHit* tmp;
  //tmp = dynamic_cast<const JPetBaseHit *>(event.getHits());

  //auto singleHit = dynamic_cast<const JPetBaseHit*>(event->getHits().at(j));
  std::vector<const JPetPhysRecoHit*> hits_tot;
  
  for(int i=0; i < event->getHits().size(); ++i){

    const JPetPhysRecoHit* singleHit = dynamic_cast<const JPetPhysRecoHit*>(event->getHits().at(i));

    //singleHit = event.getHits().at(i);
    
    //auto reco_hit = dynamic_cast<const JPetMCRecoHit*>(event->getHits().at(i));
    
    getStatistics().fillHistogram("tot_allHits", singleHit->getEnergy()/1000000 );
    
    if ( singleHit->getEnergy() < 4500000 ){
      
      hits_tot.push_back( singleHit );
      //hits_tot.push_back( i);

      
           
    }
  }
  
  return hits_tot;
  hits_tot.clear();
}

/*
std::vector<int> OPSAnalysis::scatTest(std::vector<const JPetPhysRecoHit*> hits){

  //double tDiff;                                                                                                                                                                                                
  double  c = 29.979246;   // cm/ns                                                                                                                                                                              
  std::vector<double> deltaScat;
  std::vector<int> hits_accepted, hits_sTest;
  std::vector<int> indices_rejected, indices_accepted, index_ref;

  for(auto i = hit_after_tot.begin(); i != hit_after_tot.end(); ++i)
  //for(auto i = 0; i < hit_after_tot.size(); ++i)                                                                                                                                                               
  {
    for(auto j=std::next(i); j!= hit_after_tot.end(); ++j)
      {
        if (time[*i] < time[*j]) {
          tDiff = (time[*j] - time[*i] );
        } else {
          tDiff = (time[*i] - time[*j] );
        }

        //cout<<time[*i]<<" "<<time[*j]<<endl;                                                                                                                                                                   
        distance = sqrt(pow(position[*j].X()-position[*i].X(), 2) + pow(position[*j].Y()- position[*i].Y(), 2) + pow(position[*j].Z()- position[*i].Z(), 2));

        deltaScat.push_back( fabs(distance - tDiff*c));

        h[6]->Fill(fabs(distance - tDiff*c ));
	if ( fabs(distance - tDiff*c ) < kScatterTestCut ){

          if (time[*i] < time[*j]) {

            hits_accepted.push_back(*i);
            indices_rejected.push_back(*j);
          }
          else {
            hits_accepted.push_back(*j);
            indices_rejected.push_back(*i);
          }
        }
        else {
          hits_accepted.push_back( *i );
          hits_accepted.push_back( *j );
        }
      }

    index_ref.push_back(*i);
  }

  std::sort(indices_rejected.begin(), indices_rejected.end() );
  for ( index_t a = 0; a < indices_rejected.size(); a++){

    int index_rej = indices_rejected.at(a);

    auto new_index = std::remove( index_ref.begin(), index_ref.end(), index_rej);

    index_ref.erase( new_index, index_ref.end());
  }

  for (index_t b = 0; b < index_ref.size(); b++){

    hits_sTest.push_back ( index_ref.at(b) );

    //cout<<"Out: "<<hits_sTest.at(b)<<endl;                                                                                                                                                                   

  }

  return hits_sTest;

  hits_sTest.clear();
  hits_accepted.clear();
  indices_accepted.clear();
  indices_rejected.clear();
  index_ref.clear();

  //cout<<endl;                                                                                                                                                                                                  
}
*/


std::vector<const JPetPhysRecoHit*> EventCategorizer::scatTest( std::vector<const JPetPhysRecoHit*> hits ) {

  double tDiff;
  double  c = 29.979246;   // cm/ns
  std::vector<double> deltaScatt;
  std::vector<const JPetPhysRecoHit*> hits_sTest, hits_accepted;
  std::vector<int> indices_rejected, indices_accepted, index_ref;

  getStatistics().fillHistogram("hitMult_beforeSTest", hits.size() );

  for(int i=0; i < hits.size(); ++i){
      for(int j=i+1; j < hits.size(); ++j){

	  if (hits.at(i)->getTime() < hits.at(j)->getTime()) {
	    tDiff = (hits.at(j)->getTime() - hits.at(i)->getTime())/1000.0;
	  } else {
	    tDiff = (hits.at(i)->getTime() - hits.at(j)->getTime())/1000.0;
	  }
	
	  float distance = sqrt(pow(hits.at(j)->getPosX()- hits.at(i)->getPosX(), 2) + pow(hits.at(j)->getPosY()- hits.at(i)->getPosY(), 2) + pow(hits.at(j)->getPosZ()- hits.at(i)->getPosZ(), 2));

	  deltaScatt.push_back(  fabs(distance - tDiff*c ));
	  getStatistics().fillHistogram("scattTest_moreHits", fabs(distance - tDiff*c ));
                                                                                          
	  getStatistics().fillHistogram("tDiff_scattTest", tDiff );

	  // 23
	  if ( fabs(distance - tDiff*c ) < 18 ){
	    
	    if (hits.at(i)->getTime() < hits.at(j)->getTime()) {
	      
	      hits_accepted.push_back( hits.at(i) );
	      indices_rejected.push_back(j);
	    }
	    else {
	      hits_accepted.push_back( hits.at(j) );
	      indices_rejected.push_back(i);
	    }
	  }
	  else {
	    hits_accepted.push_back( hits.at(i));
	    hits_accepted.push_back( hits.at(j));

	    //indices_accepted.push_back(i);
	    //indices_accepted.push_back(j);
	  }
	  //getStatistics().fillHistogram("tDiff_afterSTest", tDiff );	  
      }
      
      index_ref.push_back(i);
  }

  std::sort(indices_rejected.begin(), indices_rejected.end() );

  for (int a = 0; a < indices_rejected.size(); a++){
    
    int index_rej = indices_rejected.at(a);
      
    auto new_index = std::remove( index_ref.begin(), index_ref.end(), index_rej);

    index_ref.erase( new_index, index_ref.end());
   
  }
       
  for (int b = 0; b < index_ref.size(); b++){
          
    hits_sTest.push_back ( hits.at( index_ref.at(b) ));
    
  }

  getStatistics().fillHistogram("hitMult_afterSTest", hits_sTest.size() );
 
  return hits_sTest;
  
  hits_sTest.clear();
  hits_accepted.clear();
  indices_accepted.clear();
  indices_rejected.clear();
  index_ref.clear();
}




std::pair<double, double> EventCategorizer::hitDistance_tDiff ( vector<const JPetPhysRecoHit*> hits) {

  std::pair<double,double> scatterInfo;
  double tDiff;
  std::vector<double> distance, hitTDiff;
  std::vector< pair <double,int>> deltaScatt;
  double  c = 29.979246;   // cm/ns
  
  for(int i=0; i < hits.size(); ++i){
      for(int j=i+1; j < hits.size(); ++j){

	if ( hits.at(i)->getTime() < hits.at(j)->getTime()) {
	  tDiff = (hits.at(j)->getTime() - hits.at(i)->getTime())/1000.0;
	} else {
	  tDiff = ( hits.at(i)->getTime() - hits.at(j)->getTime())/1000.0;
	}
	
        distance.push_back( sqrt( pow( hits.at(j)->getPosX() - hits.at(i)->getPosX(),2) + pow( hits.at(j)->getPosY() - hits.at(i)->getPosY(),2) + pow( hits.at(j)->getPosZ() - hits.at(i)->getPosZ(),2) ));
	
	hitTDiff.push_back(tDiff);	
	
      }
  }

  for (int k =0; k<3; k++){

    deltaScatt.push_back( {fabs(distance.at(k) - hitTDiff.at(k)*c ),k} );    
  }

  std::sort ( deltaScatt.begin(), deltaScatt.end());

  scatterInfo = {distance.at(deltaScatt.at(0).second), hitTDiff.at(deltaScatt.at(0).second)};
  return scatterInfo;
  
}



std::vector<double> EventCategorizer::TOT( std::vector<double> eng){ // EvtType& event_type ){
  
  double param[4], sumTOT = 0;
  std::vector<double> tot;
  
  param[0] =  29.84;  // ns
  param[1] =  2.446; // ns/keV
  param[2] =  310; // keV
  param[3] =  2.41; // ns/keV^2
  
  for ( int i = 0; i < 3; i++){
    tot.push_back( param[0] +  (param[1] - param[0])/ (1 + pow(eng.at(i)/param[2], param[3]))); 							      
  }
  
    return tot;
    tot.clear(); 
}

double EventCategorizer:: calculateAngle(const JPetPhysRecoHit* hit1,  const JPetPhysRecoHit* hit2 ){
  
  double scalerProd = hit1->getPosX()*hit2->getPosX() + hit1->getPosY()*hit2->getPosY() ;
  double magnitude = sqrt((pow(hit1->getPosX(),2) + pow(hit1->getPosY(),2) )*(pow(hit2->getPosX(),2) + pow(hit2->getPosY(),2)));

  return acos(scalerProd/magnitude)*180/3.14159265;
}


double EventCategorizer:: calculatedLOR( vector< const JPetPhysRecoHit*> hits, EvtType event_type ){

  std::vector<Double_t> d_LOR;
  double dLOR_min;
  for(int i=0; i < hits.size(); ++i){
    for(int j=i+1; j < hits.size(); ++j){
	TVector3 vtx2g = calculateAnnihilationPoint( hits.at(i), hits.at(j));
        d_LOR.push_back( vtx2g.Mag());
      }
  }
  
  std::sort (d_LOR.begin(), d_LOR.end());
  dLOR_min = d_LOR[0];
  return dLOR_min;
    
}



std::vector<double> EventCategorizer:: calcLOR_recsAnnhPt( vector<const JPetPhysRecoHit*> hits, TVector3 decayPt ){

  std::vector<Double_t> distance;
  double distance_min;
  for(int i=0; i < hits.size(); ++i){
    for(int j=i+1; j < hits.size(); ++j){
	TVector3 vtx2g = calculateAnnihilationPoint( hits.at(i), hits.at(j));
        distance.push_back ((vtx2g - decayPt).Mag());
        
      }
  }
  
  std::sort (distance.begin(), distance.end());
  //distance_min = distance.at(0);
  return distance;

}


std::vector<double> EventCategorizer::distance_hits (vector<const JPetPhysRecoHit*> hits) {

  std::vector<double> distance;
  
  for(int i=0; i < hits.size(); ++i){
    for(int j=i+1; j < hits.size(); ++j){
	
        distance.push_back( sqrt( pow( hits.at(j)->getPosX() - hits.at(i)->getPosX(),2) + pow( hits.at(j)->getPosY() - hits.at(i)->getPosY(),2) + pow( hits.at(j)->getPosZ() - hits.at(i)->getPosZ(),2) ));
	
      }
  }

  std::sort (distance.begin(), distance.end());
  
  return distance;
}


std::vector<double> EventCategorizer::hitsDis_diff ( vector<const JPetPhysRecoHit*> hits ) {

  std::vector<double> distance_diff;
  
  for(int i=0; i < hits.size(); ++i){
    for(int j=i+1; j < hits.size(); ++j){
 
	distance_diff.push_back( fabs (distance_hits(hits).at(i) - distance_hits(hits).at(j)) );
      }
  }

  std::sort (distance_diff.begin(), distance_diff.end());
  
  return distance_diff;
}

double EventCategorizer:: calculateAngle_3D(const JPetPhysRecoHit* hit1,  const JPetPhysRecoHit* hit2 ){
  
  double scalerProd = hit1->getPosX()*hit2->getPosY() + hit1->getPosY()*hit2->getPosY() ;
  double magnitude = sqrt((pow(hit1->getPosX(),2) + pow(hit1->getPosY(),2) )*(pow(hit2->getPosX(),2) + pow(hit2->getPosY(),2))*
			  (pow(hit1->getPosZ(),2) + pow(hit1->getPosZ(),2) ));

  return acos(scalerProd/magnitude)*180/3.14159265;
}

double EventCategorizer::decayTime_3hits (std::vector<const JPetPhysRecoHit*> hits, EvtType& event_type){

  double time;
  double vel = 29.9792458 ;
 
  fReconstructor->setGamma(0, hits.at(0));
  fReconstructor->setGamma(1, hits.at(1));
  fReconstructor->setGamma(2, hits.at(2));
  
  int error = 0;
  TVector3 sol[2];
  double t[2];
  error = fReconstructor->getSolution(sol[1], t[1], 1);
  time = t[1];

  return time;

}

/*
double EventCategorizer::lifeTime ( const JPetEvent *event, EvtType& event_type, std::vector<double> prompt_emTime, int size){
  
  double timeDiff;
  
  for(int j = 0; j < decayTime_3hits(event,event_type).size(); j++){
    for(int i = 0; i < size; i++){
      timeDiff = decayTime_3hits(event,event_type).at(j) - prompt_emTime.at(i);
      return timeDiff/1000;
      
    }
	  
  }
}
*/


TVector3 EventCategorizer::calculateAnnihilationPoint(const JPetPhysRecoHit* hit1, const JPetPhysRecoHit* hit2)
{
  TVector3 middleOfLOR = 0.5 * (hit1->getPos() + hit2->getPos());
  TVector3 versorOnLOR = (hit2->getPos() - hit1->getPos()).Unit();

  double tof = calculateTOF(hit1, hit2);
  double shift = 0.5 * tof * kLightVelocity_cm_ps;

  TVector3 annihilationPoint;
  annihilationPoint.SetX(middleOfLOR.X() + shift * versorOnLOR.X());
  annihilationPoint.SetY(middleOfLOR.Y() + shift * versorOnLOR.Y());
  annihilationPoint.SetZ(middleOfLOR.Z() + shift * versorOnLOR.Z());
  return annihilationPoint;
}

double EventCategorizer::calculateTOF(const JPetPhysRecoHit* hitA, const JPetPhysRecoHit* hitB)
{
  return calculateTOF(hitA->getTime(), hitB->getTime());
}

double EventCategorizer::calculateTOF(double time1, double time2) { return (time1 - time2); }



bool EventCategorizer::ellipse_cut_values( double angle_sum, double dlor_min){
  double X, Y;
  double x,y;
  X = angle_sum, Y = dlor_min;
  double center_x = 217, center_y = 17;
  double u = 13, v = 7;
  double angle = (M_PI*30)/180; 

  coordinates ellipse_boundary;
  double ellip_cut;

  x =  pow((((X - center_x) * cos(angle) + (Y - center_y) * sin(angle)) / u), 2);
  y = pow((((X - center_x) * sin(angle) - (Y - center_y) * cos(angle)) / v), 2);
  
  if ((x+y) <= 1 ){
    ellipse_boundary.first = X;
    ellipse_boundary.second = Y;
    return true;
  }
  else{
    return false;
  }
  
}


void EventCategorizer::lifeTime ( std::vector<double> promptTime, std::vector<std::pair<double, EvtType>> hitTime){
  
  double timeDiff;

  if (promptTime.size() != 0){
    getStatistics().fillHistogram( "decayTime_prompt_size", promptTime.size());
  
    for (auto& i: EventTypeNames){
   
      for (int j = 0; j < promptTime.size(); j++){
      
	getStatistics().fillHistogram( "decayTime_prompt", promptTime.at(j));
	    
	for (int k = 0; k < hitTime.size(); k++){
			   
	  if (hitTime.at(k).second  ==  i.first) {
	 
	    timeDiff = hitTime.at(k).first - promptTime.at(j);
	  

	    getStatistics().fillHistogram(("dt_"+i.second).c_str(), timeDiff/1000 );
	    getStatistics().fillHistogram(("dt_ns_"+i.second).c_str(), timeDiff );
	    getStatistics().fillHistogram(("decayTime_"+ i.second).c_str(), hitTime.at(k).first );
     
 
	    }
      
      }
    }
    }
    
  }
}


bool EventCategorizer::terminate() {
  INFO("Event analysis completed.");
  return true;
}
	
void EventCategorizer::initialiseHistograms()
{

  if (fIsMC == true){
  	EventTypeNames = EventTypeNames_MC;
    }
  else {
    EventTypeNames = EventTypeNames_data;
  }
   
  std::vector<std::pair<unsigned, std::string>> binLabels;
  std::vector<int> z_position = {1, 11, 21, 22};

  getStatistics().createHistogramWithAxes(new TH1D( "scattTest_moreHits", "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D( "tDiff_scattTest", " ", 11, -0.5, 10.5), "tDiff [ns]", "Counts");

  getStatistics().createHistogramWithAxes(new TH1D( "hitMult_all", "Number of hits", 25, -0.5, 24.5), "Number of hits", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D( "hitMult_3greater", "Number of hits", 25, -0.5, 24.5), "Number of hits", "Counts");
  
  getStatistics().createHistogramWithAxes(new TH1D( "hitMult_beforeSTest", "Number of hits", 25, -0.5, 24.5), "Number of hits", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D( "hitMult_afterSTest", "Number of hits", 25, -0.5, 24.5), "Number of hits", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D( "hitMult_sTest", "Number of hits", 25, -0.5, 24.5), "Number of hits", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D( "hitMult_tot", "Number of hits", 25, -0.5, 24.5), "Number of hits", "Counts");

  getStatistics().createHistogramWithAxes(new TH1D( "tot_allHits", "", 40, -0.5, 19.5), "Eng [keV]", "Counts");

  getStatistics().createHistogramWithAxes(new TH1D( "ops_events_beforeCut", "Number of hits", 3, 0.5, 3.5), "Number of hits", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D( "ops_events_afterCut", "Number of hits", 3, 0.5, 3.5), "Number of hits", "Counts");

  
  for (auto& i: EventTypeNames){
    
    getStatistics().createHistogramWithAxes(new TH1D( ("HitMult_"+ i.second).c_str(), "Number of Hits", 25, -0.5, 24.5), "Number of hits", "Counts");
    binLabels.push_back(std::make_pair(i.first, i.second));

    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_beforeTri_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_afterTri_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_afterEng_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_afterEngDep_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_after2DAngle_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    
    getStatistics().createHistogramWithAxes(new TH1D( ("AnnhPointZComp_"+ i.second).c_str(), " Annihilation Point (Z) ", 120, -60, 60), " AnnihilationPoint_Z [cm]", "Counts");
    
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXZComp_"+ i.second).c_str(), "Annihilation Point XZ ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Z [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXYComp_"+ i.second).c_str(), " Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
   
    getStatistics().createHistogramWithAxes(new TH2D( ("hit_xy_"+ i.second).c_str(), " Hit position", 120, -60, 60, 120, -60, 60), "  X [cm]", "Y [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hit_xz_"+ i.second).c_str(), " Hit position ", 120, -60, 60, 120, -60, 60), " X [cm]", "Z [cm]");
    getStatistics().createHistogramWithAxes(new TH1D( ("hit_x_"+ i.second).c_str(), " Hit position ", 120, -60, 60), "X [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("hit_y_"+ i.second).c_str(), " Hit position ", 120, -60, 60), "Y [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("hit_z_"+ i.second).c_str(), " Hit position ", 120, -60, 60), "Z [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("hit_t_"+ i.second).c_str(), " Hit position ", 210, -10.5, 199.5), "t [ns]", "Counts");

    getStatistics().createHistogramWithAxes(new TH2D( ("mchit_xy_"+ i.second).c_str(), " Hit position", 120, -60, 60, 120, -60, 60), "  X [cm]", "Y [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("mchit_xz_"+ i.second).c_str(), " Hit position ", 120, -60, 60, 120, -60, 60), " X [cm]", "Z [cm]");
    //getStatistics().createHistogramWithAxes(new TH1D( ("mchit_x_"+ i.second).c_str(), " Hit position ", 120, -60, 60), "X [cm]", "Counts");
    //getStatistics().createHistogramWithAxes(new TH1D( ("mchit_y_"+ i.second).c_str(), " Hit position ", 120, -60, 60), "Y [cm]", "Counts");
    //getStatistics().createHistogramWithAxes(new TH1D( ("mchit_z_"+ i.second).c_str(), " Hit position ", 120, -60, 60), "Z [cm]", "Counts");
    //getStatistics().createHistogramWithAxes(new TH1D( ("mchit_t_"+ i.second).c_str(), " Hit position ", 210, -10.5, 199.5), "t [ns]", "Counts");

   
    //getStatistics().createHistogramWithAxes(new TH1D( ("res_x_"+ i.second).c_str(), " ", 80, -20, 20), "X (rec-gen) [cm]", "Counts");
    //getStatistics().createHistogramWithAxes(new TH1D( ("res_y_"+ i.second).c_str(), " ", 80, -20, 20), "Y (rec-gen) [cm]", "Counts");
    //getStatistics().createHistogramWithAxes(new TH1D( ("res_z_"+ i.second).c_str(), " ", 80, -20, 20), "Z (rec-gen) [cm]", "Counts");
    //getStatistics().createHistogramWithAxes(new TH1D( ("res_t_"+ i.second).c_str(), " ", 210, -10.5, 199.5), "t (rec-gen) [ns]", "Counts");


   

    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_beforeTri_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_afterEngDep_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_afterSTest_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_afterTri_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_afterEng_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_after2DAngle_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");

    getStatistics().createHistogramWithAxes(new TH2D( ("hitDis_tDiff_min_beforeTri_"+ i.second).c_str(), "", 160, -10.5, 149.5, 210, -10.5, 199.5), "hit distance [cm]", "hit time diff [ns]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitDis_tDiff_min_afterEngDep_"+ i.second).c_str(), "", 160, -10.5, 149.5, 210, -10.5, 199.5), "hit distance [cm]", "hit time diff [ns]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitDis_tDiff_min_afterSTest_"+ i.second).c_str(), "", 160, -10.5, 149.5, 210, -10.5, 199.5), "hit distance [cm]", "hit time diff [ns]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitDis_tDiff_min_after2DAngle_"+ i.second).c_str(), "", 160, -10.5, 149.5, 210, -10.5, 199.5), "hit distance [cm]", "hit time diff [ns]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitDis_tDiff_min_afterTri_"+ i.second).c_str(), "", 160, -10.5, 149.5, 210, -10.5, 199.5), "hit distance [cm]", "hit time diff [ns]");
    //getStatistics().createHistogramWithAxes(new TH2D( ("hitDis_tDiff_min_afterEllipCut_"+ i.second).c_str(), "", 160, -10.5, 149.5, 210, -10.5, 199.5), "hit distance [cm]", "hit time diff [ns]");
    
    //getStatistics().createHistogramWithAxes( new TH1F(("annh_radius_"+ i.second).c_str(), " ", 120, -60, 60), " radius [cm]", "Counts");
    //getStatistics().createHistogram( new TH1F(("annh_radius_sph_jacobian_"+ i.second).c_str(), "Radius of Annh point (weighted);" " radius [cm]; Counts", 120, -60, 60));

    getStatistics().createHistogramWithAxes( new TH2D (("theta12_beforeTri_"+ i.second).c_str(), "Angle12_vs_Angle23", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{23} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("theta12_afterSTest_"+ i.second).c_str(), "Angle12_vs_Angle23 ", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{23} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("theta12_afterEngDep_"+ i.second).c_str(), "Angle12_vs_Angle23 ", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{23} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("theta12_after2DAngle_"+ i.second).c_str(), "Angle12_vs_Angle23 ", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{23} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("theta12_afterTri_"+ i.second).c_str(), "Angle12_vs_Angle23", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{23} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("theta12_afterEng_"+ i.second).c_str(), "Angle12_vs_Angle23", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{23} [deg]");
    // getStatistics().createHistogramWithAxes( new TH2D (("theta12_afterEllipCut_"+ i.second).c_str(), "Angle12_vs_Angle23 ", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{23} [deg]");

    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_angles_beforeTri_"+ i.second).c_str()," ", 300, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{2} - #theta_{1} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_angles_afterEngDep_"+ i.second).c_str()," ", 300, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{2} - #theta_{1} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_angles_afterSTest_"+ i.second).c_str()," ", 300, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{2} - #theta_{1} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_angles_after2DAngle_"+ i.second).c_str()," ", 300, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{2} - #theta_{1} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_angles_afterTri_"+ i.second).c_str()," ", 300, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{2} - #theta_{1} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_angles_afterEng_"+ i.second).c_str()," ", 300, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{2} - #theta_{1} [deg]");
    //getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_angles_afterEllipCut_"+ i.second).c_str()," ", 300, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{2} - #theta_{1} [deg]");

    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_beforeTri_"+ i.second).c_str()," ", 300, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_afterEng_"+ i.second).c_str()," ", 300, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_afterEngDep_"+ i.second).c_str()," ", 300, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_afterSTest_"+ i.second).c_str()," ", 300, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_after2DAngle_"+ i.second).c_str()," ", 300, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_afterTri_"+ i.second).c_str()," ", 300, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_afterEllipCut_"+ i.second).c_str()," ", 300, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");

    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_afterEngDep_"+ i.second).c_str()," ", 300, -0.5, 299.5, 300, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_beforeTri_"+ i.second).c_str()," ", 300, -0.5, 299.5, 300, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_afterSTest_"+ i.second).c_str()," ", 300, -0.5, 299.5, 300, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_after2DAngle_"+ i.second).c_str()," ", 300, -0.5, 299.5, 300, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_afterTri_"+ i.second).c_str()," ", 300, -0.5, 299.5, 300, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_afterEng_"+ i.second).c_str()," ", 300, -0.5, 299.5, 300, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    //getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_afterEllipCut_"+ i.second).c_str()," ", 300, -0.5, 299.5, 300, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    
    getStatistics().createHistogramWithAxes( new TH2D (("E1E2_beforeTri_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E1E2_afterEngDep_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E1E2_afterSTest_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E1E2_after2DAngle_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E1E2_afterTri_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E1E2_afterEng_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
    //getStatistics().createHistogramWithAxes( new TH2D (("E1E2_afterEllipCut_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
   
      
    getStatistics().createHistogramWithAxes( new TH1D(("S_k1_"+ i.second).c_str(), "S.k1 ", 100, -1, 1), "S.k1", "Counts");
    getStatistics().createHistogramWithAxes( new TH1D(("S_k1_k2_"+ i.second).c_str(), "S.(k1*k2) ", 100, -1, 1), "S.(k1*k2) ", "Counts");

    std::vector<int> z_position = {1, 11, 21, 22};
    for (int n = 0; n < z_position.size(); n++){
       getStatistics().createHistogramWithAxes(new TH2D( Form(("AnnhPointXY_Zto%d_"+ i.second).c_str(),z_position.at(n)), "Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");

     }
         
    getStatistics().createHistogramWithAxes( new TH1D(("error_tri_"+ i.second).c_str(), "error_tri", 100, -0.5, 9.5), "error ", "Counts");
    
    getStatistics().createHistogramWithAxes(new TH1D( ("dt_"+ i.second).c_str(), "", 100, -50.5, 49.5), "#delta t [micro s]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("dt_ns_"+ i.second).c_str(), "", 1000, -500.5, 499.5), "#delta t [ns]", "Counts");
    getStatistics().createHistogramWithAxes( new TH1D(("decayTime_"+ i.second).c_str(), " ",1001, -500.5, 500.5), " o-ps decay time", "Counts"); 
					     
  }
  
  getStatistics().createHistogramWithAxes(new TH2D( "HitMult_ScinID", "Number of Hits v/s scinID", 200, -0.5, 199.5,  11, -0.5, 10.5), "Scin_ID", "Number of Hits");
    
  getStatistics().createHistogramWithAxes(new TH1D( "HitMult", "Number of hits", 23, 0.5, 23.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes( new TH1D("error_tri", "error_tri", 100, -0.5, 9.5), "error ", "Counts");
 
  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_beforeTri", "Number of hits", 23, 0.5, 23.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_beforeTri", getStatistics().AxisLabel::kXaxis, binLabels);
					     
  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterEng", "Number of hits", 23, 0.5, 23.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_afterEng", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterEngDep", "Number of hits", 23, 0.5, 23.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_afterEngDep", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterSTest", "Number of hits", 23, 0.5, 23.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_afterSTest", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_after2DAngle", "Number of hits", 23, 0.5, 23.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_after2DAngle", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterTri", "Number of hits", 23, 0.5, 23.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_afterTri", getStatistics().AxisLabel::kXaxis, binLabels);

  // getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterEllipCut", "Number of hits", 23, 0.5, 23.5), "Event Type", "Counts");
  //getStatistics().setHistogramBinLabel("HitMult_afterEllipCut", getStatistics().AxisLabel::kXaxis, binLabels);

  
}

