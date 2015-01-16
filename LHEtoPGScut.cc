//
// compile: g++ -m32 `root-config --libs` -I`root-config --incdir` -o LHEtoPGScut LHEtoPGScut.cc
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include "TLorentzVector.h"

using namespace std;

int nParticlesPrint;
int NDUMMY;

double pi=3.141592;

enum Charge {PLUS=0,MINUS=1};
Charge charge(int i) { return (i>0) ? PLUS :  MINUS; };



class SpectrometerCut // for basic geometries, no rotation
{
public:
  SpectrometerCut(double xmin,double xmax, double ymax, double pmin, double pmax);
  bool pass(vector<TLorentzVector> momenta, vector<Charge> signs);
  static int debug() { DEBUG=true;}

private:
  double xmin_, xmax_, ymax_, pmin_, pmax_;
  static bool DEBUG;
};

class SpectrometerCutMainz // for asymmetric geometries
{
public:
  SpectrometerCutMainz(Charge qRef, int nrotMax);
  std::string setSpec(Charge q, double xmin,double xmax, double ymax, double pmin, double pmax);
  vector<double > rotatePass(vector<TLorentzVector> momenta, vector<Charge> signs) const;
  bool pass(vector<TLorentzVector> momenta, vector<Charge> signs, Charge qref) const;
  //obsolete:  void doRotate(bool dr) { if(dr==false) ncopy_=0;  };
  static int sgn(Charge q) { if(q==PLUS) return 1; else return -1; }; 
  static int debug(){ DEBUG=true;}

  int ncopy() const { return ncopy_;}
  double xmin(Charge q) const  { return xmin_[int(q)];}
  double xmax(Charge q) const  { return xmax_[int(q)];}
  double ymax(Charge q) const  { return ymax_[int(q)];}
  double pmin(Charge q) const  { return pmin_[int(q)];}
  double pmax(Charge q) const  { return pmax_[int(q)];}


  Charge qRef_;
  int ncopy_;
  double phi_max_;

  vector<double> xmin_, xmax_, ymax_, pmin_, pmax_;
  static bool DEBUG;
};

bool SpectrometerCut::DEBUG=false;
bool SpectrometerCutMainz::DEBUG=false;



void getP(istream& ins, TLorentzVector& part)
{
  double px,py,pz,pt,mass;
  ins >> px >> py >> pz >> pt >> mass;
  
  part=TLorentzVector(px, py, pz, pt);
  //  cout << "reading particle:  x " << part.Px() << " y " <<  part.Py() << " z " <<  part.Pz() <<  " E " <<  part.E() << " mass " << part.M() << endl;
}

void printP(ostream& outs, TLorentzVector& part, string typeString, int pid,double Ebeam, int ipart)
{

        outs << ipart << setw(7) << typeString
             << " " << setw(10) << setprecision(5) << part.Eta()
             << " " << setw(10) << setprecision(5) << part.Phi()
             << " " << setw(10) << setprecision(5) << part.Pt()
             << " " << setw(10) << setprecision(5) << part.M()
             << " " << setw(6) << setprecision(2) << (pid>0? -1 : 1) << ".00";
        outs << " " << setw(6) << setprecision(2) << Ebeam;
        for(int i=0; i<NDUMMY-1; ++i) {
          outs << setw(6) << "0.00";
        }
        outs << endl;

}

class Event
{
public:
  //  ~Event() { cout << " -- end event -- " << endl;}

  void addP(TLorentzVector& part, string typeString, int pid, double Ebeam);
  void print(ostream& outs, double rotAngle);

private:
  vector<TLorentzVector> particles_;
  vector<string> typeString_;
  vector<int> pid_;
  double Ebeam_;
  
};
  
void Event::print(ostream& outs, double rotAngle)
{
  for(int i=0; i<particles_.size(); ++i) {
    TLorentzVector rotParticle(particles_[i]);
    rotParticle.RotateZ(rotAngle);

    printP(outs, rotParticle, typeString_[i], pid_[i], Ebeam_, i+1);
  }

  //  cout << " ---- good event printed ---- " << endl;
}

void Event::addP(TLorentzVector& part, string typeString, int pid, double Ebeam)
{
  particles_.push_back(part);
  typeString_.push_back(typeString);
  pid_.push_back(pid);
  Ebeam_=Ebeam;
  
  /*if(true) {if(abs(pid)==11)
      cout << setw(4) << pid << ": p=" << setw(6) << part.P() << " theta=" << setw(6) << part.Theta() << " phi=" << setw(6) << part.Phi() << endl;
      }*/
}

bool readEvent(istream& ins, ostream& passouts, int& nev, int& nevpass, SpectrometerCutMainz& cutfunc)
{
  string nextline;
  string ignoreRestOfLine;

  do {
    getline(ins, nextline);
    //    cout << "NEXT : " << nextline << endl;
    if(ins.eof()) return false; // failed to find an event
  } while (nextline.find("<event>") == string::npos);// continue while test string isn't found.
  
  //  cout << "AN EVENT: " << endl;
  // cout << "0" << setw(5) << nevpass+1 << setw(7) <<  "9999" << endl; 



  // assume there's an event coming next.
  int npart;
  ins >> npart;
  getline(ins, ignoreRestOfLine);

  //  cout << " NPART = " << npart << " (rest of line ignored ) " << endl;
  
  TLorentzVector initialNucMomentum;
  TLorentzVector finalNucMomentum;
  TLorentzVector APrimeMomentum;
  TLorentzVector initialElectronMomentum;
  vector<TLorentzVector> outgoingMomenta;
  vector<Charge> outgoingCharges;
  double EBeam;
  Event ev;
  TVector3 labBoost;
  nParticlesPrint=0;
  NDUMMY=4;
  for(int i=0; i< npart; ++i)
    {
      int pid, stat, idummy; 
      ins >> pid >> stat >> idummy >> idummy >> idummy >> idummy;

      if(pid==999) {
        // vertex isn't counted in npart.
        ++npart;
      }


      if((pid==611 || pid==11) && stat==-1) { // it's the incoming "electron"
        getP(ins, initialElectronMomentum);
      }
      if(pid==-623 && stat==-1) { // it's the incoming "nucleon"
        getP(ins, initialNucMomentum);
        labBoost= -initialNucMomentum.BoostVector();
        initialNucMomentum.Boost(labBoost);
        initialElectronMomentum.Boost(labBoost);
        EBeam=initialElectronMomentum.E();
      }

      if(pid==-623 && stat==1) { // it's the outgoing "nucleon"
        ++nParticlesPrint;
        getP(ins, finalNucMomentum);
        finalNucMomentum.Boost(labBoost);
        ev.addP(finalNucMomentum,"4", pid,EBeam);

        ++nParticlesPrint;
        TLorentzVector q = finalNucMomentum-initialNucMomentum;
        //        cout << " q px=" << q.Px() << " py=" << q.Py() << " pz=" << q.Pz();
        ev.addP(q,"5", pid,EBeam);

      }

      if(pid==622) { // it's the Aprime
        ++nParticlesPrint;
        getP(ins, APrimeMomentum);
        APrimeMomentum.Boost(labBoost);
        ev.addP(APrimeMomentum,"0", pid,EBeam);
      }
      
      if(abs(pid)==611 && stat==1) { // it's an outgoing beam electron 
        ++nParticlesPrint;
        TLorentzVector eMomentum;
        getP(ins,eMomentum);
        eMomentum.Boost(labBoost);
        ev.addP(eMomentum,"2",pid,EBeam);   // reversed to agree with the conventions in new code

        outgoingMomenta.push_back(eMomentum);
        outgoingCharges.push_back(charge(-pid));
      }
      

      if(abs(pid)==11 && stat==1) { // it's an electron 
        ++nParticlesPrint;
        TLorentzVector eMomentum;
        getP(ins,eMomentum);
        eMomentum.Boost(labBoost);
        ev.addP(eMomentum,"2",pid,EBeam);


        outgoingMomenta.push_back(eMomentum);
        outgoingCharges.push_back(charge(-pid));

      }

      if(abs(pid)==15) { // it's an electron (actually a tau)
        ++nParticlesPrint;
        TLorentzVector eMomentum;
        getP(ins,eMomentum);
        eMomentum.Boost(labBoost);
        ev.addP(eMomentum,"2",pid,EBeam);   // reversed to agree with the conventions in new code
      } 

      getline(ins, ignoreRestOfLine); // ignore the stuff in the rest of the line.
    }

  //    cout << EBeam << "   " << APrimeMomentum.M() << "  " << APrimeMomentum.E()/EBeam << " " 
  //         << APrimeMomentum.Pt()/APrimeMomentum.Pz()   << "  " << -(finalNucMomentum-initialNucMomentum).M2() << endl;

  //  cout << " i am happy" << endl;

  // if everything succeeded, let calling function know that you're happy.
  vector<double> goodRotAngles=cutfunc.rotatePass(outgoingMomenta,outgoingCharges);
		
  for(int i=0; i<goodRotAngles.size(); ++i) {
      //      cout << outs.str();
      passouts << "0" << setw(5) << nevpass+1 << setw(7) <<  "9999" << endl; 
      ev.print(passouts, goodRotAngles[i]);  //// THIS NEEDS TO BE COMPLETELY REWRITTEN -- CAN'T WRITE EVENT UNTIL WE KNOW WHAT TO PHASE IT BY.
      ++nevpass;
    }
	
  ++nev;

  return true;
}

bool readEventBasic(istream& ins, ostream& passouts, int& nev, int& nevpass, SpectrometerCut& cutfuncBasic)
{
  string nextline;
  string ignoreRestOfLine;

  do {
    getline(ins, nextline);
    //    cout << "NEXT : " << nextline << endl;
    if(ins.eof()) return false; // failed to find an event
  } while (nextline.find("<event>") == string::npos);// continue while test string isn't found.
  
  //  cout << "AN EVENT: " << endl;
  // cout << "0" << setw(5) << nevpass+1 << setw(7) <<  "9999" << endl; 



  // assume there's an event coming next.
  int npart;
  ins >> npart;
  getline(ins, ignoreRestOfLine);

  //  cout << " NPART = " << npart << " (rest of line ignored ) " << endl;
  
  TLorentzVector initialNucMomentum;
  TLorentzVector finalNucMomentum;
  TLorentzVector APrimeMomentum;
  TLorentzVector initialElectronMomentum;
  vector<TLorentzVector> outgoingMomenta;
  vector<Charge> outgoingCharges;
  double EBeam;
  Event ev;
  TVector3 labBoost;
  nParticlesPrint=0;
  NDUMMY=4;
  for(int i=0; i< npart; ++i)
    {
      int pid, stat, idummy; 
      ins >> pid >> stat >> idummy >> idummy >> idummy >> idummy;

      if(pid==999) {
        // vertex isn't counted in npart.
        ++npart;
      }


      if((pid==611 || pid==11) && stat==-1) { // it's the incoming "electron"
        getP(ins, initialElectronMomentum);
      }
      if(pid==-623 && stat==-1) { // it's the incoming "nucleon"
        getP(ins, initialNucMomentum);
        labBoost= -initialNucMomentum.BoostVector();
        initialNucMomentum.Boost(labBoost);
        initialElectronMomentum.Boost(labBoost);
        EBeam=initialElectronMomentum.E();
      }

      if(pid==-623 && stat==1) { // it's the outgoing "nucleon"
        ++nParticlesPrint;
        getP(ins, finalNucMomentum);
        finalNucMomentum.Boost(labBoost);
        ev.addP(finalNucMomentum,"4", pid,EBeam);

        ++nParticlesPrint;
        TLorentzVector q = finalNucMomentum-initialNucMomentum;
        //        cout << " q px=" << q.Px() << " py=" << q.Py() << " pz=" << q.Pz();
        ev.addP(q,"5", pid,EBeam);

      }

      if(pid==622) { // it's the Aprime
        ++nParticlesPrint;
        getP(ins, APrimeMomentum);
        APrimeMomentum.Boost(labBoost);
        ev.addP(APrimeMomentum,"0", pid,EBeam);
      }
      
      if(abs(pid)==611 && stat==1) { // it's an outgoing beam electron 
        ++nParticlesPrint;
        TLorentzVector eMomentum;
        getP(ins,eMomentum);
        eMomentum.Boost(labBoost);
        ev.addP(eMomentum,"2",pid,EBeam);   // reversed to agree with the conventions in new code

        outgoingMomenta.push_back(eMomentum);
        outgoingCharges.push_back(charge(-pid));
      }
      

      if(abs(pid)==11 && stat==1) { // it's an electron 
        ++nParticlesPrint;
        TLorentzVector eMomentum;
        getP(ins,eMomentum);
        eMomentum.Boost(labBoost);
        ev.addP(eMomentum,"2",pid,EBeam);


        outgoingMomenta.push_back(eMomentum);
        outgoingCharges.push_back(charge(-pid));

      }

      if(abs(pid)==15) { // it's an electron (actually a tau)
        ++nParticlesPrint;
        TLorentzVector eMomentum;
        getP(ins,eMomentum);
        eMomentum.Boost(labBoost);
        ev.addP(eMomentum,"2",pid,EBeam);   // reversed to agree with the conventions in new code
      } 

      getline(ins, ignoreRestOfLine); // ignore the stuff in the rest of the line.
    }

  //    cout << EBeam << "   " << APrimeMomentum.M() << "  " << APrimeMomentum.E()/EBeam << " " 
  //         << APrimeMomentum.Pt()/APrimeMomentum.Pz()   << "  " << -(finalNucMomentum-initialNucMomentum).M2() << endl;

  //  cout << " i am happy" << endl;

  // if everything succeeded, let calling function know that you're happy.
  bool pass = cutfuncBasic.pass(outgoingMomenta,outgoingCharges);
		
  if(pass) {
      passouts << "0" << setw(5) << nevpass+1 << setw(7) <<  "9999" << endl; 
      ev.print(passouts, 0.);  //// THIS NEEDS TO BE COMPLETELY REWRITTEN -- CAN'T WRITE EVENT UNTIL WE KNOW WHAT TO PHASE IT BY.
      ++nevpass;
    }
	
  ++nev;

  return true;
}

void printUsage()
{

  cerr << "LHEtoPGScut -c (5degHRS|6degHRS|12p5degHRS) -p0 pcenter [-r nrot|0|-1] [-basic] [-debug] <inputfile> <outputfile>" << endl;

  //FOR DETAILED LOOP  cerr << "LHEtoPGS [-x0 xmin- xmin+] [-x1 xmax- xmax+] [-y ymax- ymax+] [-p0 pmin- pmin+] [-p1 pmax- pmax+] <inputfile> <outputfile>" << endl;
  return;
}

double atod(const string str)
{
  istringstream s(str);
  double temp;

  s >> temp;
  return temp;
}

int main(int argc, char* const argv[] )
{

  int numOpt=0;
  string inputFileName;
  string outputFileName;
  string config="UNKNOWN";
  double p0=100.;
  
  bool useHPGS=false, debug=false;

  Charge qRef=PLUS;
  int nrot=0; // 0: no rotation; <0 : will compute largest possible number of rotations assuming circular symmetry of generated sample; >0: fixed number of rotations.

  bool doBasic=false; // uses the basic code, i.e. SpectrometerCut not more general SpectrometerCutMainz.

  /// SIMPLIFIED OPTIONS LOOP, LOADING ANGULAR CONFIGS
  while(argc-numOpt>1 && argv[1+numOpt][0]=='-') {
    numOpt++;
    if(strcmp(argv[numOpt],"-p0")==0) {
      p0=atod(argv[1+numOpt]);
      numOpt++;
    }
    else if(strcmp(argv[numOpt],"-c")==0) {
      config = argv[1+numOpt];
      numOpt++;
    }
    else if(strcmp(argv[numOpt],"-r")==0) {
      nrot=atod(argv[1+numOpt]);
      numOpt++;
    }
    /*    THIS IS PRESENTLY BROKEN -- AS WRITTEN THE CODE WILL DOUBLE-COUNT EVENTS WHEN QREF=-.  DO NOT EVER USE IT UNTIL CODE IS FIXED -- THERE IS NO REASON WE EVER NEED QREF=-, IN ANY CASE. --NT 2012-07-17*/
	  else if(strcmp(argv[numOpt],"-qref")==0) {
      qRef=charge(atod(argv[1+numOpt]));
      numOpt++;
      } */   
    else if(strcmp(argv[numOpt],"-basic")==0) {
      doBasic=true;
    }    
    else if(strcmp(argv[numOpt],"-debug")==0) {
      SpectrometerCut::debug();
	SpectrometerCutMainz::debug();
    }    
    else  {
      cerr << "Unknown option " << argv[numOpt] << endl;
      printUsage();
      return 1;
    }
  } //End options loop 

  if(argc==3+numOpt) {
    inputFileName=argv[1+numOpt];
    outputFileName=argv[2+numOpt];
  }
  else {
    printUsage();
    return 1;
  }

  /// SYMMETRIC DEFAULT ANGLES: "STANDARD" CONFIGURATIONS:
  double x0A, p0A,  dpA, dxA, dyA; // A=positive
  double x0B, p0B,  dpB, dxB, dyB; // B=negative


  if(config=="5degHRS") {
    x0A=0.089; p0A=p0;  dpA=0.045*p0A; dxA=0.017; dyA=0.037; // A=positive
    x0B=0.089; p0B=p0;  dpB=0.045*p0B; dxB=0.017; dyB=0.037; // B=negative
  }
  else if(config=="6degHRS") {
    x0A=0.105; p0A=p0;  dpA=0.045*p0; dxA=0.017; dyA=0.037; // A=positive
    x0B=0.105; p0B=p0;  dpB=0.045*p0; dxB=0.017; dyB=0.037; // B=negative
  }
  else if(config=="12p5degHRS") {
    x0A=0.218; p0A=p0;  dpA=0.045*p0; dxA=0.027; dyA=0.055; // A=positive
    x0B=0.218; p0B=p0;  dpB=0.045*p0; dxB=0.027; dyB=0.055; // B=negative
  }
  else {
    cerr << " I know nothing about the configuration " << config << " -- EXITING! " << endl;
    assert(false);
  }

  /*  FOLLOWING COMMENTS CORRESPOND TO OLD DEFAULT SETTINGS!
  // OLD default spectrometer params: 
  //USING SETUP 1 (which appears to be used in fig. 6
  //double x0A=0.089, p0A=1.131,  dpA=0.04*p0A, dxA=0.017, dyA=0.037; // A=positive
  //double x0B=0.089, p0B=1.131, dpB=0.04*p0B, dxB=0.017, dyB=0.037; // B=negative

  //  true spec. settings
  //  double x0A=22.8*pi/180., p0A=0.338,  dpA=0.1*p0A,   dxA=0.10, dyA=0.07; // A=positive
  //  double x0B=15.2*pi/180., p0B=0.4699, dpB=0.075*p0B, dxB=0.02, dyB=0.07; // B=negative

  //double xmin_p=0.0, xmax_p=100.0, ymax_p=100.0, pmin_p=0.0, pmax_p=100.0;
  //double xmin_m=0.0, xmax_m=100.0, ymax_m=100.0, pmin_m=0.0, pmax_m=100.0;

  DETAILED OPTIONS LOOP -- WE WANT TO SIMPLIFY THIS!
  while(argc-numOpt>1 && argv[1+numOpt][0]=='-') {
    numOpt++;
    if(strcmp(argv[numOpt],"-p0")==0) {
      pmin_m=atod(argv[1+numOpt]);
      pmin_p=atod(argv[2+numOpt]);
      numOpt+=2;
    }
    else if(strcmp(argv[numOpt],"-p1")==0) {
      pmax_m=atod(argv[1+numOpt]);
      pmax_p=atod(argv[2+numOpt]);
      numOpt+=2;
    }
    else if(strcmp(argv[numOpt],"-y")==0) {
      ymax_m=atod(argv[1+numOpt]);
      ymax_p=atod(argv[2+numOpt]);
      numOpt+=2;
    }
    else if(strcmp(argv[numOpt],"-x0")==0) {
      xmin_m=atod(argv[1+numOpt]);
      xmin_p=atod(argv[2+numOpt]);
      numOpt+=2;
    }
    else if(strcmp(argv[numOpt],"-x1")==0) {
      xmax_m=atod(argv[1+numOpt]);
      xmax_p=atod(argv[2+numOpt]);
      numOpt+=2;
    }
    else if(strcmp(argv[numOpt],"-r")==0) {
      int i=atoi(argv[1+numOpt]);
      if(i==0){
	doRot=false;
      }
      else qRef=charge(i);
      numOpt++;
    }
    
    else  {
      cerr << "Unknown option " << argv[numOpt] << endl;
      printUsage();
      return 1;
    }
  } //End options loop 
  */

  double xmin_p=x0A-dxA, xmax_p=x0A+dxA, ymax_p=dyA, pmin_p=p0A-dpA, pmax_p=p0A+dpA;
  double xmin_m=x0B-dxB, xmax_m=x0B+dxB, ymax_m=dyB, pmin_m=p0B-dpB, pmax_m=p0B+dpB;

  // set up both a "basic"style spectrometer and a Mainz-style fancy spectrometer so we can check both.
  SpectrometerCutMainz mainzSpec(qRef, nrot);
  mainzSpec.setSpec(MINUS,xmin_m,xmax_m,ymax_m,pmin_m,pmax_m);
  mainzSpec.setSpec(PLUS,xmin_p,xmax_p,ymax_p,pmin_p,pmax_p);

  SpectrometerCut basicSpec(xmin_m,xmax_m,ymax_m,pmin_m,pmax_m);

  //  mainzSpec.doRotate(doRot);

  ifstream ins(inputFileName.c_str());
  ofstream outs(outputFileName.c_str());
  
  if(!ins.good()){
    cerr << " ERROR: COULD NOT OPEN " << inputFileName << endl;
    return 0;
  }

  
  if(!outs.good()){
    cerr << " ERROR: COULD NOT OPEN " << outputFileName << endl;
    return 0;
  }

  // record provenance
  outs << "# PROVENANCE: " ;
  for(int i=0; i<argc; ++i)
    {
      outs << argv[i] << "  ";
    }
  outs << endl;

  cout << "  xmin_p=" << xmin_p << "  xmax_p=" << xmax_p << "  ymax_p=" << ymax_p << "  pmin_p=" << pmin_p << "  pmax_p=" << pmax_p << endl;
  cout << "  xmin_m=" << xmin_m << "  xmax_m=" << xmax_m << "  ymax_m=" << ymax_m << "  pmin_m=" << pmin_m << "  pmax_m=" << pmax_m << endl;
  cout << "  config=" << config << "  rotation setting=" << nrot << "  qref=" << qRef  << "  doBasic=" << doBasic << endl;

  outs << "#   xmin_p=" << xmin_p << "  xmax_p=" << xmax_p << "  ymax_p=" << ymax_p << "  pmin_p=" << pmin_p << "  pmax_p=" << pmax_p << endl;
  outs << "#   xmin_m=" << xmin_m << "  xmax_m=" << xmax_m << "  ymax_m=" << ymax_m << "  pmin_m=" << pmin_m << "  pmax_m=" << pmax_m << endl;
  outs << "#   config=" << config << "  rotation setting=" << nrot << "  qref=" << qRef << "  doBasic=" << doBasic << endl;


  // end provenance.

  // the run loop

  int nev=0;
  int nevpass=0;

  outs << "  #  typ      eta    phi      pt    jmas  ntrk  btag   had/em  dum1  dum2" << endl;

  if(doBasic){
    while( nev < 5000000 && readEventBasic(ins, outs, nev, nevpass, basicSpec)){
      if(nev%5000==0) cout << " read " << nev << " events (" << nevpass << " pass: " << 100.0*nevpass/nev << "%)" << endl;
    }
  } else {
    while( nev < 5000000 && readEvent(ins, outs, nev, nevpass, mainzSpec)){
      if(nev%5000==0) cout << " read " << nev << " events (" << nevpass << " pass: " << 100.0*nevpass/nev << "%)" << endl;
    }
  }


  // record statistics
  outs << "# STATISTICS: " << nevpass << " EVENTS WRITTEN OUT OF " << nev << " EVENTS READ, " << mainzSpec.ncopy() << " copies." << endl;
  outs << "# STATISTICS:      TRUE PASS PERCENTAGE " << 100.0*nevpass/nev / max(mainzSpec.ncopy(),1) << "%" << endl;
  cout << "# STATISTICS: " << nevpass << " EVENTS WRITTEN OUT OF " << nev << " EVENTS READ, " << mainzSpec.ncopy() << " copies." << endl;
  cout << "# STATISTICS:      TRUE PASS PERCENTAGE " << 100.0*nevpass/nev / max(mainzSpec.ncopy(),1) << "%" << endl;
  //ADD "TRUE ACCEPTANCE PER
  // end statistics
  

}


//////////////////////////////////////////////
//////////     USED CLASS FUNCTIONS //////////
//enum Charge {PLUS=0,MINUS=1};

/*class SpectrometerCutMainz
{
public:
  SpectrometerCutMainz();
  void setSpec(Charge q, double xmin,double xmax, double ymax, double pmin, double pmax);
  bool rotatePass(vector<TLorentzVector> momenta, vector<int> signs, Charge chargeToRotate);
  bool pass(vector<TLorentzVector> momenta, vector<int> signs);

  private:
  double xmin_[2], xmax_[2], ymax_[2], pmin_[2], pmax_[2];
  const static bool DEBUG;
};

*/
 
SpectrometerCutMainz::SpectrometerCutMainz(Charge qref, int ncopy=-1):
  pmax_(2,0.),  // Ensures that nothing will work until initialized
  pmin_(2,10.),
  xmax_(2,0.),
  xmin_(2,0.),
  ymax_(2,100.),
  qRef_(qref),
  ncopy_(ncopy), // added 2012-07-17 : user control over max number of rotations
  phi_max_(pi) {}

string SpectrometerCutMainz::setSpec(Charge q, double xmin,double xmax, double ymax, double pmin, double pmax)
{
  int qi(q);
  xmin_[qi]=xmin;
  xmax_[qi]=xmax;
  assert(xmax>xmin);
  ymax_[qi]=ymax;
  pmin_[qi]=pmin;
  pmax_[qi]=pmax;
  assert(pmin<pmax);
  ostringstream report;

  if(q==qRef_) {
    phi_max_ = (xmin==0? 2*pi : atan(ymax/xmin));
    int ncopyMax = floor(2*pi/(2*phi_max_));
    if(ncopy()<0 || ncopy() > ncopyMax ) ncopy_=ncopyMax; // added 2012-07-17 : ncopy<0 => automatic calculation; never use more than ncopyMax copies.
    report << " NCopy used = " << ncopy() << endl;
  }

  return report.str();
}

vector<double> SpectrometerCutMainz::rotatePass(vector<TLorentzVector> trueMomenta, vector<Charge> signs) const
{
  /* Procedure:  modified 2012-07-17
     1) Interpret x_min/max as theta_min/max, determine central theta0
     2) use inner theta to derive phi_max = arctan(y_max/theta_min).  
     
                            /|   ^
                         /   |   ymax
                      /phi)  |   |
                    ---------|   v
                    <-theta0->
     3) Define ncopy = floor(2pi/(2 phi_max)) = # of full copies of 2*phi_max regions you can fit in an annulus.   [up to here is done in setSpec routine]
     4) Define phishift=2pi/ncopy to be the amount of each rotations; redefine phi range so that it runs from -phishift/2 to 2*pi - phishift/2 (i.e. map the range [2*pi-ps/2,2pi] to [-ps/2, 0])
     5) Loop over particles of desired charge, let p be its momentum
        a] Find nshift = floor[phi(p) / phishift + 0.5)]
	   This is 0 if phi is between -phishift/2 and phishift/2 , etc.  With this definition, the whole acceptance is within the n=0, and its mirror image wtihin the n=1.
	b] If nshift < ncopy, check whether event passes after rotation 
           by nshift*phishift.  
     **Note: using theta_min in the determination of phi_max introduces
        'dead space' between copies of the spectrometer, but allows us to
        use rectangular acceptance
  */
  vector<double> passRots(0);

  if(ncopy_==0) {
    if(pass(trueMomenta, signs, qRef_))
      passRots.push_back(0.);  
  }
  else {
    for(int i=0; i<signs.size(); ++i) {
      if(signs[i]==qRef_) { // try to rotate this guy into 'nominal spectrometer'
	double shift=2*pi/ncopy();  // changed 2012-07-17
	double phi=trueMomenta[i].Phi(); 
	if(phi>2*pi-shift/2) phi-=2*pi; // changed 2012-07-17
	int nshift = floor(phi/shift + 0.5); // changed 2012-07-17
	if(nshift < ncopy_) { // rotate by nshift
	  vector<TLorentzVector> rotMomenta(trueMomenta);
	  double rotAngle=-nshift*shift;
	  for(int j=0; j<rotMomenta.size(); ++j) {	 
	    rotMomenta[j].RotateZ(rotAngle);
	  }
	  assert(fabs(rotMomenta[i].Phi())<shift/2);  // the rotated momentum i should be w.in shift/2 whenever n_shift < n_copy.
	  if(pass(rotMomenta, signs, qRef_))
	    passRots.push_back(rotAngle);
	}
      }
    }
  }
  return passRots;
}

bool SpectrometerCutMainz::pass(vector<TLorentzVector> momenta, vector<Charge> signs, Charge qref) const
{
  int pass[2]={0,0};
  
  int npart=momenta.size();

  if(DEBUG) cout << " hi npart = " << npart << endl;

  if(npart<2) return false;
  
  for(int i=0; i<npart; ++i)
    {
      double thx=atan(momenta[i].Px()/momenta[i].Pz());
      double thy=atan(momenta[i].Py()/momenta[i].Pz());
      double P = momenta[i].P();
      Charge q = signs[i];
      int sgnX=sgn(q)*sgn(qRef_);

      if(DEBUG) cout << "(" << thx << "  " ;
      if(DEBUG) cout << thy << ")  " ;



      if(P>pmin(q) && P < pmax(q) && abs(thy) < ymax(q)) {
        if(DEBUG) cout << "ok" << endl;
        if(sgnX*thx>xmin(q) && sgnX*thx<xmax(q)) {       
          pass[int(q)]++;
          if(DEBUG) cout << " pass[" << q<< "]=" << pass[int(q)] << endl;
        }
	/*	else {
	  cout << " thetaX fail (" << sgn(q) << ") " << thx << endl;
	  }*/
      }
      /*      else if(P>0.7){
	if(P<pmin(q))
	  cout << "   low  P   " << P << "  ("<< sgn(q) << ")" << endl;
	if(P>pmax(q))
	  cout << "   high P   " << P << "  ("<< sgn(q) << ")" << endl;
	if(abs(thy) > ymax(q))
	  cout << "   thy fail " << thy << " ("<< sgn(q) << ")" << endl;
	  }*/
    }
  /*
  if(pass[PLUS] ==0) cout << " NO PLUS" << endl;
    
  if(pass[MINUS] ==0) cout << " NO MINUS" << endl;*/

  if(DEBUG) cout << endl;
  if(DEBUG) if(pass[PLUS]>0 || pass[MINUS] > 0) cout << "passP= " << pass[PLUS] << "passM= " << pass[MINUS] << endl;

  if(pass[int(PLUS)]*pass[(MINUS)]>=1) {
    if(DEBUG) cout << "PASS!";
    return true;
  }
  else return false;
}




//////////////////////////////////////////////
//////////  BASIC SPECTROMETER CLASS FUNCTIONS -- SYMMETRIC SPECS, NO ROTATIONS//////////

 
SpectrometerCut::SpectrometerCut(double xmin,double xmax, double ymax, double pmin, double pmax):
  xmin_(xmin),
  xmax_(xmax),
  ymax_(ymax),
  pmin_(pmin),
  pmax_(pmax) {return;}

bool SpectrometerCut::pass(vector<TLorentzVector> momenta, vector<Charge> signs)
{
  // revised July 17, 2012: to make sure that events where an e+ and e- are both within geometric acceptance of a spectrometer are *not* rejected, but that opposite signs must land in different spectrometers.  This is done by assigning each spectrometer a charge that it looks for, but also including "flip" solutions.  Note that similar behavior is (I think) already accounted for in the SpectrometerCutMainz algorithm, which uses rotational symmetry -- but this code is appropriate for cases where we can only assume reflection symmetry about the y axis (should get equivalent results from including 2 rotated configurations). 
  // revised July 17, 2012: added pmax argument
  // look for pairs that are good, OR reflections that are good
  //   (for non-overlapping spectrometers and only one e+, there will never be 2 solutions)

  int passP=0;
  int passM=0;
  int passEV=0;

  int passPflip=0;
  int passMflip=0;
  int passEVflip=0;

  
  int npart=momenta.size();

  if(DEBUG) cout << " hi npart = " << npart << endl;

  if(npart<2) return false;
  


  for(int i=0; i<npart; ++i)
    {
      double thx=atan(momenta[i].Px()/momenta[i].Pz());
      double thy=atan(momenta[i].Py()/momenta[i].Pz());
      double P = momenta[i].P();
      Charge q =signs[i];

      if(DEBUG) cout << "(" << thx << "  " ;
      if(DEBUG) cout << thy << ")  " ;
      

      if(P>pmin_ && P < pmax_ && abs(thy) < ymax_) {
        if(DEBUG) cout << "p and y ok" << endl;
	// look for good particles: x+>0, x-<0
        if(q==PLUS && thx>xmin_ && thx<xmax_) {       
          ++passP;
          if(DEBUG) cout << " passP " << passP << endl;
        }
        if(q==MINUS && thx< -xmin_ && thx> - xmax_) {
          ++passM;
          if(DEBUG) cout << " passM " << passM << endl;
        }
	// look for particles that are good upon reflection: x+<0, x->0
        if(q==MINUS && thx>xmin_ && thx<xmax_) {       
          ++passMflip;
          if(DEBUG) cout << " passPflip " << passP << endl;
        }
        if(q==PLUS && thx< -xmin_ && thx> - xmax_) {
          ++passPflip;
          if(DEBUG) cout << " passMflip " << passM << endl;
        }
      }
    }

  if(DEBUG) cout << endl;
  if(DEBUG) if(passP>0 || passP>0 || passP>0 || passMflip>0) 
	      cout << "P:" << passP << " M:" << passM << "  /   Pf:" << passPflip << " Mf:" << passMflip << endl;
  if(passP == 1 && passM == 1) {
    passEV=1;
    if(DEBUG) cout << "passEV" << endl;
  }
  if(passPflip == 1 && passMflip == 1) {
    passEVflip=1;
    if(DEBUG) cout << "passEVflip" << endl;
  }

  if(passEV && passEVflip) { 
    cout << "PASSED BOTH WAYS!! IMPOSSIBLE!!" << endl;
    assert(false);
  }
  if(passEV || passEVflip)
    return true;
  else return false;
}

